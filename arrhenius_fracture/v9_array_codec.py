"""Atomic, hash-verified production array generations for v9."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil
import uuid

import numpy as np

from .v9_transactional import canonical_json


ARRAY_SCHEMA_VERSION = "V9_ARRAY_GENERATION_1"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


class ArrayGenerationError(RuntimeError):
    pass


class AtomicArrayGenerationStore:
    def __init__(self, root: Path):
        self.root = Path(root)

    def write(self, arrays, metadata, summary) -> str:
        self.root.mkdir(parents=True, exist_ok=True)
        generation = "generation_" + uuid.uuid4().hex
        temporary = self.root / ("." + generation + ".tmp")
        final = self.root / generation
        temporary.mkdir()
        try:
            normalized = {}
            schema = {}
            for name, value in sorted(arrays.items()):
                if not name or "/" in name or name == "metadata_json":
                    raise ArrayGenerationError(f"invalid array name: {name!r}")
                array = np.asarray(value)
                if array.dtype.hasobject:
                    raise ArrayGenerationError(f"object/pickle array forbidden: {name}")
                normalized[name] = np.ascontiguousarray(array)
                schema[name] = {"dtype": str(array.dtype), "shape": list(array.shape)}
            arrays_path = temporary / "state_arrays.npz"
            with arrays_path.open("wb") as stream:
                np.savez_compressed(stream, **normalized)
                stream.flush()
                os.fsync(stream.fileno())
            metadata_payload = metadata | {"schema_version": ARRAY_SCHEMA_VERSION, "array_schema": schema}
            self._write_file(temporary / "state_metadata.json", canonical_json(metadata_payload) + "\n")
            self._write_file(temporary / "summary.json", canonical_json(summary) + "\n")
            components = {}
            for name in ("state_arrays.npz", "state_metadata.json", "summary.json"):
                path = temporary / name
                components[name] = {"sha256": _sha256(path), "bytes": path.stat().st_size}
            manifest = {
                "schema_version": ARRAY_SCHEMA_VERSION,
                "generation": generation,
                "components": components,
                "array_schema": schema,
            }
            self._write_file(temporary / "manifest.json", canonical_json(manifest) + "\n")
            _fsync_directory(temporary)
            os.replace(temporary, final)
            _fsync_directory(self.root)
            self._write_file(self.root / ".active.tmp", canonical_json({"generation": generation}) + "\n")
            os.replace(self.root / ".active.tmp", self.root / "ACTIVE.json")
            _fsync_directory(self.root)
            return generation
        except Exception:
            if temporary.exists():
                shutil.rmtree(temporary)
            raise

    def load(self, generation: str | None = None):
        if generation is None:
            active = self.root / "ACTIVE.json"
            if not active.is_file():
                raise ArrayGenerationError("no atomically activated generation")
            generation = json.loads(active.read_text())["generation"]
        if not isinstance(generation, str) or not generation.startswith("generation_") or not generation[11:].isalnum():
            raise ArrayGenerationError("invalid generation identity")
        directory = self.root / generation
        manifest_path = directory / "manifest.json"
        if not manifest_path.is_file():
            raise ArrayGenerationError("incomplete generation: manifest missing")
        manifest = json.loads(manifest_path.read_text())
        if manifest.get("schema_version") != ARRAY_SCHEMA_VERSION or manifest.get("generation") != generation:
            raise ArrayGenerationError("generation schema/identity mismatch")
        for name, record in manifest["components"].items():
            if name not in {"state_arrays.npz", "state_metadata.json", "summary.json"}:
                raise ArrayGenerationError(f"unknown generation component: {name}")
            path = directory / name
            if not path.is_file() or _sha256(path) != record["sha256"] or path.stat().st_size != record["bytes"]:
                raise ArrayGenerationError(f"component hash/size mismatch: {name}")
        metadata = json.loads((directory / "state_metadata.json").read_text())
        if metadata.get("array_schema") != manifest.get("array_schema"):
            raise ArrayGenerationError("array schema metadata mismatch")
        with np.load(directory / "state_arrays.npz", allow_pickle=False) as archive:
            if set(archive.files) != set(manifest["array_schema"]):
                raise ArrayGenerationError("array key schema mismatch")
            arrays = {name: np.asarray(archive[name]).copy() for name in archive.files}
        for name, array in arrays.items():
            expected = manifest["array_schema"][name]
            if str(array.dtype) != expected["dtype"] or list(array.shape) != expected["shape"]:
                raise ArrayGenerationError(f"array dtype/shape mismatch: {name}")
        return arrays, metadata, json.loads((directory / "summary.json").read_text()), manifest

    @staticmethod
    def _write_file(path: Path, payload: str) -> None:
        with path.open("w") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
