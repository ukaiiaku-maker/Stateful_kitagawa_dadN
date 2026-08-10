"""Fail-closed bridge to the canonical audited v9.13 four-class rows."""
from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path


SCHEMA = "V9_CANONICAL_FOUR_CLASS_REGISTRY_BRIDGE_1"
EXPECTED = {
    "v913_paper_peak01_0242980_persistent_sites": ("peak", "v913_zeroD_sobol_0242980"),
    "v913_paper_dbtt01_0202500_persistent_sites": ("DBTT", "v913_zeroD_sobol_0202500"),
    "v913_paper_weakT01_0129902_persistent_sites": ("weakT", "v913_zeroD_sobol_0129902"),
    "v913_paper_ceramic01_0077080_persistent_sites": ("ceramic", "v913_zeroD_sobol_0077080"),
}
FORBIDDEN_ALTERNATES = {
    "v913_paper_weakT01_0257068_persistent_sites",
    "v913_paper_ceramic01_0189364_persistent_sites",
}
EXPECTED_HASHES = {
    "parameter_registry_v9111.py": "6946759c5786344a7111b9e25ce8b2513e5335d78fc590ba860fc730bb7d83a7",
    "v10_2_27_v913_four_class_paper_registry.csv": "4ba723c80abcfdd7101cee0afaa9b4104eebf2a8e847892528128664829c7760",
    "v10_2_27_v913_four_class_paper_selection.json": "3267ac05595cee0800a81a6300e80cb0e71b216067c35cc7a50d82ac4c70e4a2",
}
ACTIVE_PARAMETER_FINGERPRINT = "eeb992cb8f956351833ff148d7611c8fc9e210f2a02eb44a00caa85558001fc5"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _paths(source_root: str | Path) -> dict[str, Path]:
    root = Path(source_root).expanduser().resolve()
    materials = root / "arrhenius_fracture" / "data" / "materials"
    return {
        "parameter_registry_v9111.py": root / "arrhenius_fracture" / "parameter_registry_v9111.py",
        "v10_2_27_v913_four_class_paper_registry.csv": materials / "v10_2_27_v913_four_class_paper_registry.csv",
        "v10_2_27_v913_four_class_paper_selection.json": materials / "v10_2_27_v913_four_class_paper_selection.json",
    }


def select_canonical_option(option_id: str, source_root: str | Path):
    if option_id in FORBIDDEN_ALTERNATES or option_id not in EXPECTED:
        raise ValueError(f"not a canonical stable-crack-birth four-class option: {option_id!r}")
    paths = _paths(source_root)
    actual = {name: _sha256(path) for name, path in paths.items()}
    if actual != EXPECTED_HASHES:
        raise RuntimeError(f"audited four-class source hash mismatch: {actual!r}")
    selector_path = paths["parameter_registry_v9111.py"]
    spec = importlib.util.spec_from_file_location("_audited_parameter_registry_v9111", selector_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load audited selector: {selector_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    class_label, candidate_id = EXPECTED[option_id]
    registry_name = "v10_2_27_v913_four_class_paper_registry.csv"
    selection_name = "v10_2_27_v913_four_class_paper_selection.json"
    selection = json.loads(paths[selection_name].read_text())
    if selection.get("source_active_parameter_fingerprint_sha256") != ACTIVE_PARAMETER_FINGERPRINT:
        raise RuntimeError("canonical active-parameter fingerprint mismatch")
    if selection.get("canonical_option_order") != list(EXPECTED):
        raise RuntimeError("canonical final four-class option order mismatch")
    if set(selection.get("canonical_option_order", ())) & FORBIDDEN_ALTERNATES:
        raise RuntimeError("rejected transfer candidate entered canonical mapping")
    selected = module.select_option(option_id, paths[registry_name], canonical_stage3_only=False)
    if selected.candidate_id != candidate_id:
        raise RuntimeError(f"canonical candidate mismatch for {option_id}")
    candidates = selection.get("primary_candidates", [])
    if sum(row.get("option_key") == option_id for row in candidates) != 1:
        raise RuntimeError(f"selection manifest does not uniquely authorize {option_id}")
    return selected, {
        "schema": SCHEMA,
        "class_label": class_label,
        "option_id": option_id,
        "candidate_id": candidate_id,
        "source_repository": str(Path(source_root).expanduser().resolve()),
        "selector_path": str(selector_path),
        "registry_path": str(paths[registry_name]),
        "selection_path": str(paths[selection_name]),
        "source_hashes": actual,
        "active_parameter_fingerprint_sha256": ACTIVE_PARAMETER_FINGERPRINT,
        "active_parameter_row_sha256": hashlib.sha256(json.dumps(
            {key: selected.row[key] for key in selection["active_parameter_fields"]},
            sort_keys=True, separators=(",", ":"),
        ).encode()).hexdigest(),
        "exact_registry_row": dict(selected.row),
    }
