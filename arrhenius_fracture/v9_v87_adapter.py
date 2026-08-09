"""Fail-closed conversion of authoritative v8.7 checkpoints to v9 state."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Mapping

import numpy as np

from .sn_pd2d_stateful_v8_7_generalized_features import (
    MODEL_ID, SOURCE_SHA256, _CHECKPOINT_VERSION,
)
from .v9_solver import BASELINE_MECHANISM, V9SolverRequest
from .v9_transactional import PersistentClock, V9StateCapsule, canonical_json


ADAPTER_VERSION = "V87_CHECKPOINT_TO_V9_CAPSULE_1"

BASE_ARRAYS = {
    "mesh_nodes", "root_xy", "ep_gp", "rho_gp", "epsp_acc_gp", "u", "last_residual",
}
PD_POINT_ARRAYS = {
    "available", "embryo", "stable", "inactive", "candidate_sites", "available_sites",
    "embryo_sites", "stable_sites", "inactive_sites", "born_sites_cumulative",
    "healed_sites_cumulative", "birth_cumulative_hazard", "delivery_memory", "completion",
    "growth", "crack_normal_c2", "crack_normal_s2", "crack_orientation_weight",
    "primary_seed_rejected", "healed_cumulative", "born_cumulative",
}
PD_SITE_ARRAYS = {
    "site_node_index", "site_status", "site_birth_threshold", "site_birth_cycle", "site_stable_cycle",
}
PD_BOND_ARRAYS = {
    "bond_damage", "active_front_bonds", "front_backbone_bonds", "front_wake_bonds", "front_process_bonds",
}
PD_VECTOR_ARRAYS = {"active_front_contact_xy", "active_front_tip_xy"}
PD_PATH_ARRAYS = {"active_front_path_xy"}
REQUIRED_PD_SCALARS = {
    "primary_seed_node", "primary_seed_stall_updates", "primary_seed_reselections",
    "primary_seed_last_progress", "primary_seed_selected_cycles", "active_front",
    "active_front_normal_c2", "active_front_normal_s2", "active_front_length_m",
    "front_candidate_mode", "front_eligible_preferred", "front_eligible_fallback",
    "front_stall_updates", "front_last_advance_cycles", "front_max_link_rate_per_cycle",
    "diffuse_bad_updates", "cycles_first_embryo", "cycles_first_stable",
    "cycles_first_expected_embryo", "cycles_first_expected_stable", "cycles_first_softening",
    "cycles_root_connected", "cycles_two_horizon_crack", "cycles_front_capture",
    "cycles_primary_seed_reselected", "cycles_precapture_stalled", "cycles_front_stalled",
    "cycles_diffuse_abort", "cycles_connected",
}
EXPECTED_NPZ_KEYS = BASE_ARRAYS | {"metadata_json"} | {
    f"pd__{name}" for name in PD_POINT_ARRAYS | PD_SITE_ARRAYS | PD_BOND_ARRAYS | PD_VECTOR_ARRAYS | PD_PATH_ARRAYS
}
REQUIRED_METADATA = {
    "checkpoint_version", "model_id", "source_sha256", "signature", "next_block", "cycles",
    "Wp_total", "pd_scalars", "candidate_rng_state", "event_rng_state", "rows",
}


class V87CheckpointValidationError(RuntimeError):
    pass


@dataclass(frozen=True)
class AdaptedV87Checkpoint:
    capsule: V9StateCapsule
    arrays: Mapping[str, np.ndarray]
    metadata: Mapping[str, object]
    checkpoint_sha256: str
    historical_fields: tuple[str, ...]
    deterministic_reconstructed_structure: tuple[str, ...]
    new_v9_controller_state: tuple[str, ...]


def _rng_digest(state: object) -> str:
    return hashlib.sha256(canonical_json(state).encode()).hexdigest()


def _validate_rng(name: str, state: object) -> None:
    if not isinstance(state, dict) or not isinstance(state.get("bit_generator"), str):
        raise V87CheckpointValidationError(f"missing or invalid {name}")
    if "state" not in state:
        raise V87CheckpointValidationError(f"missing {name} internal state")


def _validate_shapes(arrays: Mapping[str, np.ndarray]) -> None:
    mesh = arrays["mesh_nodes"]
    if mesh.ndim != 2 or mesh.shape[1] != 2 or not np.all(np.isfinite(mesh)):
        raise V87CheckpointValidationError("inconsistent geometry: mesh_nodes must be finite (n,2)")
    if arrays["root_xy"].shape != (2,) or not np.all(np.isfinite(arrays["root_xy"])):
        raise V87CheckpointValidationError("inconsistent geometry: root_xy")
    if arrays["ep_gp"].ndim != 2 or arrays["ep_gp"].shape[0] != 3:
        raise V87CheckpointValidationError("ep_gp must have shape (3,n_element)")
    ne = arrays["ep_gp"].shape[1]
    for name in ("rho_gp", "epsp_acc_gp"):
        if arrays[name].shape != (ne,):
            raise V87CheckpointValidationError(f"shape mismatch: {name}")
    nn = len(mesh)
    if arrays["u"].shape != (2 * nn,) or arrays["last_residual"].shape != (nn,):
        raise V87CheckpointValidationError("shape mismatch: displacement/residual versus mesh")
    npoint = len(arrays["pd__delivery_memory"])
    for name in PD_POINT_ARRAYS:
        if arrays[f"pd__{name}"].shape != (npoint,):
            raise V87CheckpointValidationError(f"shape mismatch: pd__{name}")
    nsite = len(arrays["pd__site_node_index"])
    for name in PD_SITE_ARRAYS:
        if arrays[f"pd__{name}"].shape != (nsite,):
            raise V87CheckpointValidationError(f"shape mismatch: pd__{name}")
    nodes = arrays["pd__site_node_index"]
    if np.any(nodes < 0) or np.any(nodes >= npoint):
        raise V87CheckpointValidationError("site_node_index outside PD point range")
    status = arrays["pd__site_status"]
    if np.any(status > 3):
        raise V87CheckpointValidationError("unknown candidate site status")
    thresholds = arrays["pd__site_birth_threshold"]
    if np.any(~np.isfinite(thresholds)) or np.any(thresholds <= 0.0):
        raise V87CheckpointValidationError("missing/invalid persistent thresholds")
    if np.any(~np.isfinite(arrays["pd__birth_cumulative_hazard"])) or np.any(arrays["pd__birth_cumulative_hazard"] < 0.0):
        raise V87CheckpointValidationError("invalid cumulative per-site-node hazard")
    nbond = len(arrays["pd__bond_damage"])
    for name in PD_BOND_ARRAYS:
        if arrays[f"pd__{name}"].shape != (nbond,):
            raise V87CheckpointValidationError(f"shape mismatch: pd__{name}")
    for name in PD_VECTOR_ARRAYS:
        if arrays[f"pd__{name}"].shape != (2,):
            raise V87CheckpointValidationError(f"shape mismatch: pd__{name}")
    if arrays["pd__active_front_path_xy"].ndim != 2 or arrays["pd__active_front_path_xy"].shape[1] != 2:
        raise V87CheckpointValidationError("shape mismatch: pd__active_front_path_xy")


def adapt_v87_checkpoint(
    checkpoint: Path,
    *,
    immutable_signature: Mapping[str, object],
    expected_cycle: float | None = None,
    expected_next_block: int | None = None,
) -> AdaptedV87Checkpoint:
    checkpoint = Path(checkpoint)
    if not checkpoint.is_file():
        raise V87CheckpointValidationError(f"checkpoint missing: {checkpoint}")
    digest = hashlib.sha256(checkpoint.read_bytes()).hexdigest()
    try:
        with np.load(checkpoint, allow_pickle=False) as archive:
            keys = set(archive.files)
            missing = EXPECTED_NPZ_KEYS - keys
            unknown = keys - EXPECTED_NPZ_KEYS
            if missing:
                raise V87CheckpointValidationError("missing required arrays: " + ", ".join(sorted(missing)))
            if unknown:
                raise V87CheckpointValidationError("unknown checkpoint fields: " + ", ".join(sorted(unknown)))
            metadata = json.loads(str(archive["metadata_json"].item()))
            arrays = {name: np.asarray(archive[name]).copy() for name in keys if name != "metadata_json"}
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
        raise V87CheckpointValidationError(f"invalid checkpoint container: {exc}") from exc

    missing_metadata = REQUIRED_METADATA - set(metadata)
    unknown_metadata = set(metadata) - REQUIRED_METADATA
    if missing_metadata:
        raise V87CheckpointValidationError("missing metadata: " + ", ".join(sorted(missing_metadata)))
    if unknown_metadata:
        raise V87CheckpointValidationError("unknown metadata fields: " + ", ".join(sorted(unknown_metadata)))
    if metadata["checkpoint_version"] != _CHECKPOINT_VERSION:
        raise V87CheckpointValidationError("wrong checkpoint version")
    if metadata["model_id"] != MODEL_ID:
        raise V87CheckpointValidationError("wrong model")
    if metadata["source_sha256"] != SOURCE_SHA256:
        raise V87CheckpointValidationError("wrong v8.7 source hashes")
    if metadata["signature"] != dict(immutable_signature):
        raise V87CheckpointValidationError("complete run signature mismatch")
    if expected_cycle is not None and float(metadata["cycles"]) != float(expected_cycle):
        raise V87CheckpointValidationError("inconsistent cycle coordinate")
    if expected_next_block is not None and int(metadata["next_block"]) != int(expected_next_block):
        raise V87CheckpointValidationError("inconsistent next_block")
    if not math.isfinite(float(metadata["cycles"])) or float(metadata["cycles"]) < 0.0:
        raise V87CheckpointValidationError("invalid cycle coordinate")
    pd_scalars = metadata["pd_scalars"]
    if not isinstance(pd_scalars, dict) or set(pd_scalars) != REQUIRED_PD_SCALARS:
        missing = REQUIRED_PD_SCALARS - set(pd_scalars or {})
        unknown = set(pd_scalars or {}) - REQUIRED_PD_SCALARS
        raise V87CheckpointValidationError(f"PD scalar schema mismatch; missing={sorted(missing)} unknown={sorted(unknown)}")
    _validate_rng("candidate_rng_state", metadata["candidate_rng_state"])
    _validate_rng("event_rng_state", metadata["event_rng_state"])
    _validate_shapes(arrays)

    site_nodes = arrays["pd__site_node_index"].astype(np.int64, copy=False)
    statuses = arrays["pd__site_status"].astype(np.uint8, copy=False)
    thresholds = arrays["pd__site_birth_threshold"].astype(float, copy=False)
    node_hazard = arrays["pd__birth_cumulative_hazard"].astype(float, copy=False)
    clocks = tuple(
        PersistentClock(
            identity=f"v87-site-{index}", event_kind="birth", threshold=float(threshold),
            cumulative_hazard=float(node_hazard[node]), active=int(status) == 0,
        )
        for index, (node, status, threshold) in enumerate(zip(site_nodes, statuses, thresholds))
    )
    request = V9SolverRequest(
        mechanism=BASELINE_MECHANISM,
        temperature_K=float(immutable_signature["T"]),
        frequency_Hz=float(immutable_signature["frequency_Hz"]),
        source_hashes=tuple(sorted(SOURCE_SHA256.items())),
        configuration_version=ADAPTER_VERSION,
    )
    rng_states = {
        "candidate": metadata["candidate_rng_state"],
        "event": metadata["event_rng_state"],
    }
    rng_digests = {name: _rng_digest(value) for name, value in rng_states.items()}
    array_schema = {name: {"dtype": str(value.dtype), "shape": list(value.shape)} for name, value in sorted(arrays.items())}
    capsule = V9StateCapsule(
        cycle=float(metadata["cycles"]),
        fem_plastic_state_json=canonical_json({"array_keys": sorted(BASE_ARRAYS), "Wp_total": metadata["Wp_total"]}),
        rho_backstress_shielding_json=canonical_json({"array_keys": ["rho_gp"], "source": "historical_v87"}),
        delivery_memory_json=canonical_json({"array_key": "pd__delivery_memory"}),
        completion_json=canonical_json({"array_key": "pd__completion"}),
        candidate_sites_json=canonical_json({"array_keys": sorted(f"pd__{name}" for name in PD_SITE_ARRAYS | {"candidate_sites", "available_sites", "embryo_sites", "stable_sites", "inactive_sites"})}),
        clocks=clocks,
        embryo_stable_transition_json=canonical_json({"pd_scalars": pd_scalars, "array_keys": sorted(f"pd__{name}" for name in PD_POINT_ARRAYS)}),
        pd_damage_front_json=canonical_json({"array_keys": sorted(f"pd__{name}" for name in PD_BOND_ARRAYS | PD_VECTOR_ARRAYS | PD_PATH_ARRAYS), "pd_scalars": pd_scalars}),
        rng_states_json=canonical_json(rng_states),
        rng_digests_json=canonical_json(rng_digests),
        geometry_identity_state_json=canonical_json({"mesh_nodes": "historical_exact", "root_xy": "historical_exact", "array_schema": array_schema}),
        source_hashes_json=canonical_json(SOURCE_SHA256),
        request_hash=request.request_hash,
        solver_version=ADAPTER_VERSION,
        aggregate_cumulative_hazard=float(sum(clock.cumulative_hazard for clock in clocks)),
        controller_next_block_cycles=0.0,
        cap_floor_diagnostics=("new_v9_controller_state_initialized_at_conversion",),
    )
    for value in arrays.values():
        value.setflags(write=False)
    return AdaptedV87Checkpoint(
        capsule=capsule, arrays=arrays, metadata=metadata, checkpoint_sha256=digest,
        historical_fields=tuple(sorted(keys | {f"metadata.{name}" for name in REQUIRED_METADATA})),
        deterministic_reconstructed_structure=("mesh connectivity from immutable request", "feature surface indices", "PD neighborhoods/bonds", "geometry caches via rebuild_mesh_geometry and patch.update_geometry"),
        new_v9_controller_state=("controller_next_block_cycles", "commit_count", "committed_cycle_measure", "v9 generation identity"),
    )


def write_v87_roundtrip_checkpoint(adapted: AdaptedV87Checkpoint, path: Path) -> None:
    """Write exact historical fields back to a v8.7-compatible NPZ atomically.

    This is for adapter-equivalence proofs only. New v9 controller fields are
    deliberately excluded and no historical physical/stochastic value changes.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {name: np.asarray(value) for name, value in adapted.arrays.items()}
    payload["metadata_json"] = np.asarray(json.dumps(dict(adapted.metadata), allow_nan=True))
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("wb") as stream:
        np.savez_compressed(stream, **payload)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)
