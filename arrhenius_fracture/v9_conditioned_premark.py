"""Atomic conversion of a conditioned first-passage boundary into a PD branch."""
from __future__ import annotations

import copy
import hashlib
import json
import math
from dataclasses import dataclass

import numpy as np

from .v9_pd_shared_root_marked_cleavage import normalized_available_site_marks


SCHEMA = "V9_CONDITIONED_PREMARK_ATOMIC_BRANCH_1"


@dataclass(frozen=True)
class ConditionedBranchStreams:
    branch_id: str
    mark_stream_id: str
    transition_stream_id: str
    renewal_stream_id: str

    @classmethod
    def from_branch_id(cls, branch_id):
        value = str(branch_id)
        return cls(value, value, value, value)


def _stream_seed(namespace, branch_id, label):
    payload = f"{namespace}\0{branch_id}\0{label}".encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little")


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_verified_conditioned_capsule(directory):
    """Verify capsule, referenced files, and active generation identities."""
    from pathlib import Path
    directory = Path(directory)
    capsule_path = directory / "conditioned_premark_capsule.json"
    manifest_path = directory / "manifest.json"
    capsule = json.loads(capsule_path.read_text())
    manifest = json.loads(manifest_path.read_text())
    if _sha256(capsule_path) != manifest.get("capsule_sha256"):
        raise RuntimeError("conditioned capsule SHA-256 mismatch")
    files = {Path(name): expected for name, expected in manifest.get("source_files", {}).items()}
    for path, expected in files.items():
        if not path.is_file() or _sha256(path) != expected:
            raise RuntimeError(f"conditioned source hash mismatch: {path}")
    checkpoints = [path for path in files if path.name == "checkpoint_latest.npz"]
    if len(checkpoints) != 2:
        raise RuntimeError("conditioned manifest must reference source and replay checkpoints")
    replay = next((path for path in checkpoints if "replay" in str(path)), None)
    source = next((path for path in checkpoints if path != replay), None)
    if replay is None or source is None:
        raise RuntimeError("conditioned source/replay checkpoint roles are ambiguous")
    for role, path in (("source", source), ("replay", replay)):
        active = json.loads((path.parent / "v9_generations" / "ACTIVE.json").read_text())
        if active.get("generation") != capsule[f"{role}_generation"]:
            raise RuntimeError(f"conditioned {role} generation identity mismatch")
    return capsule, manifest, source, replay


def commit_conditioned_attempt(clock, patch, pd_state, capsule, branch):
    """Return a committed branch without mutating any supplied object.

    The conditioned action is not a physical threshold draw.  Branch-specific
    mark, transition, and renewal streams begin only at this transaction.
    """
    if capsule.get("protocol") != "conditioned_attempt_premark":
        raise RuntimeError("conditioned branch requires the dedicated premark protocol")
    required = ("conditioned_first_attempt_pending", "topology_continuation_permitted",
                "exact_phase_localization_present")
    if not all(capsule.get(k) is True for k in required):
        raise RuntimeError("conditioned capsule is incomplete or topology-forbidden")
    if capsule.get("physical_threshold_draw") is not False or capsule.get("mark_rng_consumed") is not False:
        raise RuntimeError("conditioned capsule has invalid threshold/mark provenance")
    streams = (ConditionedBranchStreams.from_branch_id(branch)
               if isinstance(branch, str) else branch)
    if not isinstance(streams, ConditionedBranchStreams):
        raise TypeError("branch must be an ID or ConditionedBranchStreams")
    loc = capsule["localization"]
    target = float(capsule["conditioned_action"])
    trial_clock = clock.copy()
    trial_state = copy.deepcopy(pd_state)
    if trial_clock.attempt_count != 0 or not math.isclose(
            trial_clock.global_cumulative_action, target, rel_tol=0.0, abs_tol=2e-15):
        raise RuntimeError("restored clock is not the pristine conditioned boundary")
    candidate_hash = hashlib.sha256(
        np.ascontiguousarray(trial_state.candidate_sites).view(np.uint8)
    ).hexdigest()
    geometry_hash = hashlib.sha256(
        np.ascontiguousarray(patch.xy).view(np.uint8)
    ).hexdigest()
    if candidate_hash != capsule["candidate_population_sha256"]:
        raise RuntimeError("conditioned candidate population hash mismatch")
    if geometry_hash != capsule["geometry_sha256"]:
        raise RuntimeError("conditioned geometry hash mismatch")
    # The conditioned crossing substitutes for the original hazard draw.  Only
    # the *next* threshold is sampled, from the branch renewal stream.
    trial_clock.global_threshold_action = target
    ns = str(capsule["branch_seed_namespace"])
    seeds = {
        "spatial_mark": _stream_seed(ns, streams.mark_stream_id, "spatial_mark"),
        "embryo_transition": _stream_seed(ns, streams.transition_stream_id, "embryo_transition"),
        "renewal_hazard": _stream_seed(ns, streams.renewal_stream_id, "renewal_hazard"),
    }
    trial_clock._mark_rng = np.random.default_rng(seeds["spatial_mark"])
    trial_clock._hazard_rng = np.random.default_rng(seeds["renewal_hazard"])
    patch_trial_rng = np.random.default_rng(seeds["embryo_transition"])

    available = np.asarray(trial_state.site_status, np.uint8) == 0
    ids, probabilities, _ = normalized_available_site_marks(
        trial_state.site_node_index, available,
        np.asarray(loc["local_cleavage_log_propensity_s"], float),
        patch.initiation_weight,
    )
    expected_ids = np.asarray(loc["available_site_ids"], np.int64)
    expected_p = np.asarray(loc["normalized_mark_probability"], float)
    if not np.array_equal(ids, expected_ids) or not np.allclose(
            probabilities, expected_p, rtol=2e-13, atol=2e-15):
        raise RuntimeError("exact-phase mark distribution does not replay")
    event = trial_clock.select_mark(
        trial_state.site_node_index, available,
        np.asarray(loc["local_cleavage_log_propensity_s"], float),
        patch.initiation_weight,
        cycle=float(loc["crossing_cycle"]), phase_index=int(loc["phase_index"]),
        local_fields={"phase_fraction": float(loc["phase_fraction"]),
                      "conditioned_action": target},
    )
    selected_site = int(event["site_id"])
    patch.create_external_marked_embryo(
        trial_state, selected_site, float(loc["crossing_cycle"]), rng=patch_trial_rng
    )
    transition_threshold = float(trial_state.site_transition_threshold[selected_site])
    transition_outcome = float(trial_state.site_transition_outcome_uniform[selected_site])
    node = int(event["pd_node_id"])
    event["local_fields"].update({
        "effective_opening_stress_Pa": float(loc["local_opening_stress_Pa"][node]),
        "local_backstress_Pa": float(loc["local_backstress_Pa"][node]),
        "local_state_shift_eV": float(loc["local_state_shift_eV"][node]),
        "local_cleavage_log_propensity_s": float(loc["local_cleavage_log_propensity_s"][node]),
    })
    return {
        "schema": SCHEMA, "branch_id": streams.branch_id,
        "stream_ids": {"mark_stream_id": streams.mark_stream_id,
                       "transition_stream_id": streams.transition_stream_id,
                       "renewal_stream_id": streams.renewal_stream_id},
        "stream_seeds": seeds,
        "clock": trial_clock, "pd_state": trial_state, "event": event,
        "transition_rng_state": copy.deepcopy(patch_trial_rng.bit_generator.state),
        "selected_site_transition_threshold": transition_threshold,
        "selected_site_transition_cumulative_action": 0.0,
        "selected_site_transition_outcome_uniform": transition_outcome,
        "mark_entropy_nats": float(loc["mark_entropy_nats"]),
        "source_analysis_checkpoint": capsule["source_analysis_checkpoint"],
        "source_replay_checkpoint": capsule["source_replay_checkpoint"],
    }


def json_branch_summary(result):
    """Serializable provenance summary; physical arrays remain in checkpoints."""
    return json.loads(json.dumps({
        "schema": result["schema"], "branch_id": result["branch_id"],
        "stream_ids": result["stream_ids"], "stream_seeds": result["stream_seeds"],
        "event": result["event"],
        "transition_rng_state": result["transition_rng_state"],
        "selected_site_transition_threshold": result["selected_site_transition_threshold"],
        "selected_site_transition_cumulative_action": result["selected_site_transition_cumulative_action"],
        "selected_site_transition_outcome_uniform": result["selected_site_transition_outcome_uniform"],
        "mark_entropy_nats": result["mark_entropy_nats"],
        "source_analysis_checkpoint": result["source_analysis_checkpoint"],
        "source_replay_checkpoint": result["source_replay_checkpoint"],
        "renewed_global_threshold_action": result["clock"].global_threshold_action,
        "attempt_count": result["clock"].attempt_count,
    }))
