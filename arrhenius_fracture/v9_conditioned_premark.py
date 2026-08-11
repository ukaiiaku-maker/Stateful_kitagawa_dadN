"""Atomic conversion of a conditioned first-passage boundary into a PD branch."""
from __future__ import annotations

import copy
import hashlib
import json
import math

import numpy as np

from .v9_pd_shared_root_marked_cleavage import normalized_available_site_marks


SCHEMA = "V9_CONDITIONED_PREMARK_ATOMIC_BRANCH_1"


def _stream_seed(namespace, branch_id, label):
    payload = f"{namespace}\0{branch_id}\0{label}".encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little")


def commit_conditioned_attempt(clock, patch, pd_state, capsule, branch_id):
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
    loc = capsule["localization"]
    target = float(capsule["conditioned_action"])
    trial_clock = clock.copy()
    trial_state = copy.deepcopy(pd_state)
    if trial_clock.attempt_count != 0 or not math.isclose(
            trial_clock.global_cumulative_action, target, rel_tol=0.0, abs_tol=2e-15):
        raise RuntimeError("restored clock is not the pristine conditioned boundary")
    # The conditioned crossing substitutes for the original hazard draw.  Only
    # the *next* threshold is sampled, from the branch renewal stream.
    trial_clock.global_threshold_action = target
    ns = str(capsule["branch_seed_namespace"])
    seeds = {name: _stream_seed(ns, branch_id, name)
             for name in ("spatial_mark", "embryo_transition", "renewal_hazard")}
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
    patch.create_marked_embryo(trial_state, event["site_id"], float(loc["crossing_cycle"]))
    node = int(event["pd_node_id"])
    event["local_fields"].update({
        "effective_opening_stress_Pa": float(loc["local_opening_stress_Pa"][node]),
        "local_backstress_Pa": float(loc["local_backstress_Pa"][node]),
        "local_state_shift_eV": float(loc["local_state_shift_eV"][node]),
        "local_cleavage_log_propensity_s": float(loc["local_cleavage_log_propensity_s"][node]),
    })
    return {
        "schema": SCHEMA, "branch_id": str(branch_id), "stream_seeds": seeds,
        "clock": trial_clock, "pd_state": trial_state, "event": event,
        "transition_rng_state": copy.deepcopy(patch_trial_rng.bit_generator.state),
        "mark_entropy_nats": float(loc["mark_entropy_nats"]),
        "source_analysis_checkpoint": capsule["source_analysis_checkpoint"],
        "source_replay_checkpoint": capsule["source_replay_checkpoint"],
    }


def json_branch_summary(result):
    """Serializable provenance summary; physical arrays remain in checkpoints."""
    return json.loads(json.dumps({
        "schema": result["schema"], "branch_id": result["branch_id"],
        "stream_seeds": result["stream_seeds"], "event": result["event"],
        "transition_rng_state": result["transition_rng_state"],
        "mark_entropy_nats": result["mark_entropy_nats"],
        "source_analysis_checkpoint": result["source_analysis_checkpoint"],
        "source_replay_checkpoint": result["source_replay_checkpoint"],
        "renewed_global_threshold_action": result["clock"].global_threshold_action,
        "attempt_count": result["clock"].attempt_count,
    }))
