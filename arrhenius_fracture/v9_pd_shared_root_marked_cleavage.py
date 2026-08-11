"""Shared-root signed-MPZ/global-cleavage marked PD process.

The total attempt intensity is owned solely by the authoritative root-local
signed MPZ.  Realized PD sites provide a normalized location mark and never
multiply that intensity.
"""
from __future__ import annotations

import copy
import math

import numpy as np
from scipy.special import logsumexp

from .v9_canonical_four_class_birth import canonical_effective_cleavage_log_rate
from .v9_four_class_signed_mpz import SignedMPZPreBirthState


MODEL_ID = "v9_four_class_stateful_PD_shared_root_MPZ_marked_cleavage_v1"
CAPSULE_SCHEMA = "V9_SHARED_ROOT_MPZ_MARKED_CLEAVAGE_CAPSULE_1"


def normalized_available_site_marks(site_node_index, site_available,
                                    node_log_propensity,
                                    node_initiation_weight):
    """Return site ids, normalized probabilities, and unnormalized log weights.

    A node's weight is divided among its available colocated sites.  Thus
    duplicating identical sites changes only mark resolution, never the global
    attempt intensity or the total probability assigned to that node.
    """
    nodes = np.asarray(site_node_index, dtype=np.int64)
    available = np.asarray(site_available, dtype=bool)
    log_local = np.asarray(node_log_propensity, dtype=float)
    initiation = np.asarray(node_initiation_weight, dtype=float)
    if nodes.shape != available.shape:
        raise ValueError("site node and availability ledgers must align")
    if log_local.shape != initiation.shape:
        raise ValueError("node propensity and initiation arrays must align")
    ids = np.flatnonzero(available)
    if not len(ids):
        raise RuntimeError("positive global hazard has no available PD mark support")
    selected_nodes = nodes[ids]
    if np.any(selected_nodes < 0) or np.any(selected_nodes >= len(log_local)):
        raise RuntimeError("site ledger references an invalid PD node")
    counts = np.bincount(selected_nodes, minlength=len(log_local))
    weights_ok = np.isfinite(log_local[selected_nodes]) & (
        np.isfinite(initiation[selected_nodes]) & (initiation[selected_nodes] > 0.0)
    )
    logw = np.full(len(ids), -math.inf)
    logw[weights_ok] = (
        log_local[selected_nodes[weights_ok]]
        + np.log(initiation[selected_nodes[weights_ok]])
        - np.log(counts[selected_nodes[weights_ok]])
    )
    norm = float(logsumexp(logw))
    if not math.isfinite(norm):
        raise RuntimeError("positive global hazard has no finite PD mark distribution")
    probabilities = np.exp(logw - norm)
    probabilities /= float(np.sum(probabilities))
    if not np.isclose(np.sum(probabilities), 1.0, rtol=0.0, atol=2e-15):
        raise RuntimeError("PD spatial mark failed probability conservation")
    return ids, probabilities, logw


class SharedRootMarkedCleavageState:
    """Atomic persistent global clock plus authoritative root-local MPZ."""

    def __init__(self, option_id, source_root, *, shear_modulus_Pa, poisson,
                 burgers_m, initial_tip_radius_m, hazard_seed, mark_seed,
                 engine_id=0, minimum_threshold=1e-12, m_hits=3.0,
                 tau_c_s=1e-6, mpz=None):
        self.mpz = mpz if mpz is not None else SignedMPZPreBirthState(
            option_id, source_root, shear_modulus_Pa=shear_modulus_Pa,
            poisson=poisson, burgers_m=burgers_m,
            initial_tip_radius_m=initial_tip_radius_m,
        )
        self.option_id = str(option_id)
        self.initial_tip_radius_m = float(initial_tip_radius_m)
        self.m_hits = float(m_hits)
        self.tau_c_s = float(tau_c_s)
        self.minimum_threshold = max(float(minimum_threshold), 1e-300)
        self.hazard_seed, self.mark_seed, self.engine_id = (
            int(hazard_seed), int(mark_seed), int(engine_id)
        )
        self._hazard_rng = np.random.default_rng(
            np.random.SeedSequence([self.hazard_seed, self.engine_id])
        )
        self._mark_rng = np.random.default_rng(
            np.random.SeedSequence([self.mark_seed, self.engine_id])
        )
        self.global_threshold_action = self._draw_threshold()
        self.global_cumulative_action = 0.0
        self.log_global_cumulative_action = -math.inf
        self.attempt_count = 0
        self.time_s = 0.0
        self.last_attempt = None
        self.audit = copy.deepcopy(getattr(self.mpz, "audit", {})) | {
            "model_id": MODEL_ID,
            "global_cleavage_clock_count": 1,
            "signed_mpz_state_count": 1,
            "legacy_K2_delivery_gate": False,
            "candidate_density_multiplies_global_hazard": False,
            "mark_rng_separate_from_cleavage_rng": True,
            "primary_PD_endpoint": "first_front_capture_after_stable_spatial_seed",
        }

    def _draw_threshold(self):
        return max(float(self._hazard_rng.exponential(1.0)), self.minimum_threshold)

    def copy(self):
        return copy.deepcopy(self)

    def capsule(self):
        return {
            "schema": CAPSULE_SCHEMA, "audit": copy.deepcopy(self.audit),
            "mpz": self.mpz.capsule(),
            "hazard_rng_state": copy.deepcopy(self._hazard_rng.bit_generator.state),
            "mark_rng_state": copy.deepcopy(self._mark_rng.bit_generator.state),
            "global_threshold_action": self.global_threshold_action,
            "global_cumulative_action": self.global_cumulative_action,
            "log_global_cumulative_action": self.log_global_cumulative_action,
            "attempt_count": self.attempt_count, "time_s": self.time_s,
            "last_attempt": copy.deepcopy(self.last_attempt),
        }

    def restore_capsule(self, capsule):
        if capsule.get("schema") != CAPSULE_SCHEMA or capsule.get("audit") != self.audit:
            raise RuntimeError("marked-cleavage capsule provenance/schema mismatch")
        self.mpz.restore_capsule(capsule["mpz"])
        self._hazard_rng.bit_generator.state = copy.deepcopy(capsule["hazard_rng_state"])
        self._mark_rng.bit_generator.state = copy.deepcopy(capsule["mark_rng_state"])
        for key in ("global_threshold_action", "global_cumulative_action",
                    "log_global_cumulative_action", "time_s"):
            setattr(self, key, float(capsule[key]))
        self.attempt_count = int(capsule["attempt_count"])
        self.last_attempt = copy.deepcopy(capsule["last_attempt"])

    def effective_opening_stress_Pa(self, nominal_root_stress_Pa):
        nominal = max(float(nominal_root_stress_Pa), 0.0)
        K = nominal * math.sqrt(2 * math.pi * self.initial_tip_radius_m)
        summary = self.mpz.summary()
        K = max(K - summary["signed_active_K_shield_Pa_sqrt_m"], 0.0)
        return K / math.sqrt(2 * math.pi * max(summary["tip_radius_m"], 1e-30))

    def cleavage_log_rate_s(self, nominal_root_stress_Pa, T_K):
        sigma = self.effective_opening_stress_Pa(nominal_root_stress_Pa)
        raw = self.mpz.cleavage_log_rate_s(sigma, T_K)
        return canonical_effective_cleavage_log_rate(raw, self.m_hits, self.tau_c_s)

    def add_log_action(self, log_increment):
        self.log_global_cumulative_action = float(np.logaddexp(
            self.log_global_cumulative_action, float(log_increment)
        ))
        self.global_cumulative_action = (
            math.exp(self.log_global_cumulative_action)
            if self.log_global_cumulative_action > -745 else 0.0
        )

    def threshold_crossed(self):
        return self.log_global_cumulative_action >= math.log(self.global_threshold_action)

    def select_mark(self, site_node_index, site_available, node_log_propensity,
                    node_initiation_weight, *, cycle, phase_index,
                    local_fields=None):
        if not self.threshold_crossed():
            raise RuntimeError("cannot draw a PD mark before global first passage")
        ids, p, logw = normalized_available_site_marks(
            site_node_index, site_available, node_log_propensity,
            node_initiation_weight,
        )
        chosen_offset = int(self._mark_rng.choice(len(ids), p=p))
        site = int(ids[chosen_offset])
        node = int(np.asarray(site_node_index, dtype=np.int64)[site])
        event = {
            "attempt_index": self.attempt_count + 1,
            "cycle": float(cycle), "phase_index": int(phase_index),
            "threshold_action": self.global_threshold_action,
            "cumulative_action": self.global_cumulative_action,
            "site_id": site, "pd_node_id": node,
            "mark_probability": float(p[chosen_offset]),
            "mark_log_weight": float(logw[chosen_offset]),
            "local_fields": copy.deepcopy(local_fields or {}),
        }
        self.attempt_count += 1
        self.last_attempt = event
        # Renewal increments are placed ahead of the cumulative action; no
        # action is discarded and a rejected/healed embryo consumes its Xi.
        self.global_threshold_action = (
            self.global_cumulative_action + self._draw_threshold()
        )
        return event

