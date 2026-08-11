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
        clone = copy.copy(self)
        clone.mpz = self.mpz.copy()
        clone.audit = copy.deepcopy(self.audit)
        clone.last_attempt = copy.deepcopy(self.last_attempt)
        clone._hazard_rng = np.random.default_rng()
        clone._hazard_rng.bit_generator.state = copy.deepcopy(
            self._hazard_rng.bit_generator.state
        )
        clone._mark_rng = np.random.default_rng()
        clone._mark_rng.bit_generator.state = copy.deepcopy(
            self._mark_rng.bit_generator.state
        )
        return clone

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

    def _advance_phase_block_exact(self, cycles, frequency_Hz, T_K, root_tensors):
        """Advance the authoritative MPZ/action using only root FEM tensors."""
        tensors = np.asarray(root_tensors, float)
        if tensors.ndim != 3 or tensors.shape[1:] != (2, 2) or len(tensors) < 1:
            raise ValueError("root phase tensors must have shape (n_phase,2,2)")
        cycles = max(float(cycles), 0.0)
        frequency = float(frequency_Hz)
        if frequency <= 0.0:
            raise ValueError("frequency must be positive")
        drives = [self.mpz.resolve_root_tensor(tensor) for tensor in tensors]
        opening = float(np.mean([d["opening_stress_Pa"] for d in drives]))
        signed = np.mean(np.stack([d["tau_signed_Pa"] for d in drives]), axis=0)
        dt = cycles / frequency
        self.mpz.advance(0.5 * dt, T_K, opening, signed)
        logs = np.asarray([
            self.cleavage_log_rate_s(d["opening_stress_Pa"], T_K) for d in drives
        ], float)
        log_average = float(logsumexp(logs) - math.log(len(logs)))
        log_increment = log_average + math.log(dt) if dt > 0.0 else -math.inf
        self.add_log_action(log_increment)
        self.mpz.advance(0.5 * dt, T_K, opening, signed)
        self.time_s += dt
        phase_action_fraction = np.exp(logs - float(logsumexp(logs)))
        return {
            "log_action_increment": log_increment,
            "action_increment": math.exp(log_increment) if log_increment > -745 else 0.0,
            "phase_log_rate_s": logs,
            "phase_action_fraction": phase_action_fraction,
            "phase_average_opening_stress_Pa": opening,
            "phase_average_signed_shear_Pa": signed,
        }

    def propose_phase_block(self, cycles, frequency_Hz, T_K, root_tensors,
                            *, localization_relative_tolerance=1e-12):
        """Side-effect-free ordered localization of at most the first crossing.

        Complete cycles are advanced first.  The crossing cycle is then
        evaluated in phase order and the final crossing is linearized only
        inside its constant-rate phase interval.
        """
        requested = max(float(cycles), 0.0)
        start = self.copy()
        trial = start.copy()
        detail = trial._advance_phase_block_exact(
            requested, frequency_Hz, T_K, root_tensors
        )
        if not trial.threshold_crossed():
            return {"state": trial, "cycles_consumed": requested,
                    "cycles_unused": 0.0, "crossed": False, "detail": detail}
        # Find the largest complete-cycle boundary strictly before crossing.
        max_complete = int(math.floor(requested))
        lo, hi = 0, max_complete
        while lo < hi:
            mid = (lo + hi + 1) // 2
            probe = start.copy()
            if mid:
                probe._advance_phase_block_exact(mid, frequency_Hz, T_K, root_tensors)
            if probe.threshold_crossed(): hi = mid - 1
            else: lo = mid
        complete = lo
        final = start.copy()
        if complete:
            final._advance_phase_block_exact(
                float(complete), frequency_Hz, T_K, root_tensors
            )

        tensors = np.asarray(root_tensors, float)
        drives = [final.mpz.resolve_root_tensor(tensor) for tensor in tensors]
        opening = float(np.mean([d["opening_stress_Pa"] for d in drives]))
        signed = np.mean(np.stack([d["tau_signed_Pa"] for d in drives]), axis=0)
        period = 1.0 / float(frequency_Hz)
        # Qualified canonical ordering: plastic half-step -> phase-ordered
        # midpoint cleavage -> plastic half-step.  At a crossing the latter is
        # not executed because the endpoint is inside the cleavage stage.
        final.mpz.advance(0.5 * period, T_K, opening, signed)
        logs = np.asarray([
            final.cleavage_log_rate_s(d["opening_stress_Pa"], T_K)
            for d in drives
        ], float)
        phase_dt = period / len(logs)
        phase_actions = np.exp(np.clip(logs + math.log(phase_dt), -745.0, 709.0))
        residual = max(final.global_threshold_action - final.global_cumulative_action, 0.0)
        cumulative = np.cumsum(phase_actions)
        phase = int(np.searchsorted(cumulative, residual, side="left"))
        if phase >= len(logs):
            raise RuntimeError("ordered phase localization failed to bracket crossing")
        before_phase = float(cumulative[phase - 1]) if phase else 0.0
        within = float(np.clip(
            (residual - before_phase) / max(phase_actions[phase], 1e-300),
            0.0, 1.0,
        ))
        consumed = float(complete) + (float(phase) + within) / len(logs)
        final.global_cumulative_action = final.global_threshold_action
        final.log_global_cumulative_action = math.log(final.global_threshold_action)
        final.time_s += (float(phase) + within) * phase_dt
        detail = {
            "log_action_increment": math.log(max(
                final.global_threshold_action - start.global_cumulative_action, 1e-300
            )),
            "action_increment": final.global_threshold_action - start.global_cumulative_action,
            "phase_log_rate_s": logs,
            "phase_action_per_cycle": phase_actions,
            "phase_average_opening_stress_Pa": opening,
            "phase_average_signed_shear_Pa": signed,
            "complete_cycles_before_crossing": complete,
            "crossing_phase_fraction": within,
            "localization_mode": "ordered_complete_cycles_final_cycle_phase_within_phase",
        }
        return {"state": final, "cycles_consumed": consumed,
                "cycles_unused": requested - consumed, "crossed": True,
                "phase_index": phase, "phase_fraction": within,
                "detail": detail}

    def propose_mpz_only_block(self, cycles, frequency_Hz, T_K, root_tensors):
        """Advance deterministic signed MPZ while the attempt clock is paused/stopped."""
        trial = self.copy()
        action = trial.global_cumulative_action
        log_action = trial.log_global_cumulative_action
        detail = trial._advance_phase_block_exact(
            cycles, frequency_Hz, T_K, root_tensors
        )
        trial.global_cumulative_action = action
        trial.log_global_cumulative_action = log_action
        detail["action_increment"] = 0.0
        detail["log_action_increment"] = -math.inf
        detail["clock_advanced"] = False
        return {"state": trial, "cycles_consumed": float(cycles),
                "cycles_unused": 0.0, "crossed": False, "detail": detail}

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
