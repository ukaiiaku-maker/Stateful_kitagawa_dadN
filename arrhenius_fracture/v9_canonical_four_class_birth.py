"""Canonical four-class stable-crack-birth endpoint for v9.

This path is deliberately separate from the legacy completion-gated PD site
model.  It stops at the first audited sharp-front cleavage first passage and
never executes the associated crack-advance reward.
"""
from __future__ import annotations

import copy
import math

import numpy as np
from scipy.special import gammainc, gammaln

from .v9_four_class_signed_mpz import SignedMPZPreBirthState


MODEL_ID = "v9_canonical_four_class_cleavage_stable_birth_v1"
ENDPOINT = "canonical_cleavage_stable_crack_birth"
HAZARD_SCHEMA = "v10.1.7.2_stochastic_cleavage_threshold"


def canonical_effective_cleavage_rate(raw_rate_s: float, m_hits: float = 3.0,
                                      tau_c_s: float = 1.0e-6) -> float:
    """Exact ``UnifiedMPZFrontEngine.lambda_cleave`` renewal transform."""
    raw = max(float(raw_rate_s), 0.0)
    m = max(float(m_hits), 1.0)
    if m <= 1.0 + 1.0e-12:
        return raw
    tau = max(float(tau_c_s), 1.0e-30)
    return float(gammainc(m, min(raw * tau, 1.0e12)) / tau)


def canonical_effective_cleavage_log_rate(log_raw_rate_s: float,
                                          m_hits: float = 3.0,
                                          tau_c_s: float = 1.0e-6) -> float:
    """Log of the identical renewal law without interpreting underflow as zero."""
    m = max(float(m_hits), 1.0)
    if m <= 1.0 + 1.0e-12:
        return float(log_raw_rate_s)
    tau = max(float(tau_c_s), 1.0e-30)
    log_x = float(log_raw_rate_s) + math.log(tau)
    if log_x < math.log(1.0e-4):
        # P(m,x)=x^m/Gamma(m+1)*(1-m*x/(m+1)+O(x^2)).
        correction = 0.0 if log_x < -745.0 else math.log1p(
            -m * math.exp(log_x) / (m + 1.0)
        )
        return m * log_x - float(gammaln(m + 1.0)) + correction - math.log(tau)
    x = min(math.exp(min(log_x, math.log(1.0e12))), 1.0e12)
    probability = float(gammainc(m, x))
    return math.log(probability) - math.log(tau)


class CanonicalFourClassBirthState:
    """Restartable signed-MPZ state plus the one canonical cleavage clock."""

    def __init__(self, option_id, source_root, *, shear_modulus_Pa, poisson,
                 burgers_m, initial_tip_radius_m, hazard_seed, engine_id=0,
                 minimum_threshold=1.0e-12, m_hits=3.0, tau_c_s=1.0e-6):
        self.mpz = SignedMPZPreBirthState(
            option_id, source_root, shear_modulus_Pa=shear_modulus_Pa,
            poisson=poisson, burgers_m=burgers_m,
            initial_tip_radius_m=initial_tip_radius_m,
        )
        self.option_id = str(option_id)
        self.initial_tip_radius_m = float(initial_tip_radius_m)
        self.m_hits = float(m_hits)
        self.tau_c_s = float(tau_c_s)
        self.hazard_seed = int(hazard_seed)
        self.engine_id = int(engine_id)
        self.minimum_threshold = max(float(minimum_threshold), 1.0e-300)
        self._hazard_rng = np.random.default_rng(
            np.random.SeedSequence([self.hazard_seed, self.engine_id])
        )
        self.hazard_threshold_action = self._draw_threshold()
        self.cumulative_cleavage_hazard = 0.0
        self.log_cumulative_cleavage_hazard = -math.inf
        self.time_s = 0.0
        self.stable_crack_birth_time_s = None
        self.audit = copy.deepcopy(self.mpz.audit) | {
            "model_id": MODEL_ID,
            "fatigue_endpoint": ENDPOINT,
            "stable_birth_definition": "first_canonical_cleavage_first_passage",
            "cleavage_hazard_schema": HAZARD_SCHEMA,
            "cleavage_hazard_mode": "exponential",
            "cleavage_hazard_seed": self.hazard_seed,
            "cleavage_engine_id": self.engine_id,
            "cleavage_multihit_m": self.m_hits,
            "cleavage_multihit_tau_s": self.tau_c_s,
            "pd_candidate_site_multiplier": False,
            "k2_completion_gate": False,
            "post_cleavage_stabilization_healing": False,
            "post_birth_growth": False,
        }

    @property
    def fired(self):
        return self.stable_crack_birth_time_s is not None

    def _draw_threshold(self):
        return max(float(self._hazard_rng.exponential(1.0)), self.minimum_threshold)

    def copy(self):
        clone = copy.copy(self)
        clone.mpz = copy.copy(self.mpz)
        clone.mpz.state = self.mpz.state.copy()
        clone.mpz.audit = copy.deepcopy(self.mpz.audit)
        clone.audit = copy.deepcopy(self.audit)
        clone._hazard_rng = np.random.default_rng()
        clone._hazard_rng.bit_generator.state = copy.deepcopy(
            self._hazard_rng.bit_generator.state
        )
        return clone

    def capsule(self):
        return {
            "schema": "V9_CANONICAL_FOUR_CLASS_BIRTH_CAPSULE_1",
            "audit": copy.deepcopy(self.audit),
            "mpz": self.mpz.capsule(),
            "hazard_rng_state": copy.deepcopy(self._hazard_rng.bit_generator.state),
            "hazard_threshold_action": self.hazard_threshold_action,
            "cumulative_cleavage_hazard": self.cumulative_cleavage_hazard,
            "log_cumulative_cleavage_hazard": self.log_cumulative_cleavage_hazard,
            "time_s": self.time_s,
            "stable_crack_birth_time_s": self.stable_crack_birth_time_s,
        }

    def restore_capsule(self, capsule):
        if capsule.get("schema") != "V9_CANONICAL_FOUR_CLASS_BIRTH_CAPSULE_1":
            raise RuntimeError("unsupported canonical-birth capsule schema")
        if capsule.get("audit") != self.audit:
            raise RuntimeError("canonical-birth capsule provenance/request mismatch")
        self.mpz.restore_capsule(capsule["mpz"])
        self._hazard_rng.bit_generator.state = copy.deepcopy(capsule["hazard_rng_state"])
        self.hazard_threshold_action = float(capsule["hazard_threshold_action"])
        self.cumulative_cleavage_hazard = float(capsule["cumulative_cleavage_hazard"])
        self.log_cumulative_cleavage_hazard = float(
            capsule["log_cumulative_cleavage_hazard"]
        )
        self.time_s = float(capsule["time_s"])
        value = capsule["stable_crack_birth_time_s"]
        self.stable_crack_birth_time_s = None if value is None else float(value)

    def effective_opening_stress_Pa(self, nominal_root_stress_Pa):
        """Apply the canonical K shielding/blunted-radius tip transformation."""
        nominal = max(float(nominal_root_stress_Pa), 0.0)
        K_nominal = nominal * math.sqrt(2.0 * math.pi * self.initial_tip_radius_m)
        summary = self.mpz.summary()
        K_eff = max(K_nominal - summary["signed_active_K_shield_Pa_sqrt_m"], 0.0)
        return K_eff / math.sqrt(2.0 * math.pi * max(summary["tip_radius_m"], 1.0e-30))

    def cleavage_rates(self, nominal_root_stress_Pa, T_K):
        sigma = self.effective_opening_stress_Pa(nominal_root_stress_Pa)
        log_raw = self.mpz.cleavage_log_rate_s(sigma, T_K)
        raw = math.exp(log_raw) if log_raw > math.log(np.finfo(float).tiny) else 0.0
        effective = canonical_effective_cleavage_rate(raw, self.m_hits, self.tau_c_s)
        log_effective = canonical_effective_cleavage_log_rate(
            log_raw, self.m_hits, self.tau_c_s
        )
        return {"sigma_cleave_eff_Pa": sigma, "cleavage_log_rate_raw_s": log_raw,
                "cleavage_rate_raw_s": raw, "cleavage_rate_effective_s": effective,
                "cleavage_log_rate_effective_s": log_effective}

    def _advance_exact_interval(self, dt_s, T_K, nominal_root_stress_Pa,
                                signed_shear_Pa):
        """Canonical Strang ordering without a post-first-passage reward."""
        dt = max(float(dt_s), 0.0)
        if dt == 0.0:
            rates = self.cleavage_rates(nominal_root_stress_Pa, T_K)
            return rates | {"hazard_increment": 0.0}
        half = 0.5 * dt
        self.mpz.advance(half, T_K, nominal_root_stress_Pa, signed_shear_Pa)
        rates = self.cleavage_rates(nominal_root_stress_Pa, T_K)
        log_dH = rates["cleavage_log_rate_effective_s"] + math.log(dt)
        dH = math.exp(log_dH) if log_dH > math.log(np.finfo(float).tiny) else 0.0
        self.log_cumulative_cleavage_hazard = float(np.logaddexp(
            self.log_cumulative_cleavage_hazard, log_dH
        ))
        self.cumulative_cleavage_hazard = (
            math.exp(self.log_cumulative_cleavage_hazard)
            if self.log_cumulative_cleavage_hazard > math.log(np.finfo(float).tiny)
            else 0.0
        )
        self.mpz.advance(half, T_K, nominal_root_stress_Pa, signed_shear_Pa)
        self.time_s += dt
        return rates | {"hazard_increment": dH, "log_hazard_increment": log_dH}

    def advance(self, dt_s, T_K, nominal_root_stress_Pa, signed_shear_Pa,
                *, localization_relative_tolerance=1.0e-12):
        """Advance transactionally, stopping exactly at the first crossing."""
        if self.fired:
            return self.diagnostics() | {"dt_consumed_s": 0.0, "dt_unused_s": float(dt_s)}
        requested = max(float(dt_s), 0.0)
        start = self.copy()
        result = self._advance_exact_interval(
            requested, T_K, nominal_root_stress_Pa, signed_shear_Pa
        )
        if self.log_cumulative_cleavage_hazard < math.log(self.hazard_threshold_action):
            return self.diagnostics() | result | {"dt_consumed_s": requested, "dt_unused_s": 0.0}

        lo, hi = 0.0, requested
        tolerance = max(requested * float(localization_relative_tolerance), 1.0e-15)
        while hi - lo > tolerance:
            mid = 0.5 * (lo + hi)
            trial = start.copy()
            trial._advance_exact_interval(mid, T_K, nominal_root_stress_Pa, signed_shear_Pa)
            if trial.log_cumulative_cleavage_hazard >= math.log(trial.hazard_threshold_action):
                hi = mid
            else:
                lo = mid
        final = start.copy()
        result = final._advance_exact_interval(hi, T_K, nominal_root_stress_Pa, signed_shear_Pa)
        final.cumulative_cleavage_hazard = final.hazard_threshold_action
        final.log_cumulative_cleavage_hazard = math.log(final.hazard_threshold_action)
        final.stable_crack_birth_time_s = final.time_s
        self.__dict__.clear()
        self.__dict__.update(final.__dict__)
        return self.diagnostics() | result | {
            "dt_consumed_s": hi, "dt_unused_s": requested - hi,
            "first_passage_localized": True,
        }

    def _advance_phase_block_exact(self, cycles, frequency_Hz, T_K, root_tensors):
        tensors = np.asarray(root_tensors, float)
        if tensors.ndim != 3 or tensors.shape[1:] != (2, 2) or len(tensors) < 1:
            raise ValueError("root phase history must have shape (n_phase,2,2)")
        cycles = max(float(cycles), 0.0)
        frequency = float(frequency_Hz)
        if frequency <= 0.0:
            raise ValueError("frequency must be positive")
        dt = cycles / frequency
        drives = [self.mpz.resolve_root_tensor(tensor) for tensor in tensors]
        opening = float(np.mean([drive["opening_stress_Pa"] for drive in drives]))
        signed = np.mean(np.stack([drive["tau_signed_Pa"] for drive in drives]), axis=0)
        self.mpz.advance(0.5 * dt, T_K, opening, signed)
        phase_rates = [
            self.cleavage_rates(drive["opening_stress_Pa"], T_K)
            for drive in drives
        ]
        logs = np.asarray([rate["cleavage_log_rate_effective_s"] for rate in phase_rates])
        pivot = float(np.max(logs))
        log_average = pivot + math.log(float(np.mean(np.exp(logs - pivot))))
        log_dH = log_average + math.log(dt) if dt > 0.0 else -math.inf
        self.log_cumulative_cleavage_hazard = float(np.logaddexp(
            self.log_cumulative_cleavage_hazard, log_dH
        ))
        self.cumulative_cleavage_hazard = (
            math.exp(self.log_cumulative_cleavage_hazard)
            if self.log_cumulative_cleavage_hazard > math.log(np.finfo(float).tiny)
            else 0.0
        )
        self.mpz.advance(0.5 * dt, T_K, opening, signed)
        self.time_s += dt
        return {
            "phase_count": len(tensors),
            "phase_average_opening_stress_Pa": opening,
            "phase_average_signed_shear_Pa": signed,
            "cleavage_log_rate_cycle_average_s": log_average,
            "log_hazard_increment": log_dH,
            "hazard_increment": math.exp(log_dH) if log_dH > -745.0 else 0.0,
        }

    def advance_fem_phase_block(self, cycles, frequency_Hz, T_K, root_tensors,
                                *, localization_relative_tolerance=1.0e-12):
        """Advance one accepted full-FEM cyclic block transactionally."""
        if self.fired:
            return self.diagnostics() | {"cycles_consumed": 0.0, "cycles_unused": float(cycles)}
        requested = max(float(cycles), 0.0)
        start = self.copy()
        result = self._advance_phase_block_exact(
            requested, frequency_Hz, T_K, root_tensors
        )
        log_threshold = math.log(self.hazard_threshold_action)
        if self.log_cumulative_cleavage_hazard < log_threshold:
            return self.diagnostics() | result | {
                "cycles_consumed": requested, "cycles_unused": 0.0,
            }
        lo, hi = 0.0, requested
        tolerance = max(requested * float(localization_relative_tolerance), 1.0e-12)
        while hi - lo > tolerance:
            mid = 0.5 * (lo + hi)
            trial = start.copy()
            trial._advance_phase_block_exact(mid, frequency_Hz, T_K, root_tensors)
            if trial.log_cumulative_cleavage_hazard >= log_threshold:
                hi = mid
            else:
                lo = mid
        final = start.copy()
        result = final._advance_phase_block_exact(hi, frequency_Hz, T_K, root_tensors)
        final.cumulative_cleavage_hazard = final.hazard_threshold_action
        final.log_cumulative_cleavage_hazard = log_threshold
        final.stable_crack_birth_time_s = final.time_s
        self.__dict__.clear()
        self.__dict__.update(final.__dict__)
        return self.diagnostics() | result | {
            "cycles_consumed": hi, "cycles_unused": requested - hi,
            "first_passage_localized": True,
        }

    def diagnostics(self):
        return {
            "model_id": MODEL_ID,
            "fatigue_endpoint": ENDPOINT,
            "cumulative_cleavage_hazard": self.cumulative_cleavage_hazard,
            "log_cumulative_cleavage_hazard": self.log_cumulative_cleavage_hazard,
            "hazard_threshold_action": self.hazard_threshold_action,
            "stable_crack_birth": self.fired,
            "stable_crack_birth_time_s": self.stable_crack_birth_time_s,
            "stable_crack_birth_survival": math.exp(-self.cumulative_cleavage_hazard),
            "survival_semantics": "stable_crack_birth_survival_exp_minus_H_cleave",
            "time_s": self.time_s,
        }
