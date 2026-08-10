"""Renewing canonical cleavage attempts with terminal energy-admissible birth."""
from __future__ import annotations

import copy
import math

from .v9_canonical_four_class_birth import CanonicalFourClassBirthState


MODEL_ID = "v9_canonical_four_class_energy_gated_stable_birth_v2"
ENDPOINT = "energy_admissible_stable_crack_birth"


class EnergyGatedStableBirthState(CanonicalFourClassBirthState):
    """The v10.2.30 threshold stream, including consumed rejected attempts."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.hazard_action_origin = 0.0
        self.pending_attempt = None
        self.cleavage_attempt_count = 0
        self.nonpropagating_attempt_count = 0
        self.attempt_history = []
        self.audit.update({
            "model_id": MODEL_ID, "fatigue_endpoint": ENDPOINT,
            "stable_birth_definition": "first_energy_admissible_canonical_cleavage_attempt",
            "post_first_passage_energy_gate": "v10.2.30_mesh_consistent_fixed_opening",
            "threshold_scaled_event_proposal": True,
            "rejected_attempt_consumes_threshold": True,
            "cleavage_hazard_gated_before_first_passage": False,
            "post_birth_growth": False,
        })

    @property
    def fired(self):
        return self.stable_crack_birth_time_s is not None

    def capsule(self):
        base = super().capsule()
        base["schema"] = "V9_ENERGY_GATED_STABLE_BIRTH_CAPSULE_2"
        base.update({
            "hazard_action_origin": self.hazard_action_origin,
            "pending_attempt": copy.deepcopy(self.pending_attempt),
            "cleavage_attempt_count": self.cleavage_attempt_count,
            "nonpropagating_attempt_count": self.nonpropagating_attempt_count,
            "attempt_history": copy.deepcopy(self.attempt_history),
        })
        return base

    def restore_capsule(self, capsule):
        if capsule.get("schema") != "V9_ENERGY_GATED_STABLE_BIRTH_CAPSULE_2":
            raise RuntimeError("unsupported energy-gated stable-birth capsule schema")
        base = copy.deepcopy(capsule); base["schema"] = "V9_CANONICAL_FOUR_CLASS_BIRTH_CAPSULE_1"
        super().restore_capsule(base)
        self.hazard_action_origin = float(capsule["hazard_action_origin"])
        self.pending_attempt = copy.deepcopy(capsule["pending_attempt"])
        self.cleavage_attempt_count = int(capsule["cleavage_attempt_count"])
        self.nonpropagating_attempt_count = int(capsule["nonpropagating_attempt_count"])
        self.attempt_history = copy.deepcopy(capsule["attempt_history"])

    def advance_fem_phase_block(self, cycles, frequency_Hz, T_K, root_tensors,
                                *, localization_relative_tolerance=1e-12):
        if self.fired or self.pending_attempt is not None:
            return self.diagnostics() | {"cycles_consumed": 0.0, "cycles_unused": float(cycles)}
        requested = max(float(cycles), 0.0); start = self.copy()
        result = self._advance_phase_block_exact(requested, frequency_Hz, T_K, root_tensors)
        target = self.hazard_action_origin + self.hazard_threshold_action
        if self.cumulative_cleavage_hazard < target:
            return self.diagnostics() | result | {"cycles_consumed": requested, "cycles_unused": 0.0}
        lo, hi = 0.0, requested
        tolerance = max(requested * float(localization_relative_tolerance), 1e-12)
        while hi - lo > tolerance:
            mid = 0.5 * (lo + hi); trial = start.copy()
            trial._advance_phase_block_exact(mid, frequency_Hz, T_K, root_tensors)
            if trial.cumulative_cleavage_hazard >= target: hi = mid
            else: lo = mid
        final = start.copy(); result = final._advance_phase_block_exact(hi, frequency_Hz, T_K, root_tensors)
        final.cumulative_cleavage_hazard = target
        final.log_cumulative_cleavage_hazard = math.log(target)
        final.pending_attempt = {
            "attempt_index": final.cleavage_attempt_count,
            "threshold_action": final.hazard_threshold_action,
            "cumulative_cleavage_hazard": target,
            "time_s": final.time_s,
        }
        self.__dict__.clear(); self.__dict__.update(final.__dict__)
        return self.diagnostics() | result | {"cycles_consumed": hi,
            "cycles_unused": requested-hi, "first_passage_localized": True,
            "cleavage_attempt_pending_energy_gate": True}

    def resolve_pending_attempt(self, gate_result):
        if self.pending_attempt is None:
            raise RuntimeError("no localized cleavage attempt awaits the energy gate")
        gate = {k: copy.deepcopy(v) for k, v in gate_result.items()
                if k not in {"equilibrated_displacement", "trial_rows"}}
        admitted = float(gate_result.get("committed_event_length_m", 0.0)) > 0.0
        record = copy.deepcopy(self.pending_attempt) | gate | {
            "classification": "stable_crack_birth" if admitted else "nonpropagating_cleavage_attempt"
        }
        self.attempt_history.append(record); self.cleavage_attempt_count += 1
        if admitted:
            self.stable_crack_birth_time_s = self.time_s
        else:
            self.nonpropagating_attempt_count += 1
            self.hazard_action_origin = self.cumulative_cleavage_hazard
            self.hazard_threshold_action = self._draw_threshold()
        self.pending_attempt = None
        return record

    def diagnostics(self):
        no_attempt_survival = math.exp(-self.cumulative_cleavage_hazard)
        return {
            "model_id": MODEL_ID, "fatigue_endpoint": ENDPOINT,
            "cumulative_cleavage_hazard": self.cumulative_cleavage_hazard,
            "log_cumulative_cleavage_hazard": self.log_cumulative_cleavage_hazard,
            "hazard_threshold_action": self.hazard_threshold_action,
            "hazard_action_origin": self.hazard_action_origin,
            "cleavage_attempt_count": self.cleavage_attempt_count,
            "nonpropagating_attempt_count": self.nonpropagating_attempt_count,
            "cleavage_attempt_pending_energy_gate": self.pending_attempt is not None,
            "stable_crack_birth": self.fired,
            "stable_crack_birth_time_s": self.stable_crack_birth_time_s,
            "no_cleavage_attempt_survival_exp_minus_H": no_attempt_survival,
            "stable_birth_survival": None,
            "survival_semantics": "exp_minus_H_is_no_cleavage_attempt_survival_only",
            "time_s": self.time_s,
        }

