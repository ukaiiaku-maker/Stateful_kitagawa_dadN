"""Spatial Stateful-PD adapter for :mod:`v9_pd_high_cycle`.

The driver supplies ``cycle_evaluator`` because it owns the accepted FEM phase
ordering.  The evaluator receives a private adapter clone and must return phase
data plus the endpoint active arrays.  This class owns all state separation and
commit/protection rules, so an evaluator cannot accidentally advance clocks.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import fields
import math
from typing import Any, Callable

import numpy as np

from .v9_pd_high_cycle import ActiveState, CycleEvaluation, ProtectedSignatures, _digest


ACTIVE_NAMES = ("ep_gp", "rho_gp", "epsp_acc_gp", "u", "log_delivery_memory")
LEDGER_NAMES = (
    "born_cumulative", "healed_cumulative", "born_sites_cumulative",
    "healed_sites_cumulative", "log_birth_cumulative_hazard",
)
STOCHASTIC_NAMES = (
    "site_birth_threshold", "birth_cumulative_hazard",
    "site_transition_threshold", "site_transition_cumulative_hazard",
    "site_transition_outcome_uniform",
)
TOPOLOGY_NAMES = (
    "site_status", "bond_damage", "primary_seed_node", "active_front",
    "active_front_bonds", "front_backbone_bonds", "front_wake_bonds",
    "front_process_bonds", "active_front_path_xy",
)


class SpatialPDDormantAdapter:
    def __init__(self, *, patch, pd_state, mesh, ep_gp, rho_gp, epsp_acc_gp, u,
                 cycles: float, plastic_work: float,
                 cycle_evaluator: Callable[["SpatialPDDormantAdapter"], dict[str, Any]]):
        self.patch = patch
        self.pd_state = pd_state
        self.mesh = mesh
        self.ep_gp = np.asarray(ep_gp, float).copy()
        self.rho_gp = np.asarray(rho_gp, float).copy()
        self.epsp_acc_gp = np.asarray(epsp_acc_gp, float).copy()
        self.u = np.asarray(u, float).copy()
        self.cycles = float(cycles)
        self.plastic_work = float(plastic_work)
        self.cycle_evaluator = cycle_evaluator

    def dormant_eligibility(self):
        s = self.pd_state
        checks = (
            (not np.any(np.asarray(s.site_status) == 1), "active_embryo"),
            (not np.any(np.asarray(s.site_status) == 2), "stabilized_site"),
            (not np.any(np.asarray(s.bond_damage) > 0.0), "bond_damage"),
            (not bool(s.active_front), "captured_front"),
            (int(s.primary_seed_node) < 0, "primary_seed"),
        )
        for passed, reason in checks:
            if not passed: return False, reason
        return True, "dormant_fixed_topology"

    def _active_arrays(self):
        return (self.ep_gp, self.rho_gp, self.epsp_acc_gp, self.u,
                np.asarray(self.pd_state.log_delivery_memory, float))

    def active_state(self):
        arrays = self._active_arrays()
        spec = tuple((name, a.shape, str(a.dtype)) for name, a in zip(ACTIVE_NAMES, arrays))
        vector = np.concatenate([a.ravel() for a in arrays]) if arrays else np.empty(0)
        # log memory uses -inf for exact zero. Map it reversibly to the smallest
        # finite working coordinate; no physical rate floor is introduced.
        vector = np.where(np.isneginf(vector), -1.0e300, vector)
        return ActiveState(vector, spec)

    def restore_active_state(self, snapshot, vector):
        x = np.asarray(vector, float)
        offset, restored = 0, []
        for _, shape, _ in snapshot.specification:
            size = int(np.prod(shape, dtype=int))
            value = x[offset:offset + size].reshape(shape).copy(); offset += size
            restored.append(value)
        if offset != x.size or len(restored) != 5:
            raise ValueError("active state layout mismatch")
        self.ep_gp, self.rho_gp, self.epsp_acc_gp, self.u, logmem = restored
        self.pd_state.log_delivery_memory = np.where(logmem <= -5.0e299, -math.inf, logmem)

    def protected_signatures(self):
        s = self.pd_state
        ledgers = {name: getattr(s, name) for name in LEDGER_NAMES}
        ledgers["plastic_work"] = self.plastic_work
        stochastic = {name: getattr(s, name) for name in STOCHASTIC_NAMES}
        stochastic["candidate_rng"] = self.patch._candidate_rng.bit_generator.state
        stochastic["event_rng"] = self.patch._event_rng.bit_generator.state
        topology = {name: getattr(s, name) for name in TOPOLOGY_NAMES}
        topology["mesh_nodes"] = np.asarray(self.mesh.nodes)
        topology["bonds"] = np.asarray(self.patch.bonds)
        return ProtectedSignatures(_digest(ledgers), _digest(stochastic), _digest(topology))

    def exact_private_cycle(self):
        clone = deepcopy(self)
        before = clone.protected_signatures()
        start = clone.active_state()
        payload = clone.cycle_evaluator(clone)
        end = clone.active_state()
        after = clone.protected_signatures()
        if before != after or clone.cycles != self.cycles:
            raise RuntimeError("cycle evaluator mutated protected state or physical cycles")
        log_action = np.asarray(payload["log_birth_action"], float)
        return CycleEvaluation(start, end, log_action,
            dict(payload.get("ledger_increments", {})),
            np.asarray(payload.get("phase", []), float),
            np.asarray(payload.get("phase_log_birth_rate", []), float),
            dict(payload.get("diagnostics", {})),
            str(payload.get("transition_signature", "dormant")), before.topology)

    def commit_private_cycle(self, evaluation):
        self.restore_active_state(evaluation.state_end, evaluation.state_end.vector)
        for name, increment in evaluation.ledger_increments.items():
            if name == "plastic_work": self.plastic_work += float(increment)
            elif hasattr(self.pd_state, name):
                setattr(self.pd_state, name, getattr(self.pd_state, name) + increment)

    def remaining_birth_actions(self):
        s = self.pd_state
        nodes = np.asarray(s.site_node_index, dtype=np.int64)
        available = np.asarray(s.site_status) == 0
        residual = np.full(len(nodes), math.inf)
        residual[available] = np.maximum(
            np.asarray(s.site_birth_threshold)[available]
            - np.asarray(s.birth_cumulative_hazard)[nodes[available]], 0.0)
        return residual

    def commit_birth_action(self, increment, cycles):
        inc = np.asarray(increment, float)
        s = self.pd_state
        if inc.shape == np.asarray(s.birth_cumulative_hazard).shape:
            node_inc = inc
        elif inc.shape == np.asarray(s.site_birth_threshold).shape:
            node_inc = np.zeros_like(s.birth_cumulative_hazard, dtype=float)
            np.maximum.at(node_inc, np.asarray(s.site_node_index, int), inc)
        else:
            raise ValueError("birth action shape does not match PD nodes or sites")
        s.birth_cumulative_hazard = np.asarray(s.birth_cumulative_hazard) + node_inc
        log_inc = np.where(node_inc > 0.0, np.log(node_inc), -math.inf)
        s.log_birth_cumulative_hazard = np.logaddexp(s.log_birth_cumulative_hazard, log_inc)

    def physical_cycles(self): return self.cycles
    def set_physical_cycles(self, cycles): self.cycles = float(cycles)


__all__ = ["SpatialPDDormantAdapter", "ACTIVE_NAMES", "LEDGER_NAMES",
           "STOCHASTIC_NAMES", "TOPOLOGY_NAMES"]
