"""Event-to-event high-cycle machinery for dormant Stateful-PD v9 states.

This module intentionally contains no sharp-front MPZ assumptions.  A physical
adapter supplies an exact, phase-resolved one-cycle map and the four disjoint
state inventories.  Periodic/projective trials operate on active continuous
coordinates only; ledgers, clocks/RNG, and topology are protected invariants.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import math
from typing import Any, Callable, Protocol
from pathlib import Path
import os

import numpy as np


MODEL_ID = "v9_stateful_pd_dormant_high_cycle_event_engine_v1"


def _digest(value: Any) -> str:
    """Stable digest for nested scalar/array state, including RNG dictionaries."""
    h = hashlib.sha256()

    def add(item: Any) -> None:
        if isinstance(item, np.ndarray):
            a = np.ascontiguousarray(item)
            h.update(b"array\0" + str(a.dtype).encode() + b"\0")
            h.update(json.dumps(a.shape).encode() + b"\0" + a.view(np.uint8).tobytes())
        elif isinstance(item, dict):
            h.update(b"dict\0")
            for key in sorted(item, key=str):
                add(str(key)); add(item[key])
        elif isinstance(item, (list, tuple)):
            h.update(b"seq\0")
            for child in item: add(child)
        elif isinstance(item, np.generic):
            add(item.item())
        elif isinstance(item, float):
            h.update(np.float64(item).tobytes())
        else:
            h.update(repr(item).encode() + b"\0")

    add(value)
    return h.hexdigest()


@dataclass(frozen=True)
class ActiveState:
    vector: np.ndarray
    specification: tuple[tuple[str, tuple[int, ...], str], ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "vector", np.asarray(self.vector, float).copy())
        if not np.all(np.isfinite(self.vector)):
            raise ValueError("active high-cycle coordinates must be finite")


@dataclass(frozen=True)
class ProtectedSignatures:
    ledgers: str
    stochastic: str
    topology: str


@dataclass
class CycleEvaluation:
    state_start: ActiveState
    state_end: ActiveState
    log_birth_action: np.ndarray
    ledger_increments: dict[str, float]
    phase: np.ndarray
    phase_log_birth_rate: np.ndarray
    diagnostics: dict[str, Any]
    transition_signature: str
    topology_signature: str
    private_invariants_preserved: bool = True

    @property
    def birth_action(self) -> np.ndarray:
        out = np.zeros_like(self.log_birth_action, dtype=float)
        finite = np.isfinite(self.log_birth_action)
        out[finite] = np.exp(np.minimum(self.log_birth_action[finite], 709.0))
        return out


class DormantPDAdapter(Protocol):
    def dormant_eligibility(self) -> tuple[bool, str]: ...
    def active_state(self) -> ActiveState: ...
    def restore_active_state(self, snapshot: ActiveState, vector: np.ndarray) -> None: ...
    def protected_signatures(self) -> ProtectedSignatures: ...
    def exact_private_cycle(self) -> CycleEvaluation: ...
    def commit_private_cycle(self, evaluation: CycleEvaluation) -> None: ...
    def remaining_birth_actions(self) -> np.ndarray: ...
    def commit_birth_action(self, increment: np.ndarray, cycles: float) -> None: ...
    def physical_cycles(self) -> float: ...
    def set_physical_cycles(self, cycles: float) -> None: ...


@dataclass(frozen=True)
class HighCycleConfig:
    periodic_relative_tolerance: float = 1e-9
    periodic_admission_distance: float = 1e-8
    periodic_max_iterations: int = 24
    projective_state_tolerance: float = 2e-5
    projective_log_hazard_tolerance: float = 2e-4
    projective_initial_cycles: int = 16
    projective_max_cycles: int = 10**9
    exact_retry_cycles: int = 8
    event_guard_cycles: float = 2.0
    minimum_positive_coordinate: float = -math.inf


@dataclass
class ModeRecord:
    mode: str
    cycles: float
    exact_map_evaluations: int
    accepted: bool
    detail: dict[str, Any] = field(default_factory=dict)


@dataclass
class AdvanceResult:
    cycles_consumed: float
    event_guard_reached: bool
    exact_map_evaluations: int
    accepted_projected_cycles: float
    modes: list[ModeRecord]
    cache_invalidated: bool = False


def relative_distance(a: np.ndarray, b: np.ndarray) -> float:
    aa, bb = np.asarray(a, float), np.asarray(b, float)
    scale = np.maximum.reduce([np.abs(aa), np.abs(bb), np.ones_like(aa)])
    return float(np.max(np.abs(aa - bb) / scale)) if aa.size else 0.0


def _assert_private(adapter: DormantPDAdapter, before: ProtectedSignatures) -> None:
    after = adapter.protected_signatures()
    if after != before:
        raise RuntimeError("private high-cycle operation mutated ledger, stochastic, or topology state")


def private_cycle(adapter: DormantPDAdapter) -> CycleEvaluation:
    eligible, reason = adapter.dormant_eligibility()
    if not eligible:
        raise RuntimeError(f"dormant high-cycle mode is ineligible: {reason}")
    protected = adapter.protected_signatures()
    cycles = adapter.physical_cycles()
    result = adapter.exact_private_cycle()
    _assert_private(adapter, protected)
    if adapter.physical_cycles() != cycles:
        raise RuntimeError("private one-cycle map consumed physical cycles")
    if result.topology_signature != protected.topology:
        raise RuntimeError("private one-cycle result has a changed topology signature")
    return result


def solve_periodic_state(adapter: DormantPDAdapter, cfg: HighCycleConfig) -> tuple[ActiveState, float, int, int]:
    """Private fixed-point iteration; never adopts a remote fixed point."""
    initial = adapter.active_state()
    protected = adapter.protected_signatures()
    cycles = adapter.physical_cycles()
    x = initial.vector.copy()
    maps = 0
    residual = math.inf
    try:
        for iteration in range(1, cfg.periodic_max_iterations + 1):
            adapter.restore_active_state(initial, x)
            ev = private_cycle(adapter); maps += 1
            y = ev.state_end.vector
            residual = relative_distance(x, y)
            if residual <= cfg.periodic_relative_tolerance:
                return ev.state_end, residual, iteration, maps
            x = y
        return ActiveState(x, initial.specification), residual, cfg.periodic_max_iterations, maps
    finally:
        adapter.restore_active_state(initial, initial.vector)
        adapter.set_physical_cycles(cycles)
        _assert_private(adapter, protected)


def _cycles_before_guard(remaining_action: np.ndarray, per_cycle: np.ndarray, guard: float) -> float:
    valid = (per_cycle > 0.0) & np.isfinite(remaining_action)
    if not np.any(valid):
        return math.inf
    waits = np.maximum(remaining_action[valid], 0.0) / per_cycle[valid]
    return max(float(np.min(waits)) - guard, 0.0)


def _hazard_error(predicted: np.ndarray, exact: np.ndarray) -> float:
    p, e = np.asarray(predicted, float), np.asarray(exact, float)
    mask = np.isfinite(p) | np.isfinite(e)
    if not np.any(mask): return 0.0
    return float(np.max(np.abs(np.where(np.isfinite(p), p, -1000.0)[mask] - np.where(np.isfinite(e), e, -1000.0)[mask])))


class DormantPDHighCycleEngine:
    def __init__(self, adapter: DormantPDAdapter, config: HighCycleConfig | None = None):
        self.adapter = adapter
        self.config = config or HighCycleConfig()
        self.mode_history: list[ModeRecord] = []
        self.exact_map_evaluations = 0
        self.accepted_projected_cycles = 0.0
        self._cache_topology = adapter.protected_signatures().topology

    def invalidate_if_needed(self) -> bool:
        eligible, _ = self.adapter.dormant_eligibility()
        topology = self.adapter.protected_signatures().topology
        changed = (not eligible) or topology != self._cache_topology
        if changed:
            self._cache_topology = topology
            self.mode_history.append(ModeRecord("cache_invalidation", 0.0, 0, True))
        return changed

    def _stationary(self, requested: float, cycle: CycleEvaluation, distance: float) -> float:
        q = cycle.birth_action
        allowed = min(requested, _cycles_before_guard(
            self.adapter.remaining_birth_actions(), q, self.config.event_guard_cycles
        ))
        whole = float(max(math.floor(allowed), 0))
        if whole:
            self.adapter.commit_birth_action(q * whole, whole)
            self.adapter.set_physical_cycles(self.adapter.physical_cycles() + whole)
        self.mode_history.append(ModeRecord("stationary", whole, 0, True, {
            "admission_distance": distance,
            "maximum_action_per_cycle": float(np.max(q)) if q.size else 0.0,
        }))
        return whole

    def _private_at(self, template: ActiveState, vector: np.ndarray) -> CycleEvaluation:
        base = self.adapter.active_state()
        protected = self.adapter.protected_signatures()
        cycles = self.adapter.physical_cycles()
        try:
            self.adapter.restore_active_state(template, vector)
            ev = private_cycle(self.adapter)
            self.exact_map_evaluations += 1
            return ev
        finally:
            self.adapter.restore_active_state(base, base.vector)
            self.adapter.set_physical_cycles(cycles)
            _assert_private(self.adapter, protected)

    def _projective_trial(self, horizon: int) -> tuple[bool, dict[str, Any]]:
        start = self.adapter.active_state()
        first = private_cycle(self.adapter); self.exact_map_evaluations += 1
        drift = first.state_end.vector - start.vector
        midpoint = max(horizon // 2, 1)
        xmid = start.vector + midpoint * drift
        xend = start.vector + horizon * drift
        if np.any(~np.isfinite(xend)) or np.any(xend < self.config.minimum_positive_coordinate):
            return False, {"reason": "physical_bounds"}
        mid = self._private_at(start, xmid)
        end = self._private_at(start, xend)
        predicted_mid_next = xmid + drift
        predicted_end_next = xend + drift
        state_error = max(relative_distance(predicted_mid_next, mid.state_end.vector),
                          relative_distance(predicted_end_next, end.state_end.vector))
        hazard_error = max(_hazard_error(first.log_birth_action, mid.log_birth_action),
                           _hazard_error(first.log_birth_action, end.log_birth_action))
        signature_ok = first.transition_signature == mid.transition_signature == end.transition_signature
        accepted = state_error <= self.config.projective_state_tolerance and hazard_error <= self.config.projective_log_hazard_tolerance and signature_ok
        q0, qm, q1 = first.birth_action, mid.birth_action, end.birth_action
        integrated = horizon * (q0 + 4.0 * qm + q1) / 6.0
        upper = horizon * np.maximum.reduce([q0, qm, q1])
        event_safe = np.all(upper < self.adapter.remaining_birth_actions())
        accepted = bool(accepted and event_safe)
        return accepted, {"start": start, "end_vector": xend, "action": integrated,
                          "state_error": state_error, "hazard_error": hazard_error,
                          "transition_preserved": signature_ok, "event_safe": bool(event_safe)}

    def advance(self, cycles_requested: float) -> AdvanceResult:
        requested = max(float(cycles_requested), 0.0)
        consumed = 0.0
        invalid = self.invalidate_if_needed()
        if invalid:
            return AdvanceResult(0.0, False, self.exact_map_evaluations,
                                 self.accepted_projected_cycles, self.mode_history, True)
        while consumed < requested:
            ev = private_cycle(self.adapter); self.exact_map_evaluations += 1
            periodic, residual, iterations, maps = solve_periodic_state(self.adapter, self.config)
            self.exact_map_evaluations += maps
            distance = relative_distance(self.adapter.active_state().vector, periodic.vector)
            self.mode_history.append(ModeRecord("periodic_search", 0.0, maps,
                residual <= self.config.periodic_relative_tolerance,
                {"residual": residual, "iterations": iterations, "distance": distance}))
            remaining = requested - consumed
            if residual <= self.config.periodic_relative_tolerance and distance <= self.config.periodic_admission_distance:
                advanced = self._stationary(remaining, ev, distance)
                consumed += advanced
                if advanced >= remaining: break
                self.mode_history.append(ModeRecord("event_guard", 0.0, 0, True))
                break

            proposal = min(int(remaining), self.config.projective_initial_cycles)
            accepted = False
            while proposal >= 2:
                accepted, trial = self._projective_trial(proposal)
                if accepted:
                    self.adapter.restore_active_state(trial["start"], trial["end_vector"])
                    self.adapter.commit_birth_action(trial["action"], proposal)
                    self.adapter.set_physical_cycles(self.adapter.physical_cycles() + proposal)
                    consumed += proposal
                    self.accepted_projected_cycles += proposal
                    self.mode_history.append(ModeRecord("projective", proposal, 3, True,
                        {k: v for k, v in trial.items() if k not in {"start", "end_vector", "action"}}))
                    accepted = True
                    break
                self.mode_history.append(ModeRecord("projective_reject", 0.0, 3, False,
                    {k: v for k, v in trial.items() if k not in {"start", "end_vector", "action"}}))
                proposal //= 2
            if accepted: continue

            burst = min(self.config.exact_retry_cycles, int(math.floor(remaining)))
            if burst <= 0: break
            for _ in range(burst):
                cycle = private_cycle(self.adapter); self.exact_map_evaluations += 1
                if _cycles_before_guard(self.adapter.remaining_birth_actions(), cycle.birth_action,
                                        self.config.event_guard_cycles) <= 0.0:
                    self.mode_history.append(ModeRecord("event_guard", 0.0, 1, True))
                    return AdvanceResult(consumed, True, self.exact_map_evaluations,
                                         self.accepted_projected_cycles, self.mode_history)
                self.adapter.commit_private_cycle(cycle)
                self.adapter.commit_birth_action(cycle.birth_action, 1.0)
                self.adapter.set_physical_cycles(self.adapter.physical_cycles() + 1.0)
                consumed += 1.0
            self.mode_history.append(ModeRecord("exact_burst", burst, burst, True))
        return AdvanceResult(consumed, any(m.mode == "event_guard" for m in self.mode_history),
                             self.exact_map_evaluations, self.accepted_projected_cycles,
                             self.mode_history)

    def write_atomic_mode_checkpoint(self, path: str | Path) -> None:
        """Write controller/mode state atomically at an accepted boundary.

        Complete physical arrays and clocks remain in the driver's atomic v9
        checkpoint.  This sidecar is deliberately diagnostic/controller state
        only and is safe to reconstruct after a crash.
        """
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "model_id": MODEL_ID,
            "physical_cycles": self.adapter.physical_cycles(),
            "topology_signature": self.adapter.protected_signatures().topology,
            "exact_map_evaluations": self.exact_map_evaluations,
            "accepted_projected_cycles": self.accepted_projected_cycles,
            "mode_history": [
                {"mode": m.mode, "cycles": m.cycles,
                 "exact_map_evaluations": m.exact_map_evaluations,
                 "accepted": m.accepted, "detail": m.detail}
                for m in self.mode_history
            ],
        }
        tmp = target.with_name(target.name + ".tmp")
        tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        os.replace(tmp, target)


def state_inventory() -> dict[str, tuple[str, ...]]:
    """Authoritative first-milestone inventory for the v9 spatial-PD adapter."""
    return {
        "active_continuous": ("fem.ep_gp", "fem.rho_gp", "fem.epsp_acc_gp", "fem.u",
                              "pd.log_delivery_memory"),
        "monotone_ledgers": ("fem.plastic_work", "pd.born_cumulative",
                             "pd.healed_cumulative", "pd.born_sites_cumulative",
                             "pd.healed_sites_cumulative", "pd.log_birth_cumulative_hazard"),
        "persistent_stochastic": ("pd.site_birth_threshold", "pd.birth_cumulative_hazard",
                                  "pd.site_transition_threshold", "pd.site_transition_cumulative_hazard",
                                  "pd.site_transition_outcome_uniform", "candidate_rng", "event_rng"),
        "discrete_topology": ("pd.site_status", "pd.bond_damage", "pd.active_front",
                              "pd.front_masks", "pd.front_path", "mesh.nodes", "bond_connectivity"),
    }


__all__ = ["MODEL_ID", "ActiveState", "ProtectedSignatures", "CycleEvaluation",
           "HighCycleConfig", "ModeRecord", "AdvanceResult", "DormantPDHighCycleEngine",
           "private_cycle", "solve_periodic_state", "state_inventory", "_digest"]
