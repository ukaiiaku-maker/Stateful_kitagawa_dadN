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
    ledger_increments: dict[str, Any]
    phase: np.ndarray
    phase_log_birth_rate: np.ndarray
    diagnostics: dict[str, Any]
    transition_signature: str
    topology_signature: str
    private_invariants_preserved: bool = True
    # Multiplying a phase rate by this value converts its time basis to one
    # physical cycle.  Canonical cleavage rates are per second, whereas older
    # synthetic adapters expose rates per cycle and retain the default 1.0.
    phase_rate_seconds_per_cycle: float = 1.0

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
    def commit_ledger_increments(self, increments: dict[str, Any]) -> None: ...
    def remaining_birth_actions(self) -> np.ndarray: ...
    def commit_birth_action(self, increment: np.ndarray, cycles: float) -> None: ...
    def commit_log_birth_action(self, log_increment: np.ndarray, cycles: float) -> None: ...
    def physical_cycles(self) -> float: ...
    def set_physical_cycles(self, cycles: float) -> None: ...
    def active_residual(self, a: ActiveState, b: ActiveState) -> tuple[float, dict[str, float]]: ...


@dataclass(frozen=True)
class HighCycleConfig:
    periodic_relative_tolerance: float = 1e-9
    periodic_admission_distance: float = 1e-8
    periodic_max_iterations: int = 24
    projective_state_tolerance: float = 2e-5
    projective_log_hazard_tolerance: float = 2e-4
    projective_initial_cycles: int = 16
    projective_max_cycles: int = 10**9
    projective_growth_factor: float = 4.0
    projective_curvature_tolerance: float = 2e-5
    exact_retry_cycles: int = 8
    event_guard_cycles: float = 2.0
    minimum_positive_coordinate: float = -math.inf
    max_exact_map_evaluations: int = 128
    minimum_projected_cycles_per_exact_map: float = 0.0
    private_window_training: bool = True
    private_window_initial_cycles: int = 64
    private_window_max_cycles: int = 10**8
    private_window_state_tolerance: float = 2e-7
    private_window_log_hazard_tolerance: float = 2e-4


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


def active_distance(adapter: DormantPDAdapter, a: ActiveState, b: ActiveState) -> tuple[float, dict[str, float]]:
    method = getattr(adapter, "active_residual", None)
    if method is not None:
        return method(a, b)
    value = relative_distance(a.vector, b.vector)
    return value, {"active_vector": value}


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
            residual, _ = active_distance(adapter, ActiveState(x, initial.specification), ev.state_end)
            if residual <= cfg.periodic_relative_tolerance:
                return ev.state_end, residual, iteration, maps
            x = y
        return ActiveState(x, initial.specification), residual, cfg.periodic_max_iterations, maps
    finally:
        adapter.restore_active_state(initial, initial.vector)
        adapter.set_physical_cycles(cycles)
        _assert_private(adapter, protected)


def _log_remaining_action(remaining_action: np.ndarray) -> np.ndarray:
    remaining = np.asarray(remaining_action, float)
    return np.where(remaining > 0.0, np.log(remaining), -math.inf)


def _cycles_before_guard_log(remaining_action: np.ndarray, log_per_cycle: np.ndarray, guard: float) -> float:
    log_q = np.asarray(log_per_cycle, float)
    log_r = _log_remaining_action(remaining_action)
    valid = np.isfinite(log_q) & np.isfinite(log_r)
    if not np.any(valid):
        return math.inf
    log_wait = float(np.min(log_r[valid] - log_q[valid]))
    wait = math.inf if log_wait > math.log(np.finfo(float).max) else math.exp(log_wait)
    return max(wait - guard, 0.0)


def _log_weighted_sum(log_values: tuple[np.ndarray, ...], weights: tuple[float, ...]) -> np.ndarray:
    terms = [np.asarray(value, float) + math.log(weight) for value, weight in zip(log_values, weights)]
    out = terms[0]
    for term in terms[1:]: out = np.logaddexp(out, term)
    return out


def _commit_log_action(adapter: DormantPDAdapter, log_increment: np.ndarray, cycles: float) -> None:
    method = getattr(adapter, "commit_log_birth_action", None)
    if method is not None:
        method(np.asarray(log_increment, float), cycles)
        return
    linear = np.where(np.asarray(log_increment) >= math.log(np.nextafter(0.0, 1.0)),
                      np.exp(np.asarray(log_increment)), 0.0)
    adapter.commit_birth_action(linear, cycles)


def _hazard_error(predicted: np.ndarray, exact: np.ndarray) -> float:
    p, e = np.asarray(predicted, float), np.asarray(exact, float)
    mask = np.isfinite(p) | np.isfinite(e)
    if not np.any(mask): return 0.0
    return float(np.max(np.abs(np.where(np.isfinite(p), p, -1000.0)[mask] - np.where(np.isfinite(e), e, -1000.0)[mask])))


def _linear_log_prediction(start: np.ndarray, second: np.ndarray, cycles: int) -> np.ndarray:
    """Predict a log rate from the first two trained exact-cycle samples."""
    a, b = np.asarray(start, float), np.asarray(second, float)
    out = a.copy()
    finite = np.isfinite(a) & np.isfinite(b)
    out[finite] = a[finite] + int(cycles) * (b[finite] - a[finite])
    # A rate absent in both samples remains absent.  A newly appearing rate is
    # deliberately not extrapolated and will fail exact validation if relevant.
    out[~finite & np.isneginf(a) & np.isneginf(b)] = -math.inf
    return out


def _log_geometric_cycle_sum(log_start: np.ndarray, log_second: np.ndarray, cycles: int) -> np.ndarray:
    """Log of sum(q_0 ... q_{cycles-1}) for a trained geometric rate."""
    q0, q1 = np.asarray(log_start, float), np.asarray(log_second, float)
    out = np.full(np.broadcast_shapes(q0.shape, q1.shape), -math.inf, dtype=float)
    finite = np.isfinite(q0) & np.isfinite(q1)
    if not np.any(finite) or cycles <= 0:
        return out
    slope = q1[finite] - q0[finite]
    value = np.empty_like(slope)
    near = np.abs(slope) < 1e-12
    value[near] = q0[finite][near] + math.log(cycles)
    pos = (~near) & (slope > 0.0)
    value[pos] = (q0[finite][pos] + (cycles - 1) * slope[pos]
                  + np.log(-np.expm1(-cycles * slope[pos]))
                  - np.log(-np.expm1(-slope[pos])))
    neg = (~near) & (slope < 0.0)
    value[neg] = (q0[finite][neg]
                  + np.log(-np.expm1(cycles * slope[neg]))
                  - np.log(-np.expm1(slope[neg])))
    out[finite] = value
    return out


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
        log_q = np.asarray(cycle.log_birth_action, float)
        allowed = min(requested, _cycles_before_guard_log(
            self.adapter.remaining_birth_actions(), log_q, self.config.event_guard_cycles
        ))
        whole = float(max(math.floor(allowed), 0))
        if whole:
            _commit_log_action(self.adapter, log_q + math.log(whole), whole)
            for name, rate in cycle.ledger_increments.items():
                if np.any(np.asarray(rate, float) < 0.0):
                    raise RuntimeError(f"monotone ledger {name} has negative stationary rate")
            self.adapter.commit_ledger_increments(
                {name: np.asarray(rate) * whole for name, rate in cycle.ledger_increments.items()}
            )
            self.adapter.set_physical_cycles(self.adapter.physical_cycles() + whole)
        self.mode_history.append(ModeRecord("stationary", whole, 0, True, {
            "admission_distance": distance,
            "maximum_log_action_per_cycle": float(np.max(log_q)) if log_q.size else -math.inf,
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

    def _private_window_at(self, template: ActiveState, vector: np.ndarray,
                           cycles: float, *, start_cycles: float) -> CycleEvaluation:
        method = getattr(self.adapter, "exact_private_window", None)
        if method is None:
            raise RuntimeError("adapter has no exact private-window evaluator")
        base = self.adapter.active_state()
        protected = self.adapter.protected_signatures()
        physical_cycles = self.adapter.physical_cycles()
        try:
            self.adapter.restore_active_state(template, vector)
            self.adapter.set_physical_cycles(float(start_cycles))
            evaluation = method(float(cycles))
            self.exact_map_evaluations += 1
            return evaluation
        finally:
            self.adapter.restore_active_state(base, base.vector)
            self.adapter.set_physical_cycles(physical_cycles)
            _assert_private(self.adapter, protected)

    def _private_window_trial(self, horizon: int) -> tuple[bool, dict[str, Any]]:
        """Qualify one exact macro-window against an independently partitioned map."""
        if horizon < 2:
            return False, {"reason": "window_too_short"}
        start = self.adapter.active_state()
        start_cycles = self.adapter.physical_cycles()
        left_cycles = horizon // 2
        right_cycles = horizon - left_cycles
        full = self._private_window_at(
            start, start.vector, horizon, start_cycles=start_cycles
        )
        left = self._private_window_at(
            start, start.vector, left_cycles, start_cycles=start_cycles
        )
        right = self._private_window_at(
            start, left.state_end.vector, right_cycles,
            start_cycles=start_cycles + left_cycles,
        )
        state_error, state_error_by_field = active_distance(
            self.adapter, full.state_end, right.state_end
        )
        partition_log_action = np.logaddexp(left.log_birth_action, right.log_birth_action)
        hazard_error = _hazard_error(full.log_birth_action, partition_log_action)
        names = set(full.ledger_increments) | set(left.ledger_increments) | set(right.ledger_increments)
        ledgers = dict(full.ledger_increments)
        ledger_error = 0.0
        ledger_error_by_name = {}
        for name in names:
            whole = np.asarray(full.ledger_increments.get(name, 0.0), float)
            split = (np.asarray(left.ledger_increments.get(name, 0.0), float)
                     + np.asarray(right.ledger_increments.get(name, 0.0), float))
            # These are cumulative population/work increments.  In a dormant
            # tail their values can be far below one, where a relative error
            # divided by the tiny increment is meaningless.  Normalize by one
            # physical ledger unit while retaining relative control above it.
            scale = max(float(np.max(np.abs(whole))), float(np.max(np.abs(split))), 1.0)
            error = float(np.max(np.abs(whole - split))) / scale
            ledger_error_by_name[name] = error
            ledger_error = max(ledger_error, error)
        transition_ok = (full.transition_signature == left.transition_signature
                         == right.transition_signature)
        remaining_log = _log_remaining_action(self.adapter.remaining_birth_actions())
        phase_rates = [np.asarray(ev.phase_log_birth_rate, float).ravel()
                       for ev in (full, left, right) if np.asarray(ev.phase_log_birth_rate).size]
        max_log_rate = max((float(np.max(rate)) for rate in phase_rates), default=-math.inf)
        time_scales = [float(ev.phase_rate_seconds_per_cycle) for ev in (full, left, right)]
        if any(not np.isfinite(scale) or scale <= 0.0 for scale in time_scales):
            return False, {"reason": "invalid_phase_rate_seconds_per_cycle"}
        max_time_scale = max(time_scales)
        # The action extrapolation has been checked against the split exact
        # window, but a conservative crossing guard must also cover the
        # measured validation residual.  hazard_error is already logarithmic;
        # log1p(state_error) makes the normalized state mismatch an additional
        # nonnegative multiplicative allowance.
        validation_log_margin = hazard_error + math.log1p(state_error)
        guard_log = (-math.inf if not np.isfinite(max_log_rate) else
                     max_log_rate
                     + math.log(max(self.config.event_guard_cycles * max_time_scale, 1e-300))
                     + validation_log_margin)
        guarded_action = np.logaddexp(full.log_birth_action, guard_log)
        event_safe = bool(np.all(guarded_action < remaining_log))
        efficiency = horizon / 3.0
        accepted = bool(
            state_error <= self.config.private_window_state_tolerance
            and hazard_error <= self.config.private_window_log_hazard_tolerance
            and ledger_error <= self.config.private_window_log_hazard_tolerance
            and transition_ok and event_safe
            and efficiency >= self.config.minimum_projected_cycles_per_exact_map
        )
        reason = "accepted" if accepted else (
            "event_guard" if not event_safe else
            "insufficient_window_efficiency" if efficiency < self.config.minimum_projected_cycles_per_exact_map else
            "partition_mismatch"
        )
        return accepted, {
            "reason": reason, "start": start, "end_vector": full.state_end.vector,
            "log_action": full.log_birth_action, "ledgers": ledgers,
            "state_error": state_error, "state_error_by_field": state_error_by_field,
            "log_hazard_error": hazard_error, "ledger_error": ledger_error,
            "phase_rate_seconds_per_cycle": max_time_scale,
            "validation_log_margin": validation_log_margin,
            "guard_log_action_upper_bound": guard_log,
            "ledger_error_by_name": ledger_error_by_name,
            "transition_preserved": transition_ok, "event_safe": event_safe,
            "projected_cycles_per_exact_map": efficiency,
            "partition_cycles": [left_cycles, right_cycles],
        }

    def _projective_trial(self, horizon: int) -> tuple[bool, dict[str, Any]]:
        start = self.adapter.active_state()
        first = private_cycle(self.adapter); self.exact_map_evaluations += 1
        drift = first.state_end.vector - start.vector
        second = self._private_at(start, first.state_end.vector)
        second_drift = second.state_end.vector - first.state_end.vector
        curvature_metric = getattr(self.adapter, "projective_curvature", None)
        curvature = (float(curvature_metric(start, first.state_end, second.state_end))
                     if curvature_metric is not None else relative_distance(drift, second_drift))
        midpoint = max(horizon // 2, 1)
        projector = getattr(self.adapter, "project_active_state", None)
        if projector is None:
            project = lambda h: start.vector + h * drift
        else:
            project = lambda h: np.asarray(
                projector(start, first.state_end, second.state_end, h), float
            )
        xmid = project(midpoint)
        xend = project(horizon)
        if np.any(~np.isfinite(xend)) or np.any(xend < self.config.minimum_positive_coordinate):
            return False, {"reason": "physical_bounds"}
        try:
            mid = self._private_at(start, xmid)
            end = self._private_at(start, xend)
        except ValueError as exc:
            return False, {"reason": "projected_constitutive_domain", "detail": str(exc)}
        predicted_mid_next = project(midpoint + 1)
        predicted_end_next = project(horizon + 1)
        state_error = max(relative_distance(predicted_mid_next, mid.state_end.vector),
                          relative_distance(predicted_end_next, end.state_end.vector))
        predicted_mid_log = _linear_log_prediction(
            first.log_birth_action, second.log_birth_action, midpoint
        )
        predicted_end_log = _linear_log_prediction(
            first.log_birth_action, second.log_birth_action, horizon
        )
        hazard_error = max(_hazard_error(predicted_mid_log, mid.log_birth_action),
                           _hazard_error(predicted_end_log, end.log_birth_action))
        signature_ok = first.transition_signature == mid.transition_signature == end.transition_signature
        accepted = (state_error <= self.config.projective_state_tolerance
                    and hazard_error <= self.config.projective_log_hazard_tolerance
                    and curvature <= self.config.projective_curvature_tolerance
                    and signature_ok)
        log_integrated = _log_geometric_cycle_sum(
            first.log_birth_action, second.log_birth_action, horizon
        )
        log_upper = math.log(horizon) + np.maximum.reduce(
            [first.log_birth_action, mid.log_birth_action, end.log_birth_action]
        )
        event_safe = np.all(log_upper < _log_remaining_action(self.adapter.remaining_birth_actions()))
        ledger_names = set(first.ledger_increments) | set(mid.ledger_increments) | set(end.ledger_increments)
        ledgers = {}
        ledger_variation = 0.0
        ledger_prediction_error = 0.0
        ledger_prediction_error_by_name = {}
        for name in ledger_names:
            rates = np.stack([np.asarray(first.ledger_increments.get(name, 0.0),float),
                              np.asarray(mid.ledger_increments.get(name, 0.0),float),
                              np.asarray(end.ledger_increments.get(name, 0.0),float)])
            if np.any(rates < 0.0):
                return False, {"reason": f"negative_monotone_ledger:{name}"}
            r0 = np.asarray(first.ledger_increments.get(name, 0.0), float)
            r1 = np.asarray(second.ledger_increments.get(name, 0.0), float)
            positive_pair = (r0 > 0.0) & (r1 > 0.0)
            log_r0 = np.full_like(r0, -math.inf, dtype=float)
            log_r1 = np.full_like(r1, -math.inf, dtype=float)
            np.log(r0, out=log_r0, where=positive_pair)
            np.log(r1, out=log_r1, where=positive_pair)
            log_ledger_sum = _log_geometric_cycle_sum(
                log_r0, log_r1, horizon,
            )
            geometric_sum = np.where(np.isfinite(log_ledger_sum), np.exp(np.minimum(log_ledger_sum, 709.0)), 0.0)
            # Exactly zero rates stay zero. A newly appearing/disappearing rate
            # uses sampled Simpson integration and must still pass validation.
            unresolved = (r0 == 0.0) ^ (r1 == 0.0)
            if np.any(unresolved):
                geometric_sum = np.asarray(geometric_sum)
                geometric_sum[unresolved] = horizon * (
                    rates[0][unresolved] + 4.0*rates[1][unresolved] + rates[2][unresolved]
                ) / 6.0
            ledgers[name] = geometric_sum
            temporal_range = np.ptp(rates, axis=0)
            ledger_variation = max(
                ledger_variation,
                float(np.max(temporal_range) / max(float(np.max(np.abs(rates))), 1e-300)),
            )
            ratio = np.ones_like(r0, dtype=float)
            np.divide(r1, r0, out=ratio, where=r0 > 0.0)
            predicted_mid_rate = np.where(positive_pair, r0 * np.power(ratio, midpoint), r0 + midpoint * (r1 - r0))
            predicted_end_rate = np.where(positive_pair, r0 * np.power(ratio, horizon), r0 + horizon * (r1 - r0))
            scale = max(float(np.max(np.abs(rates))), float(np.max(np.abs(r0))),
                        float(np.max(np.abs(r1))), 1e-300)
            name_error = max(
                float(np.max(np.abs(predicted_mid_rate - rates[1]))) / scale,
                float(np.max(np.abs(predicted_end_rate - rates[2]))) / scale,
            )
            ledger_prediction_error_by_name[name] = name_error
            ledger_prediction_error = max(ledger_prediction_error, name_error)
        conservative = getattr(self.adapter, "conservative_population_ledgers", None)
        if conservative is not None:
            try:
                ledgers, constrained_names = conservative(start, xend, ledgers)
            except ValueError as exc:
                return False, {"reason": "projected_population_ledger_domain", "detail": str(exc)}
            for name in constrained_names:
                ledger_prediction_error_by_name[name] = 0.0
            ledger_prediction_error = max(ledger_prediction_error_by_name.values(), default=0.0)
        accepted = bool(accepted and ledger_prediction_error <= self.config.projective_log_hazard_tolerance)
        accepted = bool(accepted and event_safe)
        return accepted, {"start": start, "end_vector": xend, "log_action": log_integrated,
                          "ledgers": ledgers, "ledger_rate_variation": ledger_variation,
                          "ledger_prediction_error": ledger_prediction_error,
                          "ledger_prediction_error_by_name": ledger_prediction_error_by_name,
                          "state_error": state_error, "hazard_error": hazard_error,
                          "curvature": curvature,
                          "transition_preserved": signature_ok, "event_safe": bool(event_safe)}

    def advance(self, cycles_requested: float) -> AdvanceResult:
        requested = max(float(cycles_requested), 0.0)
        consumed = 0.0
        next_projective = max(int(self.config.projective_initial_cycles), 2)
        invalid = self.invalidate_if_needed()
        if invalid:
            return AdvanceResult(0.0, False, self.exact_map_evaluations,
                                 self.accepted_projected_cycles, self.mode_history, True)
        while consumed < requested:
            if self.exact_map_evaluations >= self.config.max_exact_map_evaluations:
                self.mode_history.append(ModeRecord(
                    "efficiency_budget", 0.0, 0, False,
                    {"exact_map_evaluations": self.exact_map_evaluations,
                     "accepted_projected_cycles": self.accepted_projected_cycles},
                ))
                break
            remaining = requested - consumed
            window_method = getattr(self.adapter, "exact_private_window", None)
            if self.config.private_window_training and window_method is not None and remaining >= 2:
                proposal = min(int(remaining), int(self.config.private_window_max_cycles))
                proposal = max(min(proposal, max(self.config.private_window_initial_cycles, 2)), 2)
                accepted_window = False
                while proposal >= 2 and self.exact_map_evaluations + 3 <= self.config.max_exact_map_evaluations:
                    accepted_window, window = self._private_window_trial(proposal)
                    public_detail = {k: v for k, v in window.items()
                                     if k not in {"start", "end_vector", "log_action", "ledgers"}}
                    self.mode_history.append(ModeRecord(
                        "exact_private_window" if accepted_window else "exact_private_window_reject",
                        proposal if accepted_window else 0.0, 3, accepted_window, public_detail,
                    ))
                    if accepted_window:
                        self.adapter.restore_active_state(window["start"], window["end_vector"])
                        _commit_log_action(self.adapter, window["log_action"], proposal)
                        self.adapter.commit_ledger_increments(window["ledgers"])
                        self.adapter.set_physical_cycles(self.adapter.physical_cycles() + proposal)
                        consumed += proposal
                        self.accepted_projected_cycles += proposal
                        break
                    proposal //= 2
                if consumed >= requested:
                    break
                if accepted_window:
                    continue
            ev = private_cycle(self.adapter); self.exact_map_evaluations += 1
            current_residual, current_fields = active_distance(self.adapter, ev.state_start, ev.state_end)
            if current_residual <= self.config.periodic_admission_distance:
                periodic, residual, iterations, maps = solve_periodic_state(self.adapter, self.config)
                self.exact_map_evaluations += maps
                distance, field_distance = active_distance(self.adapter, self.adapter.active_state(), periodic)
                verified = self._private_at(periodic, periodic.vector)
                verified2 = self._private_at(periodic, verified.state_end.vector)
                verify_residual, verify_fields = active_distance(self.adapter, periodic, verified.state_end)
                verify_hazard = _hazard_error(verified.log_birth_action, verified2.log_birth_action)
                ledger_names = set(verified.ledger_increments) | set(verified2.ledger_increments)
                verify_ledger = 0.0
                for n in ledger_names:
                    a=np.asarray(verified.ledger_increments.get(n,0.0),float); b=np.asarray(verified2.ledger_increments.get(n,0.0),float)
                    verify_ledger=max(verify_ledger,float(np.max(np.abs(a-b)/np.maximum.reduce([np.abs(a),np.abs(b),np.full_like(a,1e-300)]))))
                stationary_verified = (verify_residual <= self.config.periodic_relative_tolerance
                                       and verify_hazard <= self.config.projective_log_hazard_tolerance
                                       and verify_ledger <= self.config.projective_log_hazard_tolerance
                                       and verified.transition_signature == verified2.transition_signature)
                self.mode_history.append(ModeRecord("periodic_search", 0.0, maps,
                    residual <= self.config.periodic_relative_tolerance,
                    {"residual": residual, "iterations": iterations, "distance": distance,
                     "field_distance": field_distance,"verification_residual":verify_residual,
                     "verification_field_residual":verify_fields,"verification_log_hazard_error":verify_hazard,
                     "verification_ledger_error":verify_ledger,"stationary_verified":stationary_verified}))
                if residual <= self.config.periodic_relative_tolerance and stationary_verified and distance <= self.config.periodic_admission_distance:
                    advanced = self._stationary(remaining, verified, distance)
                    consumed += advanced
                    if advanced >= remaining: break
                    self.mode_history.append(ModeRecord("event_guard", 0.0, 0, True))
                    break
            else:
                self.mode_history.append(ModeRecord("periodic_precheck",0.0,1,False,
                    {"current_residual":current_residual,"field_residual":current_fields}))

            proposal = min(int(remaining), next_projective, int(self.config.projective_max_cycles))
            accepted = False
            budget_exhausted = False
            while proposal >= 2:
                # One projective trial evaluates four exact maps: start,
                # second training cycle, projected midpoint, and endpoint.
                # Enforce the hard budget before starting a side-effect-free
                # trial so proposal halving can never overshoot it.
                if self.exact_map_evaluations + 4 > self.config.max_exact_map_evaluations:
                    self.mode_history.append(ModeRecord(
                        "efficiency_budget", 0.0, 0, False,
                        {"exact_map_evaluations": self.exact_map_evaluations,
                         "required_next_maps": 4,
                         "accepted_projected_cycles": self.accepted_projected_cycles},
                    ))
                    budget_exhausted = True
                    break
                accepted, trial = self._projective_trial(proposal)
                efficiency = proposal / 4.0
                if accepted and efficiency < self.config.minimum_projected_cycles_per_exact_map:
                    accepted = False
                    trial = dict(trial) | {
                        "reason": "insufficient_projective_efficiency",
                        "projected_cycles_per_exact_map": efficiency,
                        "minimum_projected_cycles_per_exact_map": self.config.minimum_projected_cycles_per_exact_map,
                    }
                if accepted:
                    self.adapter.restore_active_state(trial["start"], trial["end_vector"])
                    _commit_log_action(self.adapter, trial["log_action"], proposal)
                    self.adapter.commit_ledger_increments(trial["ledgers"])
                    self.adapter.set_physical_cycles(self.adapter.physical_cycles() + proposal)
                    consumed += proposal
                    self.accepted_projected_cycles += proposal
                    next_projective = min(int(max(proposal + 1, proposal * self.config.projective_growth_factor)),
                                          int(self.config.projective_max_cycles))
                    self.mode_history.append(ModeRecord("projective", proposal, 4, True,
                        {k: v for k, v in trial.items() if k not in {"start", "end_vector", "log_action", "ledgers"}}))
                    accepted = True
                    break
                self.mode_history.append(ModeRecord("projective_reject", 0.0, 4, False,
                    {k: v for k, v in trial.items() if k not in {"start", "end_vector", "log_action", "ledgers"}}))
                proposal //= 2
                next_projective = max(proposal, 2)
            if accepted: continue
            if budget_exhausted: break

            burst = min(self.config.exact_retry_cycles, int(math.floor(remaining)))
            if burst <= 0: break
            for _ in range(burst):
                cycle = private_cycle(self.adapter); self.exact_map_evaluations += 1
                if _cycles_before_guard_log(self.adapter.remaining_birth_actions(), cycle.log_birth_action,
                                        self.config.event_guard_cycles) <= 0.0:
                    self.mode_history.append(ModeRecord("event_guard", 0.0, 1, True))
                    return AdvanceResult(consumed, True, self.exact_map_evaluations,
                                         self.accepted_projected_cycles, self.mode_history)
                self.adapter.commit_private_cycle(cycle)
                _commit_log_action(self.adapter, cycle.log_birth_action, 1.0)
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
    """Authoritative executable inventory for the v9 spatial-PD adapter."""
    return {
        "active_continuous": ("fem.ep_gp", "fem.rho_gp", "fem.epsp_acc_gp", "fem.u",
                              "pd.log_delivery_memory", "pd.available", "pd.embryo",
                              "pd.stable", "pd.inactive", "pd.completion"),
        "warm_start_only": ("fem.u",),
        "reconstructed_cycle_diagnostics": ("pd.completion",),
        "monotone_ledgers": ("fem.plastic_work", "pd.born_cumulative",
                             "pd.healed_cumulative", "pd.born_sites_cumulative",
                             "pd.healed_sites_cumulative", "pd.log_birth_cumulative_hazard",
                             "adapter.born_expectation", "adapter.healed_expectation"),
        "persistent_stochastic": ("pd.site_birth_threshold", "pd.birth_cumulative_hazard",
                                  "pd.site_transition_threshold", "pd.site_transition_cumulative_hazard",
                                  "pd.site_transition_outcome_uniform", "candidate_rng", "event_rng"),
        "discrete_topology": ("pd.site_status", "pd.bond_damage", "pd.primary_seed_node",
                              "pd.active_front", "pd.active_front_bonds", "pd.front_backbone_bonds",
                              "pd.front_wake_bonds", "pd.front_process_bonds",
                              "pd.active_front_path_xy", "mesh.nodes", "bond_connectivity"),
    }


__all__ = ["MODEL_ID", "ActiveState", "ProtectedSignatures", "CycleEvaluation",
           "HighCycleConfig", "ModeRecord", "AdvanceResult", "DormantPDHighCycleEngine",
           "private_cycle", "solve_periodic_state", "state_inventory", "_digest"]
