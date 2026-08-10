from __future__ import annotations

from copy import deepcopy
import math

import numpy as np

from arrhenius_fracture.v9_pd_high_cycle import (
    ActiveState, CycleEvaluation, DormantPDHighCycleEngine, HighCycleConfig,
    ProtectedSignatures, _digest, private_cycle,
)


class SyntheticDormantPD:
    """Small protocol implementation with exact linear/contractive dynamics."""
    def __init__(self, *, rate=1e-14, threshold=2.0, contraction=0.0, drift=0.0):
        self.x = np.array([0.0])
        self.rate = float(rate)
        self.log_rate = math.log(rate) if rate > 0.0 else -math.inf
        self.threshold = np.array([float(threshold)])
        self.action = np.zeros(1)
        self.log_action = np.full(1, -math.inf)
        self.cycles = 0.0
        self.ledger = 0.0
        self.rng = np.random.default_rng(19)
        self.status = np.zeros(1, dtype=np.uint8)
        self.topology = np.array([0, 1], dtype=np.int64)
        self.contraction = float(contraction)
        self.drift = float(drift)

    def dormant_eligibility(self):
        return (not np.any(self.status), "dormant" if not np.any(self.status) else "event")
    def active_state(self): return ActiveState(self.x, (("x", (1,), "float64"),))
    def restore_active_state(self, snapshot, vector): self.x = np.asarray(vector, float).copy()
    def active_residual(self, a, b):
        value=float(np.max(np.abs(a.vector-b.vector)))
        return value,{"x":value}
    def protected_signatures(self):
        return ProtectedSignatures(_digest(self.ledger), _digest((self.threshold, self.action, self.rng.bit_generator.state)), _digest((self.status, self.topology)))
    def exact_private_cycle(self):
        start = self.active_state()
        end = ActiveState(self.contraction * self.x + self.drift, start.specification)
        return CycleEvaluation(start, end, np.array([self.log_rate]),
                               {"ledger": 1.0}, np.array([0.0, 1.0]),
                               np.array([[self.log_rate]]),
                               {}, "dormant", self.protected_signatures().topology)
    def commit_private_cycle(self, ev): self.x = ev.state_end.vector.copy(); self.commit_ledger_increments(ev.ledger_increments)
    def commit_ledger_increments(self, increments): self.ledger += float(increments.get("ledger", 0.0))
    def remaining_birth_actions(self): return self.threshold - self.action
    def commit_birth_action(self, increment, cycles): self.action += np.asarray(increment)
    def commit_log_birth_action(self, increment, cycles):
        self.log_action = np.logaddexp(self.log_action, np.asarray(increment, float))
        self.action = np.where(self.log_action >= math.log(np.nextafter(0.0, 1.0)), np.exp(self.log_action), 0.0)
    def physical_cycles(self): return self.cycles
    def set_physical_cycles(self, cycles): self.cycles = float(cycles)


def cfg(**kwargs):
    base = dict(periodic_relative_tolerance=1e-12, periodic_admission_distance=1e-12,
                periodic_max_iterations=8, projective_initial_cycles=16,
                exact_retry_cycles=2, event_guard_cycles=2.0)
    base.update(kwargs)
    return HighCycleConfig(**base)


def test_private_cycle_preserves_clocks_rng_ledgers_topology_and_time():
    model = SyntheticDormantPD(rate=0.2, contraction=0.5)
    before = model.protected_signatures(); state = model.active_state(); cycles = model.cycles
    result = private_cycle(model)
    assert result.birth_action[0] == 0.2
    assert model.protected_signatures() == before
    assert np.array_equal(model.active_state().vector, state.vector)
    assert model.cycles == cycles


def test_stationary_synthetic_1e12_and_1e14_right_censors():
    for horizon in (1e12, 1e14):
        model = SyntheticDormantPD(rate=1e-20, threshold=10.0, contraction=0.0)
        result = DormantPDHighCycleEngine(model, cfg()).advance(horizon)
        assert result.cycles_consumed == horizon
        assert result.event_guard_reached is False
        assert abs(model.action[0] - horizon * 1e-20) <= 1e-15
        assert model.ledger == horizon
        assert any(row.mode == "stationary" for row in result.modes)


def test_large_cycle_first_passage_returns_exact_event_guard():
    target = 10**12 + 0.25
    rate = 1e-12
    model = SyntheticDormantPD(rate=rate, threshold=target * rate, contraction=0.0)
    result = DormantPDHighCycleEngine(model, cfg(event_guard_cycles=2.0)).advance(1e14)
    assert result.event_guard_reached
    assert model.cycles == math.floor(target - 2.0)
    assert model.action[0] < model.threshold[0]


def test_partition_equivalence_stationary():
    one = SyntheticDormantPD(rate=2e-15, threshold=10.0)
    split = deepcopy(one)
    DormantPDHighCycleEngine(one, cfg()).advance(1e12)
    for part in (1e8, 9.999e11):
        DormantPDHighCycleEngine(split, cfg()).advance(part)
    assert one.cycles == split.cycles
    assert np.allclose(one.action, split.action, rtol=2e-15)
    assert np.array_equal(one.x, split.x)


def test_projective_linear_drift_is_accepted_without_false_stationarity():
    model = SyntheticDormantPD(rate=1e-30, threshold=10.0, contraction=1.0, drift=1e-6)
    engine = DormantPDHighCycleEngine(model, cfg(periodic_max_iterations=2,
        projective_state_tolerance=1e-12, projective_log_hazard_tolerance=1e-12))
    result = engine.advance(128)
    assert result.cycles_consumed == 128
    assert result.accepted_projected_cycles > 0
    assert abs(model.x[0] - 128e-6) < 1e-14
    assert any(row.mode == "projective" for row in result.modes)


def test_cache_invalidates_on_discrete_change():
    model = SyntheticDormantPD()
    engine = DormantPDHighCycleEngine(model, cfg())
    model.status[0] = 1
    result = engine.advance(100)
    assert result.cache_invalidated
    assert result.cycles_consumed == 0.0


def test_restart_equivalence_periodic_and_projective(tmp_path):
    for kwargs in ({"contraction": 0.0}, {"contraction": 1.0, "drift": 1e-6}):
        continuous = SyntheticDormantPD(rate=1e-20, threshold=10.0, **kwargs)
        restarted = deepcopy(continuous)
        config = cfg(periodic_max_iterations=2, projective_state_tolerance=1e-12,
                     projective_log_hazard_tolerance=1e-12)
        DormantPDHighCycleEngine(continuous, config).advance(128)
        first = DormantPDHighCycleEngine(restarted, config)
        first.advance(64)
        first.write_atomic_mode_checkpoint(tmp_path / "mode.json")
        DormantPDHighCycleEngine(restarted, config).advance(64)
        assert continuous.cycles == restarted.cycles
        assert np.allclose(continuous.action, restarted.action, rtol=2e-15)
        assert np.allclose(continuous.x, restarted.x, rtol=2e-15, atol=1e-18)
        assert (tmp_path / "mode.json").is_file()


def test_log_action_below_float_range_survives_large_formal_skip():
    model = SyntheticDormantPD(rate=1.0, threshold=1.0)
    model.rate = 0.0
    model.log_rate = -1000.0  # deliberately below linear representability
    result = DormantPDHighCycleEngine(model, cfg()).advance(1e300)
    assert result.cycles_consumed == 1e300
    assert abs(model.log_action[0] - (-1000.0 + math.log(1e300))) < 1e-12
    assert model.ledger == 1e300
