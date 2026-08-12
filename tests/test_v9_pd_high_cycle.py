from __future__ import annotations

from copy import deepcopy
import math
from types import SimpleNamespace

import numpy as np

from arrhenius_fracture.v9_pd_high_cycle import (
    ActiveState, CycleEvaluation, DormantPDHighCycleEngine, HighCycleConfig,
    ProtectedSignatures, _digest, private_cycle,
)
from arrhenius_fracture.v9_pd_high_cycle_adapter import SpatialPDDormantAdapter, SharedRootSpatialPDDormantAdapter
from arrhenius_fracture.v9_pd_shared_root_marked_cleavage import SharedRootMarkedCleavageState


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


def test_shared_root_adapter_carries_mpz_and_protects_threshold_and_rng():
    class MPZ:
        _CAPSULE_ARRAYS=("mobile",)
        audit={}
        def __init__(self):self.mobile=np.array([2.]);self.time=0.
        def capsule(self):return {"arrays":{"mobile":self.mobile.copy(),"diagnostic":np.array([0.])},"scalars":{"time_s":self.time,"signed_last_source_activations":0.},"schema":"test"}
        def restore_capsule(self,c):self.mobile=np.asarray(c["arrays"]["mobile"]).copy();self.time=float(c["scalars"]["time_s"])
        def copy(self):x=MPZ();x.mobile=self.mobile.copy();x.time=self.time;return x
        def summary(self):return {"signed_active_K_shield_Pa_sqrt_m":0.,"tip_radius_m":1e-4}
    clock=SharedRootMarkedCleavageState("x",".",shear_modulus_Pa=1.,poisson=.3,burgers_m=1e-10,initial_tip_radius_m=1e-4,hazard_seed=2,mark_seed=3,mpz=MPZ(),m_hits=1.)
    n=2;s=SimpleNamespace(log_delivery_memory=np.zeros(n),delivery_memory=np.ones(n),available=np.ones(n),embryo=np.zeros(n),stable=np.zeros(n),inactive=np.zeros(n),completion=np.zeros(n),born_cumulative=np.zeros(n),healed_cumulative=np.zeros(n),born_sites_cumulative=np.zeros(n),healed_sites_cumulative=np.zeros(n),log_birth_cumulative_hazard=np.full(n,-math.inf),site_birth_threshold=np.ones(n),birth_cumulative_hazard=np.zeros(n),site_transition_threshold=np.ones(n),site_transition_cumulative_hazard=np.zeros(n),site_transition_outcome_uniform=np.zeros(n),site_status=np.zeros(n,dtype=np.uint8),bond_damage=np.zeros(1),primary_seed_node=-1,active_front=False,active_front_bonds=np.zeros(1,dtype=bool),front_backbone_bonds=np.zeros(1,dtype=bool),front_wake_bonds=np.zeros(1,dtype=bool),front_process_bonds=np.zeros(1,dtype=bool),active_front_path_xy=np.empty((0,2)))
    patch=SimpleNamespace(_candidate_rng=np.random.default_rng(1),_event_rng=np.random.default_rng(2),bonds=np.array([[0,1]]))
    mesh=SimpleNamespace(nodes=np.zeros((2,2)))
    adapter=SharedRootSpatialPDDormantAdapter(shared_clock=clock,patch=patch,pd_state=s,mesh=mesh,ep_gp=np.zeros((3,1)),rho_gp=np.full(1,1e12),epsp_acc_gp=np.zeros(1),u=np.zeros(4),cycles=0.,plastic_work=0.,cycle_evaluator=lambda _: {})
    snap=adapter.active_state();threshold=adapter.shared_clock.global_threshold_action;hrng=deepcopy(adapter.shared_clock._hazard_rng.bit_generator.state);mrng=deepcopy(adapter.shared_clock._mark_rng.bit_generator.state)
    adapter.restore_active_state(snap,snap.vector);adapter.commit_log_birth_action(np.array([math.log(.1)]),10.)
    assert math.isclose(adapter.shared_clock.global_cumulative_action,.1,rel_tol=2e-16)
    assert adapter.shared_clock.global_threshold_action==threshold
    assert adapter.shared_clock._hazard_rng.bit_generator.state==hrng
    assert adapter.shared_clock._mark_rng.bit_generator.state==mrng


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


def test_exact_private_window_training_partition_guard_and_restart():
    class Windowed(SyntheticDormantPD):
        def exact_private_window(self, dN):
            n = float(dN)
            start = self.active_state()
            # Include an exactly integrable cycle-coordinate-dependent term so
            # the split right window must begin at N + dN_left.  Evaluating it
            # at the original N would fail the partition guard.
            n0 = self.physical_cycles()
            coordinate_increment = 1e-14 * ((n0 + n) ** 2 - n0 ** 2)
            end = ActiveState(
                start.vector + n * self.drift + coordinate_increment,
                start.specification,
            )
            return CycleEvaluation(
                start, end, np.array([self.log_rate + math.log(n)]),
                {"ledger": n}, np.array([0.0, 1.0]),
                np.array([[self.log_rate], [self.log_rate]]), {}, "dormant",
                self.protected_signatures().topology,
            )

    config = cfg(
        periodic_admission_distance=0.0, private_window_initial_cycles=96,
        private_window_state_tolerance=1e-13,
        private_window_log_hazard_tolerance=1e-13,
        minimum_projected_cycles_per_exact_map=16.0,
    )
    continuous = Windowed(rate=1e-8, threshold=10.0, contraction=1.0, drift=2e-6)
    restarted = deepcopy(continuous)
    result = DormantPDHighCycleEngine(continuous, config).advance(192)
    assert result.cycles_consumed == 192
    assert result.accepted_projected_cycles == 192
    assert all(row.detail["projected_cycles_per_exact_map"] >= 16.0
               for row in result.modes if row.mode == "exact_private_window")
    DormantPDHighCycleEngine(restarted, config).advance(96)
    DormantPDHighCycleEngine(restarted, config).advance(96)
    np.testing.assert_allclose(continuous.x, restarted.x, rtol=0.0, atol=1e-18)
    np.testing.assert_allclose(continuous.action, restarted.action, rtol=2e-15)
    assert continuous.ledger == restarted.ledger == 192.0

    guarded = Windowed(rate=1e-2, threshold=0.5, contraction=1.0, drift=0.0)
    guarded_result = DormantPDHighCycleEngine(guarded, config).advance(96)
    assert not any(row.mode == "exact_private_window" for row in guarded_result.modes)
    assert guarded.action[0] < guarded.threshold[0]


def test_rejected_private_window_refines_without_mutating_physical_state():
    class RejectLargeWindows(SyntheticDormantPD):
        def exact_private_window(self, dN):
            if dN > 4:
                raise RuntimeError("synthetic physical transaction rejection")
            n = float(dN)
            start = self.active_state()
            end = ActiveState(start.vector + n * self.drift, start.specification)
            return CycleEvaluation(
                start, end, np.array([self.log_rate + math.log(n)]),
                {"ledger": n}, np.array([0.0, 1.0]),
                np.array([[self.log_rate], [self.log_rate]]), {}, "dormant",
                self.protected_signatures().topology,
            )

    model = RejectLargeWindows(rate=1e-8, threshold=10.0,
                               contraction=1.0, drift=2e-6)
    result = DormantPDHighCycleEngine(model, cfg(
        periodic_admission_distance=0.0,
        private_window_initial_cycles=16,
        minimum_projected_cycles_per_exact_map=0.0,
    )).advance(16)
    assert result.cycles_consumed == 16
    assert any(row.mode == "exact_private_window_reject"
               and row.detail.get("reason") ==
               "private_window_physical_transaction_reject"
               for row in result.modes)
    assert model.cycles == 16
    assert model.ledger == 16
    np.testing.assert_allclose(model.x, [32e-6], rtol=0.0, atol=1e-18)


def test_rejected_one_cycle_admission_returns_to_direct_path_unchanged():
    class RejectPrivateCycle(SyntheticDormantPD):
        def exact_private_cycle(self):
            raise RuntimeError("embedded FEM step needs subdivision")

    model = RejectPrivateCycle(rate=1e-8, threshold=10.0,
                               contraction=1.0, drift=2e-6)
    before = model.protected_signatures()
    result = DormantPDHighCycleEngine(model, cfg(
        private_window_training=False,
    )).advance(100)
    assert result.cycles_consumed == 0.0
    assert result.accepted_projected_cycles == 0.0
    assert model.cycles == 0.0
    assert model.protected_signatures() == before
    assert result.modes[-1].mode == "exact_private_cycle_reject"
    assert result.modes[-1].detail["reason"] == \
        "private_cycle_physical_transaction_reject"


def test_rejected_projective_probe_returns_to_direct_path_unchanged():
    class RejectProjectedProbe(SyntheticDormantPD):
        calls = 0

        def exact_private_cycle(self):
            self.calls += 1
            if self.calls >= 3:
                raise RuntimeError("projected FEM state needs subdivision")
            return super().exact_private_cycle()

    model = RejectProjectedProbe(rate=1e-8, threshold=10.0,
                                 contraction=1.0, drift=2e-6)
    before = model.protected_signatures()
    result = DormantPDHighCycleEngine(model, cfg(
        private_window_training=False,
        periodic_admission_distance=0.0,
        exact_retry_cycles=0,
    )).advance(100)
    assert result.cycles_consumed == 0.0
    assert result.accepted_projected_cycles == 0.0
    assert model.cycles == 0.0
    assert model.protected_signatures() == before
    assert any(row.detail.get("reason") ==
               "private_cycle_physical_transaction_reject"
               for row in result.modes)


def test_private_window_event_guard_converts_per_second_rate_at_nondefault_frequency():
    class PerSecondWindow(SyntheticDormantPD):
        def __init__(self, frequency_hz):
            super().__init__(rate=1e-30, threshold=0.015, contraction=1.0)
            self.frequency_hz = float(frequency_hz)
        def exact_private_window(self, dN):
            n=float(dN); start=self.active_state(); end=ActiveState(start.vector,start.specification)
            action=n/self.frequency_hz
            return CycleEvaluation(
                start,end,np.array([math.log(action)]),{"ledger":n},
                np.array([0.,1.]),np.array([[0.],[0.]]),{},"dormant",
                self.protected_signatures().topology,
                phase_rate_seconds_per_cycle=1.0/self.frequency_hz,
            )

    config=cfg(periodic_admission_distance=0.0,private_window_initial_cycles=10,
               private_window_max_cycles=10,minimum_projected_cycles_per_exact_map=0.0)
    # Ten cycles carry action 0.01 at 1 kHz.  A two-cycle guard is 0.002,
    # therefore this is safe.  Treating the phase rate as per-cycle would
    # incorrectly reject it by orders of magnitude.
    model=PerSecondWindow(1000.0)
    result=DormantPDHighCycleEngine(model,config).advance(10)
    assert result.accepted_projected_cycles == 10
    row=next(row for row in result.modes if row.mode=="exact_private_window")
    assert row.detail["phase_rate_seconds_per_cycle"] == 1e-3
    assert math.isclose(math.exp(row.detail["guard_log_action_upper_bound"]),0.002)


def test_private_window_event_guard_frequency_conversion_is_not_1000hz_cancellation():
    for frequency_hz in (20.0, 2500.0):
        class Window(SyntheticDormantPD):
            def exact_private_window(self,dN):
                n=float(dN);start=self.active_state();end=ActiveState(start.vector,start.specification)
                return CycleEvaluation(start,end,np.array([math.log(n/frequency_hz)]),
                    {"ledger":n},np.array([0.,1.]),np.array([[0.],[0.]]),{},"dormant",
                    self.protected_signatures().topology,
                    phase_rate_seconds_per_cycle=1.0/frequency_hz)
        model=Window(rate=1e-30,threshold=100.,contraction=1.)
        config=cfg(periodic_admission_distance=0.0,private_window_initial_cycles=10,
                   private_window_max_cycles=10,minimum_projected_cycles_per_exact_map=0.0)
        result=DormantPDHighCycleEngine(model,config).advance(10)
        row=next(row for row in result.modes if row.mode=="exact_private_window")
        assert math.isclose(math.exp(row.detail["guard_log_action_upper_bound"]),2.0/frequency_hz)


def test_log_action_below_float_range_survives_large_formal_skip():
    model = SyntheticDormantPD(rate=1.0, threshold=1.0)
    model.rate = 0.0
    model.log_rate = -1000.0  # deliberately below linear representability
    result = DormantPDHighCycleEngine(model, cfg()).advance(1e300)
    assert result.cycles_consumed == 1e300
    assert abs(model.log_action[0] - (-1000.0 + math.log(1e300))) < 1e-12
    assert model.ledger == 1e300


def test_large_cycle_sub_ulp_event_boundary_is_next_representable_float():
    cycles=1.8058541488375247e6
    residual_wait=5.56e-11
    assert cycles+residual_wait==cycles
    step=np.nextafter(cycles,math.inf)-cycles
    assert step>residual_wait and cycles+step>cycles


def test_spatially_nonuniform_but_temporally_constant_ledger_does_not_reject_projection():
    class SpatialLedger(SyntheticDormantPD):
        def __init__(self):
            super().__init__(rate=1e-30,threshold=10.0,contraction=1.0,drift=1e-6)
            self.spatial_ledger=np.zeros(2)
        def exact_private_cycle(self):
            ev=super().exact_private_cycle(); ev.ledger_increments={"spatial":np.array([1.0,10.0])}; return ev
        def commit_ledger_increments(self,increments):
            if "spatial" in increments: self.spatial_ledger+=np.asarray(increments["spatial"])
        def protected_signatures(self):
            base=super().protected_signatures(); return ProtectedSignatures(_digest(self.spatial_ledger),base.stochastic,base.topology)
    model=SpatialLedger(); engine=DormantPDHighCycleEngine(model,cfg(periodic_max_iterations=2,
        projective_state_tolerance=1e-12,projective_log_hazard_tolerance=1e-12))
    result=engine.advance(16)
    assert result.accepted_projected_cycles==16
    np.testing.assert_allclose(model.spatial_ledger,[16.0,160.0])


def test_affine_population_projection_handles_zero_to_positive_source_without_overflow():
    spec = (("embryo", (3,), "float64"),)
    start = ActiveState(np.array([0.0, 0.2, 0.7]), spec)
    # Exact affine recurrences: x' = 0.5*x + [0.1, 0.0, 0.0].
    first = ActiveState(np.array([0.1, 0.1, 0.35]), spec)
    second = ActiveState(np.array([0.15, 0.05, 0.175]), spec)
    adapter = SpatialPDDormantAdapter.__new__(SpatialPDDormantAdapter)
    projected = adapter.project_active_state(start, first, second, 16)
    expected = np.array([0.2 * (1.0 - 0.5**16), 0.2 * 0.5**16, 0.7 * 0.5**16])
    np.testing.assert_allclose(projected, expected, rtol=1e-13, atol=1e-15)
    assert np.all((projected >= 0.0) & (projected <= 1.0))


def test_population_ledgers_close_from_validated_available_and_inactive_endpoints():
    names = ("available", "embryo", "stable", "inactive", "completion")
    spec = tuple((name, (1,), "float64") for name in names)
    start = ActiveState(np.array([1.0, 0.0, 0.0, 0.0, 0.0]), spec)
    end = np.array([0.91, 0.02, 0.03, 0.04, 0.0])
    adapter = SpatialPDDormantAdapter.__new__(SpatialPDDormantAdapter)
    adapter.patch = SimpleNamespace(cfg=SimpleNamespace(heal_return_fraction=0.2))
    ledgers, constrained = adapter.conservative_population_ledgers(start, end, {})
    np.testing.assert_allclose(ledgers["healed_cumulative"], [0.05])
    np.testing.assert_allclose(ledgers["born_cumulative"], [0.10])
    assert constrained == {"healed_cumulative", "born_cumulative"}


def test_repeated_projective_halving_obeys_hard_exact_map_budget():
    model = SyntheticDormantPD(rate=1e-20, threshold=10.0, contraction=0.5, drift=1e-3)
    engine = DormantPDHighCycleEngine(model, cfg(
        periodic_admission_distance=0.0,
        projective_initial_cycles=64,
        projective_state_tolerance=0.0,
        projective_log_hazard_tolerance=0.0,
        projective_curvature_tolerance=0.0,
        exact_retry_cycles=0,
        max_exact_map_evaluations=5,
    ))
    result = engine.advance(1000)
    # One periodic precheck plus exactly one four-map rejected proposal. The
    # next halved proposal is refused before any private map can run.
    assert result.exact_map_evaluations == 5
    assert sum(m.mode == "projective_reject" for m in result.modes) == 1
    assert result.modes[-1].mode == "efficiency_budget"
    assert result.modes[-1].detail["required_next_maps"] == 4


def test_projective_efficiency_floor_rejects_locally_valid_tiny_skip():
    model = SyntheticDormantPD(rate=1e-20, threshold=10.0, contraction=1.0, drift=1e-6)
    engine = DormantPDHighCycleEngine(model, cfg(
        periodic_admission_distance=0.0,
        projective_initial_cycles=16,
        exact_retry_cycles=0,
        max_exact_map_evaluations=5,
        minimum_projected_cycles_per_exact_map=8.0,
    ))
    result = engine.advance(1000)
    reject = next(m for m in result.modes if m.mode == "projective_reject")
    assert reject.detail["reason"] == "insufficient_projective_efficiency"
    assert reject.detail["projected_cycles_per_exact_map"] == 4.0
    assert result.accepted_projected_cycles == 0.0
