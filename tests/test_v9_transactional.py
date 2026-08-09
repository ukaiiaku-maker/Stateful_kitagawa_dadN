import json
from pathlib import Path
import tempfile
import unittest
from dataclasses import replace

from arrhenius_fracture.v9_transactional import (
    AnalyticHazardModel,
    AtomicGenerationStore,
    EngineConfig,
    PersistentClock,
    TransactionalEngine,
    empty_capsule,
)
from arrhenius_fracture.v9_solver import V9SolverRequest, build_physical_capsule


class V9TransactionalTests(unittest.TestCase):
    def run_model(self, initial_block, growth, end=1000.0):
        clock = PersistentClock("site-0", "birth", threshold=1.0e30)
        initial = empty_capsule(clocks=(clock,))
        model = AnalyticHazardModel(((2.0e-4, 3.0e-7),), state_rate=0.25)
        engine = TransactionalEngine(model, EngineConfig(initial_block_cycles=initial_block, quiet_growth_factor=growth))
        return engine.advance(initial, end)

    def test_block_partition_invariance_over_more_than_tenfold_schedules(self):
        a = self.run_model(0.1, 2.0)
        b = self.run_model(10.0, 10.0)
        self.assertAlmostEqual(a.cycle, b.cycle, places=12)
        self.assertAlmostEqual(a.aggregate_cumulative_hazard, b.aggregate_cumulative_hazard, places=12)
        self.assertAlmostEqual(a.clocks[0].cumulative_hazard, b.clocks[0].cumulative_hazard, places=12)
        self.assertAlmostEqual(json.loads(a.fem_plastic_state_json)["x"], json.loads(b.fem_plastic_state_json)["x"], places=12)

    def test_exact_constant_hazard_threshold_crossing(self):
        initial = empty_capsule(clocks=(PersistentClock("s", "birth", 5.0),))
        engine = TransactionalEngine(AnalyticHazardModel(((2.0, 0.0),)), EngineConfig(initial_block_cycles=10.0))
        final = engine.advance(initial, 10.0)
        self.assertEqual(final.event_sequence[0][0], "birth")
        self.assertAlmostEqual(final.event_sequence[0][2], 2.5, places=13)
        self.assertEqual(final.clocks[0].cumulative_hazard, 5.0)

    def test_exact_time_varying_hazard_threshold_crossing(self):
        initial = empty_capsule(clocks=(PersistentClock("s", "birth", 9.0),))
        engine = TransactionalEngine(AnalyticHazardModel(((0.0, 2.0),)), EngineConfig(initial_block_cycles=10.0))
        final = engine.advance(initial, 10.0)
        self.assertAlmostEqual(final.event_sequence[0][2], 3.0, places=13)

    def test_earliest_event_across_birth_embryo_and_progression_clocks(self):
        clocks = (
            PersistentClock("birth-0", "birth", 5.0),
            PersistentClock("embryo-0", "embryo_transition", 2.0),
            PersistentClock("front-0", "progression", 8.0),
        )
        initial = empty_capsule(clocks=clocks)
        model = AnalyticHazardModel(((1.0, 0.0), (1.0, 0.0), (1.0, 0.0)))
        final = TransactionalEngine(model, EngineConfig(initial_block_cycles=10.0)).advance(initial, 10.0)
        self.assertEqual(final.event_sequence[0][:2], ("embryo_transition", "embryo-0"))
        self.assertAlmostEqual(final.event_sequence[0][2], 2.0, places=13)
        self.assertEqual({event[0] for event in final.event_sequence}, {"birth", "embryo_transition", "progression"})

    def test_error_control_rejects_and_subdivides_before_commit(self):
        base = AnalyticHazardModel(((1.0e-3, 0.0),), state_rate=1.0)

        class ErrorWrappedModel:
            def propose(self, state, dN):
                return replace(base.propose(state, dN), error_estimate=dN * dN)

            def crossing_fraction(self, state, clock_index, target_increment, dN):
                return base.crossing_fraction(state, clock_index, target_increment, dN)

        initial = empty_capsule(clocks=(PersistentClock("s", "birth", 1.0e9),))
        config = EngineConfig(initial_block_cycles=8.0, quiet_growth_factor=2.0, error_tolerance=1.0)
        final = TransactionalEngine(ErrorWrappedModel(), config).advance(initial, 5.0)
        self.assertAlmostEqual(final.cycle, 5.0, places=13)
        self.assertAlmostEqual(final.committed_cycle_measure, 5.0, places=13)
        self.assertGreaterEqual(final.commit_count, 5)

    def test_restart_equivalence_at_accepted_boundary(self):
        initial = empty_capsule(clocks=(PersistentClock("s", "birth", 80.0),), request_hash="restart-test")
        model = AnalyticHazardModel(((0.5, 0.01),), state_rate=0.125)
        config = EngineConfig(initial_block_cycles=1.0, quiet_growth_factor=10.0)
        uninterrupted = TransactionalEngine(model, config).advance(initial, 111.0)
        partial = TransactionalEngine(model, config).advance(initial, 11.0)
        with tempfile.TemporaryDirectory() as directory:
            store = AtomicGenerationStore(Path(directory))
            store.write_generation(partial, {"status": "checkpoint"})
            restored = store.load_latest()
            resumed = TransactionalEngine(model, config).advance(restored, 111.0)
        self.assertEqual(resumed.event_sequence, uninterrupted.event_sequence)
        self.assertEqual(resumed.clocks, uninterrupted.clocks)
        self.assertEqual(resumed.fem_plastic_state_json, uninterrupted.fem_plastic_state_json)
        self.assertAlmostEqual(resumed.aggregate_cumulative_hazard, uninterrupted.aggregate_cumulative_hazard, places=13)
        self.assertEqual(resumed.rng_states_json, uninterrupted.rng_states_json)

    def test_no_double_commit(self):
        final = self.run_model(1.0, 3.0, end=1234.5)
        self.assertAlmostEqual(final.committed_cycle_measure, final.cycle, places=12)
        self.assertGreater(final.commit_count, 0)

    def test_cap_activation_is_reported(self):
        initial = empty_capsule(clocks=(PersistentClock("s", "birth", 1.0e20),))
        model = AnalyticHazardModel(((1.0e-6, 0.0),), state_rate=1.0, state_cap=5.0)
        final = TransactionalEngine(model, EngineConfig(initial_block_cycles=10.0)).advance(initial, 20.0)
        self.assertIn("synthetic_state_cap", final.cap_floor_diagnostics)
        self.assertEqual(json.loads(final.fem_plastic_state_json)["x"], 5.0)
        uncapped = TransactionalEngine(
            AnalyticHazardModel(((1.0e-6, 0.0),), state_rate=1.0),
            EngineConfig(initial_block_cycles=10.0),
        ).advance(initial, 20.0)
        self.assertAlmostEqual(final.aggregate_cumulative_hazard, uncapped.aggregate_cumulative_hazard, places=15)

    def test_capsule_contains_required_identity_and_state_sections(self):
        state = empty_capsule(clocks=(PersistentClock("s", "birth", 1.0),))
        payload = state.to_dict()
        for key in (
            "fem_plastic_state_json", "rho_backstress_shielding_json", "delivery_memory_json",
            "completion_json", "candidate_sites_json", "clocks", "embryo_stable_transition_json",
            "pd_damage_front_json", "rng_states_json", "rng_digests_json",
            "geometry_identity_state_json", "source_hashes_json", "request_hash", "solver_version",
        ):
            self.assertIn(key, payload)

    def test_physical_capsule_refuses_missing_authoritative_state(self):
        request = V9SolverRequest(source_hashes=(("v8.7", "abc"),))
        with self.assertRaisesRegex(ValueError, "rng_states"):
            build_physical_capsule(
                request,
                cycle=1.0,
                fem_plastic_state={}, rho_backstress_shielding={}, delivery_memory={},
                completion={}, candidate_sites={}, clocks=(), embryo_stable_transition={},
                pd_damage_front={}, rng_states=None, rng_digests={}, geometry_identity_state={},
            )


if __name__ == "__main__":
    unittest.main()
