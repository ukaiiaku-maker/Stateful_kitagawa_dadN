import copy
import pickle
import unittest
from pathlib import Path

import numpy as np

from arrhenius_fracture.config import ElasticProperties
from arrhenius_fracture.sn_feature_geometry_v8_7 import BluntNotchGeometry, make_blunt_edge_notch_mesh
from arrhenius_fracture.sn_intact_fem import plane_strain_D, solve_symmetric_tension, assemble_intact_mechanics
from arrhenius_fracture.v9_energy_gated_stable_birth import EnergyGatedStableBirthState
from arrhenius_fracture.v9_four_class_registry import EXPECTED
from arrhenius_fracture.v9_stable_birth_energy_gate import (
    StableBirthEnergyGateConfig, evaluate_stable_birth_energy_gate,
    threshold_scaled_event_length,
)


SOURCE = Path("/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1")


class StableBirthEnergyGateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        geom = BluntNotchGeometry(2e-3, 4e-3, 1.5e-4, 3e-4, feature_type="ellipse")
        cls.mesh, cls.bnd, root = make_blunt_edge_notch_mesh(
            geom, nx=10, ny=18, jitter=0.0, root_h_fine=2.0e-5, seed=1)
        cls.root = np.asarray(root, float)
        cls.mat = ElasticProperties(); cls.D = plane_strain_D(cls.mat)
        ep = np.zeros((3, cls.mesh.ne)); u0 = np.zeros(cls.mesh.ndof)
        K, R, *_ = assemble_intact_mechanics(cls.mesh, u0, ep, cls.D, cls.mat)
        cls.u, _ = solve_symmetric_tension(K, R, u0, cls.bnd, cls.mesh, 2e-6, -2e-6)
        cls.ep = ep

    def test_threshold_scaled_reward_is_same_mean_preserving_law(self):
        cfg = StableBirthEnergyGateConfig()
        low, f_low = threshold_scaled_event_length(0.01, cfg)
        high, f_high = threshold_scaled_event_length(10.0, cfg)
        self.assertEqual(low, cfg.base_checkpoint_m * f_low)
        self.assertEqual(high, cfg.base_checkpoint_m * f_high)
        self.assertLess(f_low, f_high)

    def gate(self, barrier):
        return evaluate_stable_birth_energy_gate(
            mesh=self.mesh, boundaries=self.bnd, displacement=self.u,
            ep_gp=self.ep, Dmat=self.D, root_xy=self.root,
            event_direction=np.array([1.0, 0.0]), event_K_Pa_sqrt_m=20e6,
            cleavage_barrier_J=barrier, cooperative_hits=3.0,
            burgers_m=self.mat.b, threshold_action=1.0,
            plane_strain_modulus_Pa=self.mat.Eprime,
            config=StableBirthEnergyGateConfig(base_checkpoint_m=2.0e-4))

    def test_gate_has_rejected_and_admitted_branches_without_mutation(self):
        before = self.u.copy()
        rejected = self.gate(1e-12)
        admitted = self.gate(0.0)
        self.assertEqual(rejected["committed_event_length_m"], 0.0)
        self.assertGreater(admitted["committed_event_length_m"], 0.0)
        np.testing.assert_array_equal(self.u, before)
        self.assertFalse(rejected["hazard_gated_before_first_passage"])

    def test_phase_cache_preserves_exact_energy_gate_result(self):
        kwargs = dict(
            mesh=self.mesh, boundaries=self.bnd, displacement=self.u,
            ep_gp=self.ep, Dmat=self.D, root_xy=self.root,
            event_direction=np.array([1.0, 0.0]), event_K_Pa_sqrt_m=20e6,
            cleavage_barrier_J=0.0, cooperative_hits=3.0,
            burgers_m=self.mat.b, threshold_action=1.0,
            plane_strain_modulus_Pa=self.mat.Eprime,
            config=StableBirthEnergyGateConfig(base_checkpoint_m=2.0e-4),
        )
        uncached = evaluate_stable_birth_energy_gate(**kwargs)
        cache = {}
        first = evaluate_stable_birth_energy_gate(**kwargs, evaluation_cache=cache)
        second = evaluate_stable_birth_energy_gate(**kwargs, evaluation_cache=cache)
        self.assertEqual(pickle.dumps(first), pickle.dumps(uncached))
        self.assertEqual(pickle.dumps(second), pickle.dumps(uncached))
        self.assertTrue(cache)

    def test_rejected_attempt_consumes_xi_and_restart_preserves_next_draw(self):
        option = next(iter(EXPECTED))
        state = EnergyGatedStableBirthState(
            option, SOURCE, shear_modulus_Pa=160.15625e9, poisson=.28,
            burgers_m=2.74e-10, initial_tip_radius_m=1e-6, hazard_seed=7123)
        first = state.hazard_threshold_action
        state.pending_attempt = {"attempt_index": 0, "threshold_action": first,
                                 "cumulative_cleavage_hazard": first, "time_s": 1.0}
        state.cumulative_cleavage_hazard = first; state.time_s = 1.0
        state.resolve_pending_attempt({"committed_event_length_m": 0.0,
                                       "arrest_reason": "no_mesh_resolved_admissible_increment"})
        self.assertEqual(state.nonpropagating_attempt_count, 1)
        self.assertNotEqual(state.hazard_threshold_action, first)
        capsule = state.capsule(); restored = EnergyGatedStableBirthState(
            option, SOURCE, shear_modulus_Pa=160.15625e9, poisson=.28,
            burgers_m=2.74e-10, initial_tip_radius_m=1e-6, hazard_seed=7123)
        restored.restore_capsule(copy.deepcopy(capsule))
        self.assertEqual(restored.capsule()["hazard_rng_state"], capsule["hazard_rng_state"])
        self.assertEqual(restored.hazard_threshold_action, state.hazard_threshold_action)
        self.assertIsNone(restored.diagnostics()["stable_birth_survival"])

    def test_eleven_way_partition_and_atomic_restart_preserve_renewal_state(self):
        option = next(iter(EXPECTED))
        def make():
            return EnergyGatedStableBirthState(
                option, SOURCE, shear_modulus_Pa=160.15625e9, poisson=.28,
                burgers_m=2.74e-10, initial_tip_radius_m=1e-6, hazard_seed=8119)
        tensors = np.array([[[2e9, .3e9], [.3e9, 3e9]],
                            [[1e9, -.2e9], [-.2e9, 1.5e9]]])
        whole, split = make(), make()
        whole.hazard_threshold_action = split.hazard_threshold_action = 1e100
        whole.advance_fem_phase_block(0.11, 1000.0, 300.0, tensors)
        for _ in range(11):
            split.advance_fem_phase_block(0.01, 1000.0, 300.0, tensors)
        self.assertAlmostEqual(whole.cumulative_cleavage_hazard,
                               split.cumulative_cleavage_hazard,
                               delta=max(split.cumulative_cleavage_hazard*4e-4, 1e-300))
        capsule = split.capsule(); restored = make(); restored.restore_capsule(capsule)
        self.assertEqual(restored.capsule()["hazard_rng_state"], capsule["hazard_rng_state"])
        self.assertEqual(restored.nonpropagating_attempt_count, split.nonpropagating_attempt_count)

    def test_admitted_birth_stops_without_post_birth_rng_draw(self):
        option = next(iter(EXPECTED))
        state = EnergyGatedStableBirthState(
            option, SOURCE, shear_modulus_Pa=160.15625e9, poisson=.28,
            burgers_m=2.74e-10, initial_tip_radius_m=1e-6, hazard_seed=991)
        rng_before = copy.deepcopy(state._hazard_rng.bit_generator.state)
        state.pending_attempt = {"attempt_index": 0, "threshold_action": state.hazard_threshold_action,
                                 "cumulative_cleavage_hazard": state.hazard_threshold_action,
                                 "time_s": 2.0}
        state.time_s = 2.0
        state.resolve_pending_attempt({"committed_event_length_m": 1e-6,
                                       "arrest_reason": "hazard_derived_energy_arrest"})
        self.assertTrue(state.fired)
        self.assertEqual(rng_before, state._hazard_rng.bit_generator.state)
        self.assertEqual(state.nonpropagating_attempt_count, 0)


if __name__ == "__main__":
    unittest.main()
