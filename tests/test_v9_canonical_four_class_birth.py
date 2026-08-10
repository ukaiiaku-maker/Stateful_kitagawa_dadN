import copy
import unittest
from pathlib import Path

import numpy as np
from scipy.special import gammainc

from arrhenius_fracture.v9_canonical_four_class_birth import (
    ENDPOINT, CanonicalFourClassBirthState, canonical_effective_cleavage_log_rate,
    canonical_effective_cleavage_rate,
)
from arrhenius_fracture.v9_four_class_registry import EXPECTED


SOURCE = Path("/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1")


class CanonicalFourClassBirthTests(unittest.TestCase):
    def make(self, option=None, seed=1720):
        return CanonicalFourClassBirthState(
            option or next(iter(EXPECTED)), SOURCE,
            shear_modulus_Pa=160.15625e9, poisson=0.28,
            burgers_m=2.74e-10, initial_tip_radius_m=1.0e-6,
            hazard_seed=seed,
        )

    def test_exact_audited_multihit_rate_all_four_rows(self):
        for option in EXPECTED:
            state = self.make(option)
            for temperature in (300.0, 900.0, 1300.0):
                rates = state.cleavage_rates(3.0e9, temperature)
                raw = rates["cleavage_rate_raw_s"]
                expected = float(gammainc(3.0, min(raw * 1.0e-6, 1.0e12)) / 1.0e-6)
                self.assertEqual(rates["cleavage_rate_effective_s"], expected)
                self.assertEqual(canonical_effective_cleavage_rate(raw), expected)

    def test_threshold_is_exact_authoritative_seedsequence_draw(self):
        state = self.make(seed=2420)
        rng = np.random.default_rng(np.random.SeedSequence([2420, 0]))
        self.assertEqual(state.hazard_threshold_action, max(float(rng.exponential(1.0)), 1e-12))
        self.assertEqual(state.audit["fatigue_endpoint"], ENDPOINT)
        self.assertFalse(state.audit["pd_candidate_site_multiplier"])
        self.assertFalse(state.audit["k2_completion_gate"])

    def test_ultrasmall_canonical_rate_remains_physical_in_log_domain(self):
        log_rate = canonical_effective_cleavage_log_rate(-1000.0)
        self.assertTrue(np.isfinite(log_rate))
        self.assertLess(log_rate, -2000.0)
        self.assertEqual(canonical_effective_cleavage_rate(0.0), 0.0)

    def test_audited_root_tensor_projection_preserves_signed_channels(self):
        state = self.make()
        tensor = np.array([[2.0e9, 0.4e9], [0.4e9, 3.0e9]])
        drive = state.mpz.resolve_root_tensor(tensor)
        direct = state.mpz.modules["anisotropic_emission_v10174"].resolve_channel_drives(
            tensor, [tensor, tensor], 0.0,
        )
        np.testing.assert_array_equal(drive["tau_signed_Pa"], direct["tau_signed_Pa"])
        self.assertEqual(drive["opening_stress_Pa"], 3.0e9)

    def test_phase_block_restart_is_exact_and_has_no_post_birth_state(self):
        tensors = np.array([
            [[2.0e9, 0.3e9], [0.3e9, 3.0e9]],
            [[1.0e9, -0.2e9], [-0.2e9, 1.5e9]],
        ])
        a = self.make(); b = self.make()
        a.hazard_threshold_action = b.hazard_threshold_action = 1e100
        a.advance_fem_phase_block(0.4, 1000.0, 300.0, tensors)
        capsule = a.capsule()
        resumed = self.make(); resumed.restore_capsule(capsule)
        a.advance_fem_phase_block(0.6, 1000.0, 300.0, tensors)
        resumed.advance_fem_phase_block(0.6, 1000.0, 300.0, tensors)
        self.assertEqual(a.capsule()["hazard_rng_state"], resumed.capsule()["hazard_rng_state"])
        self.assertEqual(a.log_cumulative_cleavage_hazard, resumed.log_cumulative_cleavage_hazard)
        for key in a.mpz._CAPSULE_ARRAYS:
            np.testing.assert_array_equal(getattr(a.mpz.state, key), getattr(resumed.mpz.state, key))
        self.assertNotIn("pd", a.diagnostics())

    def test_atomic_restart_preserves_signed_and_stochastic_state(self):
        history = [(2e-7, 3.0e9, [1.4e9, -1.1e9]),
                   (3e-7, 3.2e9, [-1.2e9, 1.5e9]),
                   (1e-7, 2.8e9, [1.0e9, -0.9e9])]
        uninterrupted = self.make()
        interrupted = self.make()
        for dt, stress, shear in history:
            uninterrupted.advance(dt, 300.0, stress, shear)
        interrupted.advance(history[0][0], 300.0, history[0][1], history[0][2])
        capsule = interrupted.capsule()
        resumed = self.make()
        resumed.restore_capsule(copy.deepcopy(capsule))
        for dt, stress, shear in history[1:]:
            resumed.advance(dt, 300.0, stress, shear)
        self.assertEqual(uninterrupted.capsule()["hazard_rng_state"], resumed.capsule()["hazard_rng_state"])
        self.assertEqual(uninterrupted.hazard_threshold_action, resumed.hazard_threshold_action)
        self.assertEqual(uninterrupted.cumulative_cleavage_hazard, resumed.cumulative_cleavage_hazard)
        for key in uninterrupted.mpz._CAPSULE_ARRAYS:
            np.testing.assert_array_equal(
                getattr(uninterrupted.mpz.state, key), getattr(resumed.mpz.state, key)
            )

    def test_first_passage_localizes_and_does_not_draw_post_birth_threshold(self):
        state = self.make()
        state.hazard_threshold_action = 1.0e-12
        rng_before = copy.deepcopy(state._hazard_rng.bit_generator.state)
        result = state.advance(1.0, 1300.0, 20.0e9, [0.0, 0.0])
        self.assertTrue(result["stable_crack_birth"])
        self.assertLess(result["dt_consumed_s"], 1.0)
        self.assertGreater(result["dt_unused_s"], 0.0)
        self.assertEqual(state.cumulative_cleavage_hazard, state.hazard_threshold_action)
        self.assertEqual(rng_before, state._hazard_rng.bit_generator.state)
        self.assertEqual(state.stable_crack_birth_time_s, state.time_s)

    def test_ten_way_partition_converges_for_full_signed_hazard_path(self):
        whole = self.make()
        split = self.make()
        whole.hazard_threshold_action = split.hazard_threshold_action = 1.0e100
        whole.advance(1e-6, 300.0, 3.0e9, [1.4e9, -1.1e9])
        for _ in range(10):
            split.advance(1e-7, 300.0, 3.0e9, [1.4e9, -1.1e9])
        self.assertAlmostEqual(
            whole.cumulative_cleavage_hazard,
            split.cumulative_cleavage_hazard,
            delta=max(split.cumulative_cleavage_hazard * 3e-4, 1e-300),
        )
        for key in ("mobile_positive", "mobile_negative", "retained_positive", "retained_negative"):
            np.testing.assert_allclose(
                whole.mpz.summary()[key], split.mpz.summary()[key], rtol=3e-4, atol=1e-18
            )


if __name__ == "__main__":
    unittest.main()
