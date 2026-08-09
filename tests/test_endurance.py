import unittest

import numpy as np

from arrhenius_fracture.endurance import analyze_endurance_tail, classify_physical_handoff_endurance


class EnduranceTailTests(unittest.TestCase):
    def setUp(self):
        self.cycles = np.geomspace(10.0, 1.0e8, 240)

    def classify(self, hazard):
        return analyze_endurance_tail(self.cycles, hazard, tail_model_is_exact=True).classification

    def test_constant_positive_hazard_has_no_endurance(self):
        self.assertEqual(self.classify(np.full_like(self.cycles, 2.0e-8)), "no_endurance")

    def test_exponential_decay_supports_endurance(self):
        cycles = np.linspace(1.0, 1000.0, 240)
        result = analyze_endurance_tail(cycles, 3.0e-4 * np.exp(-cycles / 170.0), tail_model_is_exact=True)
        self.assertEqual(result.classification, "endurance_supported")
        self.assertEqual(result.tail_model, "exponential")
        self.assertGreater(result.physics_survival_asymptote_lower_bound, 0.0)

    def test_power_half_has_no_endurance(self):
        self.assertEqual(self.classify(1.0e-3 * self.cycles ** -0.5), "no_endurance")

    def test_power_one_has_no_endurance(self):
        self.assertEqual(self.classify(1.0e-3 * self.cycles ** -1.0), "no_endurance")

    def test_power_one_and_half_supports_endurance(self):
        result = analyze_endurance_tail(self.cycles, 1.0e-3 * self.cycles ** -1.5, tail_model_is_exact=True)
        self.assertEqual(result.classification, "endurance_supported")
        self.assertIsNotNone(result.empirical_extrapolated_remaining_integrated_hazard)

    def test_empirical_integrable_fit_requires_physics_bound(self):
        result = analyze_endurance_tail(self.cycles, 1.0e-3 * self.cycles ** -1.5)
        self.assertEqual(result.empirical_tail_evidence, "integrable_tail_observed")
        self.assertEqual(result.physics_asymptotic_classification, "undetermined")

    def test_exponential_with_hidden_positive_floor_is_not_supported(self):
        cycles = np.linspace(1.0, 1000.0, 240)
        hazard = 3.0e-4 * np.exp(-cycles / 170.0) + 1.0e-20
        result = analyze_endurance_tail(cycles, hazard)
        self.assertNotEqual(result.physics_asymptotic_classification, "endurance_supported")

    def test_power_law_with_hidden_positive_floor_is_not_supported(self):
        hazard = 1.0e-3 * self.cycles ** -1.5 + 1.0e-30
        result = analyze_endurance_tail(self.cycles, hazard)
        self.assertNotEqual(result.physics_asymptotic_classification, "endurance_supported")

    def test_late_crossover_to_divergent_power_is_detected(self):
        cycles = np.geomspace(10.0, 1.0e12, 360)
        crossover = 1.0e8
        early = 1.0e-3 * cycles ** -1.5
        continuity = 1.0e-3 * crossover ** (-1.5 + 0.8)
        hazard = np.where(cycles <= crossover, early, continuity * cycles ** -0.8)
        result = analyze_endurance_tail(cycles, hazard, tail_model_is_exact=True)
        self.assertNotEqual(result.classification, "endurance_supported")

    def test_short_integrable_looking_dynamic_range_is_undetermined(self):
        cycles = np.linspace(1.0, 2.0, 60)
        hazard = np.exp(-0.1 * cycles)
        result = analyze_endurance_tail(cycles, hazard)
        self.assertEqual(result.classification, "undetermined")
        self.assertEqual(result.reason_code, "insufficient_tail_dynamic_range")

    def test_physics_bound_can_support_endurance(self):
        result = analyze_endurance_tail(
            self.cycles,
            1.0e-3 * self.cycles ** -1.5,
            physics_remaining_integrated_hazard_upper_bound=0.01,
        )
        self.assertEqual(result.physics_asymptotic_classification, "endurance_supported")
        endpoints = classify_physical_handoff_endurance(result)
        self.assertEqual(endpoints.physical_handoff_endurance_classification, "endurance_supported")

    def test_divergent_birth_does_not_prove_physical_handoff(self):
        birth = analyze_endurance_tail(
            self.cycles,
            np.full_like(self.cycles, 2.0e-8),
            tail_model_is_exact=True,
        )
        self.assertEqual(birth.classification, "no_endurance")
        endpoints = classify_physical_handoff_endurance(birth)
        self.assertEqual(endpoints.physical_handoff_endurance_classification, "undetermined")

    def test_short_noisy_tail_is_undetermined(self):
        rng = np.random.default_rng(20260809)
        cycles = np.geomspace(10.0, 1.0e4, 18)
        hazard = 1.0e-4 * cycles ** -1.05 * np.exp(rng.normal(0.0, 0.8, len(cycles)))
        self.assertEqual(analyze_endurance_tail(cycles, hazard).classification, "undetermined")

    def test_zero_hazard_cannot_manufacture_endurance(self):
        hazard = self.cycles ** -1.5
        hazard[-1] = 0.0
        result = analyze_endurance_tail(self.cycles, hazard)
        self.assertEqual(result.classification, "undetermined")
        self.assertEqual(result.reason_code, "invalid_or_clipped_hazard")


if __name__ == "__main__":
    unittest.main()
