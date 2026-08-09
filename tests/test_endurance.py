import unittest

import numpy as np

from arrhenius_fracture.endurance import analyze_endurance_tail


class EnduranceTailTests(unittest.TestCase):
    def setUp(self):
        self.cycles = np.geomspace(10.0, 1.0e8, 240)

    def classify(self, hazard):
        return analyze_endurance_tail(self.cycles, hazard).classification

    def test_constant_positive_hazard_has_no_endurance(self):
        self.assertEqual(self.classify(np.full_like(self.cycles, 2.0e-8)), "no_endurance")

    def test_exponential_decay_supports_endurance(self):
        cycles = np.linspace(1.0, 1000.0, 240)
        result = analyze_endurance_tail(cycles, 3.0e-4 * np.exp(-cycles / 170.0))
        self.assertEqual(result.classification, "endurance_supported")
        self.assertEqual(result.tail_model, "exponential")
        self.assertGreater(result.survival_asymptote, 0.0)

    def test_power_half_has_no_endurance(self):
        self.assertEqual(self.classify(1.0e-3 * self.cycles ** -0.5), "no_endurance")

    def test_power_one_has_no_endurance(self):
        self.assertEqual(self.classify(1.0e-3 * self.cycles ** -1.0), "no_endurance")

    def test_power_one_and_half_supports_endurance(self):
        result = analyze_endurance_tail(self.cycles, 1.0e-3 * self.cycles ** -1.5)
        self.assertEqual(result.classification, "endurance_supported")
        self.assertIsNotNone(result.residual_integrated_hazard)

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
