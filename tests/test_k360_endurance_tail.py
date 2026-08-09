import json
from pathlib import Path
import unittest

from scripts.analyze_k360_endurance_tail import analyze


class K360EnduranceTailRegressionTests(unittest.TestCase):
    def test_committed_metrics_are_reproduced(self):
        root = Path(__file__).resolve().parents[1]
        result = analyze(root / "reference/k360_v2_8_tail/K360_CENSOR_DERIVED_TAIL.csv")
        expected = json.loads((root / "reference/k360_v2_8_tail/K360_CENSOR_TAIL_AUDIT.json").read_text())
        self.assertEqual(result["row_count"], 190)
        self.assertAlmostEqual(result["cumulative_expected_births_end"], expected["cumulative_expected_births_end"], places=14)
        for actual, reference in zip(result["nested_power_law_fits"], expected["nested_power_law_fits"]):
            self.assertAlmostEqual(actual["exponent"], reference["exponent"], places=12)
            self.assertAlmostEqual(actual["r_squared"], reference["r_squared"], places=14)
            self.assertAlmostEqual(actual["empirical_extrapolated_remaining_integrated_hazard"], reference["residual_integrated_hazard"], places=14)
        self.assertAlmostEqual(result["naive_fit_based_survival_asymptote"], expected["naive_fit_based_survival_asymptote"], places=14)
        self.assertEqual(result["empirical_tail_evidence"], "integrable_tail_observed")
        self.assertEqual(result["physics_asymptotic_classification"], "undetermined")


if __name__ == "__main__":
    unittest.main()
