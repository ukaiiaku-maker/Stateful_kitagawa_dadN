import unittest

from scripts.update_v9_sn_campaign import canonical_hash, classify, select_next_stress


class V9SNCampaignTests(unittest.TestCase):
    def test_request_hash_is_order_independent(self):
        self.assertEqual(canonical_hash({"a": 1, "b": 2}), canonical_hash({"b": 2, "a": 1}))

    def test_censor_is_not_promoted_to_endurance(self):
        summary = {"status": "right_censored", "geometry_saturated": False}
        self.assertEqual(classify(summary, None), "right_censored")

    def test_invalid_handoff_fails_closed(self):
        summary = {"status": "physical_handoff", "pd_handoff_pass_final": False}
        self.assertEqual(classify(summary, None), "invalid")

    def test_selector_interpolates_one_condition_into_log_life_gap(self):
        conditions = [
            {"classification": "physical_handoff", "sigma_a_MPa": 800.0, "cycles": 1e4},
            {"classification": "physical_handoff", "sigma_a_MPa": 700.0, "cycles": 1e8},
            {"classification": "right_censored", "sigma_a_MPa": 650.0, "cycles": 1e9},
        ]
        result = select_next_stress(conditions)
        self.assertIsNotNone(result)
        self.assertEqual(result["sigma_a_MPa"], 750.0)
        self.assertEqual(result["target_log10_cycles"], 6.0)


if __name__ == "__main__":
    unittest.main()
