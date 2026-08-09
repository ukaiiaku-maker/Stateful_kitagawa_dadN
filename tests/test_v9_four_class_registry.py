import unittest
from pathlib import Path

from arrhenius_fracture.v9_four_class_registry import (
    EXPECTED, FORBIDDEN_ALTERNATES, select_canonical_option,
)


SOURCE = Path("/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1")


class V9FourClassRegistryTests(unittest.TestCase):
    def test_exact_four_options_resolve_through_audited_selector(self):
        for option_id, (class_label, candidate_id, _) in EXPECTED.items():
            selected, audit = select_canonical_option(option_id, SOURCE)
            self.assertEqual(selected.option_key, option_id)
            self.assertEqual(selected.candidate_id, candidate_id)
            self.assertEqual(audit["class_label"], class_label)
            self.assertEqual(audit["exact_registry_row"], selected.row)

    def test_later_bulk_plasticity_alternates_fail_closed(self):
        for option_id in FORBIDDEN_ALTERNATES:
            with self.assertRaises(ValueError):
                select_canonical_option(option_id, SOURCE)


if __name__ == "__main__":
    unittest.main()
