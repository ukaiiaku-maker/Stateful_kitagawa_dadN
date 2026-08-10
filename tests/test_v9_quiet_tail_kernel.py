import unittest
import csv
import tempfile
from pathlib import Path

import numpy as np

from arrhenius_fracture.v9_quiet_tail_kernel import stationary_marked_renewal
from scripts.run_v9_quiet_tail_kernel import atomic_merge_csv


class QuietTailKernelTests(unittest.TestCase):
    def test_atomic_csv_keys_normalize_equivalent_numeric_spellings(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "result.csv"
            atomic_merge_csv(path, ["class", "stress", "value"],
                             [{"class": "A", "stress": "329", "value": 1}],
                             ("class", "stress"))
            atomic_merge_csv(path, ["class", "stress", "value"],
                             [{"class": "A", "stress": 329.0, "value": 2}],
                             ("class", "stress"))
            with path.open(newline="") as stream:
                rows = list(csv.DictReader(stream))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["value"], "2")

    def test_all_admitted_marks_reduce_to_poisson_attempt_survival(self):
        rows = [
            {"phase_index": index, "phase": phase, "Xi": xi,
             "admitted_length_m": 1.0e-6}
            for index, phase in enumerate((0.0, 0.5)) for xi in (0.1, 1.0, 4.0)
        ]
        result = stationary_marked_renewal(
            0.7, 1.0e-4, np.array([5.0e-5, 1.0e-4]), rows)
        self.assertEqual(result["stationary_attempt_admission_probability"], 1.0)
        self.assertAlmostEqual(result["stable_birth_survival"], np.exp(-0.7))

    def test_all_rejected_marks_preserve_stable_birth_survival(self):
        rows = [
            {"phase_index": index, "phase": phase, "Xi": xi,
             "admitted_length_m": 0.0}
            for index, phase in enumerate((0.0, 0.5)) for xi in (0.1, 1.0, 4.0)
        ]
        result = stationary_marked_renewal(
            12.0, 1.0e-4, np.array([5.0e-5, 1.0e-4]), rows)
        self.assertEqual(result["stationary_attempt_admission_probability"], 0.0)
        self.assertEqual(result["stable_birth_survival"], 1.0)


if __name__ == "__main__":
    unittest.main()
