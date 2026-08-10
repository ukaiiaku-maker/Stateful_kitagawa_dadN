import unittest
import csv
import tempfile
from pathlib import Path

import numpy as np

from arrhenius_fracture.v9_quiet_tail_kernel import stationary_marked_renewal
from scripts.run_v9_quiet_tail_kernel import atomic_merge_csv
from scripts.rebuild_v9_endurance_descent import reconstruct_direct_checkpoint_survival


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

    def test_direct_checkpoint_H_replaces_stationary_projection(self):
        states = [
            {"material_class": "weak-T", "option_id": "row", "sigma_a_MPa": 329,
             "N": 1e8, "H_cleave": 7.8e-6},
            {"material_class": "weak-T", "option_id": "row", "sigma_a_MPa": 329,
             "N": 1e10, "H_cleave": 7.7e-4},
        ]
        gates = [
            {"material_class": "weak-T", "sigma_a_MPa": 329, "N": N,
             "admitted_length_m": 1e-6}
            for N in (1e8, 1e10)
        ]
        rows = reconstruct_direct_checkpoint_survival(states, gates)
        self.assertEqual(rows[-1]["H_cleave"], 7.7e-4)
        self.assertAlmostEqual(rows[-1]["S_stable_birth"], np.exp(-7.7e-4))
        self.assertEqual(rows[-1]["estimator"],
                         "direct_checkpoint_H_plus_evaluated_gate_admission")

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
