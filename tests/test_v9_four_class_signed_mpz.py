import unittest
from pathlib import Path
import inspect

import numpy as np

from arrhenius_fracture.v9_four_class_registry import EXPECTED
from arrhenius_fracture.v9_four_class_signed_mpz import SignedMPZPreBirthState


SOURCE = Path("/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1")


class SignedMPZPortTests(unittest.TestCase):
    def make(self, option):
        return SignedMPZPreBirthState(
            option, SOURCE, shear_modulus_Pa=160.15625e9, poisson=0.28,
            burgers_m=2.74e-10, initial_tip_radius_m=1.0e-6,
        )

    def test_all_four_rows_execute_independent_signed_closure(self):
        for option in EXPECTED:
            port = self.make(option)
            result = port.advance(1e-3, 300.0, 3.0e9, [1.4e9, -1.1e9])
            self.assertEqual(result["mobile_positive"].shape[0], 2)
            self.assertEqual(result["mobile_negative"].shape[0], 2)
            self.assertGreaterEqual(result["peierls_rate_s"], 0.0)
            self.assertGreaterEqual(result["taylor_completion_rate_s"], 0.0)
            self.assertIn("constitutive_source_hashes", port.audit)

    def test_copy_partition_keeps_authoritative_signed_coordinates(self):
        option = next(iter(EXPECTED))
        whole = self.make(option)
        split = self.make(option)
        whole.advance(1e-6, 300.0, 3.0e9, [1.4e9, -1.1e9])
        for _ in range(10):
            split.advance(1e-7, 300.0, 3.0e9, [1.4e9, -1.1e9])
        # Backward-Euler source emission is convergent, not algebraically
        # partition invariant.  This is the numerical error gate for the v9
        # transaction controller, while every signed coordinate is retained.
        for key in ("mobile_positive", "mobile_negative", "retained_positive", "retained_negative"):
            np.testing.assert_allclose(whole.summary()[key], split.summary()[key], rtol=3e-4, atol=1e-18)

    def test_audited_persistent_emission_has_no_discrete_site_rng_ledger(self):
        port = self.make(next(iter(EXPECTED)))
        source = inspect.getsource(
            port.modules["persistent_site_source_v10221"]._persistent_emit
        )
        self.assertIn("solve_backstress_limited_activations", source)
        self.assertIn("Legacy capacity arrays remain full", source)
        for forbidden in ("rng", "threshold", "exponential(", "site_id"):
            self.assertNotIn(forbidden, source)


if __name__ == "__main__":
    unittest.main()
