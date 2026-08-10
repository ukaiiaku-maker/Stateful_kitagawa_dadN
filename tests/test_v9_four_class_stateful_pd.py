import math
from types import SimpleNamespace
import unittest

import numpy as np

from arrhenius_fracture.sn_arrhenius_chain import AuditedExpFloorBarrier
from arrhenius_fracture.sn_feature_geometry_v8_7 import BluntNotchGeometry
from arrhenius_fracture.sn_pd2d_stateful_v9_transactional import (
    ScratchExpFloorBarrier, build_parser as build_pd_parser,
)
from scripts.run_v9_four_class_stateful_pd import build_parser as build_four_class_parser


class FourClassStatefulPDTests(unittest.TestCase):
    def test_historical_blunt_notch_contract(self):
        geometry = BluntNotchGeometry()
        self.assertEqual(geometry.depth_a, 150e-6)
        self.assertEqual(geometry.half_height_b, 300e-6)
        self.assertAlmostEqual(geometry.root_radius, 600e-6)

    def test_four_class_production_image_default_is_none(self):
        args = build_four_class_parser().parse_args([
            "--material-class", "Peak", "--sigma-a-MPa", "500",
            "--source-root", "/tmp/source",
        ])
        self.assertEqual(args.pd_image_policy, "none")
        legacy = build_pd_parser().parse_args([])
        self.assertEqual(legacy.pd_image_policy, "selected")

    def test_audited_exp_floor_matches_definition(self):
        barrier = AuditedExpFloorBarrier(
            G00_eV=2.4, gT_eV_per_K=0.002, sigc0_Pa=5e9,
            sT_Pa_per_K=-1e6, alpha=0.7, exponent=1.2,
            floor_fraction=0.03, rate_prefactor=1e12,
        )
        stress = np.array([0.0, 1e9, 4e9])
        T = 300.0
        dT = T - 481.33
        G0 = 2.4 + 0.002*dT
        sigc = 5e9 - 1e6*dT
        floor = max(1e-4, 0.03*G0)
        expected = floor + (G0-floor)*np.exp(-0.7*(stress/sigc)**1.2)
        np.testing.assert_allclose(barrier.deltaG_eV(stress, T), expected, rtol=0, atol=0)

    def test_pd_cleavage_audited_linear_surface(self):
        args = SimpleNamespace(
            crack_G00_eV=2.4, crack_sigc0_GPa=5.0,
            crack_exp_a=0.7, crack_exp_n=1.2, crack_floor_frac=0.03,
            crack_T_mode="audited_linear", crack_Tref_K=481.33,
            crack_mu_dlnmu_dT_per_K=0.0, crack_G0_mu_power=1.0,
            crack_sigc_mu_power=1.0, S_crack_kB=0.0, nu0_crack=1e12,
            crack_gT_eV_per_K=0.002, crack_sT_GPa_per_K=-0.001,
        )
        pd_barrier = ScratchExpFloorBarrier(args)
        direct = AuditedExpFloorBarrier(
            G00_eV=2.4, gT_eV_per_K=0.002, sigc0_Pa=5e9,
            sT_Pa_per_K=-1e6, alpha=0.7, exponent=1.2,
            floor_fraction=0.03, rate_prefactor=1e12,
        )
        stress = np.array([0.0, 1e9, 4e9])
        np.testing.assert_allclose(
            pd_barrier.deltaG_eV(stress, 300.0),
            direct.deltaG_eV(stress, 300.0), rtol=2e-15, atol=0,
        )


if __name__ == "__main__":
    unittest.main()
