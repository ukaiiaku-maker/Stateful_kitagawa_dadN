import unittest

import numpy as np

from arrhenius_fracture.config import ElasticProperties
from arrhenius_fracture.sn_arrhenius_chain import build_chain_from_namespace
from arrhenius_fracture.sn_feature_geometry_v8_7 import (
    BluntNotchGeometry, make_blunt_edge_notch_mesh,
)
from arrhenius_fracture.sn_intact_fem import (
    affine_stress_control_displacements, cycle_stress_histories, plane_strain_D,
)
from arrhenius_fracture.sn_pd2d_stateful_v8_7_generalized_features import (
    apply_representative_fatigue_model, build_parser,
)
from arrhenius_fracture.v9_fem_transaction import (
    EmbeddedFEMTransaction, FEMPhysicalState,
)
from arrhenius_fracture.v9_cached_fem import CachedIntactFEM


class V9FEMTransactionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        args = build_parser().parse_args([
            "--resolution-profile", "custom", "--nx", "10", "--ny", "18",
            "--jitter", "0", "--plastic-n-phase", "8",
            "--feature-type", "ellipse", "--notch-depth-m", "1.5e-4",
            "--notch-half-height-m", "3e-4",
        ])
        apply_representative_fatigue_model(args)
        geom = BluntNotchGeometry(
            args.Lx, args.Ly, args.notch_depth_m, args.notch_half_height_m,
            feature_type=args.feature_type,
        )
        mesh, bnd, _ = make_blunt_edge_notch_mesh(
            geom, nx=args.nx, ny=args.ny, jitter=args.jitter,
            root_h_fine=args.root_h_fine, seed=args.seed,
        )
        mat = ElasticProperties(E=args.E_GPa * 1e9, nu=args.nu, b=args.b_m, Tm=args.Tm_K)
        chain = build_chain_from_namespace(args, mat.b)
        sigma_max = 2.0 * 700e6 / (1.0 - args.R)
        cls.model = EmbeddedFEMTransaction(
            mesh=mesh, boundaries=bnd, material=mat, Dmat=plane_strain_D(mat),
            plastic_chain=chain, args=args, sigma_max_Pa=sigma_max,
            sigma_min_Pa=args.R * sigma_max, relative_tolerance=1e-4,
        )
        cls.cached_model = EmbeddedFEMTransaction(
            mesh=mesh, boundaries=bnd, material=mat, Dmat=plane_strain_D(mat),
            plastic_chain=chain, args=args, sigma_max_Pa=sigma_max,
            sigma_min_Pa=args.R * sigma_max, relative_tolerance=1e-4,
            cached_fem=CachedIntactFEM(mesh, bnd, mat, plane_strain_D(mat)),
        )
        cls.mesh, cls.bnd, cls.mat, cls.args = mesh, bnd, mat, args
        cls.initial = FEMPhysicalState(
            np.zeros((3, mesh.ne)), np.full(mesh.ne, args.rho0),
            np.zeros(mesh.ne), np.zeros(mesh.ndof), 0.0,
        )

    def test_proposal_is_side_effect_free_and_embedded_error_decreases(self):
        ep_before = self.initial.ep_gp.copy()
        coarse = self.model.propose(self.initial, 1.0)
        fine = self.model.propose(self.initial, 0.05)
        np.testing.assert_array_equal(self.initial.ep_gp, ep_before)
        self.assertTrue(np.isfinite(coarse.normalized_error))
        self.assertLess(fine.normalized_error, coarse.normalized_error)

    def test_absolute_strain_scale_does_not_overresolve_inactive_component(self):
        a = np.asarray([0.0])
        b = np.asarray([1.0e-13])
        controlled = self.model._scaled_error(a, b, 1.0e-8, 2.0e-4)
        machine_zero = self.model._scaled_error(a, b, 1.0e-18, 2.0e-4)
        self.assertLess(controlled, 0.11)
        self.assertGreater(machine_zero, 1.0e3)

    def test_physical_absolute_strain_tolerance_matches_tighter_reference(self):
        prior = self.cached_model.ep_atol
        try:
            self.cached_model.ep_atol = 1.0e-8
            production = self.cached_model.advance(
                self.initial, cycle_start=0.0, cycle_end=2.0, initial_block_dN=2.0
            )
            self.cached_model.ep_atol = 1.0e-9
            tighter = self.cached_model.advance(
                self.initial, cycle_start=0.0, cycle_end=2.0, initial_block_dN=2.0
            )
        finally:
            self.cached_model.ep_atol = prior
        np.testing.assert_allclose(production.state.ep_gp, tighter.state.ep_gp, rtol=3e-4, atol=2e-8)
        np.testing.assert_allclose(production.state.rho_gp, tighter.state.rho_gp, rtol=2e-7, atol=2.0)

    def test_cached_real_fem_is_physics_equivalent(self):
        native = self.model.propose(self.initial, 0.05)
        cached = self.cached_model.propose(self.initial, 0.05)
        np.testing.assert_allclose(cached.state.ep_gp, native.state.ep_gp, rtol=2e-12, atol=1e-20)
        np.testing.assert_allclose(cached.state.rho_gp, native.state.rho_gp, rtol=2e-14, atol=0.02)
        np.testing.assert_allclose(cached.state.u, native.state.u, rtol=2e-12, atol=1e-20)

    def test_affine_cached_phase_history_matches_repeated_solves(self):
        Dmat = plane_strain_D(self.mat)
        Umax, Umin, u_zero, _, _ = affine_stress_control_displacements(
            self.mesh, self.bnd, self.mat, Dmat, self.initial.ep_gp,
            self.model.sigma_max_Pa, self.model.sigma_min_Pa, self.initial.u,
        )
        native = cycle_stress_histories(
            self.mesh, self.bnd, self.mat, Dmat, self.initial.ep_gp,
            Umax, Umin, 16, u_zero,
        )
        cached = self.cached_model.cached_fem.stress_histories(
            self.initial.ep_gp, Umax, Umin, 16, u_zero,
        )
        for key in ("sigma_node", "seq_node", "s1_node", "Ftop"):
            np.testing.assert_allclose(cached[key], native[key], rtol=3e-10, atol=1e-4)
        # Tensile-energy gating is discontinuous at s1=0; roundoff can toggle
        # a few essentially unloaded nodes without affecting the rate drivers.
        np.testing.assert_allclose(cached["psi_node"], native["psi_node"], rtol=1e-3, atol=2e3)
        np.testing.assert_allclose(cached["u_end"], native["u_end"], rtol=3e-12, atol=2e-20)

    def test_accepted_heun_state_converges_under_partition(self):
        whole = self.model.propose(self.initial, 0.02).state
        half = self.model.propose(self.initial, 0.01).state
        split = self.model.propose(half, 0.01).state
        np.testing.assert_allclose(whole.ep_gp, split.ep_gp, rtol=3e-4, atol=1e-19)
        np.testing.assert_allclose(whole.rho_gp, split.rho_gp, rtol=3e-8, atol=1.0)

    def test_adaptive_real_trajectory_is_initial_partition_invariant(self):
        a = self.model.advance(
            self.initial, cycle_start=0.0, cycle_end=0.2, initial_block_dN=0.2
        )
        b = self.model.advance(
            self.initial, cycle_start=0.0, cycle_end=0.2, initial_block_dN=0.01
        )
        np.testing.assert_allclose(a.state.ep_gp, b.state.ep_gp, rtol=2e-4, atol=2e-18)
        np.testing.assert_allclose(a.state.rho_gp, b.state.rho_gp, rtol=2e-7, atol=2.0)
        np.testing.assert_allclose(a.state.epsp_acc_gp, b.state.epsp_acc_gp, rtol=2e-4, atol=2e-18)

    def test_adaptive_real_trajectory_restart_is_equivalent(self):
        uninterrupted = self.model.advance(
            self.initial, cycle_start=0.0, cycle_end=0.2, initial_block_dN=0.03
        )
        first = self.model.advance(
            self.initial, cycle_start=0.0, cycle_end=0.08, initial_block_dN=0.03
        )
        resumed = self.model.advance(
            first.state, cycle_start=first.cycle, cycle_end=0.2,
            initial_block_dN=first.next_block_dN,
        )
        np.testing.assert_allclose(uninterrupted.state.ep_gp, resumed.state.ep_gp, rtol=2e-4, atol=2e-18)
        np.testing.assert_allclose(uninterrupted.state.rho_gp, resumed.state.rho_gp, rtol=2e-7, atol=2.0)


if __name__ == "__main__":
    unittest.main()
