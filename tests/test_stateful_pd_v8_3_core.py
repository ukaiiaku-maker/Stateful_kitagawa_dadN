from __future__ import annotations

import tempfile
import unittest
import math
from pathlib import Path
from types import SimpleNamespace
import numpy as np

from arrhenius_fracture.config import ElasticProperties
from arrhenius_fracture.sn_geometry import (
    BluntNotchGeometry,
    identify_feature_surface_nodes,
    local_root_xy,
    make_blunt_edge_notch_mesh,
)
from arrhenius_fracture.stateful_peridynamics_v8_3 import StatefulPDConfig, StatefulPDPatch
from arrhenius_fracture.stateful_peridynamics_v8_7_local_front_spacing import (
    StatefulPDPatch as StatefulPDPatchV87,
)
from arrhenius_fracture.sn_pd2d_stateful_v8_3 import (
    MODEL_ID,
    SOURCE_SHA256,
    apply_representative_fatigue_model,
    build_parser,
    phase_resolved_delivery_rate,
    apply_resolution_profile,
    _save_case_checkpoint,
    _load_case_checkpoint,
    _checkpoint_signature,
    _geometry_resolution_audit,
    _geometry_taper_factor,
)


class StatefulPDV83CoreTests(unittest.TestCase):
    def make_patch(self):
        geom = BluntNotchGeometry()
        mesh, _, _ = make_blunt_edge_notch_mesh(
            geom, nx=18, ny=36, jitter=0.0, root_h_fine=35e-6, seed=1
        )
        feature = identify_feature_surface_nodes(mesh, geom)
        root = local_root_xy(mesh, feature)
        cfg = StatefulPDConfig(
            patch_radius_m=0.46e-3,
            horizon_m=105e-6,
            boundary_shell_m=100e-6,
            initiation_radius_m=280e-6,
            initiation_taper_m=60e-6,
            initiation_back_extent_m=120e-6,
            handoff_mode="physical",
            handoff_min_length_horizons=0.75,
            handoff_max_width_ratio=0.50,
            handoff_max_tip_width_horizons=1.5,
            handoff_min_boundary_clearance_horizons=-1.0,
            handoff_min_connected_bonds=3,
            handoff_min_orientation_coherence=0.55,
            handoff_min_axial_coverage=0.45,
            handoff_max_axial_gap_horizons=1.0,
            handoff_min_slenderness=2.0,
            surface_connection_horizons=1.0,
            delivery_hit_count=2.0,
            delivery_memory_s=1e-3,
        )
        patch = StatefulPDPatch(
            mesh, geom, root, ElasticProperties(), cfg,
            feature_surface_global_nodes=feature,
        )
        return mesh, patch

    def force_horizontal_crack(self, patch, state, length_horizons=1.6):
        rel = patch.bond_midpoints - patch.root_xy[None, :]
        # A horizontal crack surface is represented by bonds with a mostly
        # vertical bond normal whose midpoints trace +x from the root.
        normal_align = np.abs(patch.n[:, 1])
        sel = (
            (rel[:, 0] >= 0.0)
            & (rel[:, 0] <= length_horizons * patch.cfg.horizon_m)
            & (np.abs(rel[:, 1]) <= 0.30 * patch.cfg.horizon_m)
            & (normal_align >= 0.80)
        )
        ids = np.where(sel)[0]
        state.bond_damage[ids] = 1.0
        return ids


    def test_v8_modules_do_not_import_legacy_diffuse_fracture(self):
        import ast
        import inspect
        import arrhenius_fracture.sn_pd2d_stateful_v8_3 as driver
        import arrhenius_fracture.stateful_peridynamics_v8_3 as pdmod
        imported = []
        for module in (driver, pdmod):
            tree = ast.parse(inspect.getsource(module))
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    imported.extend(alias.name.lower() for alias in node.names)
                elif isinstance(node, ast.ImportFrom):
                    imported.append((node.module or "").lower())
        joined = " ".join(imported)
        for forbidden in ("phase_field", ".at1", ".at2"):
            self.assertNotIn(forbidden, joined)

    def test_version_identifier(self):
        self.assertIn("v8_3_continuous_first_passage_fixed_geometry", MODEL_ID)

    def test_state_initialization_and_site_measure(self):
        _, patch = self.make_patch()
        state = patch.initial_state()
        self.assertGreater(len(patch.xy), 10)
        self.assertGreater(len(patch.bonds), len(patch.xy))
        self.assertTrue(np.allclose(state.available, 1.0))
        self.assertGreater(np.sum(patch.mean_candidate_sites), 1.0)
        self.assertGreater(np.sum(state.candidate_sites), 0)
        self.assertTrue(np.all(state.candidate_sites[patch.boundary] == 0))
        self.assertGreater(len(patch.feature_surface_xy), 4)

    def test_candidate_and_event_rng_streams_are_reproducible(self):
        _, p1 = self.make_patch(); s1 = p1.initial_state()
        _, p2 = self.make_patch(); s2 = p2.initial_state()
        self.assertTrue(np.array_equal(s1.candidate_sites, s2.candidate_sites))
        self.assertEqual(p1._candidate_rng.bit_generator.state["bit_generator"], "PCG64")
        self.assertEqual(p1._event_rng.bit_generator.state["bit_generator"], "PCG64")

    def test_surface_connection_uses_physical_feature(self):
        _, patch = self.make_patch(); state = patch.initial_state()
        ids = self.force_horizontal_crack(patch, state, 1.4)
        self.assertGreater(len(ids), 2)
        audit = patch.crack_handoff_audit(state, 1.5e9, 2.5e9)
        self.assertTrue(audit.root_connected)
        self.assertLessEqual(audit.surface_gap_m, patch.cfg.horizon_m)

    def test_directional_seed_gate_prefers_crack_normal_bonds(self):
        _, patch = self.make_patch(); state = patch.initial_state()
        node = int(np.argmin(np.linalg.norm(patch.xy - patch.root_xy[None, :], axis=1)))
        state.stable_sites[node] = 1
        state.growth[node] = 1.0
        # normal = y => cos(2 theta)=-1, sin(2 theta)=0
        state.crack_normal_c2[node] = -1.0
        state.crack_normal_s2[node] = 0.0
        state.crack_orientation_weight[node] = 1.0
        activity, seed, _ = patch._directional_link_activity(state)
        vertical = np.abs(patch.n[:, 1]) > 0.85
        horizontal = np.abs(patch.n[:, 0]) > 0.85
        self.assertGreater(np.max(seed[vertical]), 0.0)
        self.assertGreater(np.max(seed[vertical]), 5.0 * max(np.max(seed[horizontal]), 1e-30))
        self.assertTrue(np.all(activity >= 0.0))

    def test_principal_direction_for_y_tension(self):
        _, patch = self.make_patch()
        sig = np.zeros((4, 3, len(patch.xy)))
        sig[:, 1, :] = 2.0e9
        normal, _, _ = patch._principal_normal_from_stress_history(sig)
        self.assertTrue(np.all(np.abs(normal[:, 1]) > 0.999))


    def test_directional_update_does_not_create_radial_star(self):
        mesh, patch = self.make_patch(); state = patch.initial_state()
        available = np.where(state.available_sites > 0)[0]
        node = int(available[np.argmin(np.linalg.norm(patch.xy[available] - patch.root_xy[None, :], axis=1))])
        state.available_sites[node] -= 1
        state.stable_sites[node] += 1
        state.available[node] = max(state.available[node] - 1.0 / max(state.candidate_sites[node], 1), 0.0)
        state.stable[node] = min(state.stable[node] + 1.0 / max(state.candidate_sites[node], 1), 1.0)
        state.growth[node] = 1.0
        state.crack_normal_c2[node] = -1.0
        state.crack_normal_s2[node] = 0.0
        state.crack_orientation_weight[node] = 1.0
        class InactiveBarrier:
            rate_prefactor = 0.0
            @staticmethod
            def deltaG_eV(sigma, T): return np.ones_like(np.asarray(sigma, float))*10.0
        sigma_hist = np.zeros((4, 3, mesh.nn)); sigma_hist[:, 1, :] = 5e9
        patch.update(
            state, InactiveBarrier(), sigma_hist, np.zeros((4, mesh.nn)),
            300.0, 1000.0, 2e6, 0.0, np.zeros(mesh.nn), np.zeros(mesh.nn),
            0.0, np.zeros(mesh.nn), np.ones(len(patch.xy)), np.ones(len(patch.bonds)),
        )
        ids = np.where(state.bond_damage > 0.5)[0]
        self.assertGreater(len(ids), 0)
        _, coherence = patch._director_average(patch.n[ids], state.bond_damage[ids])
        self.assertGreater(coherence, 0.60)

    def test_starburst_is_rejected(self):
        _, patch = self.make_patch(); state = patch.initial_state()
        node = int(np.argmin(np.linalg.norm(patch.xy - patch.root_xy[None, :], axis=1)))
        for b in patch.incident[node]:
            state.bond_damage[b] = 1.0
        audit = patch.crack_handoff_audit(state, 1.5e9, 2.5e9)
        self.assertFalse(audit.handoff_pass)
        self.assertTrue((not audit.orientation_pass) or (not audit.coverage_pass) or (not audit.slenderness_pass))

    def test_slender_surface_crack_can_pass(self):
        _, patch = self.make_patch(); state = patch.initial_state()
        ids = self.force_horizontal_crack(patch, state, 1.8)
        self.assertGreaterEqual(len(ids), 3)
        audit = patch.crack_handoff_audit(state, 1.5e9, 2.5e9)
        self.assertTrue(audit.root_connected)
        self.assertGreater(audit.centerline_length_m, 0.0)
        self.assertGreater(audit.orientation_coherence, 0.55)
        self.assertGreater(audit.remote_K_MPam05, 0.0)

    def test_detached_component_is_not_root_connected(self):
        _, patch = self.make_patch(); state = patch.initial_state()
        rel = patch.bond_midpoints - patch.root_xy[None, :]
        ids = np.where(
            (rel[:, 0] > 2.0 * patch.cfg.horizon_m)
            & (rel[:, 0] < 3.0 * patch.cfg.horizon_m)
            & (np.abs(rel[:, 1]) < 0.25 * patch.cfg.horizon_m)
            & (np.abs(patch.n[:, 1]) > 0.8)
        )[0]
        state.bond_damage[ids] = 1.0
        audit = patch.crack_handoff_audit(state, 1.5e9, 2.5e9)
        self.assertFalse(audit.root_connected)

    def test_separated_delivery_and_nucleation_kernels(self):
        _, patch = self.make_patch()
        nphase, npoint = 16, len(patch.xy)
        memory0 = np.zeros(npoint)
        delivery = np.full((nphase, npoint), 500.0)
        nuc1 = np.full((nphase, npoint), 2.0)
        nuc2 = np.full((nphase, npoint), 6.0)
        a = patch._phase_resolved_delivery_nucleation(delivery, nuc1, 1000.0, memory0)
        b = patch._phase_resolved_delivery_nucleation(delivery, nuc2, 1000.0, memory0)
        self.assertTrue(np.allclose(a["memory_current_end"], b["memory_current_end"]))
        self.assertTrue(np.allclose(b["birth_hazard_current_cycle"], 3.0*a["birth_hazard_current_cycle"], rtol=1e-12, atol=1e-15))

    def test_delivery_pulse_increases_short_memory_completion(self):
        _, patch = self.make_patch()
        patch.cfg.delivery_memory_s = 1e-5
        nphase, npoint = 16, len(patch.xy)
        uniform = np.full((nphase, npoint), 100.0)
        pulse = np.zeros((nphase, npoint)); pulse[0, :] = 1600.0
        nucleation = np.full((nphase, npoint), 10.0)
        a = patch._phase_resolved_delivery_nucleation(uniform, nucleation, 1000.0, np.zeros(npoint))
        b = patch._phase_resolved_delivery_nucleation(pulse, nucleation, 1000.0, np.zeros(npoint))
        self.assertTrue(np.allclose(a["delivery_events_cycle"], b["delivery_events_cycle"]))
        self.assertGreater(np.max(b["birth_hazard_periodic_cycle"]), np.max(a["birth_hazard_periodic_cycle"]))

    def test_no_bond_softening_without_realized_stable_defect(self):
        mesh, patch = self.make_patch(); state = patch.initial_state()
        state.stable[:] = 0.2; state.growth[:] = 1.0
        class InactiveBarrier:
            rate_prefactor = 0.0
            @staticmethod
            def deltaG_eV(sigma, T): return np.ones_like(np.asarray(sigma, float))*10.0
        sigma_hist = np.zeros((4, 3, mesh.nn)); sigma_hist[:, 1, :] = 5e9
        diag = patch.update(state, InactiveBarrier(), sigma_hist, np.zeros((4, mesh.nn)), 300.0, 1000.0, 1e6, 0.0, np.zeros(mesh.nn), np.zeros(mesh.nn), 0.0, np.zeros(mesh.nn), np.ones(len(patch.xy)), np.ones(len(patch.bonds)))
        self.assertEqual(diag.realized_stable, 0)
        self.assertEqual(diag.max_bond_damage, 0.0)

    def test_representative_case64_m1_preset_and_defaults(self):
        args = build_parser().parse_args([]); apply_representative_fatigue_model(args)
        self.assertEqual(args.fatigue_model, "plastic_shielded_case64_M1")
        self.assertAlmostEqual(args.crack_G00_eV, 1.0)
        self.assertAlmostEqual(args.crack_sigc0_GPa, 3.0)
        self.assertEqual(args.delivery_source, "completed_flow")
        self.assertAlmostEqual(args.delivery_hit_count, 2.0)

    def test_completed_flow_delivery_source(self):
        args = build_parser().parse_args([])
        class DummyChain:
            @staticmethod
            def rates(seq, rho, T):
                shape = np.asarray(seq).shape
                return {"lambda_flow": np.full(shape, 7.0), "lambda_emit": np.full(shape, 11.0), "lambda_peierls": np.full(shape, 13.0), "lambda_taylor": np.full(shape, 17.0), "dot_ep": np.full(shape, 19.0)}
        rate = phase_resolved_delivery_rate(args, DummyChain(), np.zeros((4,6)), np.ones(6), 300.0)
        self.assertTrue(np.allclose(rate, 7.0))

    def test_h15_and_h10_resolution_profiles(self):
        args = build_parser().parse_args([]); apply_resolution_profile(args)
        self.assertEqual(args.resolution_profile, "h15")
        self.assertAlmostEqual(args.pd_horizon_m, 45e-6)
        args2 = build_parser().parse_args(["--resolution-profile", "h10"]); apply_resolution_profile(args2)
        self.assertAlmostEqual(args2.pd_horizon_m, 30e-6)

    def test_max_unique_gap_is_safe_for_singleton_and_duplicate_anchors(self):
        self.assertEqual(StatefulPDPatch._max_unique_gap([]), 0.0)
        self.assertEqual(StatefulPDPatch._max_unique_gap([0.0]), 0.0)
        self.assertEqual(StatefulPDPatch._max_unique_gap([0.0, 0.0, 0.0]), 0.0)
        self.assertAlmostEqual(StatefulPDPatch._max_unique_gap([0.0, 2.0, 2.0, 5.0]), 3.0)

    def test_single_surface_bond_audit_is_nonfatal_and_rejected(self):
        _, patch = self.make_patch(); state = patch.initial_state()
        gaps = patch._bond_surface_gap(np.arange(len(patch.bonds), dtype=int))
        bond_id = int(np.argmin(gaps))
        state.bond_damage[bond_id] = 1.0
        audit = patch.crack_handoff_audit(state, 1.5e9, 2.5e9)
        self.assertTrue(audit.root_connected)
        self.assertEqual(audit.connected_bonds, 1)
        self.assertTrue(np.isfinite(audit.max_axial_gap_m))
        self.assertFalse(audit.length_pass)
        self.assertFalse(audit.handoff_pass)

    def test_empty_director_average_is_safe(self):
        director, coherence = StatefulPDPatch._director_average(np.empty((0, 2)))
        self.assertTrue(np.all(np.isfinite(director)))
        self.assertEqual(coherence, 0.0)

    def test_zero_weight_director_average_is_safe(self):
        director, coherence = StatefulPDPatch._director_average(
            np.array([[1.0, 0.0], [0.0, 1.0]]), np.zeros(2)
        )
        self.assertTrue(np.all(np.isfinite(director)))
        self.assertEqual(coherence, 0.0)


    def test_front_capture_selects_one_surface_component(self):
        _, patch = self.make_patch(); state = patch.initial_state()
        ids = self.force_horizontal_crack(patch, state, 1.2)
        self.assertGreaterEqual(len(ids), patch.cfg.front_capture_min_bonds)
        captured = patch._capture_or_update_active_front(state, 123.0)
        self.assertTrue(captured)
        self.assertTrue(state.active_front)
        self.assertEqual(state.cycles_front_capture, 123.0)
        self.assertGreater(np.count_nonzero(state.active_front_bonds), 0)
        self.assertGreater(state.active_front_length_m, 0.0)

    def test_post_capture_birth_and_stabilization_are_suppressed(self):
        mesh, patch = self.make_patch(); state = patch.initial_state()
        state.active_front = True
        state.active_front_contact_xy = patch.root_xy.copy()
        state.active_front_tip_xy = patch.root_xy + np.array([patch.cfg.horizon_m, 0.0])
        state.active_front_length_m = patch.cfg.horizon_m
        state.active_front_normal_c2 = -1.0
        state.active_front_normal_s2 = 0.0
        class Barrier:
            rate_prefactor = 1e8
            @staticmethod
            def deltaG_eV(sigma, T): return np.full_like(np.asarray(sigma, float), 0.05)
        sigma = np.zeros((8, 3, mesh.nn)); sigma[:, 1, :] = 3e9
        delivery = np.full((8, mesh.nn), 1e4)
        rates = patch.preview_rates(
            state, Barrier(), sigma, delivery, 300.0, 1000.0,
            np.zeros(mesh.nn), np.zeros(mesh.nn), 0.0, np.zeros(mesh.nn),
            np.ones(len(patch.xy)), np.ones(len(patch.bonds)),
        )
        self.assertEqual(float(np.max(rates["mu_birth"])), 0.0)
        self.assertEqual(float(np.max(rates["mu_stab"])), 0.0)

    def test_active_front_linkage_is_tip_localized(self):
        _, patch = self.make_patch(); state = patch.initial_state()
        ids = self.force_horizontal_crack(patch, state, 1.0)
        self.assertTrue(patch._capture_or_update_active_front(state, 1.0))
        activity, _, _ = patch._directional_link_activity(state)
        normal, tangent, contact, _ = patch._active_front_frame(state)
        rel = patch.bond_midpoints - contact[None, :]
        axial = rel @ tangent
        trans = np.abs(rel @ normal)
        ahead = axial - state.active_front_length_m
        near_tip = (ahead >= -0.1*patch.cfg.horizon_m) & (ahead <= patch.cfg.horizon_m) & (trans <= 0.4*patch.cfg.horizon_m)
        far_wake = (ahead < -patch.cfg.horizon_m) | (trans > 1.5*patch.cfg.horizon_m)
        self.assertGreater(float(np.max(activity[near_tip])), 0.0)
        if np.any(far_wake):
            self.assertLessEqual(float(np.max(activity[far_wake])), 1e-12)

    def test_active_front_audit_ignores_offfront_damage_islands(self):
        _, patch = self.make_patch(); state = patch.initial_state()
        ids = self.force_horizontal_crack(patch, state, 1.8)
        self.assertTrue(patch._capture_or_update_active_front(state, 1.0))
        rel = patch.bond_midpoints - patch.root_xy[None, :]
        off = np.where((np.abs(rel[:, 1]) > 1.2*patch.cfg.horizon_m) & (rel[:, 0] > 0))[0][:40]
        state.bond_damage[off] = 1.0
        audit = patch.crack_handoff_audit(state, 1.5e9, 2.5e9)
        self.assertEqual(audit.connected_bonds, int(np.count_nonzero(state.active_front_bonds)))
        self.assertGreater(audit.orientation_coherence, 0.55)

    def test_diffuse_abort_is_patience_gated(self):
        _, patch = self.make_patch(); state = patch.initial_state()
        patch.cfg.diffuse_abort_min_bonds = 10
        patch.cfg.diffuse_abort_min_offfront_fraction = 0.5
        patch.cfg.diffuse_abort_patience_updates = 2
        state.bond_damage[:20] = 1.0
        noff, frac = patch._update_diffuse_abort(state, 10.0)
        self.assertGreaterEqual(noff, 10)
        self.assertGreaterEqual(frac, 0.5)
        self.assertIsNone(state.cycles_diffuse_abort)
        patch._update_diffuse_abort(state, 20.0)
        self.assertEqual(state.cycles_diffuse_abort, 20.0)

    def test_checkpoint_roundtrip_preserves_front_and_rng(self):
        mesh, patch = self.make_patch(); state = patch.initial_state()
        ids = self.force_horizontal_crack(patch, state, 1.2)
        patch._capture_or_update_active_front(state, 7.0)
        state.primary_seed_node = 3
        args = SimpleNamespace(alpha=1, beta="x", cycles_max=1e9, max_blocks=100,
                               out="unused", resume=False, skip_existing=False,
                               checkpoint_every_blocks=25, checkpoint_path="",
                               snapshot_every=0, print_every=0)
        with tempfile.TemporaryDirectory() as td:
            cp = Path(td) / "state.npz"
            _save_case_checkpoint(
                cp, args=args, case_name="no_shield", sigma_a_MPa=700.0,
                next_block=4, cycles=12.0, Wp_total=3.0, mesh=mesh,
                root_xy=patch.root_xy, ep_gp=np.zeros((3, mesh.ne)),
                rho_gp=np.ones(mesh.ne), epsp_acc_gp=np.zeros(mesh.ne),
                u=np.zeros(mesh.ndof), last_residual=np.zeros(mesh.nn),
                pd_state=state, patch=patch, rows=[{"block": 3}],
            )
            expected = patch._event_rng.random()
            mesh2, patch2 = self.make_patch()
            restored = _load_case_checkpoint(
                cp, args=args, case_name="no_shield", sigma_a_MPa=700.0,
                mesh=mesh2, patch=patch2,
            )
            s2 = restored["pd_state"]
            self.assertTrue(s2.active_front)
            self.assertTrue(np.array_equal(s2.active_front_bonds, state.active_front_bonds))
            self.assertTrue(np.array_equal(s2.site_node_index, state.site_node_index))
            self.assertTrue(np.allclose(s2.site_birth_threshold, state.site_birth_threshold))
            self.assertTrue(np.allclose(s2.birth_cumulative_hazard, state.birth_cumulative_hazard))
            self.assertEqual(s2.primary_seed_node, 3)
            self.assertAlmostEqual(restored["cycles"], 12.0)
            self.assertAlmostEqual(patch2._event_rng.random(), expected)


    def test_captured_slender_front_can_satisfy_handoff(self):
        _, patch = self.make_patch(); state = patch.initial_state()
        self.force_horizontal_crack(patch, state, 1.8)
        self.assertTrue(patch._capture_or_update_active_front(state, 1.0))
        audit = patch.crack_handoff_audit(state, 1.5e9, 2.5e9)
        self.assertTrue(audit.front_pass)
        self.assertTrue(audit.competition_pass)
        self.assertTrue(audit.handoff_pass, audit.failure_reasons)

    def test_handoff_rejects_large_offfront_damage_cloud(self):
        _, patch = self.make_patch(); state = patch.initial_state()
        self.force_horizontal_crack(patch, state, 1.8)
        self.assertTrue(patch._capture_or_update_active_front(state, 1.0))
        active = np.asarray(state.active_front_bonds, bool)
        off = np.where(~active)[0][:max(2*np.count_nonzero(active), 1)]
        state.bond_damage[off] = 1.0
        audit = patch.crack_handoff_audit(state, 1.5e9, 2.5e9)
        self.assertFalse(audit.competition_pass)
        self.assertFalse(audit.handoff_pass)
        self.assertIn("competition", audit.failure_reasons)

    def test_precapture_neighbor_growth_stays_on_primary_crack_plane(self):
        _, patch = self.make_patch(); state = patch.initial_state()
        available = np.where(state.available_sites > 0)[0]
        node = int(available[np.argmin(np.linalg.norm(
            patch.xy[available] - patch.root_xy[None, :], axis=1
        ))])
        state.stable_sites[node] = 1
        state.growth[node] = 1.0
        state.primary_seed_node = node
        state.crack_normal_c2[node] = -1.0  # normal approximately y
        state.crack_normal_s2[node] = 0.0
        state.crack_orientation_weight[node] = 1.0
        # Seed one admissible near-surface bond so the neighbor term is active.
        activity0, seed0, _ = patch._directional_link_activity(state)
        seed_ids = np.where(seed0 > 0.1 * np.max(seed0))[0]
        self.assertGreater(len(seed_ids), 0)
        state.bond_damage[int(seed_ids[0])] = 1.0
        activity, _, _ = patch._directional_link_activity(state)
        misaligned = np.abs(patch.n[:, 1]) < math.cos(
            math.radians(patch.cfg.front_orientation_tolerance_deg)
        )
        self.assertTrue(np.all(activity[misaligned] == 0.0))

    def test_checkpoint_signature_tracks_physics_not_run_extension(self):
        parser = build_parser()
        a = parser.parse_args([])
        b = parser.parse_args([])
        b.cycles_max = 2.0 * a.cycles_max
        self.assertEqual(
            _checkpoint_signature(a, "shielded", 700.0),
            _checkpoint_signature(b, "shielded", 700.0),
        )
        b.front_band_horizons = a.front_band_horizons * 1.5
        self.assertNotEqual(
            _checkpoint_signature(a, "shielded", 700.0),
            _checkpoint_signature(b, "shielded", 700.0),
        )

    def test_config_validation_rejects_unsafe_active_front_parameters(self):
        cfg = StatefulPDConfig(horizon_m=-1.0)
        with self.assertRaises(ValueError):
            StatefulPDPatch._validate_config(cfg)
        cfg = StatefulPDConfig(handoff_max_offfront_broken_fraction=1.5)
        with self.assertRaises(ValueError):
            StatefulPDPatch._validate_config(cfg)
        cfg = StatefulPDConfig(front_orientation_tolerance_deg=90.0)
        with self.assertRaises(ValueError):
            StatefulPDPatch._validate_config(cfg)

    def test_zero_site_density_is_a_valid_disabled_nucleation_mode(self):
        cfg = StatefulPDConfig(site_density_m2=0.0)
        StatefulPDPatch._validate_config(cfg)

    def test_source_provenance_hashes_are_recordable(self):
        self.assertEqual(set(SOURCE_SHA256), {"driver", "pd_module"})
        for digest in SOURCE_SHA256.values():
            self.assertEqual(len(digest), 64)
            int(digest, 16)


    def test_polyline_projection_returns_curved_arclength(self):
        points = np.array([[0.5, 0.2], [1.2, 0.5]])
        path = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]])
        dist, arc, closest = StatefulPDPatch._project_points_to_polyline(points, path)
        self.assertTrue(np.all(np.isfinite(dist)))
        self.assertTrue(np.all(np.isfinite(arc)))
        self.assertAlmostEqual(arc[0], 0.5, places=12)
        self.assertAlmostEqual(arc[1], 1.5, places=12)
        self.assertTrue(np.allclose(closest[0], [0.5, 0.0]))
        self.assertTrue(np.allclose(closest[1], [1.0, 0.5]))

    def test_active_front_uses_penalized_fallback_when_preferred_set_is_empty(self):
        _, patch = self.make_patch(); state = patch.initial_state()
        self.force_horizontal_crack(patch, state, 1.0)
        self.assertTrue(patch._capture_or_update_active_front(state, 1.0))
        patch.cfg.front_preferred_orientation_tolerance_deg = 0.001
        patch.cfg.front_preferred_band_horizons = 0.001
        patch.cfg.front_fallback_orientation_tolerance_deg = 80.0
        patch.cfg.front_fallback_band_horizons = 1.5
        activity, _, _ = patch._directional_link_activity(state)
        self.assertEqual(state.front_eligible_preferred, 0)
        self.assertGreater(state.front_eligible_fallback, 0)
        self.assertEqual(state.front_candidate_mode, 2)
        self.assertGreater(float(np.max(activity)), 0.0)

    def test_cumulative_wake_is_not_replaced_during_front_update(self):
        _, patch = self.make_patch(); state = patch.initial_state()
        self.force_horizontal_crack(patch, state, 1.0)
        self.assertTrue(patch._capture_or_update_active_front(state, 1.0))
        old = state.front_wake_bonds.copy()
        activity, _, _ = patch._directional_link_activity(state)
        candidates = np.where(activity > 0.0)[0]
        self.assertGreater(len(candidates), 0)
        b = int(candidates[np.argmax(activity[candidates])])
        state.bond_damage[b] = 1.0
        patch._capture_or_update_active_front(state, 2.0)
        self.assertTrue(np.all(state.front_wake_bonds[old]))
        self.assertTrue(state.front_wake_bonds[b])

    def test_mesh_lock_stall_is_explicit_and_patience_gated(self):
        _, patch = self.make_patch(); state = patch.initial_state()
        state.active_front = True
        state.active_front_contact_xy = patch.root_xy.copy()
        state.active_front_tip_xy = patch.root_xy + np.array([patch.cfg.horizon_m, 0.0])
        state.active_front_path_xy = np.vstack([state.active_front_contact_xy, state.active_front_tip_xy])
        state.active_front_length_m = patch.cfg.horizon_m
        state.active_front_normal_c2 = -1.0
        state.active_front_normal_s2 = 0.0
        state.front_candidate_mode = 0
        patch.cfg.front_stall_patience_updates = 2
        patch._capture_or_update_active_front(state, 10.0)
        self.assertIsNone(state.cycles_front_stalled)
        patch._capture_or_update_active_front(state, 20.0)
        self.assertEqual(state.cycles_front_stalled, 20.0)



    def test_primary_seed_reselection_releases_nonproductive_first_seed(self):
        _, patch = self.make_patch(); state = patch.initial_state()
        candidates = np.where(state.candidate_sites > 0)[0]
        self.assertGreaterEqual(len(candidates), 2)
        n0, n1 = int(candidates[0]), int(candidates[1])
        state.stable_sites[n0] = 1
        state.stable_sites[n1] = 2
        state.crack_orientation_weight[[n0, n1]] = 1.0
        state.primary_seed_node = n0
        state.primary_seed_last_progress = 0.0
        patch.cfg.primary_seed_reselection_patience_updates = 1
        smax = np.zeros(len(patch.xy)); smax[n1] = 2.0e9
        patch._update_primary_seed_selection(state, smax, 10.0)
        self.assertTrue(state.primary_seed_rejected[n0])
        self.assertEqual(state.primary_seed_reselections, 1)
        self.assertEqual(state.primary_seed_node, n1)
        self.assertEqual(state.cycles_primary_seed_reselected, 10.0)

    def test_primary_seed_slow_monotone_damage_is_not_partition_stall(self):
        _, patch = self.make_patch(); state = patch.initial_state()
        node = int(np.where(state.candidate_sites > 0)[0][0])
        state.stable_sites[node] = 1
        state.crack_orientation_weight[node] = 1.0
        state.primary_seed_node = node
        state.primary_seed_last_progress = 0.0
        patch.cfg.primary_seed_reselection_patience_updates = 2
        patch.cfg.primary_seed_progress_damage_increment = 0.02
        radius = max(patch.cfg.seed_influence_horizons * patch.cfg.horizon_m,
                     2.0 * patch.point_spacing_m)
        local = np.where(np.linalg.norm(patch.bond_midpoints - patch.xy[node], axis=1) <= radius)[0]
        self.assertGreater(len(local), 0)
        for update in range(8):
            state.bond_damage[int(local[0])] += 1.0e-3
            StatefulPDPatchV87._update_primary_seed_selection(
                patch,
                state, np.zeros(len(patch.xy)), float(update + 1)
            )
        self.assertEqual(state.primary_seed_node, node)
        self.assertEqual(state.primary_seed_reselections, 0)
        self.assertEqual(state.primary_seed_stall_updates, 0)

    def test_primary_seed_reselection_exhaustion_is_explicit(self):
        _, patch = self.make_patch(); state = patch.initial_state()
        candidates = np.where(state.candidate_sites > 0)[0]
        n0 = int(candidates[0])
        state.stable_sites[n0] = 1
        state.crack_orientation_weight[n0] = 1.0
        state.primary_seed_node = n0
        patch.cfg.primary_seed_reselection_patience_updates = 1
        patch.cfg.primary_seed_max_reselections = 1
        patch._update_primary_seed_selection(state, np.zeros(len(patch.xy)), 20.0)
        self.assertEqual(state.cycles_precapture_stalled, 20.0)

    def test_geometry_resolution_audit_and_freeze_taper(self):
        mesh, patch = self.make_patch()
        args = SimpleNamespace(
            geometry_min_valid_radius_spacing=5.0,
            geometry_min_global_area_fraction=0.05,
            geometry_saturation_start_radius_spacing=8.0,
            geometry_limit_mode="freeze",
        )
        audit = _geometry_resolution_audit(
            mesh, patch.feature_surface_global_nodes, patch.point_spacing_m,
            float(np.min(mesh.area_e)), args,
        )
        self.assertTrue(audit["pass"])
        synthetic = dict(audit)
        synthetic["root_radius_over_spacing"] = 6.5
        self.assertAlmostEqual(_geometry_taper_factor(synthetic, args), 0.5)
        args.geometry_limit_mode = "terminate"
        self.assertEqual(_geometry_taper_factor(synthetic, args), 1.0)

    def test_primary_seed_birth_competition_scale_is_bounded(self):
        cfg = StatefulPDConfig(pre_capture_birth_scale_with_primary=0.05)
        StatefulPDPatch._validate_config(cfg)
        with self.assertRaises(ValueError):
            StatefulPDPatch._validate_config(
                StatefulPDConfig(pre_capture_birth_scale_with_primary=-0.1)
            )


    def test_v83_defaults_to_fixed_geometry_and_aligned_birth_clocks(self):
        args = build_parser().parse_args([])
        self.assertFalse(args.enable_geometry_evolution)
        self.assertTrue(args.align_blocks_to_birth_clock)
        self.assertAlmostEqual(args.birth_clock_safety_factor, 1.000001)

    def test_site_thresholds_are_persistent_and_reproducible(self):
        _, p1 = self.make_patch(); s1 = p1.initial_state()
        _, p2 = self.make_patch(); s2 = p2.initial_state()
        self.assertEqual(len(s1.site_node_index), int(np.sum(s1.candidate_sites)))
        self.assertTrue(np.array_equal(s1.site_node_index, s2.site_node_index))
        self.assertTrue(np.allclose(s1.site_birth_threshold, s2.site_birth_threshold))
        self.assertTrue(np.all(s1.site_birth_threshold > 0.0))
        self.assertTrue(np.all(s1.birth_cumulative_hazard == 0.0))

    def test_next_birth_wait_responds_continuously_to_hazard(self):
        _, patch = self.make_patch(); state = patch.initial_state()
        state.site_birth_threshold[:] = 1000.0
        site = 0
        node = int(state.site_node_index[site])
        state.site_birth_threshold[site] = 1.0
        rate_fast = np.zeros(len(patch.xy)); rate_fast[node] = 0.1
        rate_slow = np.zeros(len(patch.xy)); rate_slow[node] = 0.05
        self.assertAlmostEqual(patch.next_birth_wait_cycles(state, rate_fast), 10.0)
        self.assertAlmostEqual(patch.next_birth_wait_cycles(state, rate_slow), 20.0)

    def test_first_passage_cycle_is_block_partition_invariant(self):
        _, p1 = self.make_patch(); s1 = p1.initial_state()
        _, p2 = self.make_patch(); s2 = p2.initial_state()
        for s in (s1, s2):
            s.site_birth_threshold[:] = 1000.0
            s.site_birth_threshold[0] = 0.75
        node = int(s1.site_node_index[0])
        dH = np.zeros(len(p1.xy)); dH[node] = 0.75
        ids1, _, _, cyc1 = p1._advance_discrete_birth_clocks(s1, dH, 0.0, 7.5)
        self.assertEqual(len(ids1), 1)
        dHa = np.zeros(len(p2.xy)); dHa[node] = 0.30
        dHb = np.zeros(len(p2.xy)); dHb[node] = 0.45
        ids_a, _, _, _ = p2._advance_discrete_birth_clocks(s2, dHa, 0.0, 3.0)
        ids_b, _, _, cyc2 = p2._advance_discrete_birth_clocks(s2, dHb, 3.0, 4.5)
        self.assertEqual(len(ids_a), 0)
        self.assertEqual(len(ids_b), 1)
        self.assertAlmostEqual(float(cyc1[0]), 7.5, places=10)
        self.assertAlmostEqual(float(cyc2[0]), 7.5, places=10)

    def test_birth_at_block_end_is_not_stabilized_for_full_block(self):
        _, patch = self.make_patch(); state = patch.initial_state()
        state.site_birth_threshold[:] = 1000.0
        node = int(state.site_node_index[0])
        state.site_birth_threshold[0] = 1.0
        dH = np.zeros(len(patch.xy)); dH[node] = 1.0
        ids, _, _, cycles = patch._advance_discrete_birth_clocks(state, dH, 0.0, 10.0)
        self.assertEqual(len(ids), 1)
        self.assertAlmostEqual(float(cycles[0]), 10.0)
        mu_s = np.zeros(len(patch.xy)); mu_s[node] = 1e6
        mu_h = np.zeros(len(patch.xy))
        stabilized, _, first = patch._advance_discrete_embryo_transitions(
            state, mu_s, mu_h, 0.0, 10.0
        )
        self.assertEqual(int(np.sum(stabilized)), 0)
        self.assertIsNone(first)
        self.assertEqual(state.embryo_sites[node], 1)

    def test_active_front_backstress_reduces_effective_link_traction(self):
        mesh, patch = self.make_patch(); state = patch.initial_state()
        self.force_horizontal_crack(patch, state, 1.2)
        self.assertTrue(patch._capture_or_update_active_front(state, 1.0))
        sigma = np.zeros((8, 3, mesh.nn)); sigma[:, 1, :] = 2.5e9
        delivery = np.zeros((8, mesh.nn))
        zero = np.zeros(mesh.nn)
        base = patch.preview_rates(
            state, type('B', (), {'rate_prefactor':0.0, 'deltaG_eV':staticmethod(lambda s,T: np.ones_like(s))})(),
            sigma, delivery, 300.0, 1000.0, zero, zero, 0.0, zero,
            np.ones(len(patch.xy)), np.ones(len(patch.bonds)),
        )
        back = np.full(mesh.nn, 1.0e9)
        shield = patch.preview_rates(
            state, type('B', (), {'rate_prefactor':0.0, 'deltaG_eV':staticmethod(lambda s,T: np.ones_like(s))})(),
            sigma, delivery, 300.0, 1000.0, zero, back, 0.6, zero,
            np.ones(len(patch.xy)), np.ones(len(patch.bonds)),
        )
        self.assertLess(float(np.max(shield['tn'])), float(np.max(base['tn'])))
        self.assertLess(float(np.max(shield['mu_link'])), float(np.max(base['mu_link'])))

    def test_positive_front_state_shift_reduces_link_rate(self):
        mesh, patch = self.make_patch(); state = patch.initial_state()
        self.force_horizontal_crack(patch, state, 1.2)
        self.assertTrue(patch._capture_or_update_active_front(state, 1.0))
        sigma = np.zeros((8, 3, mesh.nn)); sigma[:, 1, :] = 2.5e9
        delivery = np.zeros((8, mesh.nn)); zero = np.zeros(mesh.nn)
        barrier = type('B', (), {'rate_prefactor':0.0, 'deltaG_eV':staticmethod(lambda s,T: np.ones_like(s))})()
        base = patch.preview_rates(
            state, barrier, sigma, delivery, 300.0, 1000.0, zero, zero, 0.0, zero,
            np.ones(len(patch.xy)), np.ones(len(patch.bonds)),
        )
        shift = np.full(mesh.nn, 0.35)
        shifted = patch.preview_rates(
            state, barrier, sigma, delivery, 300.0, 1000.0, shift, zero, 0.0, zero,
            np.ones(len(patch.xy)), np.ones(len(patch.bonds)),
        )
        self.assertGreater(float(np.max(shifted['link_state_shift_z'])), 0.0)
        self.assertLess(float(np.max(shifted['mu_link'])), float(np.max(base['mu_link'])))

if __name__ == "__main__":
    unittest.main()
