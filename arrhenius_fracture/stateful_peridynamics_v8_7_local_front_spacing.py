"""Stateful local peridynamic initiation patch with local-front spacing.

The patch is deliberately separate from every legacy diffuse-fracture module.
It combines:

* physical candidate-site density and finite-memory plastic-delivery completion separated from cleavage nucleation,
* reversible embryos with competing stabilization and healing,
* stable-defect growth with a stored crack-plane director,
* directional cohesive softening and neighbor-to-neighbor crack-front propagation,
* physical scratch-surface connectivity,
* orientation-aware crack topology and a resolution-aware mechanical handoff audit.

This v8.3 revision adds persistent first-passage clocks for every realized
candidate site, keeps matched shielded/no-shield thresholds identical, and
retains shielding in the active-front linkage traction.  The production
baseline uses a fixed scratch geometry; ALE evolution remains an explicit
sensitivity option in the driver.

The implementation remains a one-way local coupling: the intact global FEM
supplies boundary motion, plastic eigenstrain, residual stress and cyclic stress
history.  Bond loss redistributes deformation inside the peridynamic patch but
is not yet fed back to the global FEM stiffness.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import math
from typing import Optional

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.spatial import cKDTree
from scipy.special import expit, gammainc

_GL_X, _GL_W = np.polynomial.legendre.leggauss(8)

from .config import EV_TO_J, KB, ElasticProperties
from .sn_geometry import BluntNotchGeometry
from .sn_intact_fem import lumped_nodal_area

KBEV = KB / EV_TO_J

_SITE_AVAILABLE = np.uint8(0)
_SITE_EMBRYO = np.uint8(1)
_SITE_STABLE = np.uint8(2)
_SITE_INACTIVE = np.uint8(3)


@dataclass
class StatefulPDConfig:
    patch_radius_m: float = 0.45e-3
    horizon_m: float = 90e-6
    boundary_shell_m: float = 100e-6
    kernel_power: float = 2.0
    residual_bond_stiffness: float = 1e-7

    # Physical initiation support. Candidate sites are excluded from the
    # displacement-coupling shell and smoothly tapered to zero outside the
    # root process zone.
    initiation_radius_m: float = 240e-6
    initiation_taper_m: float = 60e-6
    initiation_back_extent_m: float = 60e-6

    # Candidate-site population and separated delivery/nucleation closure.
    # Plastic-event delivery accumulates in a finite-memory state. Cleavage
    # nucleation remains an independent Arrhenius hazard and is gated by the
    # multi-hit completion Q(K, Lambda_delivery).
    site_density_m2: float = 5.0e10
    delivery_hit_count: float = 2.0
    delivery_memory_s: float = 1.0e-3
    birth_scale: float = 1.0
    # Discrete candidate sites carry persistent Exp(1) thresholds in cumulative
    # per-site birth hazard.  This removes block-end binomial timing and makes
    # paired shielded/no-shield first passage responsive to modest hazard shifts.
    birth_first_passage_enabled: bool = True

    # Reversible embryo transitions.
    nu_stabilize_s: float = 5.0e2
    nu_heal_s: float = 2.0e2
    stabilize_stress_Pa: float = 1.4e9
    stabilize_width_Pa: float = 2.5e8
    stabilize_plastic_gain: float = 2.0
    heal_return_fraction: float = 0.9

    # Stable-object growth and bond linkage.
    nu_grow_s: float = 2.0e-2
    grow_stress_Pa: float = 1.2e9
    grow_width_Pa: float = 2.5e8
    stable_count_scale: float = 2.0
    nu_link_s: float = 8.0e-3
    link_stress_Pa: float = 1.1e9
    link_width_Pa: float = 2.0e8
    # Directional crack growth.  A stabilized embryo stores the local maximum-
    # principal-tension direction as its crack normal.  Only bonds aligned with
    # that normal and lying near the corresponding crack plane are softened.
    link_orientation_power: float = 8.0
    link_orientation_floor: float = 0.0
    directional_band_horizons: float = 0.40
    seed_influence_horizons: float = 1.25
    front_neighbor_spacing_factor: float = 1.75
    front_orientation_tolerance_deg: float = 22.5
    neighbor_link_gain: float = 2.0

    # Active-front capture and competition.  Before capture, only the first
    # realized stable seed is permitted to drive a crack plane.  After a
    # coherent surface-connected component is captured, new births and
    # off-front stable-site growth are suppressed and linkage is restricted to
    # a narrow process zone immediately ahead of the active tip.
    front_capture_enabled: bool = True
    front_capture_min_bonds: int = 6
    front_capture_min_length_horizons: float = 0.50
    front_capture_min_orientation_coherence: float = 0.55
    front_capture_max_surface_gap_horizons: float = 1.0
    front_process_ahead_horizons: float = 1.50
    front_process_behind_horizons: float = 0.35
    front_band_horizons: float = 0.30
    front_wake_band_horizons: float = 0.25
    front_tip_seed_gain: float = 1.0

    # Adaptive active-front continuation.  Preferred candidates preserve a
    # narrow crack path; a penalized fallback avoids mesh locking when an
    # irregular point cloud contains no bond in the exact preferred direction.
    front_preferred_orientation_tolerance_deg: float = 25.0
    front_fallback_orientation_tolerance_deg: float = 42.0
    front_preferred_band_horizons: float = 0.35
    front_fallback_band_horizons: float = 0.70
    front_fallback_activity_scale: float = 0.35
    front_direction_smoothing: float = 0.50
    front_max_turn_deg: float = 18.0
    front_recent_segment_horizons: float = 1.25
    front_crack_tube_horizons: float = 0.80
    front_backbone_band_horizons: float = 0.40
    front_min_advance_spacing_factor: float = 0.25
    front_stall_enabled: bool = True
    front_stall_patience_updates: int = 12

    # Pre-capture seed competition. A stabilized site that fails to generate
    # measurable directional bond progress is released so that a better
    # surface-connected candidate may take over. This prevents a single poor
    # first event from blocking crack formation for the rest of a VHCF run.
    primary_seed_reselection_enabled: bool = True
    primary_seed_reselection_patience_updates: int = 24
    primary_seed_progress_damage_increment: float = 0.02
    primary_seed_max_reselections: int = 64
    primary_seed_surface_score_horizons: float = 2.0
    primary_seed_stress_score_weight: float = 1.0
    primary_seed_population_score_weight: float = 0.35
    primary_seed_damage_score_weight: float = 1.5
    pre_capture_birth_scale_with_primary: float = 0.05

    post_capture_birth_scale: float = 0.0
    post_capture_stabilization_scale: float = 0.0
    off_front_growth_scale: float = 0.0

    # Crack-tip localization of the local PD/FEM redistribution ratio.  The
    # intact FEM already contains the notch concentration; after capture only
    # the near-tip region may receive additional nonlocal amplification.
    front_amplification_ahead_horizons: float = 1.50
    front_amplification_behind_horizons: float = 0.25
    front_amplification_band_horizons: float = 1.00
    front_amplification_cap: float = 3.0

    # Active-front shielding. Back stress is subtracted from the opening
    # traction at every bond.  Once a front is captured, the stored-energy /
    # shielding barrier shift also enters the logistic linkage coordinate in a
    # bounded, dimensionless form.
    front_link_backstress_enabled: bool = True
    front_link_state_shift_weight: float = 1.0
    front_link_state_shift_scale_eV: float = 0.35
    front_link_state_shift_z_clip: float = 3.0

    # Diagnostic fail-safe.  A run is stopped as morphologically invalid if a
    # large fraction of broken bonds remains outside the captured front for
    # several consecutive updates.  This avoids spending the remainder of a
    # long S-N calculation on a diffuse damage cloud.
    diffuse_abort_enabled: bool = True
    diffuse_abort_min_bonds: int = 200
    diffuse_abort_min_offfront_fraction: float = 0.50
    diffuse_abort_patience_updates: int = 5

    # Surface-aware topology and physical handoff.
    surface_connection_horizons: float = 1.0
    topology_neighbor_spacing_factor: float = 1.75
    topology_orientation_tolerance_deg: float = 30.0
    handoff_min_orientation_coherence: float = 0.65
    handoff_min_axial_coverage: float = 0.70
    handoff_max_axial_gap_horizons: float = 0.75
    handoff_min_slenderness: float = 5.0
    handoff_require_active_front: bool = True
    handoff_max_offfront_broken_fraction: float = 0.25

    # Numerical event controls and crack definition.
    max_transition_probability: float = 0.08
    broken_damage: float = 0.95
    root_seed_radius_m: float = 120e-6
    # Legacy length-only criterion retained only for explicit compatibility mode.
    established_extent_m: float = 240e-6

    # Physical crack-handoff criterion. A root-connected high-damage component
    # must be resolved, slender, sharp at its leading tip, and sufficiently far
    # from the displacement-coupling shell.
    handoff_mode: str = "physical"
    handoff_min_length_m: float = 0.0
    handoff_min_length_horizons: float = 2.0
    handoff_max_width_ratio: float = 0.30
    handoff_max_tip_width_horizons: float = 1.0
    handoff_tip_window_horizons: float = 1.0
    handoff_width_root_exclusion_horizons: float = 0.5
    handoff_min_boundary_clearance_horizons: float = 3.0
    handoff_min_connected_bonds: int = 3
    handoff_edge_geometry_factor: float = 1.12
    handoff_min_remote_K_MPam05: float = 0.0
    pd_amplification_cap: float = 4.0
    amplification_damage_scale: float = 0.05
    strain_reference: float = 1e-8
    random_seed: int = 1
    softening_damage: float = 1e-3


@dataclass
class StatefulPDState:
    available: np.ndarray
    embryo: np.ndarray
    stable: np.ndarray
    inactive: np.ndarray
    candidate_sites: np.ndarray
    available_sites: np.ndarray
    embryo_sites: np.ndarray
    stable_sites: np.ndarray
    inactive_sites: np.ndarray
    born_sites_cumulative: np.ndarray
    healed_sites_cumulative: np.ndarray
    # Explicit realized-site ledger for persistent first-passage clocks.
    site_node_index: np.ndarray
    site_status: np.ndarray
    site_birth_threshold: np.ndarray
    site_birth_cycle: np.ndarray
    site_stable_cycle: np.ndarray
    birth_cumulative_hazard: np.ndarray
    delivery_memory: np.ndarray
    completion: np.ndarray
    growth: np.ndarray
    # Axial director representation of the crack normal at each stable site.
    # cos(2 theta), sin(2 theta) avoids a sign ambiguity between n and -n.
    crack_normal_c2: np.ndarray
    crack_normal_s2: np.ndarray
    crack_orientation_weight: np.ndarray
    bond_damage: np.ndarray

    # Active-front state.  The front is represented by a locked axial director,
    # a surface contact point, a leading tip, and a mask of bonds belonging to
    # the selected crack manifold.
    primary_seed_node: int
    primary_seed_rejected: np.ndarray
    primary_seed_stall_updates: int
    primary_seed_reselections: int
    primary_seed_last_progress: float
    primary_seed_selected_cycles: float
    active_front: bool
    active_front_normal_c2: float
    active_front_normal_s2: float
    active_front_contact_xy: np.ndarray
    active_front_tip_xy: np.ndarray
    active_front_length_m: float
    active_front_bonds: np.ndarray
    front_backbone_bonds: np.ndarray
    front_wake_bonds: np.ndarray
    front_process_bonds: np.ndarray
    active_front_path_xy: np.ndarray
    front_candidate_mode: int
    front_eligible_preferred: int
    front_eligible_fallback: int
    front_stall_updates: int
    front_last_advance_cycles: float
    front_max_link_rate_per_cycle: float
    diffuse_bad_updates: int

    healed_cumulative: np.ndarray
    born_cumulative: np.ndarray
    cycles_first_embryo: Optional[float] = None
    cycles_first_stable: Optional[float] = None
    cycles_first_expected_embryo: Optional[float] = None
    cycles_first_expected_stable: Optional[float] = None
    cycles_first_softening: Optional[float] = None
    cycles_root_connected: Optional[float] = None
    cycles_two_horizon_crack: Optional[float] = None
    cycles_front_capture: Optional[float] = None
    cycles_primary_seed_reselected: Optional[float] = None
    cycles_precapture_stalled: Optional[float] = None
    cycles_front_stalled: Optional[float] = None
    cycles_diffuse_abort: Optional[float] = None
    cycles_connected: Optional[float] = None


@dataclass
class PDUpdateDiagnostics:
    max_delivery_memory: float
    max_completion: float
    max_embryo: float
    max_stable: float
    max_growth: float
    max_bond_damage: float
    broken_bonds: int
    connected_extent_m: float
    connected_bonds: int
    expected_embryos: float
    expected_births_cumulative: float
    expected_stable: float
    realized_embryos: int
    realized_births_cumulative: int
    realized_stable: int
    max_effective_stress_Pa: float
    max_pd_amplification: float
    max_rate_per_cycle: float
    max_delivery_rate_s: float
    max_delivery_events_per_cycle: float
    max_nucleation_rate_s: float
    max_nucleation_hazard_per_cycle: float
    max_birth_rate_per_cycle: float
    temporal_transient_cycles: int
    expected_candidate_sites: float
    realized_candidate_sites: int
    crack_centerline_length_m: float
    crack_width_m: float
    crack_tip_width_m: float
    crack_width_ratio: float
    crack_tip_radius_eff_m: float
    crack_orientation_deg: float
    crack_boundary_clearance_m: float
    crack_remote_K_MPam05: float
    crack_local_upper_K_MPam05: float
    crack_remote_KI_MPam05: float
    crack_remote_KII_MPam05: float
    crack_orientation_coherence: float
    crack_axial_coverage: float
    crack_max_axial_gap_m: float
    crack_surface_gap_m: float
    crack_slenderness: float
    handoff_root_connected: bool
    handoff_length_pass: bool
    handoff_slenderness_pass: bool
    handoff_tip_pass: bool
    handoff_boundary_pass: bool
    handoff_K_pass: bool
    handoff_orientation_pass: bool
    handoff_coverage_pass: bool
    handoff_surface_pass: bool
    handoff_front_pass: bool
    handoff_competition_pass: bool
    handoff_pass: bool
    active_front: bool
    active_front_length_m: float
    active_front_bonds: int
    front_backbone_bonds: int
    front_wake_bonds: int
    front_process_bonds: int
    front_candidate_mode: int
    front_eligible_preferred: int
    front_eligible_fallback: int
    front_stall_updates: int
    front_stalled: bool
    front_max_link_rate_per_cycle: float
    primary_seed_node: int
    primary_seed_stall_updates: int
    primary_seed_reselections: int
    precapture_stalled: bool
    off_front_broken_bonds: int
    off_front_broken_fraction: float
    diffuse_abort: bool


@dataclass
class CrackHandoffAudit:
    root_connected: bool = False
    connected_bonds: int = 0
    forward_extent_m: float = 0.0
    centerline_length_m: float = 0.0
    width_m: float = 0.0
    tip_width_m: float = 0.0
    width_ratio: float = float("inf")
    tip_radius_eff_m: float = float("inf")
    orientation_deg: float = 0.0
    boundary_clearance_m: float = 0.0
    remote_K_MPam05: float = 0.0
    local_upper_K_MPam05: float = 0.0
    remote_KI_MPam05: float = 0.0
    remote_KII_MPam05: float = 0.0
    orientation_coherence: float = 0.0
    axial_coverage: float = 0.0
    max_axial_gap_m: float = float("inf")
    surface_gap_m: float = float("inf")
    surface_contact_x_m: float = float("nan")
    surface_contact_y_m: float = float("nan")
    slenderness: float = 0.0
    length_pass: bool = False
    slenderness_pass: bool = False
    tip_pass: bool = False
    boundary_pass: bool = False
    K_pass: bool = False
    orientation_pass: bool = False
    coverage_pass: bool = False
    surface_pass: bool = False
    front_pass: bool = False
    competition_pass: bool = False
    off_front_broken_fraction: float = 0.0
    handoff_pass: bool = False
    failure_reasons: tuple[str, ...] = ()


class StatefulPDPatch:
    @staticmethod
    def _validate_config(cfg: StatefulPDConfig) -> None:
        positive = {
            "patch_radius_m": cfg.patch_radius_m,
            "horizon_m": cfg.horizon_m,
            "delivery_hit_count": cfg.delivery_hit_count,
            "delivery_memory_s": cfg.delivery_memory_s,
            "front_capture_min_bonds": cfg.front_capture_min_bonds,
            "diffuse_abort_patience_updates": cfg.diffuse_abort_patience_updates,
            "front_stall_patience_updates": cfg.front_stall_patience_updates,
            "primary_seed_reselection_patience_updates": cfg.primary_seed_reselection_patience_updates,
            "primary_seed_max_reselections": cfg.primary_seed_max_reselections,
            "handoff_min_connected_bonds": cfg.handoff_min_connected_bonds,
            "front_link_state_shift_scale_eV": cfg.front_link_state_shift_scale_eV,
            "front_link_state_shift_z_clip": cfg.front_link_state_shift_z_clip,
        }
        bad = [name for name, value in positive.items() if not np.isfinite(value) or float(value) <= 0.0]
        if bad:
            raise ValueError("STATEFUL_PD_V8_3 config requires positive finite " + ", ".join(bad))
        nonnegative = {
            "boundary_shell_m": cfg.boundary_shell_m,
            "site_density_m2": cfg.site_density_m2,
            "initiation_radius_m": cfg.initiation_radius_m,
            "initiation_taper_m": cfg.initiation_taper_m,
            "initiation_back_extent_m": cfg.initiation_back_extent_m,
            "post_capture_birth_scale": cfg.post_capture_birth_scale,
            "post_capture_stabilization_scale": cfg.post_capture_stabilization_scale,
            "off_front_growth_scale": cfg.off_front_growth_scale,
            "front_tip_seed_gain": cfg.front_tip_seed_gain,
            "front_preferred_band_horizons": cfg.front_preferred_band_horizons,
            "front_fallback_band_horizons": cfg.front_fallback_band_horizons,
            "front_crack_tube_horizons": cfg.front_crack_tube_horizons,
            "front_backbone_band_horizons": cfg.front_backbone_band_horizons,
            "front_min_advance_spacing_factor": cfg.front_min_advance_spacing_factor,
            "primary_seed_progress_damage_increment": cfg.primary_seed_progress_damage_increment,
            "primary_seed_surface_score_horizons": cfg.primary_seed_surface_score_horizons,
            "primary_seed_stress_score_weight": cfg.primary_seed_stress_score_weight,
            "primary_seed_population_score_weight": cfg.primary_seed_population_score_weight,
            "primary_seed_damage_score_weight": cfg.primary_seed_damage_score_weight,
            "pre_capture_birth_scale_with_primary": cfg.pre_capture_birth_scale_with_primary,
        }
        bad = [name for name, value in nonnegative.items() if not np.isfinite(value) or float(value) < 0.0]
        if bad:
            raise ValueError("STATEFUL_PD_V8_3 config requires nonnegative finite " + ", ".join(bad))
        if cfg.boundary_shell_m >= cfg.patch_radius_m:
            raise ValueError("boundary_shell_m must be smaller than patch_radius_m")
        if not 0.0 < cfg.max_transition_probability < 1.0:
            raise ValueError("max_transition_probability must lie in (0, 1)")
        if not 0.0 < cfg.broken_damage <= 1.0:
            raise ValueError("broken_damage must lie in (0, 1]")
        if not 0.0 <= cfg.softening_damage <= cfg.broken_damage:
            raise ValueError("softening_damage must lie in [0, broken_damage]")
        fractions = {
            "front_capture_min_orientation_coherence": cfg.front_capture_min_orientation_coherence,
            "diffuse_abort_min_offfront_fraction": cfg.diffuse_abort_min_offfront_fraction,
            "handoff_min_orientation_coherence": cfg.handoff_min_orientation_coherence,
            "handoff_min_axial_coverage": cfg.handoff_min_axial_coverage,
            "handoff_max_offfront_broken_fraction": cfg.handoff_max_offfront_broken_fraction,
            "pre_capture_birth_scale_with_primary": cfg.pre_capture_birth_scale_with_primary,
            "post_capture_birth_scale": cfg.post_capture_birth_scale,
            "front_link_state_shift_weight": cfg.front_link_state_shift_weight,
            "post_capture_stabilization_scale": cfg.post_capture_stabilization_scale,
            "off_front_growth_scale": cfg.off_front_growth_scale,
        }
        bad = [name for name, value in fractions.items() if not np.isfinite(value) or not 0.0 <= float(value) <= 1.0]
        if bad:
            raise ValueError("STATEFUL_PD_V8_3 config requires fractions in [0, 1]: " + ", ".join(bad))
        for name, value in {
            "front_orientation_tolerance_deg": cfg.front_orientation_tolerance_deg,
            "front_preferred_orientation_tolerance_deg": cfg.front_preferred_orientation_tolerance_deg,
            "front_fallback_orientation_tolerance_deg": cfg.front_fallback_orientation_tolerance_deg,
            "front_max_turn_deg": cfg.front_max_turn_deg,
            "topology_orientation_tolerance_deg": cfg.topology_orientation_tolerance_deg,
        }.items():
            if not np.isfinite(value) or not 0.0 < float(value) < 90.0:
                raise ValueError(f"{name} must lie in (0, 90) degrees")
        if cfg.front_fallback_orientation_tolerance_deg < cfg.front_preferred_orientation_tolerance_deg:
            raise ValueError("front fallback tolerance must not be narrower than preferred tolerance")
        if not 0.0 <= cfg.front_fallback_activity_scale <= 1.0:
            raise ValueError("front_fallback_activity_scale must lie in [0, 1]")
        if not 0.0 <= cfg.front_direction_smoothing <= 1.0:
            raise ValueError("front_direction_smoothing must lie in [0, 1]")
        if cfg.handoff_mode not in {"physical", "legacy_extent"}:
            raise ValueError("handoff_mode must be 'physical' or 'legacy_extent'")
        if cfg.front_amplification_cap < 1.0 or cfg.pd_amplification_cap < 1.0:
            raise ValueError("PD amplification caps must be at least 1")

    def __init__(
        self,
        mesh,
        geom: BluntNotchGeometry,
        root_xy,
        mat: ElasticProperties,
        cfg: StatefulPDConfig,
        feature_surface_global_nodes=None,
    ):
        self._validate_config(cfg)
        self.cfg = cfg
        self.geom = geom
        self.root_xy = np.asarray(root_xy, float)
        if self.root_xy.shape != (2,) or not np.all(np.isfinite(self.root_xy)):
            raise ValueError("root_xy must be a finite two-component coordinate")
        if not np.all(np.isfinite(np.asarray(mesh.nodes, float))):
            raise ValueError("FEM mesh contains non-finite node coordinates")
        self.patch_center_xy = self.root_xy.copy()
        self.mat = mat

        dist = np.linalg.norm(mesh.nodes - self.root_xy[None, :], axis=1)
        self.global_nodes = np.where(dist <= cfg.patch_radius_m)[0]
        if len(self.global_nodes) < 8:
            raise RuntimeError("peridynamic patch contains too few points")
        self.xy = np.asarray(mesh.nodes[self.global_nodes], float).copy()
        self.area = np.asarray(lumped_nodal_area(mesh)[self.global_nodes], float)
        positive_area = self.area[np.isfinite(self.area) & (self.area > 0.0)]
        if positive_area.size == 0:
            raise RuntimeError("peridynamic patch has no positive finite nodal areas")
        area_floor = float(np.median(positive_area)) * 1.0e-3
        self.area = np.where(np.isfinite(self.area), np.maximum(self.area, area_floor), area_floor)

        # Preserve the physical scratch-surface node set.  Connectivity is
        # evaluated against this evolving polyline, not a nominal root circle.
        feature = np.asarray(
            [] if feature_surface_global_nodes is None else feature_surface_global_nodes,
            dtype=int,
        )
        if np.any(feature < 0) or np.any(feature >= len(mesh.nodes)):
            raise ValueError("feature-surface node index lies outside the FEM mesh")
        self.feature_surface_global_nodes = feature.copy()
        self._global_to_local = {int(g): i for i, g in enumerate(self.global_nodes)}
        self.feature_surface_local_nodes = np.asarray(
            [self._global_to_local[int(g)] for g in feature if int(g) in self._global_to_local],
            dtype=int,
        )
        self.feature_surface_xy = np.asarray(mesh.nodes[feature], float).copy() if len(feature) else self.root_xy[None, :].copy()

        radial = np.linalg.norm(self.xy - self.patch_center_xy[None, :], axis=1)
        shell0 = max(cfg.patch_radius_m - cfg.boundary_shell_m, 0.5 * cfg.patch_radius_m)
        self.shell_radius_m = float(shell0)
        self.boundary = radial >= shell0
        if np.count_nonzero(self.boundary) < 4:
            order = np.argsort(radial)
            self.boundary[order[-max(4, len(order) // 8):]] = True
        self._update_initiation_weight()
        self.remote_sigma_max_Pa = float("nan")
        tree_spacing = cKDTree(self.xy).query(self.xy, k=2)[0][:, 1]
        finite_spacing = tree_spacing[np.isfinite(tree_spacing) & (tree_spacing > 0.0)]
        if finite_spacing.size == 0:
            raise RuntimeError("peridynamic patch has no positive finite point spacing")
        self._set_spacing_metrics(finite_spacing)

        pairs = np.asarray(list(cKDTree(self.xy).query_pairs(cfg.horizon_m)), dtype=int)
        if pairs.size == 0:
            raise RuntimeError("peridynamic horizon produced no bonds")
        pairs = pairs.reshape(-1, 2)
        keep = self._filter_bonds_outside_initial_void(pairs)
        self.bonds = pairs[keep]
        if len(self.bonds) < len(self.xy):
            raise RuntimeError("peridynamic patch is under-connected; increase horizon")

        self._stiffness_scale = None
        self._update_bond_geometry(calibrate=True)
        self.incident = [[] for _ in range(len(self.xy))]
        for b, (i, j) in enumerate(self.bonds):
            self.incident[int(i)].append(b)
            self.incident[int(j)].append(b)
        seed_sequence = np.random.SeedSequence(int(self.cfg.random_seed))
        ss_candidate, ss_event = seed_sequence.spawn(2)
        self._candidate_rng = np.random.default_rng(ss_candidate)
        self._event_rng = np.random.default_rng(ss_event)
        self._rebuild_bond_topology()

    def _update_initiation_weight(self, update_site_measure=True):
        """Build a root-localized physical site measure.

        The coupling shell carries Dirichlet data only and cannot host
        nucleation sites. A cosine taper avoids a mesh-sensitive hard cutoff.
        """
        cfg = self.cfg
        radial = np.linalg.norm(self.xy - self.root_xy[None, :], axis=1)
        r_outer = max(float(cfg.initiation_radius_m), 0.0)
        taper = max(float(cfg.initiation_taper_m), 0.0)
        r_inner = max(r_outer - taper, 0.0)
        w = np.zeros(len(self.xy), dtype=float)
        if r_outer <= 0.0:
            w[:] = 1.0
        elif taper <= 0.0:
            w[radial < r_outer] = 1.0
        else:
            w[radial <= r_inner] = 1.0
            mid = (radial > r_inner) & (radial < r_outer)
            q = (radial[mid] - r_inner) / max(r_outer - r_inner, 1e-30)
            w[mid] = 0.5 * (1.0 + np.cos(np.pi * q))
        x_min = self.root_xy[0] - max(float(cfg.initiation_back_extent_m), 0.0)
        w[self.xy[:, 0] < x_min] = 0.0
        w[self.boundary] = 0.0
        self.initiation_weight = np.clip(w, 0.0, 1.0)
        if update_site_measure or not hasattr(self, "mean_candidate_sites"):
            self.mean_candidate_sites = (
                max(float(cfg.site_density_m2), 0.0) * self.area * self.initiation_weight
            )

    def _filter_bonds_outside_initial_void(self, pairs: np.ndarray) -> np.ndarray:
        p0 = self.xy[pairs[:, 0]]
        p1 = self.xy[pairs[:, 1]]
        keep = np.ones(len(pairs), dtype=bool)
        a = max(float(self.geom.depth_a), 1e-30)
        b = max(float(self.geom.half_height_b), 1e-30)
        for t in (0.25, 0.5, 0.75):
            q = (1.0 - t) * p0 + t * p1
            inside = (q[:, 0] >= -1e-15) & ((q[:, 0] / a) ** 2 + (q[:, 1] / b) ** 2 < 1.0 - 1e-9)
            keep &= ~inside
        return keep

    def _update_bond_geometry(self, calibrate=False):
        i, j = self.bonds[:, 0], self.bonds[:, 1]
        dx = self.xy[j] - self.xy[i]
        L = np.linalg.norm(dx, axis=1)
        if not np.all(np.isfinite(L)) or np.any(L <= 1e-15):
            raise RuntimeError("PD bond geometry contains collapsed or non-finite bonds")
        self.L = L
        self.n = dx / self.L[:, None]
        q = np.clip(self.L / max(self.cfg.horizon_m, 1e-30), 0.0, 1.0)
        kernel = np.exp(-(q ** self.cfg.kernel_power)) * (1.0 - q) ** 2
        raw_w = kernel * self.area[i] * self.area[j] / np.maximum(self.L**2, 1e-30)
        if calibrate or self._stiffness_scale is None:
            eps = 1.0e-6
            uy = eps * (self.xy[:, 1] - np.mean(self.xy[:, 1]))
            utest = np.column_stack([np.zeros(len(self.xy)), uy])
            ext = np.einsum("bi,bi->b", utest[j] - utest[i], self.n)
            raw_energy = 0.5 * np.sum(raw_w * ext**2)
            Dyyyy = self.mat.E * (1.0 - self.mat.nu) / ((1.0 + self.mat.nu) * (1.0 - 2.0 * self.mat.nu))
            target_energy = 0.5 * Dyyyy * eps**2 * np.sum(self.area)
            self._stiffness_scale = target_energy / max(raw_energy, 1e-300)
        self.k0 = self._stiffness_scale * raw_w


    def _rebuild_bond_topology(self):
        """Build geometry-aware bond neighborhoods used by crack propagation.

        A peridynamic crack is represented by a sequence of similarly oriented
        broken bonds whose *midpoints* trace a surface.  Requiring shared end
        nodes creates artificial radial stars, so propagation/connectivity use
        midpoint proximity and axial-orientation coherence instead.
        """
        i, j = self.bonds[:, 0], self.bonds[:, 1]
        self.bond_midpoints = 0.5 * (self.xy[i] + self.xy[j])
        r = max(
            float(self.cfg.front_neighbor_spacing_factor) * self.point_spacing_m,
            0.25 * float(self.cfg.horizon_m),
        )
        pairs = list(cKDTree(self.bond_midpoints).query_pairs(r))
        tol = math.cos(math.radians(float(self.cfg.front_orientation_tolerance_deg)))
        self.bond_neighbors = [[] for _ in range(len(self.bonds))]
        for b, c in pairs:
            if abs(float(np.dot(self.n[b], self.n[c]))) >= tol:
                self.bond_neighbors[b].append(c)
                self.bond_neighbors[c].append(b)

        rt = max(
            float(self.cfg.topology_neighbor_spacing_factor) * self.point_spacing_m,
            0.25 * float(self.cfg.horizon_m),
        )
        topo_pairs = list(cKDTree(self.bond_midpoints).query_pairs(rt))
        topo_tol = math.cos(
            math.radians(float(self.cfg.topology_orientation_tolerance_deg))
        )
        self.topology_neighbors = [[] for _ in range(len(self.bonds))]
        for b, c in topo_pairs:
            if abs(float(np.dot(self.n[b], self.n[c]))) >= topo_tol:
                self.topology_neighbors[b].append(c)
                self.topology_neighbors[c].append(b)

    @staticmethod
    def _distance_points_to_polyline(points, polyline):
        points = np.asarray(points, float)
        polyline = np.asarray(polyline, float)
        if len(polyline) == 0:
            return np.full(len(points), np.inf)
        if len(polyline) == 1:
            return np.linalg.norm(points - polyline[0], axis=1)
        out = np.full(len(points), np.inf)
        for a, b in zip(polyline[:-1], polyline[1:]):
            ab = b - a
            den = max(float(np.dot(ab, ab)), 1e-30)
            q = np.clip(((points - a) @ ab) / den, 0.0, 1.0)
            proj = a[None, :] + q[:, None] * ab[None, :]
            out = np.minimum(out, np.linalg.norm(points - proj, axis=1))
        return out

    @staticmethod
    def _project_points_to_polyline(points, polyline):
        """Return distance, arclength coordinate and closest point.

        The arclength coordinate is measured from the first polyline vertex.
        This supports curved active fronts without reducing them to one global
        straight-line projection.
        """
        points = np.asarray(points, float)
        polyline = np.asarray(polyline, float)
        npt = len(points)
        if len(polyline) == 0:
            return (
                np.full(npt, np.inf),
                np.zeros(npt),
                np.full((npt, 2), np.nan),
            )
        if len(polyline) == 1:
            return (
                np.linalg.norm(points - polyline[0], axis=1),
                np.zeros(npt),
                np.repeat(polyline[:1], npt, axis=0),
            )
        seg = np.diff(polyline, axis=0)
        seglen = np.linalg.norm(seg, axis=1)
        cumulative = np.r_[0.0, np.cumsum(seglen)]
        best_d = np.full(npt, np.inf)
        best_s = np.zeros(npt)
        best_q = np.full((npt, 2), np.nan)
        for k, (a, b) in enumerate(zip(polyline[:-1], polyline[1:])):
            ab = b - a
            den = max(float(np.dot(ab, ab)), 1e-30)
            q = np.clip(((points - a) @ ab) / den, 0.0, 1.0)
            proj = a[None, :] + q[:, None] * ab[None, :]
            d = np.linalg.norm(points - proj, axis=1)
            take = d < best_d
            best_d[take] = d[take]
            best_s[take] = cumulative[k] + q[take] * seglen[k]
            best_q[take] = proj[take]
        return best_d, best_s, best_q

    def _bond_surface_gap(self, bond_ids=None):
        ids = np.arange(len(self.bonds), dtype=int) if bond_ids is None else np.asarray(bond_ids, int)
        if len(ids) == 0:
            return np.empty(0)
        i, j = self.bonds[ids, 0], self.bonds[ids, 1]
        sample = np.vstack([self.xy[i], self.xy[j], self.bond_midpoints[ids]])
        d = self._distance_points_to_polyline(sample, self.feature_surface_xy)
        n = len(ids)
        return np.minimum(d[:n], np.minimum(d[n:2*n], d[2*n:]))

    @staticmethod
    def _principal_normal_from_stress_history(sigma_hist_local):
        """Return the max-principal tensile direction at the peak phase."""
        sig = np.asarray(sigma_hist_local, float)
        sx, sy, txy = sig[:, 0], sig[:, 1], sig[:, 2]
        savg = 0.5 * (sx + sy)
        rad = np.sqrt((0.5 * (sx - sy)) ** 2 + txy**2)
        s1 = savg + rad
        phase = np.argmax(s1, axis=0)
        cols = np.arange(sig.shape[2])
        sxp = sx[phase, cols]
        syp = sy[phase, cols]
        txyp = txy[phase, cols]
        theta = 0.5 * np.arctan2(2.0 * txyp, sxp - syp)
        normal = np.column_stack([np.cos(theta), np.sin(theta)])
        return normal, phase, s1[phase, cols]

    @staticmethod
    def _director_to_normal(c2, s2):
        theta = 0.5 * np.arctan2(s2, c2)
        return np.column_stack([np.cos(theta), np.sin(theta)])

    def _update_crack_directors(self, state, stabilized_sites, principal_normal):
        stabilized_sites = np.asarray(stabilized_sites, float)
        active = stabilized_sites > 0
        if not np.any(active):
            return
        theta = np.arctan2(principal_normal[:, 1], principal_normal[:, 0])
        add_c2 = stabilized_sites * np.cos(2.0 * theta)
        add_s2 = stabilized_sites * np.sin(2.0 * theta)
        old_w = state.crack_orientation_weight
        csum = old_w * state.crack_normal_c2 + add_c2
        ssum = old_w * state.crack_normal_s2 + add_s2
        new_w = old_w + stabilized_sites
        mag = np.hypot(csum, ssum)
        good = new_w > 0
        state.crack_normal_c2[good] = csum[good] / np.maximum(mag[good], 1e-30)
        state.crack_normal_s2[good] = ssum[good] / np.maximum(mag[good], 1e-30)
        state.crack_orientation_weight = new_w

    def _primary_seed_progress_metric(self, state, node):
        """Return directional damage progress local to one pre-capture seed."""
        node = int(node)
        if node < 0 or node >= len(self.xy) or len(self.bonds) == 0:
            return 0.0
        radius = max(
            self.cfg.seed_influence_horizons * self.cfg.horizon_m,
            2.0 * self.point_spacing_m,
        )
        local = np.linalg.norm(self.bond_midpoints - self.xy[node], axis=1) <= radius
        if not np.any(local):
            return 0.0
        return float(np.max(np.asarray(state.bond_damage, float)[local]))

    def _select_primary_seed(self, state, smax, cycles_new):
        """Choose the best viable stabilized site for pre-capture growth.

        The score rewards surface proximity, local tensile driving force,
        stable-site population, and any existing local bond softening. Rejected
        sites are excluded so a poor first event cannot permanently monopolize
        the crack-growth pathway.
        """
        cfg = self.cfg
        stable = np.asarray(state.stable_sites, dtype=np.int64)
        rejected = np.asarray(state.primary_seed_rejected, dtype=bool)
        candidates = np.where(
            (stable > 0)
            & (state.crack_orientation_weight > 0.0)
            & (~rejected)
        )[0]
        if candidates.size == 0:
            state.primary_seed_node = -1
            return -1

        gaps = self._distance_points_to_polyline(
            self.xy[candidates], self.feature_surface_xy
        )
        gap_scale = max(
            cfg.primary_seed_surface_score_horizons * cfg.horizon_m,
            self.point_spacing_m,
        )
        surface_score = np.exp(-np.maximum(gaps, 0.0) / gap_scale)

        smax = np.asarray(smax, float)
        local_stress = np.maximum(smax[candidates], 0.0)
        stress_scale = max(float(np.max(local_stress)), 1e-30)
        stress_score = local_stress / stress_scale

        population = np.log1p(stable[candidates].astype(float))
        pop_scale = max(float(np.max(population)), 1.0)
        population_score = population / pop_scale

        damage_score = np.asarray(
            [self._primary_seed_progress_metric(state, int(k)) for k in candidates],
            float,
        )
        score = (
            surface_score
            + cfg.primary_seed_stress_score_weight * stress_score
            + cfg.primary_seed_population_score_weight * population_score
            + cfg.primary_seed_damage_score_weight * damage_score
        )
        # Deterministic tie-break: prefer the closest surface site, then the
        # lowest local node index for reproducibility across platforms.
        order = np.lexsort((candidates, gaps, -score))
        selected = int(candidates[int(order[0])])
        state.primary_seed_node = selected
        state.primary_seed_stall_updates = 0
        state.primary_seed_last_progress = self._primary_seed_progress_metric(
            state, selected
        )
        state.primary_seed_selected_cycles = float(cycles_new)
        return selected

    def _update_primary_seed_selection(self, state, smax, cycles_new):
        """Release a nonproductive pre-capture seed and select a replacement."""
        cfg = self.cfg
        if state.active_front or state.cycles_precapture_stalled is not None:
            return

        node = int(state.primary_seed_node)
        valid = (
            0 <= node < len(self.xy)
            and state.stable_sites[node] > 0
            and not state.primary_seed_rejected[node]
        )
        if not valid:
            self._select_primary_seed(state, smax, cycles_new)
            return
        if not cfg.primary_seed_reselection_enabled:
            return

        progress = self._primary_seed_progress_metric(state, node)
        threshold = max(float(cfg.primary_seed_progress_damage_increment), 0.0)
        if progress >= float(state.primary_seed_last_progress) + threshold:
            state.primary_seed_last_progress = progress
            state.primary_seed_stall_updates = 0
            return

        state.primary_seed_stall_updates += 1
        if state.primary_seed_stall_updates < int(
            cfg.primary_seed_reselection_patience_updates
        ):
            return

        state.primary_seed_rejected[node] = True
        state.primary_seed_reselections += 1
        state.cycles_primary_seed_reselected = float(cycles_new)
        state.primary_seed_node = -1
        state.primary_seed_stall_updates = 0
        state.primary_seed_last_progress = 0.0

        if state.primary_seed_reselections >= int(cfg.primary_seed_max_reselections):
            state.cycles_precapture_stalled = float(cycles_new)
            return
        selected = self._select_primary_seed(state, smax, cycles_new)
        if selected < 0 and np.any(state.stable_sites > 0):
            # All currently stable sites have been rejected. New sites may still
            # appear in later blocks, so do not invalidate the run immediately.
            state.primary_seed_node = -1

    def _active_front_frame(self, state):
        """Return normal, tangent, contact and tip for the captured front."""
        theta = 0.5 * math.atan2(
            float(state.active_front_normal_s2),
            float(state.active_front_normal_c2),
        )
        normal = np.array([math.cos(theta), math.sin(theta)], dtype=float)
        tangent = np.array([-normal[1], normal[0]], dtype=float)
        if tangent[0] < 0.0:
            tangent *= -1.0
        normal = np.array([-tangent[1], tangent[0]], dtype=float)
        contact = np.asarray(state.active_front_contact_xy, float)
        tip = np.asarray(state.active_front_tip_xy, float)
        return normal, tangent, contact, tip

    def _component_geometry(self, ids):
        """Geometry of one damaged-bond component, robust for small sets."""
        ids = np.asarray(ids, dtype=int)
        if ids.size == 0:
            return None
        mids = self.bond_midpoints[ids]
        weights = np.maximum(np.ones(len(ids)), 1e-12)
        normal, coherence = self._director_average(self.n[ids], weights)
        tangent = np.array([-normal[1], normal[0]], dtype=float)
        if tangent[0] < 0.0:
            tangent *= -1.0
        normal = np.array([-tangent[1], tangent[0]], dtype=float)

        i_id, j_id = self.bonds[ids, 0], self.bonds[ids, 1]
        probes = np.vstack([self.xy[i_id], self.xy[j_id], mids])
        if len(self.feature_surface_xy):
            dist_s, idx_s = cKDTree(self.feature_surface_xy).query(probes, k=1)
            k = int(np.argmin(dist_s))
            contact = self.feature_surface_xy[int(idx_s[k])].copy()
            surface_gap = float(dist_s[k])
        else:
            contact = self.root_xy.copy()
            surface_gap = float(np.min(np.linalg.norm(probes - contact[None, :], axis=1)))
        axial = np.asarray((mids - contact[None, :]) @ tangent, float)
        finite = axial[np.isfinite(axial)]
        length = max(float(np.max(finite)), 0.0) if finite.size else 0.0
        tip = contact + length * tangent
        return {
            "normal": normal,
            "tangent": tangent,
            "contact": contact,
            "tip": tip,
            "length": length,
            "coherence": float(coherence),
            "surface_gap": surface_gap,
        }

    def _components_from_ids(self, ids):
        allowed = set(map(int, np.asarray(ids, dtype=int)))
        unseen = set(allowed)
        components = []
        while unseen:
            seed = unseen.pop()
            comp = {seed}
            stack = [seed]
            while stack:
                b = stack.pop()
                for c in self.topology_neighbors[b]:
                    if c in unseen:
                        unseen.remove(c)
                        comp.add(c)
                        stack.append(c)
            components.append(np.fromiter(comp, dtype=int))
        return components

    def _capture_or_update_active_front(self, state, cycles_new):
        """Capture and advance one persistent, curved active crack front.

        The selected crack wake is cumulative.  The local tangent is updated
        only after measurable tip advance, with smoothing and a bounded turn.
        This avoids both the v7 diffuse fan and the v8 fixed-direction mesh lock.
        """
        cfg = self.cfg
        if not cfg.front_capture_enabled:
            return False

        if not state.active_front:
            ids = self._root_connected_bond_ids(state)
            geom = self._component_geometry(ids)
            if geom is None:
                return False
            surface_tol = max(
                cfg.front_capture_max_surface_gap_horizons * cfg.horizon_m,
                self.point_spacing_m,
            )
            capture = (
                len(ids) >= int(cfg.front_capture_min_bonds)
                and geom["length"] >= cfg.front_capture_min_length_horizons * cfg.horizon_m
                and geom["coherence"] >= cfg.front_capture_min_orientation_coherence
                and geom["surface_gap"] <= surface_tol
            )
            if not capture:
                return False
            theta = math.atan2(geom["normal"][1], geom["normal"][0])
            state.active_front = True
            state.active_front_normal_c2 = math.cos(2.0 * theta)
            state.active_front_normal_s2 = math.sin(2.0 * theta)
            state.active_front_contact_xy = geom["contact"].copy()
            state.active_front_tip_xy = geom["tip"].copy()
            state.active_front_length_m = float(geom["length"])
            state.active_front_bonds[:] = False
            state.front_backbone_bonds[:] = False
            state.front_wake_bonds[:] = False
            state.front_process_bonds[:] = False
            state.active_front_bonds[ids] = True
            state.front_backbone_bonds[ids] = True
            state.front_wake_bonds[ids] = True
            state.active_front_path_xy = np.vstack([geom["contact"], geom["tip"]])
            state.front_stall_updates = 0
            state.front_last_advance_cycles = float(cycles_new)
            if state.cycles_front_capture is None:
                state.cycles_front_capture = float(cycles_new)
            return True

        normal, tangent, _, tip_old = self._active_front_frame(state)
        broken = state.bond_damage >= cfg.broken_damage
        selected = broken & (
            np.asarray(state.front_wake_bonds, bool)
            | np.asarray(state.front_process_bonds, bool)
        )
        if np.any(selected):
            state.front_wake_bonds |= selected

        wake_ids = np.where(np.asarray(state.front_wake_bonds, bool) & broken)[0]
        advanced = False
        min_advance = max(
            cfg.front_min_advance_spacing_factor * self.point_spacing_m,
            1.0e-12,
        )
        if len(wake_ids):
            mids = self.bond_midpoints[wake_ids]
            forward = np.asarray((mids - tip_old[None, :]) @ tangent, float)
            k = int(np.argmax(forward))
            if np.isfinite(forward[k]) and forward[k] >= min_advance:
                tip_candidate = mids[k].copy()
                recent_radius = max(
                    cfg.front_recent_segment_horizons * cfg.horizon_m,
                    2.0 * self.point_spacing_m,
                )
                recent = wake_ids[
                    np.linalg.norm(self.bond_midpoints[wake_ids] - tip_candidate[None, :], axis=1)
                    <= recent_radius
                ]
                if len(recent):
                    n_new, _ = self._director_average(
                        self.n[recent], np.maximum(state.bond_damage[recent], 1e-12)
                    )
                    t_new = np.array([-n_new[1], n_new[0]], dtype=float)
                    if float(t_new @ tangent) < 0.0:
                        t_new *= -1.0
                    angle_old = math.atan2(tangent[1], tangent[0])
                    angle_new = math.atan2(t_new[1], t_new[0])
                    turn = math.atan2(
                        math.sin(angle_new - angle_old),
                        math.cos(angle_new - angle_old),
                    )
                    max_turn = math.radians(cfg.front_max_turn_deg)
                    turn = float(np.clip(turn, -max_turn, max_turn))
                    turn *= float(np.clip(cfg.front_direction_smoothing, 0.0, 1.0))
                    angle = angle_old + turn
                    tangent = np.array([math.cos(angle), math.sin(angle)], dtype=float)
                    normal = np.array([-tangent[1], tangent[0]], dtype=float)

                path = np.asarray(state.active_front_path_xy, float)
                if path.ndim != 2 or path.shape[1] != 2 or len(path) == 0:
                    path = np.asarray([state.active_front_contact_xy, tip_old], float)
                if np.linalg.norm(tip_candidate - path[-1]) >= min_advance:
                    path = np.vstack([path, tip_candidate])
                    state.active_front_path_xy = path
                    state.active_front_tip_xy = tip_candidate
                    state.active_front_length_m = float(
                        np.sum(np.linalg.norm(np.diff(path, axis=0), axis=1))
                    )
                    theta_n = math.atan2(normal[1], normal[0])
                    state.active_front_normal_c2 = math.cos(2.0 * theta_n)
                    state.active_front_normal_s2 = math.sin(2.0 * theta_n)
                    state.front_last_advance_cycles = float(cycles_new)
                    state.front_stall_updates = 0
                    advanced = True

        path = np.asarray(state.active_front_path_xy, float)
        if path.ndim == 2 and len(path):
            dist, _, _ = self._project_points_to_polyline(self.bond_midpoints, path)
            backbone_band = max(
                cfg.front_backbone_band_horizons * cfg.horizon_m,
                self.point_spacing_m,
            )
            state.front_backbone_bonds = (
                np.asarray(state.front_wake_bonds, bool)
                & broken
                & (dist <= backbone_band)
            )
        state.active_front_bonds = np.asarray(state.front_wake_bonds, bool).copy()

        if not advanced:
            if state.front_candidate_mode == 0:
                state.front_stall_updates += 1
            else:
                # Eligible bonds exist; lack of advance can be physically slow
                # and must not be classified as mesh lock.
                state.front_stall_updates = 0
        if (
            cfg.front_stall_enabled
            and state.front_stall_updates >= int(cfg.front_stall_patience_updates)
            and state.cycles_front_stalled is None
        ):
            state.cycles_front_stalled = float(cycles_new)
        return advanced

    def _front_localization(self, state, for_points=False):
        """Near-tip amplification and cumulative crack-wake unloading."""
        cfg = self.cfg
        normal, tangent, _, tip = self._active_front_frame(state)
        pos = self.xy if for_points else self.bond_midpoints
        rel_tip = pos - tip[None, :]
        forward = rel_tip @ tangent
        transverse = np.abs(rel_tip @ normal)
        ahead_max = max(
            cfg.front_amplification_ahead_horizons * cfg.horizon_m,
            self.point_spacing_m,
        )
        behind_max = max(
            cfg.front_amplification_behind_horizons * cfg.horizon_m,
            self.point_spacing_m,
        )
        band = max(
            cfg.front_amplification_band_horizons * cfg.horizon_m,
            self.point_spacing_m,
        )
        longitudinal = np.exp(
            -0.5 * ((forward - 0.20 * cfg.horizon_m) / max(0.55 * ahead_max, 1e-30)) ** 2
        )
        longitudinal[(forward < -behind_max) | (forward > ahead_max)] = 0.0
        transverse_w = np.exp(-0.5 * (transverse / band) ** 2)
        tip_weight = np.clip(longitudinal * transverse_w, 0.0, 1.0)

        path = np.asarray(state.active_front_path_xy, float)
        if path.ndim != 2 or len(path) == 0:
            path = np.asarray([state.active_front_contact_xy, state.active_front_tip_xy], float)
        dist_path = self._distance_points_to_polyline(pos, path)
        wake = (
            (dist_path <= max(cfg.front_crack_tube_horizons * cfg.horizon_m, self.point_spacing_m))
            & (forward <= -0.05 * cfg.horizon_m)
        )
        return tip_weight, wake

    def _directional_link_activity(self, state):
        """Seed one crack, then advance it using adaptive preferred/fallback candidates."""
        cfg = self.cfg
        nb = len(self.bonds)

        if state.active_front:
            normal, tangent, _, tip = self._active_front_frame(state)
            rel = self.bond_midpoints - tip[None, :]
            forward = rel @ tangent
            transverse = np.abs(rel @ normal)
            abs_align = np.abs(self.n @ normal)
            unbroken = state.bond_damage < cfg.broken_damage

            front_support = np.zeros(nb)
            wake_set = set(map(int, np.where(np.asarray(state.front_wake_bonds, bool))[0]))
            for b, neigh in enumerate(self.topology_neighbors):
                vals = [float(state.bond_damage[c]) for c in neigh if c in wake_set]
                if vals:
                    front_support[b] = max(vals)
            near_tip = np.linalg.norm(rel, axis=1) <= max(
                1.25 * cfg.horizon_m, 2.0 * self.point_spacing_m
            )
            connected = (front_support > 0.0) | near_tip
            behind = max(cfg.front_process_behind_horizons * cfg.horizon_m, self.point_spacing_m)
            ahead = max(cfg.front_process_ahead_horizons * cfg.horizon_m, self.point_spacing_m)
            longitudinal = (forward >= -behind) & (forward <= ahead)

            preferred = (
                unbroken
                & connected
                & longitudinal
                & (transverse <= max(cfg.front_preferred_band_horizons * cfg.horizon_m, self.point_spacing_m))
                & (abs_align >= math.cos(math.radians(cfg.front_preferred_orientation_tolerance_deg)))
            )
            fallback = (
                unbroken
                & connected
                & longitudinal
                & (transverse <= max(cfg.front_fallback_band_horizons * cfg.horizon_m, self.point_spacing_m))
                & (abs_align >= math.cos(math.radians(cfg.front_fallback_orientation_tolerance_deg)))
            )
            state.front_eligible_preferred = int(np.count_nonzero(preferred))
            state.front_eligible_fallback = int(np.count_nonzero(fallback))
            if state.front_eligible_preferred > 0:
                selected = preferred
                band = max(cfg.front_preferred_band_horizons * cfg.horizon_m, self.point_spacing_m)
                scale = 1.0
                state.front_candidate_mode = 1
            elif state.front_eligible_fallback > 0:
                selected = fallback
                band = max(cfg.front_fallback_band_horizons * cfg.horizon_m, self.point_spacing_m)
                scale = float(np.clip(cfg.front_fallback_activity_scale, 0.0, 1.0))
                state.front_candidate_mode = 2
            else:
                selected = np.zeros(nb, dtype=bool)
                band = max(cfg.front_fallback_band_horizons * cfg.horizon_m, self.point_spacing_m)
                scale = 0.0
                state.front_candidate_mode = 0
            state.front_process_bonds = np.asarray(selected, bool).copy()

            align = abs_align ** max(cfg.link_orientation_power, 0.0)
            plane = np.exp(-0.5 * (transverse / max(band, 1e-30)) ** 2)
            window = np.exp(
                -0.5 * ((forward - 0.25 * cfg.horizon_m) / max(0.60 * ahead, 1e-30)) ** 2
            )
            tip_seed = cfg.front_tip_seed_gain * np.exp(
                -0.5 * (np.linalg.norm(rel, axis=1) / max(0.70 * cfg.horizon_m, self.point_spacing_m)) ** 2
            )
            drive = np.maximum(tip_seed, cfg.neighbor_link_gain * front_support)
            activity = scale * align * plane * window * np.clip(drive, 0.0, 1.0)
            activity[~selected] = 0.0
            seed_out = np.clip(scale * tip_seed * align * plane, 0.0, 1.0)
            seed_out[~selected] = 0.0
            return np.clip(activity, 0.0, 1.0), seed_out, front_support

        # Pre-capture: retain the validated single-primary-seed directional rule.
        seed_gate = np.zeros(nb)
        directional_eligibility = np.zeros(nb)
        node = int(state.primary_seed_node)
        if node < 0 or node >= len(self.xy) or state.stable_sites[node] <= 0:
            candidates = np.where(
                (state.stable_sites > 0)
                & (state.crack_orientation_weight > 0)
                & (state.growth > 0)
                & (~np.asarray(state.primary_seed_rejected, dtype=bool))
            )[0]
            if len(candidates):
                gaps = self._distance_points_to_polyline(self.xy[candidates], self.feature_surface_xy)
                node = int(candidates[int(np.argmin(gaps))])
        if node >= 0 and node < len(self.xy) and state.stable_sites[node] > 0:
            normal = self._director_to_normal(
                np.asarray([state.crack_normal_c2[node]]),
                np.asarray([state.crack_normal_s2[node]]),
            )[0]
            tangent = np.array([-normal[1], normal[0]], dtype=float)
            if len(self.feature_surface_xy):
                _, surface_idx = cKDTree(self.feature_surface_xy).query(self.xy[node], k=1)
                surface_point = self.feature_surface_xy[int(surface_idx)]
            else:
                surface_point = self.root_xy
            inward = self.xy[node] - surface_point
            if np.linalg.norm(inward) > 1e-30:
                if float(tangent @ inward) < 0.0:
                    tangent *= -1.0
            elif tangent[0] < 0.0:
                tangent *= -1.0
            normal = np.array([-tangent[1], tangent[0]], dtype=float)

            rel = self.bond_midpoints - self.xy[node]
            transverse = np.abs(rel @ normal)
            axial = rel @ tangent
            abs_align = np.abs(self.n @ normal)
            align = abs_align ** max(cfg.link_orientation_power, 0.0)
            hard_align = abs_align >= math.cos(math.radians(cfg.front_orientation_tolerance_deg))
            band = max(cfg.directional_band_horizons * cfg.horizon_m, self.point_spacing_m)
            seed_len = max(cfg.seed_influence_horizons * cfg.horizon_m, self.point_spacing_m)
            plane = np.exp(-0.5 * (transverse / band) ** 2)
            local = np.exp(-0.5 * (axial / seed_len) ** 2)
            forward_window = (
                (axial >= -0.20 * cfg.horizon_m)
                & (axial <= seed_len)
                & (transverse <= max(1.25 * band, self.point_spacing_m))
                & hard_align
            )
            directional_eligibility = align * plane * local
            directional_eligibility[~forward_window] = 0.0
            count_gate = 1.0 - math.exp(
                -float(state.stable_sites[node]) / max(cfg.stable_count_scale, 1e-30)
            )
            seed_gate = count_gate * float(state.growth[node]) * directional_eligibility

        front = np.zeros(nb)
        for b, neigh in enumerate(self.bond_neighbors):
            if neigh:
                front[b] = max(float(state.bond_damage[c]) for c in neigh)
        propagated = (
            np.clip(cfg.neighbor_link_gain * front, 0.0, 1.0)
            * np.clip(directional_eligibility, 0.0, 1.0)
        )
        activity = 1.0 - (1.0 - np.clip(seed_gate, 0.0, 1.0)) * (1.0 - propagated)
        return np.clip(activity, 0.0, 1.0), seed_gate, front

    def _update_diffuse_abort(self, state, cycles_new):
        broken = state.bond_damage >= self.cfg.broken_damage
        nbroken = int(np.count_nonzero(broken))
        if state.active_front:
            allowed = np.asarray(state.front_wake_bonds, bool) & broken
            off = broken & ~allowed
        else:
            off = broken
        noff = int(np.count_nonzero(off))
        frac = noff / max(nbroken, 1)
        bad = (
            self.cfg.diffuse_abort_enabled
            and nbroken >= int(self.cfg.diffuse_abort_min_bonds)
            and frac >= float(self.cfg.diffuse_abort_min_offfront_fraction)
        )
        state.diffuse_bad_updates = state.diffuse_bad_updates + 1 if bad else 0
        if (
            bad
            and state.diffuse_bad_updates >= int(self.cfg.diffuse_abort_patience_updates)
            and state.cycles_diffuse_abort is None
        ):
            state.cycles_diffuse_abort = float(cycles_new)
        return noff, float(frac)

    def update_geometry(self, mesh, root_xy=None):
        self.xy = np.asarray(mesh.nodes[self.global_nodes], float).copy()
        if root_xy is not None:
            self.root_xy = np.asarray(root_xy, float)
        if len(self.feature_surface_global_nodes):
            self.feature_surface_xy = np.asarray(
                mesh.nodes[self.feature_surface_global_nodes], float
            ).copy()
        if not np.all(np.isfinite(self.xy)):
            raise RuntimeError("ALE update produced non-finite PD coordinates")
        spacing = cKDTree(self.xy).query(self.xy, k=2)[0][:, 1]
        finite_spacing = spacing[np.isfinite(spacing) & (spacing > 0.0)]
        if not len(finite_spacing):
            raise RuntimeError("ALE update produced duplicate or invalid PD point spacing")
        self._set_spacing_metrics(finite_spacing)
        self._update_initiation_weight(update_site_measure=False)
        self._update_bond_geometry(calibrate=False)
        self._rebuild_bond_topology()

    def _set_spacing_metrics(self, finite_spacing):
        """Store global and crack-path-local spacing metrics.

        The base implementation exposes only the finite nearest-neighbor
        distances at both assignment sites; the full ``spacing`` array is not
        in scope in one constructor branch.  v8.7 therefore recomputes the
        local nearest-neighbor distances directly from ``self.xy``.  This is
        executed only when the patch is created or restored.
        """
        finite_spacing = np.asarray(finite_spacing, float)
        finite_spacing = finite_spacing[np.isfinite(finite_spacing) & (finite_spacing > 0.0)]
        if finite_spacing.size == 0:
            raise ValueError("STATEFUL_PD_V8_7 could not determine global point spacing")
        global_h = float(np.median(finite_spacing))

        points = np.asarray(self.xy, float)
        if points.ndim != 2 or points.shape[1] != 2 or points.shape[0] < 2:
            local_h = global_h
        else:
            fallback_root = points[int(np.argmin(points[:, 0]))]
            root_xy = np.asarray(getattr(self, "root_xy", fallback_root), float)
            cfg = getattr(self, "cfg", None)
            horizon = float(getattr(cfg, "horizon_m", 4.0 * global_h))
            root_seed = float(getattr(cfg, "root_seed_radius_m", 5.0 * global_h))
            patch_radius = float(getattr(cfg, "patch_radius_m", 10.0 * global_h))
            local_radius = max(3.0 * horizon, 2.0 * root_seed, 0.15 * patch_radius)
            distance = np.linalg.norm(points - root_xy[None, :], axis=1)
            local_points = points[np.isfinite(distance) & (distance <= local_radius)]
            if local_points.shape[0] >= 8:
                try:
                    from scipy.spatial import cKDTree
                    nearest = cKDTree(local_points).query(local_points, k=2)[0][:, 1]
                    nearest = nearest[np.isfinite(nearest) & (nearest > 0.0)]
                    local_h = float(np.median(nearest)) if nearest.size >= 8 else global_h
                except Exception:
                    local_h = global_h
            else:
                local_h = global_h

        self.point_spacing_global_m = global_h
        self.point_spacing_root_m = local_h
        self.point_spacing_front_m = local_h
        self.point_spacing_m = local_h

    def production_preflight_audit(self):
        ratio = self.cfg.horizon_m / max(self.point_spacing_m, 1e-30)
        boundary_sites = float(np.sum(self.mean_candidate_sites[self.boundary]))
        root_offset = float(np.linalg.norm(self.root_xy - self.patch_center_xy))
        active_clearance = (
            self.shell_radius_m - root_offset - float(self.cfg.initiation_radius_m)
        )
        checks = {
            "horizon_spacing": 2.5 <= ratio <= 5.0,
            "surface_nodes": len(self.feature_surface_xy) >= 5,
            "no_boundary_sites": boundary_sites <= 1e-12,
            "active_region_clearance": active_clearance >= 3.0 * self.cfg.horizon_m,
            "bond_connectivity": len(self.bonds) >= len(self.xy),
        }
        return {
            "point_spacing_m": float(self.point_spacing_m),
            "point_spacing_global_m": float(getattr(self, "point_spacing_global_m", self.point_spacing_m)),
            "point_spacing_root_m": float(getattr(self, "point_spacing_root_m", self.point_spacing_m)),
            "point_spacing_front_m": float(getattr(self, "point_spacing_front_m", self.point_spacing_m)),
            "spacing_metric_mode": "local_root_and_refined_path",
            "horizon_m": float(self.cfg.horizon_m),
            "horizon_over_spacing": float(ratio),
            "surface_node_count": int(len(self.feature_surface_xy)),
            "expected_boundary_sites": boundary_sites,
            "active_region_clearance_m": float(active_clearance),
            "checks": checks,
            "pass": bool(all(checks.values())),
            "failure_reasons": [k for k, v in checks.items() if not v],
        }

    def initial_state(self) -> StatefulPDState:
        npnt, nb = len(self.xy), len(self.bonds)
        # A spatial crack graph is a single specimen realization, not an
        # ensemble average.  Draw a fixed physical candidate-site population
        # once and then use binomial state transitions.  Mean-field fractions
        # are retained in parallel for diagnostics and calibration.
        mean_sites = np.maximum(self.mean_candidate_sites, 0.0)
        candidate_sites = self._candidate_rng.poisson(mean_sites).astype(np.int64)
        if np.sum(candidate_sites) == 0 and np.sum(mean_sites) > 0.0:
            candidate_sites[int(np.argmax(mean_sites))] = 1
        site_node_index = np.repeat(np.arange(npnt, dtype=np.int64), candidate_sites)
        site_birth_threshold = self._event_rng.exponential(1.0, size=len(site_node_index))
        return StatefulPDState(
            available=np.ones(npnt),
            embryo=np.zeros(npnt),
            stable=np.zeros(npnt),
            inactive=np.zeros(npnt),
            candidate_sites=candidate_sites.copy(),
            available_sites=candidate_sites.copy(),
            embryo_sites=np.zeros(npnt, dtype=np.int64),
            stable_sites=np.zeros(npnt, dtype=np.int64),
            inactive_sites=np.zeros(npnt, dtype=np.int64),
            born_sites_cumulative=np.zeros(npnt, dtype=np.int64),
            healed_sites_cumulative=np.zeros(npnt, dtype=np.int64),
            site_node_index=site_node_index,
            site_status=np.full(len(site_node_index), _SITE_AVAILABLE, dtype=np.uint8),
            site_birth_threshold=np.asarray(site_birth_threshold, float),
            site_birth_cycle=np.full(len(site_node_index), np.nan),
            site_stable_cycle=np.full(len(site_node_index), np.nan),
            birth_cumulative_hazard=np.zeros(npnt),
            delivery_memory=np.zeros(npnt),
            completion=np.zeros(npnt),
            growth=np.zeros(npnt),
            crack_normal_c2=np.ones(npnt),
            crack_normal_s2=np.zeros(npnt),
            crack_orientation_weight=np.zeros(npnt),
            bond_damage=np.zeros(nb),
            primary_seed_node=-1,
            primary_seed_rejected=np.zeros(npnt, dtype=bool),
            primary_seed_stall_updates=0,
            primary_seed_reselections=0,
            primary_seed_last_progress=0.0,
            primary_seed_selected_cycles=0.0,
            active_front=False,
            active_front_normal_c2=1.0,
            active_front_normal_s2=0.0,
            active_front_contact_xy=self.root_xy.copy(),
            active_front_tip_xy=self.root_xy.copy(),
            active_front_length_m=0.0,
            active_front_bonds=np.zeros(nb, dtype=bool),
            front_backbone_bonds=np.zeros(nb, dtype=bool),
            front_wake_bonds=np.zeros(nb, dtype=bool),
            front_process_bonds=np.zeros(nb, dtype=bool),
            active_front_path_xy=self.root_xy[None, :].copy(),
            front_candidate_mode=0,
            front_eligible_preferred=0,
            front_eligible_fallback=0,
            front_stall_updates=0,
            front_last_advance_cycles=0.0,
            front_max_link_rate_per_cycle=0.0,
            diffuse_bad_updates=0,
            healed_cumulative=np.zeros(npnt),
            born_cumulative=np.zeros(npnt),
        )

    @staticmethod
    def _counts_from_site_ledger(state, npnt):
        nodes = np.asarray(state.site_node_index, dtype=np.int64)
        status = np.asarray(state.site_status, dtype=np.uint8)
        if len(nodes) != len(status):
            raise RuntimeError("site ledger node/status length mismatch")
        out = []
        for code in (_SITE_AVAILABLE, _SITE_EMBRYO, _SITE_STABLE, _SITE_INACTIVE):
            mask = status == code
            out.append(np.bincount(nodes[mask], minlength=npnt).astype(np.int64))
        return tuple(out)

    def _refresh_node_counts_from_site_ledger(self, state):
        a, e, s, i = self._counts_from_site_ledger(state, len(self.xy))
        state.available_sites = a
        state.embryo_sites = e
        state.stable_sites = s
        state.inactive_sites = i

    def _sync_site_ledger(self, state):
        """Keep explicit site identities consistent with aggregate node counts.

        Production states remain synchronized naturally.  The reconstruction
        path exists for deterministic unit tests and for guarded loading of
        manually edited states; it never changes candidate-site totals.
        """
        nodes = np.asarray(state.site_node_index, dtype=np.int64)
        if len(nodes) != int(np.sum(state.candidate_sites)):
            raise RuntimeError("site ledger size does not match candidate-site population")
        current = self._counts_from_site_ledger(state, len(self.xy))
        target = tuple(np.asarray(x, dtype=np.int64) for x in (
            state.available_sites, state.embryo_sites, state.stable_sites, state.inactive_sites
        ))
        if all(np.array_equal(a, b) for a, b in zip(current, target)):
            return
        status = np.empty(len(nodes), dtype=np.uint8)
        for node in range(len(self.xy)):
            ids = np.where(nodes == node)[0]
            counts = [int(x[node]) for x in target]
            if sum(counts) != len(ids):
                raise RuntimeError("aggregate realized-site counts violate node conservation")
            k = 0
            for code, count in zip((_SITE_AVAILABLE, _SITE_EMBRYO, _SITE_STABLE, _SITE_INACTIVE), counts):
                status[ids[k:k+count]] = code
                k += count
        state.site_status = status
        # An available site's threshold must lie strictly ahead of the current
        # cumulative hazard.  This matters only for reconstructed/manual states.
        avail = status == _SITE_AVAILABLE
        h = np.asarray(state.birth_cumulative_hazard, float)[nodes]
        stale = avail & (np.asarray(state.site_birth_threshold, float) <= h)
        if np.any(stale):
            state.site_birth_threshold[stale] = h[stale] + self._event_rng.exponential(1.0, np.count_nonzero(stale))

    def next_birth_wait_cycles(self, state, birth_rate_per_cycle):
        """Return the next realized-site first-passage time under local rates."""
        if not self.cfg.birth_first_passage_enabled:
            return float("inf")
        self._sync_site_ledger(state)
        rate = np.maximum(np.asarray(birth_rate_per_cycle, float), 0.0)
        if rate.shape != (len(self.xy),):
            raise ValueError("birth_rate_per_cycle must have one value per PD point")
        nodes = np.asarray(state.site_node_index, dtype=np.int64)
        avail = np.asarray(state.site_status, dtype=np.uint8) == _SITE_AVAILABLE
        if not np.any(avail):
            return float("inf")
        residual = np.asarray(state.site_birth_threshold, float) - np.asarray(state.birth_cumulative_hazard, float)[nodes]
        site_rate = rate[nodes]
        valid = avail & np.isfinite(residual) & (site_rate > 0.0)
        if not np.any(valid):
            return float("inf")
        waits = np.maximum(residual[valid], 0.0) / site_rate[valid]
        return float(np.min(waits)) if waits.size else float("inf")

    def _advance_discrete_birth_clocks(self, state, hazard_increment, cycles_old, dN):
        """Advance persistent candidate-site clocks and interpolate crossings."""
        self._sync_site_ledger(state)
        dH = np.maximum(np.asarray(hazard_increment, float), 0.0)
        if dH.shape != (len(self.xy),):
            raise ValueError("birth hazard increment must have one value per PD point")
        H0 = np.asarray(state.birth_cumulative_hazard, float).copy()
        H1 = H0 + dH
        nodes = np.asarray(state.site_node_index, dtype=np.int64)
        status = np.asarray(state.site_status, dtype=np.uint8)
        avail = status == _SITE_AVAILABLE
        threshold = np.asarray(state.site_birth_threshold, float)
        tol = 64.0 * np.finfo(float).eps * np.maximum(1.0, np.abs(H1[nodes]))
        crossed = avail & (dH[nodes] > 0.0) & (threshold <= H1[nodes] + tol)
        ids = np.where(crossed)[0]
        fractions = np.empty(0, dtype=float)
        cycles = np.empty(0, dtype=float)
        if len(ids):
            denom = np.maximum(dH[nodes[ids]], 1e-300)
            fractions = np.clip((threshold[ids] - H0[nodes[ids]]) / denom, 0.0, 1.0)
            cycles = float(cycles_old) + fractions * float(dN)
            state.site_status[ids] = _SITE_EMBRYO
            state.site_birth_cycle[ids] = cycles
            counts = np.bincount(nodes[ids], minlength=len(self.xy)).astype(np.int64)
            state.born_sites_cumulative += counts
        else:
            counts = np.zeros(len(self.xy), dtype=np.int64)
        state.birth_cumulative_hazard = H1
        self._refresh_node_counts_from_site_ledger(state)
        return ids, counts, fractions, cycles

    def _advance_discrete_embryo_transitions(self, state, mu_stab, mu_heal, cycles_old, dN):
        """Advance realized embryo stabilization/healing with continuous waits."""
        self._sync_site_ledger(state)
        mu_s = np.maximum(np.asarray(mu_stab, float), 0.0)
        mu_h = np.maximum(np.asarray(mu_heal, float), 0.0)
        nodes = np.asarray(state.site_node_index, dtype=np.int64)
        embryo_ids = np.where(np.asarray(state.site_status, dtype=np.uint8) == _SITE_EMBRYO)[0]
        stabilized_counts = np.zeros(len(self.xy), dtype=np.int64)
        healed_counts = np.zeros(len(self.xy), dtype=np.int64)
        stable_cycles = []
        cycles_new = float(cycles_old) + float(dN)
        for site in embryo_ids:
            node = int(nodes[site])
            rate = float(mu_s[node] + mu_h[node])
            if rate <= 0.0:
                continue
            entered = float(state.site_birth_cycle[site])
            start = max(float(cycles_old), entered if np.isfinite(entered) else float(cycles_old))
            exposure = max(cycles_new - start, 0.0)
            if exposure <= 0.0:
                continue
            wait = float(self._event_rng.exponential(1.0) / rate)
            if wait > exposure:
                continue
            transition_cycle = start + wait
            if self._event_rng.random() < float(mu_s[node] / rate):
                state.site_status[site] = _SITE_STABLE
                state.site_stable_cycle[site] = transition_cycle
                stabilized_counts[node] += 1
                stable_cycles.append(transition_cycle)
            else:
                healed_counts[node] += 1
                state.healed_sites_cumulative[node] += 1
                if self._event_rng.random() < float(np.clip(self.cfg.heal_return_fraction, 0.0, 1.0)):
                    state.site_status[site] = _SITE_AVAILABLE
                    # Restart at the end of this accepted block.  This avoids an
                    # unresolved heal/rebirth cascade inside one block while
                    # preserving the correct memoryless clock thereafter.
                    H = float(state.birth_cumulative_hazard[node])
                    state.site_birth_threshold[site] = H + float(self._event_rng.exponential(1.0))
                    state.site_birth_cycle[site] = np.nan
                    state.site_stable_cycle[site] = np.nan
                else:
                    state.site_status[site] = _SITE_INACTIVE
        self._refresh_node_counts_from_site_ledger(state)
        return stabilized_counts, healed_counts, (min(stable_cycles) if stable_cycles else None)

    def _assemble_spring_system(self, bond_damage, ep_node):
        npnt = len(self.xy)
        ndof = 2 * npnt
        rows, cols, vals = [], [], []
        rhs = np.zeros(ndof)
        i_all, j_all = self.bonds[:, 0], self.bonds[:, 1]
        exx, eyy, gxy = ep_node[0], ep_node[1], ep_node[2]
        for b, (i, j) in enumerate(zip(i_all, j_all)):
            n = self.n[b]
            nn = np.outer(n, n)
            k = self.k0[b] * max(1.0 - float(bond_damage[b]), self.cfg.residual_bond_stiffness)
            dofi = np.array([2 * i, 2 * i + 1])
            dofj = np.array([2 * j, 2 * j + 1])
            for aa in range(2):
                for bb in range(2):
                    v = k * nn[aa, bb]
                    rows.extend([dofi[aa], dofi[aa], dofj[aa], dofj[aa]])
                    cols.extend([dofi[bb], dofj[bb], dofi[bb], dofj[bb]])
                    vals.extend([v, -v, -v, v])
            epn_i = n[0] ** 2 * exx[i] + n[1] ** 2 * eyy[i] + n[0] * n[1] * gxy[i]
            epn_j = n[0] ** 2 * exx[j] + n[1] ** 2 * eyy[j] + n[0] * n[1] * gxy[j]
            e0 = self.L[b] * 0.5 * (epn_i + epn_j)
            f0 = k * e0 * n
            rhs[dofi] -= f0
            rhs[dofj] += f0
        K = sparse.csr_matrix((vals, (rows, cols)), shape=(ndof, ndof))
        diag_reg = max(float(np.max(self.k0)), 1.0) * 1e-12
        K = K + sparse.eye(ndof, format="csr") * diag_reg
        return K, rhs

    def solve_local_mechanics(self, state, fem_u_global, ep_node_global):
        """Equilibrate the damaged nonlocal patch under FEM shell motion."""
        up = np.asarray(fem_u_global, float).reshape(-1, 2)[self.global_nodes]
        ep = np.asarray(ep_node_global, float)[:, self.global_nodes]
        K, rhs = self._assemble_spring_system(state.bond_damage, ep)
        prescribed_nodes = np.where(self.boundary)[0]
        prescribed = np.zeros(2 * len(self.xy), dtype=bool)
        prescribed[2 * prescribed_nodes] = True
        prescribed[2 * prescribed_nodes + 1] = True
        uvec = np.zeros(2 * len(self.xy))
        uvec[2 * prescribed_nodes] = up[prescribed_nodes, 0]
        uvec[2 * prescribed_nodes + 1] = up[prescribed_nodes, 1]
        free = ~prescribed
        rhs_free = rhs[free] - K[np.ix_(free, prescribed)] @ uvec[prescribed]
        if np.any(free):
            solved = spsolve(K[np.ix_(free, free)], rhs_free)
            if not np.all(np.isfinite(solved)):
                raise RuntimeError("local PD equilibrium solve returned non-finite displacement")
            uvec[free] = solved
        u = uvec.reshape(-1, 2)

        i, j = self.bonds[:, 0], self.bonds[:, 1]
        du = u[j] - u[i]
        ext_total = np.einsum("bi,bi->b", du, self.n)
        exx, eyy, gxy = ep[0], ep[1], ep[2]
        epn_i = self.n[:, 0] ** 2 * exx[i] + self.n[:, 1] ** 2 * eyy[i] + self.n[:, 0] * self.n[:, 1] * gxy[i]
        epn_j = self.n[:, 0] ** 2 * exx[j] + self.n[:, 1] ** 2 * eyy[j] + self.n[:, 0] * self.n[:, 1] * gxy[j]
        mech_strain = ext_total / self.L - 0.5 * (epn_i + epn_j)

        ufem = up
        fem_ext = np.einsum("bi,bi->b", ufem[j] - ufem[i], self.n) / self.L - 0.5 * (epn_i + epn_j)
        amp_raw = np.maximum(mech_strain, 0.0) / np.maximum(
            np.maximum(fem_ext, 0.0), self.cfg.strain_reference
        )
        if state.active_front:
            amp_raw = np.clip(
                amp_raw, 0.25, min(self.cfg.pd_amplification_cap, self.cfg.front_amplification_cap)
            )
            tip_weight, wake = self._front_localization(state, for_points=False)
            amp = 1.0 + tip_weight * (amp_raw - 1.0)
            # Behind the leading tip, the local cohesive solve is allowed to
            # unload but never to create a second broad tensile-amplification
            # zone along the crack wake.
            amp[wake] = np.minimum(amp_raw[wake], 1.0)
        else:
            amp_raw = np.clip(amp_raw, 0.25, self.cfg.pd_amplification_cap)
            # The global FEM already resolves the intact scratch concentration.
            # Activate the PD/FEM redistribution ratio only after cohesive damage
            # develops; this removes spurious intact-boundary amplification.
            gate = np.clip(
                np.asarray(state.bond_damage, float)
                / max(float(self.cfg.amplification_damage_scale), 1e-30),
                0.0,
                1.0,
            )
            amp = 1.0 + gate * (amp_raw - 1.0)

        acc = np.zeros(len(self.xy)); w = np.zeros(len(self.xy))
        for b, (ii, jj) in enumerate(self.bonds):
            acc[ii] += amp[b]; acc[jj] += amp[b]
            w[ii] += 1.0; w[jj] += 1.0
        point_amp = np.where(w > 0, acc / np.maximum(w, 1.0), 1.0)
        if state.active_front:
            point_tip_weight, point_wake = self._front_localization(state, for_points=True)
            point_amp = 1.0 + point_tip_weight * (point_amp - 1.0)
            point_amp[point_wake] = np.minimum(point_amp[point_wake], 1.0)
        return u, mech_strain, fem_ext, amp, point_amp

    def _point_drivers(self, sigma_hist_global, point_amp):
        sig = np.asarray(sigma_hist_global, float)[:, :, self.global_nodes]
        sx, sy, txy = sig[:, 0], sig[:, 1], sig[:, 2]
        szz = self.mat.nu * (sx + sy)
        hydro = (sx + sy + szz) / 3.0
        savg = 0.5 * (sx + sy)
        rad = np.sqrt((0.5 * (sx - sy)) ** 2 + txy**2)
        s1 = savg + rad
        s1_eff = np.maximum(s1, 0.0) * point_amp[None, :]
        hydro_eff = np.maximum(hydro, 0.0) * point_amp[None, :]
        return s1_eff, hydro_eff, sig

    def _bond_traction(self, sigma_hist_global, bond_amp, sigma_back_global=None, chi=0.0):
        sig = np.asarray(sigma_hist_global, float)[:, :, self.global_nodes]
        i, j = self.bonds[:, 0], self.bonds[:, 1]
        sb = 0.5 * (sig[:, :, i] + sig[:, :, j])
        nx, ny = self.n[:, 0], self.n[:, 1]
        tn_raw_phase = nx[None, :] ** 2 * sb[:, 0] + ny[None, :] ** 2 * sb[:, 1] + 2.0 * nx[None, :] * ny[None, :] * sb[:, 2]
        tn_raw = np.max(np.maximum(tn_raw_phase, 0.0), axis=0) * bond_amp
        if sigma_back_global is None or not self.cfg.front_link_backstress_enabled:
            return tn_raw, tn_raw.copy()
        back = np.asarray(sigma_back_global, float)[self.global_nodes]
        back_bond = 0.5 * (back[i] + back[j])
        tn_eff_phase = tn_raw_phase - float(chi) * back_bond[None, :]
        tn_eff = np.max(np.maximum(tn_eff_phase, 0.0), axis=0) * bond_amp
        return tn_eff, tn_raw

    @staticmethod
    def _advance_constant_delivery_phase(
        memory0,
        delivery_rate_s,
        nucleation_rate_s,
        dt_s,
        tau_s,
        hit_count,
    ):
        """Advance delivery memory and integrate the gated cleavage hazard.

        The delivery memory obeys

            dLambda/dt = r_delivery(t) - Lambda/tau_delivery,

        while embryo birth is controlled by the independent cleavage hazard

            h_birth(t) = h_cleave(t) Q(K, Lambda(t)).

        Separating these kernels avoids the v4 double use of the cleavage rate,
        which produced an effective high-order dependence on the same barrier.
        """
        memory0 = np.maximum(np.asarray(memory0, float), 0.0)
        delivery_rate_s = np.maximum(np.asarray(delivery_rate_s, float), 0.0)
        nucleation_rate_s = np.maximum(np.asarray(nucleation_rate_s, float), 0.0)
        dt_s = max(float(dt_s), 0.0)
        tau_s = max(float(tau_s), 1e-30)
        if dt_s <= 0.0:
            return memory0.copy(), np.zeros_like(memory0), memory0.copy()

        e_end = np.exp(-dt_s / tau_s)
        equilibrium = delivery_rate_s * tau_s
        memory1 = equilibrium + (memory0 - equilibrium) * e_end

        tq = 0.5 * dt_s * (_GL_X + 1.0)
        eq = np.exp(-tq[:, None] / tau_s)
        memory_q = equilibrium[None, :] + (memory0[None, :] - equilibrium[None, :]) * eq
        completion_q = gammainc(hit_count, np.maximum(memory_q, 0.0))
        birth_int = 0.5 * dt_s * np.sum(
            _GL_W[:, None] * nucleation_rate_s[None, :] * completion_q,
            axis=0,
        )
        return memory1, birth_int, np.maximum(memory0, memory1)

    def _run_delivery_cycle(
        self,
        memory0,
        delivery_rate_phase_s,
        nucleation_rate_phase_s,
        phase_dt_s,
        fraction=1.0,
    ):
        """Run one full or partial representative cycle."""
        memory = np.maximum(np.asarray(memory0, float), 0.0).copy()
        delivery_rate_phase_s = np.maximum(np.asarray(delivery_rate_phase_s, float), 0.0)
        nucleation_rate_phase_s = np.maximum(np.asarray(nucleation_rate_phase_s, float), 0.0)
        if delivery_rate_phase_s.shape != nucleation_rate_phase_s.shape:
            raise ValueError("delivery and nucleation phase arrays must have the same shape")
        fraction = float(np.clip(fraction, 0.0, 1.0))
        remaining = fraction * delivery_rate_phase_s.shape[0] * phase_dt_s
        birth = np.zeros_like(memory)
        max_memory = memory.copy()
        for r_del, r_nuc in zip(delivery_rate_phase_s, nucleation_rate_phase_s):
            if remaining <= 0.0:
                break
            dt = min(phase_dt_s, remaining)
            memory, inc, phase_max = self._advance_constant_delivery_phase(
                memory,
                r_del,
                r_nuc,
                dt,
                self.cfg.delivery_memory_s,
                self.cfg.delivery_hit_count,
            )
            birth += inc
            max_memory = np.maximum(max_memory, phase_max)
            remaining -= dt
        return memory, birth, max_memory

    def _phase_resolved_delivery_nucleation(
        self,
        delivery_rate_phase_s,
        nucleation_rate_phase_s,
        frequency_Hz,
        memory0,
        dN=None,
    ):
        """Resolve separated delivery memory and cleavage nucleation by phase.

        ``delivery_rate_phase_s`` comes from the Arrhenius plastic-event chain.
        ``nucleation_rate_phase_s`` comes from the EXP-floor cleavage barrier
        evaluated with the local FEM scratch stress and state shifts.
        """
        r_del = np.maximum(np.asarray(delivery_rate_phase_s, float), 0.0)
        r_nuc = np.maximum(np.asarray(nucleation_rate_phase_s, float), 0.0)
        if r_del.ndim != 2 or r_nuc.ndim != 2:
            raise ValueError("phase rates must have shape (n_phase, n_point)")
        if r_del.shape != r_nuc.shape:
            raise ValueError("delivery and nucleation phase arrays must have the same shape")

        nphase = max(int(r_del.shape[0]), 1)
        freq = max(float(frequency_Hz), 1e-30)
        period = 1.0 / freq
        phase_dt = period / nphase
        delivery_events_cycle = np.sum(r_del, axis=0) * phase_dt
        nucleation_hazard_cycle = np.sum(r_nuc, axis=0) * phase_dt

        memory0 = np.maximum(np.asarray(memory0, float), 0.0)
        zero = np.zeros_like(memory0)
        memory_zero_end, _, _ = self._run_delivery_cycle(
            zero, r_del, r_nuc, phase_dt
        )
        cycle_decay = float(np.exp(-period / max(self.cfg.delivery_memory_s, 1e-30)))
        memory_periodic = memory_zero_end / max(1.0 - cycle_decay, 1e-300)

        memory_current_end, birth_current, max_current = self._run_delivery_cycle(
            memory0, r_del, r_nuc, phase_dt
        )
        _, birth_periodic, max_periodic = self._run_delivery_cycle(
            memory_periodic, r_del, r_nuc, phase_dt
        )

        result = {
            "delivery_events_cycle": delivery_events_cycle,
            "delivery_rate_peak_s": np.max(r_del, axis=0),
            "nucleation_hazard_cycle": nucleation_hazard_cycle,
            "nucleation_rate_peak_s": np.max(r_nuc, axis=0),
            "memory_periodic_start": memory_periodic,
            "memory_current_end": memory_current_end,
            "birth_hazard_current_cycle": birth_current,
            "birth_hazard_periodic_cycle": birth_periodic,
            "memory_max_current_cycle": max_current,
            "memory_max_periodic_cycle": max_periodic,
            "completion_max_current_cycle": gammainc(
                self.cfg.delivery_hit_count, np.maximum(max_current, 0.0)
            ),
            "completion_max_periodic_cycle": gammainc(
                self.cfg.delivery_hit_count, np.maximum(max_periodic, 0.0)
            ),
        }
        if dN is None:
            return result

        ncycles = max(float(dN), 0.0)
        nfull = int(np.floor(ncycles + 1e-12))
        frac = max(ncycles - nfull, 0.0)
        memory = memory0.copy()
        birth_block = np.zeros_like(memory)
        max_memory = memory.copy()
        transient_used = 0

        max_transient = min(nfull, 512)
        for _ in range(max_transient):
            memory, inc, cyc_max = self._run_delivery_cycle(
                memory, r_del, r_nuc, phase_dt
            )
            birth_block += inc
            max_memory = np.maximum(max_memory, cyc_max)
            transient_used += 1
            if np.max(np.abs(memory - memory_periodic)) <= 1e-11 * (
                1.0 + np.max(memory_periodic)
            ):
                break

        remaining_full = nfull - transient_used
        if remaining_full > 0:
            if np.max(np.abs(memory - memory_periodic)) <= 1e-9 * (
                1.0 + np.max(memory_periodic)
            ):
                birth_block += remaining_full * birth_periodic
            else:
                A = cycle_decay
                m_start = memory
                m_mid = memory_periodic + (m_start - memory_periodic) * (
                    A ** (0.5 * remaining_full)
                )
                m_end = memory_periodic + (m_start - memory_periodic) * (
                    A ** remaining_full
                )
                _, j0, mx0 = self._run_delivery_cycle(
                    m_start, r_del, r_nuc, phase_dt
                )
                _, jm, mxm = self._run_delivery_cycle(
                    m_mid, r_del, r_nuc, phase_dt
                )
                _, j1, mx1 = self._run_delivery_cycle(
                    m_end, r_del, r_nuc, phase_dt
                )
                birth_block += remaining_full * (j0 + 4.0 * jm + j1) / 6.0
                max_memory = np.maximum(
                    max_memory, np.maximum(mx0, np.maximum(mxm, mx1))
                )
            memory = memory_periodic + (memory - memory_periodic) * (
                cycle_decay ** remaining_full
            )

        if frac > 1e-12:
            memory, inc, partial_max = self._run_delivery_cycle(
                memory, r_del, r_nuc, phase_dt, fraction=frac
            )
            birth_block += inc
            max_memory = np.maximum(max_memory, partial_max)

        result.update({
            "memory_end": memory,
            "memory_max_block": max_memory,
            "completion_max_block": gammainc(
                self.cfg.delivery_hit_count, np.maximum(max_memory, 0.0)
            ),
            "birth_hazard_block": birth_block,
            "birth_hazard_mean_per_cycle": np.divide(
                birth_block, max(ncycles, 1e-300)
            ),
            "transient_cycles_used": transient_used,
        })
        return result

    def preview_rates(
        self,
        state,
        crack_barrier,
        sigma_hist_global,
        delivery_rate_phase_global,
        T_K,
        frequency_Hz,
        state_shift_eV_global,
        sigma_back_global,
        chi,
        plastic_state_global=None,
        point_amp=None,
        bond_amp=None,
    ):
        if point_amp is None:
            point_amp = np.ones(len(self.xy))
        if bond_amp is None:
            bond_amp = np.ones(len(self.bonds))
        s1, hydro, sig_local = self._point_drivers(sigma_hist_global, point_amp)
        principal_normal, principal_phase, principal_s1 = (
            self._principal_normal_from_stress_history(sig_local)
        )
        shift = np.asarray(state_shift_eV_global)[self.global_nodes]
        back = np.asarray(sigma_back_global)[self.global_nodes]
        sig_open = np.maximum(s1 - chi * back[None, :], 0.0)
        G = np.maximum(
            crack_barrier.deltaG_eV(sig_open, T_K) + shift[None, :],
            1e-12,
        )
        lam_nucleation = crack_barrier.rate_prefactor * np.exp(
            np.clip(-G / max(KBEV * T_K, 1e-30), -700.0, 0.0)
        )

        delivery_global = np.asarray(delivery_rate_phase_global, float)
        if delivery_global.ndim != 2:
            raise ValueError(
                "delivery_rate_phase_global must have shape (n_phase, n_global_node)"
            )
        delivery_local = np.maximum(delivery_global[:, self.global_nodes], 0.0)
        if delivery_local.shape != lam_nucleation.shape:
            raise ValueError(
                "delivery-rate and cleavage-rate phase histories do not align"
            )

        temporal = self._phase_resolved_delivery_nucleation(
            delivery_local,
            lam_nucleation,
            frequency_Hz,
            state.delivery_memory,
            dN=None,
        )
        delivery_events = temporal["delivery_events_cycle"]
        delivery_rate_s = temporal["delivery_rate_peak_s"]
        nucleation_hazard = temporal["nucleation_hazard_cycle"]
        nucleation_rate_s = temporal["nucleation_rate_peak_s"]
        completion = temporal["completion_max_current_cycle"]
        completion_eq = temporal["completion_max_periodic_cycle"]
        site_weight = self.initiation_weight
        mu_birth = (
            self.cfg.birth_scale
            * temporal["birth_hazard_current_cycle"]
            * site_weight
        )
        mu_birth_eq = (
            self.cfg.birth_scale
            * temporal["birth_hazard_periodic_cycle"]
            * site_weight
        )
        mu_birth_bound = np.maximum(mu_birth, mu_birth_eq)
        if state.active_front:
            scale = float(np.clip(self.cfg.post_capture_birth_scale, 0.0, 1.0))
            mu_birth *= scale
            mu_birth_eq *= scale
            mu_birth_bound *= scale
        elif state.primary_seed_node >= 0:
            scale = float(np.clip(self.cfg.pre_capture_birth_scale_with_primary, 0.0, 1.0))
            mu_birth *= scale
            mu_birth_eq *= scale
            mu_birth_bound *= scale

        smax = np.max(s1, axis=0)
        zstab = (smax - self.cfg.stabilize_stress_Pa) / max(
            self.cfg.stabilize_width_Pa, 1e-30
        )
        if plastic_state_global is not None:
            zstab = zstab + self.cfg.stabilize_plastic_gain * np.asarray(
                plastic_state_global
            )[self.global_nodes]
        mu_stab = (
            self.cfg.nu_stabilize_s
            / max(frequency_Hz, 1e-30)
            * expit(zstab)
        )
        mu_heal = (
            self.cfg.nu_heal_s
            / max(frequency_Hz, 1e-30)
            * expit(-zstab)
        )
        if state.active_front:
            mu_stab *= float(np.clip(self.cfg.post_capture_stabilization_scale, 0.0, 1.0))
        zgrow = (smax - self.cfg.grow_stress_Pa) / max(
            self.cfg.grow_width_Pa, 1e-30
        )
        stable_activity = 1.0 - np.exp(
            -np.maximum(state.stable_sites.astype(float), 0.0)
            / max(self.cfg.stable_count_scale, 1e-30)
        )
        mu_grow = (
            self.cfg.nu_grow_s
            / max(frequency_Hz, 1e-30)
            * expit(zgrow)
            * stable_activity
        )
        if state.active_front:
            mu_grow *= float(np.clip(self.cfg.off_front_growth_scale, 0.0, 1.0))
        elif state.primary_seed_node >= 0:
            keep = np.zeros_like(mu_grow)
            if state.primary_seed_node < len(keep):
                keep[int(state.primary_seed_node)] = 1.0
            mu_grow *= keep

        tn, tn_raw = self._bond_traction(
            sigma_hist_global, bond_amp, sigma_back_global=sigma_back_global, chi=chi
        )
        zlink = (tn - self.cfg.link_stress_Pa) / max(
            self.cfg.link_width_Pa, 1e-30
        )
        i, j = self.bonds[:, 0], self.bonds[:, 1]
        bond_shift_eV = 0.5 * (shift[i] + shift[j])
        link_shift_z = np.zeros_like(zlink)
        if state.active_front and self.cfg.front_link_state_shift_weight > 0.0:
            link_shift_z = (
                self.cfg.front_link_state_shift_weight
                * bond_shift_eV
                / max(self.cfg.front_link_state_shift_scale_eV, 1e-30)
            )
            link_shift_z = np.clip(
                link_shift_z,
                -abs(self.cfg.front_link_state_shift_z_clip),
                abs(self.cfg.front_link_state_shift_z_clip),
            )
            zlink = zlink - link_shift_z
        directional_activity, seed_gate, front_support = (
            self._directional_link_activity(state)
        )
        mu_link = (
            self.cfg.nu_link_s
            / max(frequency_Hz, 1e-30)
            * expit(zlink)
            * directional_activity
        )
        effective_birth = state.available * mu_birth_bound
        effective_embryo_transition = state.embryo * (mu_stab + mu_heal)
        effective_grow = (1.0 - state.growth) * mu_grow
        effective_link = (1.0 - state.bond_damage) * mu_link
        max_rate = float(
            np.max(
                np.r_[
                    effective_birth,
                    effective_embryo_transition,
                    effective_grow,
                    effective_link,
                ]
            )
        )
        return {
            "delivery_rate_s": delivery_rate_s,
            "delivery_events_per_cycle": delivery_events,
            "nucleation_rate_s": nucleation_rate_s,
            "nucleation_hazard_per_cycle": nucleation_hazard,
            "mu_birth": mu_birth,
            "mu_birth_bound": mu_birth_bound,
            "mu_stab": mu_stab,
            "mu_heal": mu_heal,
            "mu_grow": mu_grow,
            "mu_link": mu_link,
            "smax": smax,
            "tn": tn,
            "tn_raw": tn_raw,
            "bond_state_shift_eV": bond_shift_eV,
            "link_state_shift_z": link_shift_z,
            "principal_normal": principal_normal,
            "principal_phase": principal_phase,
            "principal_s1": principal_s1,
            "directional_activity": directional_activity,
            "directional_seed_gate": seed_gate,
            "directional_front_support": front_support,
            "max_rate_per_cycle": max_rate,
            "_delivery_rate_phase_s": delivery_local,
            "_nucleation_rate_phase_s": lam_nucleation,
            "_completion_current": completion,
            "_completion_periodic": completion_eq,
        }

    def update(
        self,
        state: StatefulPDState,
        crack_barrier,
        sigma_hist_global,
        delivery_rate_phase_global,
        T_K,
        frequency_Hz,
        dN,
        cycles_old,
        state_shift_eV_global,
        sigma_back_global,
        chi,
        plastic_state_global,
        point_amp,
        bond_amp,
    ) -> PDUpdateDiagnostics:
        cfg = self.cfg
        self._sync_site_ledger(state)
        expected_births_old = float(np.sum(np.maximum(state.born_cumulative, 0.0) * self.mean_candidate_sites))
        expected_stable_old = float(np.sum(np.maximum(state.stable, 0.0) * self.mean_candidate_sites))
        realized_births_old = int(np.sum(state.born_sites_cumulative))
        realized_stable_old = int(np.sum(state.stable_sites))
        max_damage_old = float(np.max(state.bond_damage))
        rates = self.preview_rates(
            state,
            crack_barrier,
            sigma_hist_global,
            delivery_rate_phase_global,
            T_K,
            frequency_Hz,
            state_shift_eV_global,
            sigma_back_global,
            chi,
            plastic_state_global,
            point_amp,
            bond_amp,
        )
        temporal = self._phase_resolved_delivery_nucleation(
            rates["_delivery_rate_phase_s"],
            rates["_nucleation_rate_phase_s"],
            frequency_Hz,
            state.delivery_memory,
            dN=dN,
        )
        state.delivery_memory = temporal["memory_end"]
        state.completion = temporal["completion_max_block"]

        if state.active_front:
            birth_competition_scale = cfg.post_capture_birth_scale
        elif state.primary_seed_node >= 0:
            birth_competition_scale = cfg.pre_capture_birth_scale_with_primary
        else:
            birth_competition_scale = 1.0
        birth_scale_active = cfg.birth_scale * birth_competition_scale
        birth_hazard = (
            birth_scale_active * temporal['birth_hazard_block'] * self.initiation_weight
        )
        p_birth = 1.0 - np.exp(-np.clip(birth_hazard, 0.0, 700.0))
        mu_birth = birth_scale_active * temporal['birth_hazard_mean_per_cycle'] * self.initiation_weight
        rates['mu_birth'] = mu_birth
        new_embryo = state.available * p_birth
        state.available -= new_embryo
        state.embryo += new_embryo
        state.born_cumulative += new_embryo

        if cfg.birth_first_passage_enabled:
            birth_site_ids, born_sites, birth_fractions, birth_cycles = (
                self._advance_discrete_birth_clocks(
                    state, birth_hazard, cycles_old, dN
                )
            )
        else:
            p_birth_clip = np.clip(p_birth, 0.0, 1.0)
            born_sites = self._event_rng.binomial(state.available_sites, p_birth_clip).astype(np.int64)
            state.available_sites -= born_sites
            state.embryo_sites += born_sites
            state.born_sites_cumulative += born_sites
            birth_site_ids = np.empty(0, dtype=np.int64)
            birth_fractions = np.empty(0, dtype=float)
            birth_cycles = np.empty(0, dtype=float)

        mu_s, mu_h = rates['mu_stab'], rates['mu_heal']
        mu_tot = mu_s + mu_h
        p_any = 1.0 - np.exp(-np.clip(mu_tot * dN, 0.0, 700.0))
        trans = state.embryo * p_any
        stabilized = trans * mu_s / np.maximum(mu_tot, 1e-300)
        healed = trans * mu_h / np.maximum(mu_tot, 1e-300)
        state.embryo -= stabilized + healed
        state.stable += stabilized
        state.available += cfg.heal_return_fraction * healed
        state.inactive += (1.0 - cfg.heal_return_fraction) * healed
        state.healed_cumulative += healed

        if cfg.birth_first_passage_enabled:
            stabilized_sites, healed_sites, first_stable_cycle_block = (
                self._advance_discrete_embryo_transitions(
                    state, mu_s, mu_h, cycles_old, dN
                )
            )
        else:
            p_any_clip = np.clip(p_any, 0.0, 1.0)
            transitioned_sites = self._event_rng.binomial(state.embryo_sites, p_any_clip).astype(np.int64)
            p_stab_cond = np.divide(mu_s, np.maximum(mu_tot, 1e-300))
            stabilized_sites = self._event_rng.binomial(transitioned_sites, np.clip(p_stab_cond, 0.0, 1.0)).astype(np.int64)
            healed_sites = transitioned_sites - stabilized_sites
            returned_sites = self._event_rng.binomial(healed_sites, np.clip(cfg.heal_return_fraction, 0.0, 1.0)).astype(np.int64)
            state.embryo_sites -= transitioned_sites
            state.stable_sites += stabilized_sites
            state.available_sites += returned_sites
            state.inactive_sites += healed_sites - returned_sites
            state.healed_sites_cumulative += healed_sites
            first_stable_cycle_block = None
        self._update_crack_directors(
            state, stabilized_sites, rates["principal_normal"]
        )
        if state.primary_seed_node < 0 and np.any(stabilized_sites > 0):
            self._select_primary_seed(
                state, rates["smax"],
                first_stable_cycle_block if first_stable_cycle_block is not None else cycles_old + dN
            )

        realized_total = state.available_sites + state.embryo_sites + state.stable_sites + state.inactive_sites
        if np.any(realized_total != state.candidate_sites):
            raise RuntimeError('realized site-population conservation failure')

        state.available = np.maximum(state.available, 0.0)
        state.embryo = np.maximum(state.embryo, 0.0)
        state.stable = np.maximum(state.stable, 0.0)
        state.inactive = np.maximum(state.inactive, 0.0)
        occupied = state.available + state.embryo + state.stable + state.inactive
        over = occupied > 1.0
        if np.any(over):
            scale = 1.0 / occupied[over]
            state.available[over] *= scale
            state.embryo[over] *= scale
            state.stable[over] *= scale
            state.inactive[over] *= scale

        p_grow = 1.0 - np.exp(-np.clip(rates['mu_grow'] * dN, 0.0, 700.0))
        state.growth += (1.0 - state.growth) * p_grow
        state.growth = np.clip(state.growth, 0.0, 1.0)

        mu_link = rates['mu_link']
        if state.active_front and np.any(state.front_process_bonds):
            state.front_max_link_rate_per_cycle = float(
                np.max(mu_link[np.asarray(state.front_process_bonds, bool)])
            )
        else:
            state.front_max_link_rate_per_cycle = float(np.max(mu_link)) if len(mu_link) else 0.0
        p_link = 1.0 - np.exp(-np.clip(mu_link * dN, 0.0, 700.0))
        state.bond_damage += (1.0 - state.bond_damage) * p_link
        state.bond_damage = np.clip(state.bond_damage, 0.0, 1.0)

        expected_births = float(np.sum(np.maximum(state.born_cumulative, 0.0) * self.mean_candidate_sites))
        expected_embryo = float(np.sum(np.maximum(state.embryo, 0.0) * self.mean_candidate_sites))
        expected_stable = float(np.sum(np.maximum(state.stable, 0.0) * self.mean_candidate_sites))
        realized_births = int(np.sum(state.born_sites_cumulative))
        realized_embryo = int(np.sum(state.embryo_sites))
        realized_stable = int(np.sum(state.stable_sites))
        cycles_new = cycles_old + dN

        if state.cycles_first_expected_embryo is None and expected_births >= 1.0:
            frac = (1.0 - expected_births_old) / max(expected_births - expected_births_old, 1e-300)
            state.cycles_first_expected_embryo = cycles_old + float(np.clip(frac, 0.0, 1.0)) * dN
        if state.cycles_first_expected_stable is None and expected_stable >= 1.0:
            frac = (1.0 - expected_stable_old) / max(expected_stable - expected_stable_old, 1e-300)
            state.cycles_first_expected_stable = cycles_old + float(np.clip(frac, 0.0, 1.0)) * dN
        if state.cycles_first_embryo is None and realized_births_old < 1 <= realized_births:
            if len(birth_cycles):
                state.cycles_first_embryo = float(np.min(birth_cycles))
            else:
                state.cycles_first_embryo = cycles_new
        if state.cycles_first_stable is None and realized_stable_old < 1 <= realized_stable:
            state.cycles_first_stable = (
                float(first_stable_cycle_block)
                if first_stable_cycle_block is not None else cycles_new
            )
        soft_thr = float(np.clip(cfg.softening_damage, 0.0, 1.0))
        if state.cycles_first_softening is None and np.max(state.bond_damage) >= soft_thr:
            state.cycles_first_softening = cycles_new

        self._capture_or_update_active_front(state, cycles_new)
        if not state.active_front:
            self._update_primary_seed_selection(state, rates["smax"], cycles_new)
        off_front_broken, off_front_fraction = self._update_diffuse_abort(
            state, cycles_new
        )

        audit = self.crack_handoff_audit(
            state,
            remote_sigma_max_Pa=self.remote_sigma_max_Pa,
            local_effective_stress_Pa=float(np.max(rates["smax"])),
        )
        extent, nconn = audit.forward_extent_m, audit.connected_bonds
        if state.cycles_root_connected is None and audit.root_connected:
            state.cycles_root_connected = cycles_new
        if (
            state.cycles_two_horizon_crack is None
            and audit.root_connected
            and audit.centerline_length_m >= 2.0 * cfg.horizon_m
        ):
            state.cycles_two_horizon_crack = cycles_new
        if state.cycles_connected is None and audit.handoff_pass:
            state.cycles_connected = cycles_new

        self.last_rates = {
            k: np.asarray(v).copy()
            for k, v in rates.items()
            if isinstance(v, np.ndarray) and not k.startswith('_')
        }
        self.last_point_amp = np.asarray(point_amp, float).copy()
        self.last_bond_amp = np.asarray(bond_amp, float).copy()

        return PDUpdateDiagnostics(
            max_delivery_memory=float(np.max(temporal["memory_max_block"])),
            max_completion=float(np.max(temporal['completion_max_block'])),
            max_embryo=float(np.max(state.embryo)),
            max_stable=float(np.max(state.stable)),
            max_growth=float(np.max(state.growth)),
            max_bond_damage=float(np.max(state.bond_damage)),
            broken_bonds=int(np.count_nonzero(state.bond_damage >= cfg.broken_damage)),
            connected_extent_m=float(extent),
            connected_bonds=int(nconn),
            expected_embryos=expected_embryo,
            expected_births_cumulative=expected_births,
            expected_stable=expected_stable,
            realized_embryos=realized_embryo,
            realized_births_cumulative=realized_births,
            realized_stable=realized_stable,
            max_effective_stress_Pa=float(np.max(rates['smax'])),
            max_pd_amplification=float(max(np.max(point_amp), np.max(bond_amp))),
            max_rate_per_cycle=float(max(rates['max_rate_per_cycle'], np.max(mu_birth), np.max(mu_link))),
            max_delivery_rate_s=float(np.max(temporal["delivery_rate_peak_s"])),
            max_delivery_events_per_cycle=float(np.max(temporal["delivery_events_cycle"])),
            max_nucleation_rate_s=float(np.max(temporal["nucleation_rate_peak_s"])),
            max_nucleation_hazard_per_cycle=float(np.max(temporal["nucleation_hazard_cycle"])),
            max_birth_rate_per_cycle=float(np.max(mu_birth)),
            temporal_transient_cycles=int(temporal['transient_cycles_used']),
            expected_candidate_sites=float(np.sum(self.mean_candidate_sites)),
            realized_candidate_sites=int(np.sum(state.candidate_sites)),
            crack_centerline_length_m=float(audit.centerline_length_m),
            crack_width_m=float(audit.width_m),
            crack_tip_width_m=float(audit.tip_width_m),
            crack_width_ratio=float(audit.width_ratio),
            crack_tip_radius_eff_m=float(audit.tip_radius_eff_m),
            crack_orientation_deg=float(audit.orientation_deg),
            crack_boundary_clearance_m=float(audit.boundary_clearance_m),
            crack_remote_K_MPam05=float(audit.remote_K_MPam05),
            crack_local_upper_K_MPam05=float(audit.local_upper_K_MPam05),
            crack_remote_KI_MPam05=float(audit.remote_KI_MPam05),
            crack_remote_KII_MPam05=float(audit.remote_KII_MPam05),
            crack_orientation_coherence=float(audit.orientation_coherence),
            crack_axial_coverage=float(audit.axial_coverage),
            crack_max_axial_gap_m=float(audit.max_axial_gap_m),
            crack_surface_gap_m=float(audit.surface_gap_m),
            crack_slenderness=float(audit.slenderness),
            handoff_root_connected=bool(audit.root_connected),
            handoff_length_pass=bool(audit.length_pass),
            handoff_slenderness_pass=bool(audit.slenderness_pass),
            handoff_tip_pass=bool(audit.tip_pass),
            handoff_boundary_pass=bool(audit.boundary_pass),
            handoff_K_pass=bool(audit.K_pass),
            handoff_orientation_pass=bool(audit.orientation_pass),
            handoff_coverage_pass=bool(audit.coverage_pass),
            handoff_surface_pass=bool(audit.surface_pass),
            handoff_front_pass=bool(audit.front_pass),
            handoff_competition_pass=bool(audit.competition_pass),
            handoff_pass=bool(audit.handoff_pass),
            active_front=bool(state.active_front),
            active_front_length_m=float(state.active_front_length_m),
            active_front_bonds=int(np.count_nonzero(state.active_front_bonds)),
            front_backbone_bonds=int(np.count_nonzero(state.front_backbone_bonds)),
            front_wake_bonds=int(np.count_nonzero(state.front_wake_bonds)),
            front_process_bonds=int(np.count_nonzero(state.front_process_bonds)),
            front_candidate_mode=int(state.front_candidate_mode),
            front_eligible_preferred=int(state.front_eligible_preferred),
            front_eligible_fallback=int(state.front_eligible_fallback),
            front_stall_updates=int(state.front_stall_updates),
            front_stalled=bool(state.cycles_front_stalled is not None),
            front_max_link_rate_per_cycle=float(state.front_max_link_rate_per_cycle),
            primary_seed_node=int(state.primary_seed_node),
            primary_seed_stall_updates=int(state.primary_seed_stall_updates),
            primary_seed_reselections=int(state.primary_seed_reselections),
            precapture_stalled=bool(state.cycles_precapture_stalled is not None),
            off_front_broken_bonds=int(off_front_broken),
            off_front_broken_fraction=float(off_front_fraction),
            diffuse_abort=bool(state.cycles_diffuse_abort is not None),
        )

    def _active_bond_components(self, state: StatefulPDState):
        active = np.where(state.bond_damage >= self.cfg.broken_damage)[0]
        active_set = set(map(int, active))
        components = []
        unseen = set(active_set)
        while unseen:
            seed = unseen.pop()
            comp = {seed}
            stack = [seed]
            while stack:
                b = stack.pop()
                for c in self.topology_neighbors[b]:
                    if c in unseen:
                        unseen.remove(c)
                        comp.add(c)
                        stack.append(c)
            components.append(np.fromiter(comp, dtype=int))
        return components

    def _root_connected_bond_ids(self, state: StatefulPDState):
        """Return the selected active crack, or the best pre-capture component.

        After front capture, handoff diagnostics intentionally ignore unrelated
        broken-bond islands.  This prevents diffuse secondary damage from
        redefining the crack geometry or orientation.
        """
        if state.active_front:
            ids = np.where(
                np.asarray(state.front_wake_bonds, bool)
                & (state.bond_damage >= self.cfg.broken_damage)
            )[0]
            if len(ids):
                return ids
        components = self._active_bond_components(state)
        if not components:
            return np.empty(0, dtype=int)
        tol = max(
            self.cfg.surface_connection_horizons * self.cfg.horizon_m,
            self.point_spacing_m,
        )
        candidates = []
        for comp in components:
            gaps = np.asarray(self._bond_surface_gap(comp), float)
            gaps = gaps[np.isfinite(gaps)]
            gap = float(np.min(gaps)) if gaps.size else float("inf")
            if gap <= tol:
                mids = self.bond_midpoints[comp]
                forward = max(float(np.max(mids[:, 0]) - self.root_xy[0]), 0.0)
                candidates.append((forward, len(comp), -gap, comp))
        if not candidates:
            return np.empty(0, dtype=int)
        candidates.sort(key=lambda q: (q[0], q[1], q[2]), reverse=True)
        return candidates[0][3]

    @staticmethod
    def _robust_span(values):
        values = np.asarray(values, float).reshape(-1)
        values = values[np.isfinite(values)]
        if values.size <= 1:
            return 0.0
        return float(np.quantile(values, 0.95) - np.quantile(values, 0.05))

    @staticmethod
    def _director_average(vectors, weights=None):
        vectors = np.asarray(vectors, float)
        if vectors.size == 0:
            return np.array([1.0, 0.0]), 0.0
        vectors = np.atleast_2d(vectors)
        if vectors.shape[1] != 2:
            raise ValueError("director vectors must have shape (n, 2)")
        finite = np.all(np.isfinite(vectors), axis=1)
        vectors = vectors[finite]
        if vectors.size == 0:
            return np.array([1.0, 0.0]), 0.0
        theta = np.arctan2(vectors[:, 1], vectors[:, 0])
        if weights is None:
            w = np.ones(len(theta), dtype=float)
        else:
            w_all = np.asarray(weights, float).reshape(-1)
            if len(w_all) != len(finite):
                raise ValueError("director weights must match vector count")
            w = w_all[finite]
            w = np.where(np.isfinite(w), w, 0.0)
        c = float(np.sum(w * np.cos(2.0 * theta)))
        q = float(np.sum(w * np.sin(2.0 * theta)))
        den = float(np.sum(np.abs(w)))
        if den <= 1e-30:
            return np.array([1.0, 0.0]), 0.0
        coherence = float(np.clip(math.hypot(c, q) / den, 0.0, 1.0))
        angle = 0.5 * math.atan2(q, c)
        return np.array([math.cos(angle), math.sin(angle)]), coherence

    @staticmethod
    def _max_unique_gap(values) -> float:
        """Largest spacing between finite unique anchors, safe for 0/1 values.

        A single broken bond can project both endpoints and its midpoint onto the
        same axial coordinate.  In that case ``np.diff(np.unique(values))`` is
        empty.  The topology audit should reject the component through its
        length/bond-count criteria, not terminate the simulation.
        """
        values = np.asarray(values, float).reshape(-1)
        values = np.unique(values[np.isfinite(values)])
        if values.size < 2:
            return 0.0
        gaps = np.diff(values)
        return float(np.max(gaps)) if gaps.size else 0.0

    def crack_handoff_audit(
        self,
        state: StatefulPDState,
        remote_sigma_max_Pa=None,
        local_effective_stress_Pa=None,
    ) -> CrackHandoffAudit:
        cfg = self.cfg
        ids = self._root_connected_bond_ids(state)
        if len(ids) == 0:
            active = np.where(state.bond_damage >= cfg.broken_damage)[0]
            if len(active):
                gaps = np.asarray(self._bond_surface_gap(active), float)
                gaps = gaps[np.isfinite(gaps)]
                gap = float(np.min(gaps)) if gaps.size else float("inf")
            else:
                gap = float("inf")
            return CrackHandoffAudit(
                surface_gap_m=gap,
                failure_reasons=("no_surface_connected_component",),
            )

        mids = self.bond_midpoints[ids]
        damage_w = np.maximum(state.bond_damage[ids], 1e-12)
        normal_dir, coherence = self._director_average(self.n[ids], damage_w)

        if state.active_front:
            path = np.asarray(state.active_front_path_xy, float)
            if path.ndim != 2 or path.shape[1] != 2 or len(path) < 2:
                path = np.asarray(
                    [state.active_front_contact_xy, state.active_front_tip_xy], float
                )
            seg = np.diff(path, axis=0)
            seglen = np.linalg.norm(seg, axis=1)
            good = np.where(seglen > 1e-15)[0]
            if len(good):
                k = int(good[-1])
                tangent = seg[k] / seglen[k]
            else:
                tangent = np.array([-normal_dir[1], normal_dir[0]], float)
            if tangent[0] < 0.0:
                tangent *= -1.0
            normal = np.array([-tangent[1], tangent[0]], float)
            contact_xy = path[0].copy()
            tip_xy = path[-1].copy()
            length = float(np.sum(seglen))
        else:
            tangent = np.array([-normal_dir[1], normal_dir[0]], dtype=float)
            if tangent[0] < 0.0:
                tangent *= -1.0
            normal = np.array([-tangent[1], tangent[0]], dtype=float)
            i_id, j_id = self.bonds[ids, 0], self.bonds[ids, 1]
            probes = np.vstack([self.xy[i_id], self.xy[j_id], mids])
            if len(self.feature_surface_xy):
                dist_s, idx_s = cKDTree(self.feature_surface_xy).query(probes, k=1)
                contact_xy = self.feature_surface_xy[int(idx_s[int(np.argmin(dist_s))])]
            else:
                contact_xy = self.root_xy.copy()
            axial0 = np.asarray((mids - contact_xy[None, :]) @ tangent, float)
            finite = axial0[np.isfinite(axial0)]
            length = max(float(np.max(finite)), 0.0) if finite.size else 0.0
            tip_xy = contact_xy + length * tangent
            path = np.vstack([contact_xy, tip_xy])

        dist_path, axial, closest_path = self._project_points_to_polyline(mids, path)
        signed_path = np.asarray(np.einsum("ij,j->i", mids - closest_path, normal), float)
        finite_x = mids[:, 0][np.isfinite(mids[:, 0])]
        forward_extent = (
            max(float(np.max(finite_x) - contact_xy[0]), 0.0)
            if finite_x.size else 0.0
        )

        root_excl = cfg.handoff_width_root_exclusion_horizons * cfg.horizon_m
        body = axial >= min(root_excl, 0.5 * length)
        body_signed = signed_path[body] if np.any(body) else signed_path
        width = self._robust_span(body_signed)
        tip_window = max(
            cfg.handoff_tip_window_horizons * cfg.horizon_m,
            self.point_spacing_m,
        )
        tip_sel = axial >= max(length - tip_window, 0.0)
        tip_trans = np.asarray((mids[tip_sel] - tip_xy[None, :]) @ normal, float)
        tip_width = self._robust_span(tip_trans) if tip_trans.size else width
        width_ratio = width / max(length, 1e-30)
        rho_eff = max(0.5 * tip_width, 0.5 * self.point_spacing_m)
        slenderness = length / max(rho_eff, 1e-30)

        axial_pos = np.sort(axial[np.isfinite(axial) & (axial >= -0.25 * self.point_spacing_m)])
        bin_w = max(self.point_spacing_m, 0.25 * cfg.horizon_m)
        nbin = max(int(math.ceil(length / max(bin_w, 1e-30))), 1)
        occupied = np.zeros(nbin, dtype=bool)
        if len(axial_pos):
            idx = np.clip((axial_pos / max(bin_w, 1e-30)).astype(int), 0, nbin - 1)
            occupied[idx] = True
        axial_coverage = float(np.mean(occupied)) if len(occupied) else 0.0
        anchors = np.r_[0.0, axial_pos[(axial_pos >= 0.0) & (axial_pos <= length)], length]
        max_gap = self._max_unique_gap(anchors)

        surface_gaps = np.asarray(self._bond_surface_gap(ids), float)
        surface_gaps = surface_gaps[np.isfinite(surface_gaps)]
        surface_gap = float(np.min(surface_gaps)) if surface_gaps.size else float("inf")
        clearance = self.shell_radius_m - float(np.linalg.norm(tip_xy - self.patch_center_xy))
        remote = float(
            self.remote_sigma_max_Pa
            if remote_sigma_max_Pa is None else remote_sigma_max_Pa
        )
        local = float(0.0 if local_effective_stress_Pa is None else local_effective_stress_Pa)
        Y = float(cfg.handoff_edge_geometry_factor)
        base_remote = (
            0.0 if not np.isfinite(remote)
            else Y * max(remote, 0.0) * math.sqrt(math.pi * max(length, 0.0)) / 1e6
        )
        theta_t = math.atan2(tangent[1], tangent[0])
        KI = base_remote * math.cos(theta_t) ** 2
        KII = base_remote * math.sin(theta_t) * math.cos(theta_t)
        Keq = math.hypot(KI, KII)
        Klocal = Y * max(local, 0.0) * math.sqrt(math.pi * max(length, 0.0)) / 1e6

        min_length = max(
            cfg.handoff_min_length_m,
            cfg.handoff_min_length_horizons * cfg.horizon_m,
        )
        surface_tol = max(
            cfg.surface_connection_horizons * cfg.horizon_m,
            self.point_spacing_m,
        )
        length_pass = length >= min_length and len(ids) >= int(cfg.handoff_min_connected_bonds)
        slenderness_pass = (
            width_ratio <= cfg.handoff_max_width_ratio
            and slenderness >= cfg.handoff_min_slenderness
        )
        tip_pass = tip_width <= cfg.handoff_max_tip_width_horizons * cfg.horizon_m
        boundary_pass = clearance >= cfg.handoff_min_boundary_clearance_horizons * cfg.horizon_m
        K_pass = cfg.handoff_min_remote_K_MPam05 <= 0.0 or Keq >= cfg.handoff_min_remote_K_MPam05
        orientation_pass = coherence >= cfg.handoff_min_orientation_coherence
        coverage_pass = (
            axial_coverage >= cfg.handoff_min_axial_coverage
            and max_gap <= cfg.handoff_max_axial_gap_horizons * cfg.horizon_m
        )
        surface_pass = surface_gap <= surface_tol
        front_pass = (not cfg.handoff_require_active_front) or bool(state.active_front)
        broken_all = state.bond_damage >= cfg.broken_damage
        if state.active_front:
            allowed = np.asarray(state.front_wake_bonds, bool) & broken_all
        else:
            allowed = np.zeros(len(self.bonds), dtype=bool)
            allowed[ids] = True
        off_front = int(np.count_nonzero(broken_all & ~allowed))
        nbroken_all = int(np.count_nonzero(broken_all))
        off_front_fraction = off_front / max(nbroken_all, 1)
        competition_pass = off_front_fraction <= cfg.handoff_max_offfront_broken_fraction

        if cfg.handoff_mode == "legacy_extent":
            handoff = forward_extent >= cfg.established_extent_m
        elif cfg.handoff_mode == "physical":
            handoff = bool(
                length_pass and slenderness_pass and tip_pass and boundary_pass
                and K_pass and orientation_pass and coverage_pass and surface_pass
                and front_pass and competition_pass
            )
        else:
            handoff = False

        checks = {
            "length": length_pass,
            "slenderness": slenderness_pass,
            "tip": tip_pass,
            "boundary": boundary_pass,
            "K": K_pass,
            "orientation": orientation_pass,
            "coverage": coverage_pass,
            "surface": surface_pass,
            "front": front_pass,
            "competition": competition_pass,
        }
        reasons = tuple(name for name, passed in checks.items() if not passed)
        return CrackHandoffAudit(
            root_connected=True,
            connected_bonds=len(ids),
            forward_extent_m=forward_extent,
            centerline_length_m=length,
            width_m=width,
            tip_width_m=tip_width,
            width_ratio=width_ratio,
            tip_radius_eff_m=rho_eff,
            orientation_deg=float(np.degrees(theta_t)),
            boundary_clearance_m=clearance,
            remote_K_MPam05=Keq,
            local_upper_K_MPam05=Klocal,
            remote_KI_MPam05=KI,
            remote_KII_MPam05=KII,
            orientation_coherence=coherence,
            axial_coverage=axial_coverage,
            max_axial_gap_m=max_gap,
            surface_gap_m=surface_gap,
            surface_contact_x_m=float(contact_xy[0]),
            surface_contact_y_m=float(contact_xy[1]),
            slenderness=slenderness,
            length_pass=bool(length_pass),
            slenderness_pass=bool(slenderness_pass),
            tip_pass=bool(tip_pass),
            boundary_pass=bool(boundary_pass),
            K_pass=bool(K_pass),
            orientation_pass=bool(orientation_pass),
            coverage_pass=bool(coverage_pass),
            surface_pass=bool(surface_pass),
            front_pass=bool(front_pass),
            competition_pass=bool(competition_pass),
            off_front_broken_fraction=float(off_front_fraction),
            handoff_pass=bool(handoff),
            failure_reasons=reasons,
        )

    def connected_crack(self, state: StatefulPDState):
        audit = self.crack_handoff_audit(
            state, remote_sigma_max_Pa=self.remote_sigma_max_Pa
        )
        return float(audit.forward_extent_m), int(audit.connected_bonds)

    def state_fields_global(self, state: StatefulPDState, nn_global: int):
        out = {}
        for name in (
            "available", "embryo", "stable", "delivery_memory", "completion",
            "growth", "crack_normal_c2", "crack_normal_s2",
            "crack_orientation_weight",
        ):
            arr = np.full(nn_global, np.nan)
            arr[self.global_nodes] = getattr(state, name)
            out[name] = arr
        return out

    def plot_initiation_diagnostics(self, state, out_png: Path, title: str = ""):
        """Plot the fields controlling the location of the first event."""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        rates = getattr(self, "last_rates", {})
        delivery_rate = np.asarray(rates.get("delivery_rate_s", np.zeros(len(self.xy))), float)
        nucleation_rate = np.asarray(rates.get("nucleation_rate_s", np.zeros(len(self.xy))), float)
        birth_rate = np.asarray(rates.get("mu_birth", np.zeros(len(self.xy))), float)
        expected_birth_node = (
            np.maximum(state.born_cumulative, 0.0) * self.mean_candidate_sites
        )
        fields = [
            ("initiation weight", self.initiation_weight),
            ("candidate sites", state.candidate_sites.astype(float)),
            ("log10 delivery rate (s^-1)", np.log10(np.maximum(delivery_rate, 1e-300))),
            ("log10 cleavage rate (s^-1)", np.log10(np.maximum(nucleation_rate, 1e-300))),
            ("completion Q(K,Lambda)", state.completion),
            ("log10 birth intensity / cycle", np.log10(np.maximum(birth_rate, 1e-300))),
            ("expected cumulative births", expected_birth_node),
        ]
        fig, axes = plt.subplots(2, 4, figsize=(15.5, 7.4), constrained_layout=True)
        x = self.xy[:, 0] * 1e3
        y = self.xy[:, 1] * 1e3
        for ax, (label, vals) in zip(axes.flat, fields):
            sc = ax.scatter(x, y, c=vals, s=18, edgecolors="none")
            ax.scatter(
                x[self.boundary], y[self.boundary], s=7,
                facecolors="none", edgecolors="k", linewidths=0.3
            )
            ax.plot(self.root_xy[0] * 1e3, self.root_xy[1] * 1e3, "rx", ms=7)
            ax.set_aspect("equal")
            ax.set_title(label)
            ax.set_xlabel("x (mm)")
            ax.set_ylabel("y (mm)")
            fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.03)
        for ax in axes.flat[len(fields):]:
            ax.set_visible(False)
        fig.suptitle(title or "Stateful PD initiation diagnostics")
        fig.savefig(Path(out_png), dpi=220)
        plt.close(fig)

    def plot_snapshot(self, state, out_png: Path, title: str = ""):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.collections import LineCollection

        out_png = Path(out_png)
        out_png.parent.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(9.0, 5.8), constrained_layout=True)
        seg = np.stack([self.xy[self.bonds[:, 0]], self.xy[self.bonds[:, 1]]], axis=1) * 1e3
        lc = LineCollection(seg, array=state.bond_damage, linewidths=0.6, cmap="viridis")
        lc.set_clim(0.0, 1.0)
        ax.add_collection(lc)
        realized = state.stable_sites.astype(float)
        sc = ax.scatter(self.xy[:, 0] * 1e3, self.xy[:, 1] * 1e3, c=realized, s=12, cmap="magma", edgecolors="none")
        active = self.initiation_weight > 0.0
        ax.scatter(self.xy[active, 0] * 1e3, self.xy[active, 1] * 1e3, s=6, facecolors="none", edgecolors="tab:blue", linewidths=0.25)
        ax.scatter(self.xy[self.boundary, 0] * 1e3, self.xy[self.boundary, 1] * 1e3, s=5, facecolors="none", edgecolors="k", linewidths=0.3)
        if len(self.feature_surface_xy):
            ax.plot(self.feature_surface_xy[:, 0] * 1e3, self.feature_surface_xy[:, 1] * 1e3, "k-", lw=1.2, label="scratch surface")
        ax.plot(self.root_xy[0] * 1e3, self.root_xy[1] * 1e3, "rx", ms=8)

        orient_nodes = np.where((state.stable_sites > 0) & (state.crack_orientation_weight > 0))[0]
        if len(orient_nodes):
            normal = self._director_to_normal(state.crack_normal_c2[orient_nodes], state.crack_normal_s2[orient_nodes])
            tangent = np.column_stack([-normal[:, 1], normal[:, 0]])
            qlen = 0.45 * self.cfg.horizon_m * 1e3
            ax.quiver(
                self.xy[orient_nodes, 0] * 1e3,
                self.xy[orient_nodes, 1] * 1e3,
                tangent[:, 0], tangent[:, 1],
                angles="xy", scale_units="xy", scale=1.0 / max(qlen, 1e-30),
                width=0.004, label="stored crack tangent",
            )

        root_ids = self._root_connected_bond_ids(state)
        if len(root_ids):
            root_seg = seg[root_ids]
            root_lc = LineCollection(root_seg, linewidths=2.0, linestyles="solid")
            ax.add_collection(root_lc)
        if state.active_front:
            front_ids = np.where(np.asarray(state.active_front_bonds, bool))[0]
            if len(front_ids):
                front_lc = LineCollection(seg[front_ids], linewidths=2.8, linestyles="solid")
                ax.add_collection(front_lc)
            contact = np.asarray(state.active_front_contact_xy, float) * 1e3
            tip = np.asarray(state.active_front_tip_xy, float) * 1e3
            path = np.asarray(state.active_front_path_xy, float) * 1e3
            if path.ndim == 2 and len(path) >= 2:
                ax.plot(path[:, 0], path[:, 1], "-", lw=1.8, label="active crack backbone")
            ax.plot(contact[0], contact[1], marker="s", ms=6, linestyle="none", label="front contact")
            ax.plot(tip[0], tip[1], marker="*", ms=10, linestyle="none", label="active tip")
        audit = self.crack_handoff_audit(
            state,
            remote_sigma_max_Pa=self.remote_sigma_max_Pa,
            local_effective_stress_Pa=float(np.max(getattr(self, "last_rates", {}).get("smax", [0.0]))),
        )
        subtitle = (
            f"front={int(state.active_front)}, a={audit.centerline_length_m*1e6:.1f} um, "
            f"C={audit.orientation_coherence:.2f}, cov={audit.axial_coverage:.2f}, "
            f"off={audit.off_front_broken_fraction:.2f}, mode={state.front_candidate_mode}, "
            f"handoff={audit.handoff_pass}, diffuse_abort={state.cycles_diffuse_abort is not None}, "
            f"front_stalled={state.cycles_front_stalled is not None}"
        )
        if audit.failure_reasons:
            subtitle += " | fail: " + ",".join(audit.failure_reasons)
        ax.set_aspect("equal")
        ax.autoscale()
        ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)")
        ax.set_title((title or "Stateful peridynamic initiation patch") + "\n" + subtitle)
        cb1 = fig.colorbar(lc, ax=ax, fraction=0.045, pad=0.03)
        cb1.set_label("bond cohesive damage")
        cb2 = fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.10)
        cb2.set_label("realized stable defects / point")
        fig.savefig(out_png, dpi=220)
        plt.close(fig)

