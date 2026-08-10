"""Canonical elastic-FEM bulk coupled to the audited signed tip MPZ."""
from __future__ import annotations

import copy
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from .config import EV_TO_J, ElasticProperties
from .sn_feature_geometry_v8_7 import (
    BluntNotchGeometry, identify_feature_surface_nodes, local_root_radius,
    local_root_xy, make_blunt_edge_notch_mesh,
)
from .sn_intact_fem import plane_strain_D
from .v9_cached_fem import CachedIntactFEM
from .v9_canonical_four_class_fem import _tensor_history
from .v9_energy_gated_stable_birth import EnergyGatedStableBirthState, MODEL_ID
from .v9_fem_transaction import FEMPhysicalState
from .v9_stable_birth_energy_gate import evaluate_stable_birth_energy_gate


FEM_MODEL_ID = "v9_canonical_four_class_elastic_fem_signed_mpz_energy_gated_birth_v3"
CAPSULE_SCHEMA = "V9_CANONICAL_FOUR_CLASS_ELASTIC_FEM_CAPSULE_1"
CONTRACT = {
    "continuum_bulk_role": "elastic_fem_only",
    "bulk_state_evolves_in_fem": False,
    "bulk_hardening_law": "none_bulk_elastic",
    "moving_crack_tip_mpz_active": True,
    "tip_mpz_pt_model_active": True,
    "bulk_pt_model_active": False,
    "bulk_scalar_rho_used_for_signed_shielding": False,
}


class CanonicalElasticFourClassFEMCondition:
    """Fixed-geometry elastic bulk; all pre-birth evolution is tip-local MPZ."""

    def __init__(self, args, option_id, source_root, sigma_a_MPa, hazard_seed):
        self.args = copy.deepcopy(args)
        self.option_id = str(option_id)
        self.sigma_a_MPa = float(sigma_a_MPa)
        mat = ElasticProperties(
            E=float(args.E_GPa) * 1e9, nu=float(args.nu), b=float(args.b_m),
            Tm=float(args.Tm_K),
        )
        Dmat = plane_strain_D(mat)
        geom = BluntNotchGeometry(
            args.Lx, args.Ly, args.notch_depth_m, args.notch_half_height_m,
            feature_type=args.feature_type, root_radius_m=args.notch_root_radius_m,
            opening_angle_deg=args.notch_opening_angle_deg,
            path_refine_length_m=args.path_refine_length_m,
            path_refine_half_height_m=args.path_refine_half_height_m,
        )
        mesh, boundaries, _ = make_blunt_edge_notch_mesh(
            geom, nx=args.nx, ny=args.ny, jitter=args.jitter,
            root_h_fine=args.root_h_fine, seed=args.seed,
        )
        feature_nodes = identify_feature_surface_nodes(mesh, geom)
        root_xy = local_root_xy(mesh, feature_nodes)
        root_node = int(feature_nodes[np.argmin(
            np.linalg.norm(mesh.nodes[feature_nodes] - root_xy, axis=1)
        )])
        sigma_max = 2.0 * self.sigma_a_MPa * 1e6 / max(1.0 - args.R, 1e-30)
        sigma_min = args.R * sigma_max
        cached = CachedIntactFEM(mesh, boundaries, mat, Dmat)

        # These arrays satisfy the existing elastic assembler interface only.
        # They are immutable, neutral, and never enter a constitutive update.
        ep_reference = np.zeros((3, mesh.ne), dtype=float)
        rho_reference = np.full(mesh.ne, float(args.rho0), dtype=float)
        epsp_reference = np.zeros(mesh.ne, dtype=float)
        for array in (ep_reference, rho_reference, epsp_reference):
            array.setflags(write=False)
        Umax, Umin, u_zero, _F0, _ = cached.affine(
            ep_reference, sigma_max, sigma_min, np.zeros(mesh.ndof)
        )
        history = cached.stress_histories(
            ep_reference, Umax, Umin, args.hazard_n_phase, u_zero
        )
        self._elastic_history = {
            key: np.asarray(value, float).copy() for key, value in history.items()
        }
        for value in self._elastic_history.values():
            value.setflags(write=False)
        self._root_tensors = _tensor_history(history["sigma_node"], root_node)
        self._root_tensors.setflags(write=False)

        self.mesh = mesh
        self.boundaries = boundaries
        self.root_node = root_node
        self.root_xy = np.asarray(root_xy, float).copy()
        self.root_radius_initial_m = float(local_root_radius(mesh, feature_nodes))
        self.sigma_max_Pa = float(sigma_max)
        self.sigma_min_Pa = float(sigma_min)
        self.cached_fem = cached
        self.fem_transaction = None
        self.fem = FEMPhysicalState(
            ep_reference, rho_reference, epsp_reference,
            np.asarray(u_zero, float).copy(), 0.0,
        )
        shear_modulus = mat.E / (2.0 * (1.0 + mat.nu))
        self.birth = EnergyGatedStableBirthState(
            option_id, source_root, shear_modulus_Pa=shear_modulus,
            poisson=mat.nu, burgers_m=mat.b,
            initial_tip_radius_m=self.root_radius_initial_m,
            hazard_seed=hazard_seed,
        )
        self.cycles = 0.0
        self.accepted_blocks = 0
        self.rejected_blocks = 0
        self.next_block_cycles = max(float(args.min_block_cycles), 1.0)
        self.audit = {
            "model_id": FEM_MODEL_ID, "endpoint_model_id": MODEL_ID,
            "option_id": option_id, "sigma_a_MPa": self.sigma_a_MPa,
            "T_K": float(args.T), "R": float(args.R),
            "frequency_Hz": float(args.frequency_Hz),
            "root_node": root_node, "root_xy_m": self.root_xy.tolist(),
            "root_radius_initial_m": self.root_radius_initial_m,
            "elastic_phase_response_cached": True,
            "ep_gp_role": "immutable_zero_assembler_placeholder",
            "rho_gp_role": "immutable_reference_assembler_placeholder_not_constitutive",
            "fem_rho_cap_active": False,
            "post_birth_PD_active": False, "post_birth_growth_active": False,
            "crack_geometry_committed": False,
            **CONTRACT,
        } | copy.deepcopy(self.birth.audit)

    @classmethod
    def from_run_args(cls, run_args_path, option_id, source_root, sigma_a_MPa,
                      hazard_seed=1720):
        values = json.loads(Path(run_args_path).read_text(),
                            parse_constant=lambda _x: float("inf"))
        values["T"] = float(values["T"])
        return cls(SimpleNamespace(**values), option_id, source_root,
                   sigma_a_MPa, hazard_seed)

    def copy(self):
        clone = copy.copy(self)
        clone.birth = self.birth.copy()
        clone.audit = copy.deepcopy(self.audit)
        return clone

    def _root_history(self, _fem_state=None):
        return self._root_tensors.copy()

    def assert_neutral_fem_state(self):
        if (np.any(self.fem.ep_gp != 0.0) or np.any(self.fem.epsp_acc_gp != 0.0)
                or self.fem.plastic_work_J_per_m != 0.0):
            raise RuntimeError("canonical elastic FEM neutral state was mutated")
        if any(array.flags.writeable for array in (
                self.fem.ep_gp, self.fem.rho_gp, self.fem.epsp_acc_gp)):
            raise RuntimeError("canonical elastic FEM placeholders must be immutable")

    def _event_gate(self):
        drives = [self.birth.mpz.resolve_root_tensor(t) for t in self._root_tensors]
        peak = max(drives, key=lambda row: row["opening_stress_Pa"])
        sigma_eff = self.birth.effective_opening_stress_Pa(peak["opening_stress_Pa"])
        tip_radius = self.birth.mpz.summary()["tip_radius_m"]
        event_K = sigma_eff * np.sqrt(2.0 * np.pi * tip_radius)
        barrier_eV = float(np.asarray(
            self.birth.mpz.state.manifest.cleavage.values_eV(sigma_eff, self.args.T)
        ))
        return evaluate_stable_birth_energy_gate(
            mesh=self.mesh, boundaries=self.boundaries,
            displacement=np.asarray(self._elastic_history["u_max"], float),
            ep_gp=self.fem.ep_gp, Dmat=self.cached_fem.Dmat,
            root_xy=self.root_xy, event_direction=np.array([1.0, 0.0]),
            event_K_Pa_sqrt_m=event_K, cleavage_barrier_J=barrier_eV * EV_TO_J,
            cooperative_hits=self.birth.m_hits,
            burgers_m=self.birth.mpz.burgers_m,
            threshold_action=self.birth.pending_attempt["threshold_action"],
            plane_strain_modulus_Pa=self.cached_fem.material.Eprime,
        )

    def advance(self, requested_cycles):
        self.assert_neutral_fem_state()
        if self.birth.fired:
            return self.summary() | {
                "cycles_consumed": 0.0, "cycles_unused": float(requested_cycles)
            }
        requested = max(float(requested_cycles), 0.0)
        remaining = requested
        consumed_total = 0.0
        while remaining > 0.0 and not self.birth.fired:
            proposal_cycles = min(remaining, self.next_block_cycles)
            birth_start = self.birth.copy()
            try:
                endpoint = self.birth.advance_fem_phase_block(
                    proposal_cycles, self.args.frequency_Hz, self.args.T,
                    self._root_tensors,
                )
            except RuntimeError as exc:
                self.birth = birth_start
                if "failed to bracket persistent-site backstress root" not in str(exc):
                    raise
                proposal_cycles *= 0.5
                self.next_block_cycles = proposal_cycles
                self.rejected_blocks += 1
                if proposal_cycles < self.args.min_block_cycles:
                    raise RuntimeError("canonical elastic-FEM/MPZ block failed below minimum")
                continue
            consumed = float(endpoint["cycles_consumed"])
            self.cycles += consumed
            consumed_total += consumed
            remaining -= consumed
            self.accepted_blocks += 1
            self.next_block_cycles = max(
                float(self.args.min_block_cycles), proposal_cycles * 2.0
            )
            if self.birth.pending_attempt is not None:
                self.birth.resolve_pending_attempt(self._event_gate())
                if self.birth.fired:
                    break
            if consumed <= 0.0 and remaining > 0.0:
                raise RuntimeError("energy-gated endpoint made no progress")
        self.assert_neutral_fem_state()
        return self.summary() | {
            "cycles_consumed": consumed_total,
            "cycles_unused": requested - consumed_total,
        }

    def capsule(self):
        self.assert_neutral_fem_state()
        return {
            "schema": CAPSULE_SCHEMA, "audit": copy.deepcopy(self.audit),
            "cycles": self.cycles, "accepted_blocks": self.accepted_blocks,
            "rejected_blocks": self.rejected_blocks,
            "next_block_cycles": self.next_block_cycles,
            "birth": self.birth.capsule(),
        }

    def restore_capsule(self, capsule):
        if capsule.get("schema") != CAPSULE_SCHEMA or capsule.get("audit") != self.audit:
            raise RuntimeError("canonical elastic-FEM capsule provenance/request mismatch")
        self.birth.restore_capsule(capsule["birth"])
        self.cycles = float(capsule["cycles"])
        self.accepted_blocks = int(capsule["accepted_blocks"])
        self.rejected_blocks = int(capsule["rejected_blocks"])
        self.next_block_cycles = float(capsule["next_block_cycles"])
        self.assert_neutral_fem_state()

    def summary(self):
        return self.birth.diagnostics() | self.birth.mpz.summary() | {
            "model_id": FEM_MODEL_ID, "option_id": self.option_id,
            "sigma_a_MPa": self.sigma_a_MPa, "T_K": float(self.args.T),
            "cycles": self.cycles, "accepted_blocks": self.accepted_blocks,
            "rejected_blocks": self.rejected_blocks,
            "next_block_cycles": self.next_block_cycles,
            "post_birth_growth_executed": False, "pd_state_present": False,
            **CONTRACT,
        }
