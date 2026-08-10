"""Full-FEM transaction path for canonical four-class stable-crack birth."""
from __future__ import annotations

import copy
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from .config import ElasticProperties
from .sn_arrhenius_chain import build_chain_from_namespace
from .sn_feature_geometry_v8_7 import (
    BluntNotchGeometry, identify_feature_surface_nodes, local_root_radius,
    local_root_xy, make_blunt_edge_notch_mesh,
)
from .sn_intact_fem import plane_strain_D
from .v9_cached_fem import CachedIntactFEM
from .v9_fem_transaction import EmbeddedFEMTransaction, FEMPhysicalState
from .v9_canonical_four_class_birth import CanonicalFourClassBirthState, MODEL_ID


FEM_MODEL_ID = "v9_full_fem_canonical_four_class_stable_birth_v1"


def _tensor_history(sigma_node_phase, root_node):
    history = np.asarray(sigma_node_phase, float)
    values = (
        history[:, :, int(root_node)]
        if history.ndim == 3 and history.shape[1] == 3
        else history[:, int(root_node), :]
    )
    if values.shape[1] != 3:
        raise RuntimeError("FEM stress history must use [xx,yy,xy] components")
    out = np.empty((len(values), 2, 2), float)
    out[:, 0, 0] = values[:, 0]
    out[:, 1, 1] = values[:, 1]
    out[:, 0, 1] = out[:, 1, 0] = values[:, 2]
    return out


class CanonicalFourClassFEMCondition:
    """Accepted intact-FEM evolution coupled only to the canonical endpoint."""

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
        root_node = int(feature_nodes[np.argmin(np.linalg.norm(mesh.nodes[feature_nodes] - root_xy, axis=1))])
        root_radius = local_root_radius(mesh, feature_nodes)
        sigma_max = 2.0 * self.sigma_a_MPa * 1e6 / max(1.0 - args.R, 1e-30)
        sigma_min = args.R * sigma_max
        chain = build_chain_from_namespace(args, mat.b)
        cached = CachedIntactFEM(mesh, boundaries, mat, Dmat)
        transaction = EmbeddedFEMTransaction(
            mesh=mesh, boundaries=boundaries, material=mat, Dmat=Dmat,
            plastic_chain=chain, args=args, sigma_max_Pa=sigma_max,
            sigma_min_Pa=sigma_min, cached_fem=cached,
        )
        self.mesh = mesh
        self.boundaries = boundaries
        self.root_node = root_node
        self.root_xy = np.asarray(root_xy).copy()
        self.root_radius_initial_m = float(root_radius)
        self.sigma_max_Pa = float(sigma_max)
        self.sigma_min_Pa = float(sigma_min)
        self.cached_fem = cached
        self.fem_transaction = transaction
        self.fem = FEMPhysicalState(
            np.zeros((3, mesh.ne)), np.full(mesh.ne, args.rho0),
            np.zeros(mesh.ne), np.zeros(mesh.ndof), 0.0,
        )
        shear_modulus = mat.E / (2.0 * (1.0 + mat.nu))
        self.birth = CanonicalFourClassBirthState(
            option_id, source_root, shear_modulus_Pa=shear_modulus,
            poisson=mat.nu, burgers_m=mat.b,
            initial_tip_radius_m=self.root_radius_initial_m,
            hazard_seed=hazard_seed,
        )
        self.cycles = 0.0
        self.accepted_blocks = 0
        self.rejected_blocks = 0
        self.audit = {
            "model_id": FEM_MODEL_ID, "endpoint_model_id": MODEL_ID,
            "option_id": option_id, "sigma_a_MPa": self.sigma_a_MPa,
            "T_K": float(args.T), "R": float(args.R),
            "frequency_Hz": float(args.frequency_Hz),
            "root_node": root_node, "root_xy_m": self.root_xy.tolist(),
            "root_radius_initial_m": self.root_radius_initial_m,
            "post_birth_PD_active": False, "post_birth_growth_active": False,
        } | copy.deepcopy(self.birth.audit)

    @classmethod
    def from_run_args(cls, run_args_path, option_id, source_root, sigma_a_MPa,
                      hazard_seed=1720):
        values = json.loads(Path(run_args_path).read_text(), parse_constant=lambda _x: float("inf"))
        values["T"] = float(values["T"])
        return cls(SimpleNamespace(**values), option_id, source_root, sigma_a_MPa, hazard_seed)

    def copy(self):
        clone = copy.copy(self)
        clone.fem = FEMPhysicalState(
            self.fem.ep_gp.copy(), self.fem.rho_gp.copy(),
            self.fem.epsp_acc_gp.copy(), self.fem.u.copy(),
            float(self.fem.plastic_work_J_per_m),
        )
        clone.birth = self.birth.copy()
        clone.audit = copy.deepcopy(self.audit)
        return clone

    def _root_history(self, fem_state):
        Umax, Umin, u_zero, _F0, _ = self.cached_fem.affine(
            fem_state.ep_gp, self.sigma_max_Pa, self.sigma_min_Pa, fem_state.u
        )
        history = self.cached_fem.stress_histories(
            fem_state.ep_gp, Umax, Umin, self.args.hazard_n_phase, u_zero
        )
        return _tensor_history(history["sigma_node"], self.root_node)

    def advance(self, requested_cycles):
        if self.birth.fired:
            return self.summary() | {"cycles_consumed": 0.0, "cycles_unused": float(requested_cycles)}
        requested = max(float(requested_cycles), 0.0)
        proposal_cycles = requested
        while True:
            proposal = self.fem_transaction.propose(self.fem, proposal_cycles)
            if proposal.normalized_error <= 1.0:
                pre = self._root_history(self.fem)
                post = self._root_history(proposal.state)
                root_history = 0.5 * (pre + post)
                birth_start = self.birth.copy()
                try:
                    endpoint = self.birth.advance_fem_phase_block(
                        proposal_cycles, self.args.frequency_Hz, self.args.T,
                        root_history,
                    )
                    break
                except RuntimeError as exc:
                    self.birth = birth_start
                    if "failed to bracket persistent-site backstress root" not in str(exc):
                        raise
            proposal_cycles *= 0.5
            self.rejected_blocks += 1
            if proposal_cycles < self.args.min_block_cycles:
                raise RuntimeError("canonical FEM/MPZ transaction failed below minimum block")
        consumed = float(endpoint["cycles_consumed"])
        if consumed < proposal_cycles:
            localized = self.fem_transaction.propose(self.fem, consumed)
            if localized.normalized_error > 1.0:
                raise RuntimeError("localized canonical FEM proposal violates error gate")
            self.fem = localized.state
        else:
            self.fem = proposal.state
        self.cycles += consumed
        self.accepted_blocks += 1
        return self.summary() | endpoint | {
            "cycles_consumed": consumed,
            "cycles_unused": requested - consumed,
            "fem_normalized_error": proposal.normalized_error,
        }

    def capsule(self):
        return {
            "schema": "V9_CANONICAL_FOUR_CLASS_FEM_CAPSULE_1",
            "audit": copy.deepcopy(self.audit), "cycles": self.cycles,
            "accepted_blocks": self.accepted_blocks,
            "rejected_blocks": self.rejected_blocks,
            "fem": {
                "ep_gp": self.fem.ep_gp.copy(), "rho_gp": self.fem.rho_gp.copy(),
                "epsp_acc_gp": self.fem.epsp_acc_gp.copy(), "u": self.fem.u.copy(),
                "plastic_work_J_per_m": float(self.fem.plastic_work_J_per_m),
            },
            "birth": self.birth.capsule(),
        }

    def restore_capsule(self, capsule):
        if capsule.get("schema") != "V9_CANONICAL_FOUR_CLASS_FEM_CAPSULE_1" or capsule.get("audit") != self.audit:
            raise RuntimeError("canonical FEM capsule provenance/request mismatch")
        fem = capsule["fem"]
        self.fem = FEMPhysicalState(
            np.asarray(fem["ep_gp"], float).copy(), np.asarray(fem["rho_gp"], float).copy(),
            np.asarray(fem["epsp_acc_gp"], float).copy(), np.asarray(fem["u"], float).copy(),
            float(fem["plastic_work_J_per_m"]),
        )
        self.birth.restore_capsule(capsule["birth"])
        self.cycles = float(capsule["cycles"])
        self.accepted_blocks = int(capsule["accepted_blocks"])
        self.rejected_blocks = int(capsule["rejected_blocks"])

    def summary(self):
        return {
            "model_id": FEM_MODEL_ID, "option_id": self.option_id,
            "sigma_a_MPa": self.sigma_a_MPa, "T_K": float(self.args.T),
            "cycles": self.cycles, "accepted_blocks": self.accepted_blocks,
            "rejected_blocks": self.rejected_blocks,
            "post_birth_growth_executed": False, "pd_state_present": False,
        } | self.birth.diagnostics() | self.birth.mpz.summary()
