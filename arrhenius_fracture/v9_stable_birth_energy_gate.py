"""v10.2.30 post-first-passage energy gate for the v9 intact FEM path.

The trial crack exists only inside this calculation.  A positive admitted
length is the terminal stable-birth event; neither the trial damage nor its
equilibrated displacement is committed to the pre-birth trajectory.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import math

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

from .config import EV_TO_J
from .sn_intact_fem import _element_dofs


MODEL_ID = "v9_v10230_mesh_consistent_fixed_opening_energy_gate_v1"


@dataclass(frozen=True)
class StableBirthEnergyGateConfig:
    base_checkpoint_m: float = 5.0e-6
    minimum_factor: float = 0.5
    maximum_factor: float = 4.0
    trial_fraction: float = 0.1
    bisection_iterations: int = 24
    relative_energy_tolerance: float = 1.0e-8
    absolute_energy_tolerance_J_per_m: float = 1.0e-12
    residual_stiffness: float = 1.0e-6
    gamma_relative: float = 1.0

    def validate(self):
        if self.base_checkpoint_m <= 0 or not 0 < self.trial_fraction <= 1:
            raise ValueError("invalid v10.2.30 energy-gate length configuration")
        if self.minimum_factor <= 0 or self.maximum_factor < self.minimum_factor:
            raise ValueError("invalid threshold-scaled event bounds")
        if self.bisection_iterations < 0 or self.residual_stiffness <= 0:
            raise ValueError("invalid energy-gate numerical configuration")
        return self


def clipped_exponential_mean(a: float, b: float) -> float:
    a = max(float(a), 0.0); b = max(float(b), a)
    return max(a + math.exp(-a) - math.exp(-b), 1.0e-300)


def threshold_scaled_event_length(threshold_action: float, config=StableBirthEnergyGateConfig()):
    cfg = config.validate()
    clipped = min(max(float(threshold_action), cfg.minimum_factor), cfg.maximum_factor)
    factor = clipped / clipped_exponential_mean(cfg.minimum_factor, cfg.maximum_factor)
    return cfg.base_checkpoint_m * factor, factor


def hazard_resistance_J_per_m2(barrier_J, cooperative_hits, burgers_m, gamma_relative):
    return (max(float(gamma_relative), 1e-12) * max(float(cooperative_hits), 1.0)
            * max(float(barrier_J), 0.0) / max(abs(float(burgers_m)), 1e-300) ** 2)


def _damage_for_segment(mesh, p0, p1):
    seg = np.asarray(p1, float) - np.asarray(p0, float)
    length2 = float(seg @ seg)
    damage = np.zeros(mesh.nn)
    if length2 <= 0:
        return damage
    centroids = mesh.nodes[mesh.elems].mean(axis=1)
    t = np.clip(((centroids - p0[None, :]) @ seg) / length2, 0.0, 1.0)
    projection = p0[None, :] + t[:, None] * seg[None, :]
    distance2 = np.sum((centroids - projection) ** 2, axis=1)
    radius = 0.7 * np.sqrt(np.maximum(mesh.area_e, 1e-30))
    selected = distance2 <= radius * radius
    if np.any(selected):
        damage[mesh.elems[selected]] = 1.0
    return damage


def _assemble_damaged(mesh, u, ep_gp, D, damage, kappa):
    edofs = _element_dofs(mesh)
    eps = np.einsum("eij,ej->ei", mesh.B_e, np.asarray(u)[edofs]) - np.asarray(ep_gp).T
    dgp = np.mean(np.asarray(damage)[mesh.elems], axis=1)
    degradation = (1.0 - dgp) ** 2 + float(kappa)
    stress = degradation[:, None] * (eps @ D.T)
    gD = degradation[:, None, None] * D[None, :, :]
    BtgD = np.einsum("eji,ejk->eik", mesh.B_e, gD)
    Ke = np.einsum("eij,ejk->eik", BtgD, mesh.B_e) * mesh.area_e[:, None, None]
    Re = np.einsum("eji,ej->ei", mesh.B_e, stress) * mesh.area_e[:, None]
    ii = np.repeat(edofs[:, :, None], 6, axis=2)
    jj = np.repeat(edofs[:, None, :], 6, axis=1)
    K = sparse.csr_matrix((Ke.ravel(), (ii.ravel(), jj.ravel())), shape=(mesh.ndof, mesh.ndof))
    R = np.zeros(mesh.ndof); np.add.at(R, edofs.ravel(), Re.ravel())
    energy = float(np.sum(degradation * 0.5 * np.sum(eps * (eps @ D.T), axis=1) * mesh.area_e))
    return K, R, energy


def _fixed_opening(mesh, boundaries, u_initial, ep_gp, D, damage, kappa):
    K, R, _ = _assemble_damaged(mesh, u_initial, ep_gp, D, damage, kappa)
    prescribed = np.zeros(mesh.ndof, bool); up = np.zeros(mesh.ndof)
    prescribed[2 * boundaries.top_nodes + 1] = True
    prescribed[2 * boundaries.bot_nodes + 1] = True
    up[2 * boundaries.top_nodes + 1] = float(np.mean(u_initial[2 * boundaries.top_nodes + 1]))
    up[2 * boundaries.bot_nodes + 1] = float(np.mean(u_initial[2 * boundaries.bot_nodes + 1]))
    # Match v10.2.30 ``solve_dirichlet`` exactly: both bottom corners restrain
    # x, and the left-bottom y value is the imposed bottom opening.
    prescribed[2 * boundaries.left_bot] = True
    prescribed[2 * boundaries.left_bot + 1] = True
    up[2 * boundaries.left_bot + 1] = up[2 * boundaries.bot_nodes[0] + 1]
    prescribed[2 * boundaries.right_bot] = True
    free = ~prescribed; old = np.asarray(u_initial, float)
    rhs = -R[free] - K[np.ix_(free, prescribed)] @ (up[prescribed] - old[prescribed])
    u = old.copy(); u[free] = old[free] + spsolve(K[np.ix_(free, free)].tocsc(), rhs)
    u[prescribed] = up[prescribed]
    _, _, energy = _assemble_damaged(mesh, u, ep_gp, D, damage, kappa)
    return u, energy


def evaluate_stable_birth_energy_gate(*, mesh, boundaries, displacement, ep_gp, Dmat,
                                      root_xy, event_direction, event_K_Pa_sqrt_m,
                                      cleavage_barrier_J, cooperative_hits, burgers_m,
                                      threshold_action, plane_strain_modulus_Pa,
                                      config=StableBirthEnergyGateConfig()):
    """Exact v10.2.30 mesh-consistent search, without committing geometry."""
    cfg = config.validate()
    proposal, factor = threshold_scaled_event_length(threshold_action, cfg)
    direction = np.asarray(event_direction, float).reshape(2)
    direction /= max(float(np.linalg.norm(direction)), 1e-300)
    p0 = np.asarray(root_xy, float).reshape(2)
    resistance = hazard_resistance_J_per_m2(
        cleavage_barrier_J, cooperative_hits, burgers_m, cfg.gamma_relative)
    Eprime = max(float(plane_strain_modulus_Pa), 1e-300)
    zero_damage = np.zeros(mesh.nn)
    u_pre, energy_pre = _fixed_opening(
        mesh, boundaries, displacement, ep_gp, Dmat, zero_damage, cfg.residual_stiffness)
    ntrial = max(int(math.ceil(1.0 / cfg.trial_fraction)), 1)
    candidates = [proposal * i / ntrial for i in range(1, ntrial + 1)]
    rows = []; accepted = 0.0; accepted_u = u_pre; first_failed = None

    def evaluate(length):
        damage = _damage_for_segment(mesh, p0, p0 + length * direction)
        changed = bool(np.count_nonzero(damage) > 0)
        if changed:
            utrial, energy_post = _fixed_opening(
                mesh, boundaries, u_pre, ep_gp, Dmat, damage, cfg.residual_stiffness)
            released = max(energy_pre - energy_post, 0.0)
            source = "fixed_opening_re_equilibrated_energy_drop"
        else:
            utrial, energy_post = u_pre, energy_pre
            released = event_K_Pa_sqrt_m ** 2 / Eprime * length
            source = "directional_J_search_continuation_only"
        dissipated = resistance * length
        tol = max(cfg.absolute_energy_tolerance_J_per_m,
                  cfg.relative_energy_tolerance * max(abs(released), abs(dissipated), 1e-300))
        row = {"trial_length_m": length, "stored_energy_pre_J_per_m": energy_pre,
               "stored_energy_post_J_per_m": energy_post, "elastic_release_J_per_m": released,
               "hazard_dissipation_J_per_m": dissipated, "energy_residual_J_per_m": released-dissipated,
               "energy_tolerance_J_per_m": tol, "admissible": released-dissipated+tol >= 0,
               "newly_killed_nodes": int(np.count_nonzero(damage)), "topology_changed": changed,
               "energy_release_source": source, "subgrid_search_continuation_only": not changed}
        return released - dissipated + tol, utrial, row, changed

    for length in candidates:
        residual, utrial, row, changed = evaluate(length); rows.append(row)
        if residual >= 0 and changed:
            accepted, accepted_u = length, utrial
        elif changed:
            first_failed = length; break
    if accepted > 0 and first_failed is not None:
        lo, hi, ulo = accepted, first_failed, accepted_u
        for _ in range(cfg.bisection_iterations):
            mid = 0.5 * (lo + hi); residual, umid, row, changed = evaluate(mid)
            row["bisection"] = True; rows.append(row)
            if residual >= 0 and changed: lo, ulo = mid, umid
            else: hi = mid
        accepted, accepted_u = lo, ulo
    return {"energy_gate_model_id": MODEL_ID, "energy_gate_config": asdict(cfg),
            "threshold_action": float(threshold_action), "event_length_factor": factor,
            "stochastic_proposed_event_length_m": proposal,
            "energy_admissible_event_length_m": accepted,
            "committed_event_length_m": min(proposal, accepted),
            "arrest_reason": ("stochastic_proposal_reached" if accepted >= proposal*(1-1e-12)
                              else "hazard_derived_energy_arrest" if accepted > 0
                              else "no_mesh_resolved_admissible_increment"),
            "hazard_barrier_J": float(cleavage_barrier_J),
            "hazard_cooperative_hits": float(cooperative_hits),
            "hazard_burgers_vector_m": float(burgers_m),
            "orientation_gamma_relative": cfg.gamma_relative,
            "hazard_resistance_J_per_m2": resistance,
            "event_K_Pa_sqrt_m": float(event_K_Pa_sqrt_m), "trial_rows": rows,
            "equilibrated_displacement": accepted_u, "athermal_Gc_used": False,
            "independent_toughness_floor_used": False, "hazard_gated_before_first_passage": False}
