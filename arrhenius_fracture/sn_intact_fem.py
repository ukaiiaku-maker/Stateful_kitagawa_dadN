"""Intact FEM mechanics used by the stateful peridynamic S-N pilot.

This module intentionally contains no phase-field or diffuse-damage imports.
It provides the blunt-scratch elastic/plastic cycle mechanics needed to drive a
local peridynamic initiation patch while keeping the global FEM mesh intact.
"""
from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

from .mesh import TriMesh, BoundaryData
from .config import ElasticProperties
from .sn_arrhenius_chain import ArrheniusPlasticChain
from .sn_geometry import feature_tangent_normal


def plane_strain_D(mat: ElasticProperties) -> np.ndarray:
    E, nu = float(mat.E), float(mat.nu)
    c = E / ((1.0 + nu) * (1.0 - 2.0 * nu))
    return c * np.array([
        [1.0 - nu, nu, 0.0],
        [nu, 1.0 - nu, 0.0],
        [0.0, 0.0, 0.5 * (1.0 - 2.0 * nu)],
    ])


def _element_dofs(mesh: TriMesh) -> np.ndarray:
    conn = mesh.elems
    edofs = np.empty((mesh.ne, 6), dtype=int)
    edofs[:, 0] = 2 * conn[:, 0]
    edofs[:, 1] = 2 * conn[:, 0] + 1
    edofs[:, 2] = 2 * conn[:, 1]
    edofs[:, 3] = 2 * conn[:, 1] + 1
    edofs[:, 4] = 2 * conn[:, 2]
    edofs[:, 5] = 2 * conn[:, 2] + 1
    return edofs


def assemble_intact_mechanics(
    mesh: TriMesh,
    u: np.ndarray,
    ep_gp: np.ndarray,
    D: np.ndarray,
    mat: ElasticProperties,
):
    """Assemble small-strain plane-strain FEM with plastic eigenstrain.

    No damage or stiffness-degradation variable is present.
    """
    conn = mesh.elems
    edofs = _element_dofs(mesh)
    ue = u[edofs]
    eps_tot = np.einsum("eij,ej->ei", mesh.B_e, ue)
    eps_e = eps_tot - ep_gp.T
    sig = eps_e @ D.T

    sx, sy, txy = sig[:, 0], sig[:, 1], sig[:, 2]
    szz = mat.nu * (sx + sy)
    seq = np.sqrt(
        0.5 * ((sx - sy) ** 2 + (sy - szz) ** 2 + (szz - sx) ** 2)
        + 3.0 * txy**2
    )
    savg = 0.5 * (sx + sy)
    rad = np.sqrt((0.5 * (sx - sy)) ** 2 + txy**2)
    s1 = savg + rad
    psi = 0.5 * np.sum(eps_e * (eps_e @ D.T), axis=1)
    psi = np.where(s1 > 0.0, psi, 0.0)

    A = mesh.area_e
    BTD = np.einsum("eji,jk->eik", mesh.B_e, D)
    Ke = np.einsum("eij,ejk->eik", BTD, mesh.B_e) * A[:, None, None]
    Re = np.einsum("eji,ej->ei", mesh.B_e, sig) * A[:, None]

    ii = np.repeat(edofs[:, :, None], 6, axis=2)
    jj = np.repeat(edofs[:, None, :], 6, axis=1)
    K = sparse.csr_matrix(
        (Ke.ravel(), (ii.ravel(), jj.ravel())), shape=(mesh.ndof, mesh.ndof)
    )
    Rint = np.zeros(mesh.ndof)
    np.add.at(Rint, edofs.ravel(), Re.ravel())
    return K, Rint, sig.T, seq, s1, psi


def stress_state_intact(
    mesh: TriMesh,
    u: np.ndarray,
    ep_gp: np.ndarray,
    D: np.ndarray,
    mat: ElasticProperties,
):
    edofs = _element_dofs(mesh)
    ue = u[edofs]
    eps_tot = np.einsum("eij,ej->ei", mesh.B_e, ue)
    eps_e = eps_tot - ep_gp.T
    sig = eps_e @ D.T
    sx, sy, txy = sig[:, 0], sig[:, 1], sig[:, 2]
    szz = mat.nu * (sx + sy)
    seq = np.sqrt(
        0.5 * ((sx - sy) ** 2 + (sy - szz) ** 2 + (szz - sx) ** 2)
        + 3.0 * txy**2
    )
    savg = 0.5 * (sx + sy)
    rad = np.sqrt((0.5 * (sx - sy)) ** 2 + txy**2)
    s1 = savg + rad
    psi = 0.5 * np.sum(eps_e * (eps_e @ D.T), axis=1)
    psi = np.where(s1 > 0.0, psi, 0.0)
    return sig.T, seq, s1, psi


def project_gp_to_nodes(mesh: TriMesh, val_gp: np.ndarray) -> np.ndarray:
    arr = np.asarray(val_gp)
    if arr.ndim == 1:
        return _project_scalar(mesh, arr)
    out = np.zeros((arr.shape[0], mesh.nn), dtype=float)
    for k in range(arr.shape[0]):
        out[k] = _project_scalar(mesh, arr[k])
    return out


def _project_scalar(mesh: TriMesh, val_gp: np.ndarray) -> np.ndarray:
    conn = mesh.elems
    contrib = np.asarray(val_gp, float) * mesh.area_e / 3.0
    acc = np.zeros(mesh.nn)
    wacc = np.zeros(mesh.nn)
    for a in range(3):
        np.add.at(acc, conn[:, a], contrib)
        np.add.at(wacc, conn[:, a], mesh.area_e / 3.0)
    return acc / np.maximum(wacc, 1e-30)


def lumped_nodal_area(mesh: TriMesh) -> np.ndarray:
    area = np.zeros(mesh.nn)
    for a in range(3):
        np.add.at(area, mesh.elems[:, a], mesh.area_e / 3.0)
    return area


def solve_symmetric_tension(K, Rint, u, bnd: BoundaryData, mesh: TriMesh, Uy_top, Uy_bot):
    prescribed = np.zeros(mesh.ndof, dtype=bool)
    up = np.zeros(mesh.ndof)
    prescribed[2 * bnd.top_nodes + 1] = True
    up[2 * bnd.top_nodes + 1] = Uy_top
    prescribed[2 * bnd.bot_nodes + 1] = True
    up[2 * bnd.bot_nodes + 1] = Uy_bot
    x, y = mesh.nodes[:, 0], mesh.nodes[:, 1]
    anchor = int(np.argmin((x - x.max()) ** 2 + y**2))
    prescribed[2 * anchor] = True
    free = ~prescribed
    Kc = K.tocsr()
    du_p = up[prescribed] - u[prescribed]
    rhs = -Rint[free] - Kc[np.ix_(free, prescribed)] @ du_p
    un = u.copy()
    un[free] = u[free] + spsolve(Kc[np.ix_(free, free)], rhs)
    un[prescribed] = up[prescribed]
    Rfull = Rint + Kc @ (un - u)
    Ftop = float(np.sum(Rfull[2 * bnd.top_nodes + 1]))
    return un, Ftop


def affine_stress_control_displacements(
    mesh,
    bnd,
    mat,
    Dmat,
    ep_gp,
    sigma_max_Pa,
    sigma_min_Pa,
    u_guess,
):
    K, Rint, *_ = assemble_intact_mechanics(mesh, u_guess, ep_gp, Dmat, mat)
    uz, F0 = solve_symmetric_tension(K, Rint, u_guess, bnd, mesh, 0.0, 0.0)
    K0, R0, *_ = assemble_intact_mechanics(mesh, uz, ep_gp, Dmat, mat)
    Uprobe = max(mesh.hbar_tip, 1e-8) * 1e-3
    _, Fp = solve_symmetric_tension(K0, R0, uz, bnd, mesh, Uprobe, -Uprobe)
    slope = (Fp - F0) / max(Uprobe, 1e-30)
    if not np.isfinite(slope) or abs(slope) < 1e-30:
        raise RuntimeError("stress-control displacement calibration failed")
    area2d = max(float(mesh.nodes[:, 0].max() - mesh.nodes[:, 0].min()), 1e-30)
    Umax = (sigma_max_Pa * area2d - F0) / slope
    Umin = (sigma_min_Pa * area2d - F0) / slope
    return float(Umax), float(Umin), uz, float(F0), float(slope)


def equivalent_dep_from_tensor(dep):
    exx, eyy, gxy = dep[0], dep[1], dep[2]
    ezz = -(exx + eyy)
    tensor_sq = exx**2 + eyy**2 + ezz**2 + 0.5 * gxy**2
    return np.sqrt(np.maximum((2.0 / 3.0) * tensor_sq, 0.0))


def _von_mises_flow_direction(sigma_gp, nu):
    sx, sy, txy = sigma_gp[0], sigma_gp[1], sigma_gp[2]
    szz = nu * (sx + sy)
    mean = (sx + sy + szz) / 3.0
    sxx, syy, szzd = sx - mean, sy - mean, szz - mean
    seq = np.sqrt(np.maximum(1.5 * (sxx**2 + syy**2 + szzd**2 + 2.0 * txy**2), 0.0))
    norm = np.sqrt(np.maximum(sxx**2 + syy**2 + szzd**2 + 2.0 * txy**2, 1e-60))
    return seq, sxx / norm, syy / norm, txy / norm


def arrhenius_chain_plastic_update(
    ep_gp,
    rho_gp,
    sigma_gp,
    mat,
    T_K,
    dt_s,
    chain: ArrheniusPlasticChain,
    k_store,
    k_dyn,
    rho_floor,
    rho_cap,
    max_dep_phase,
    max_rho_rel_phase,
):
    seq, nxx, nyy, nxy = _von_mises_flow_direction(sigma_gp, mat.nu)
    rho = np.clip(np.asarray(rho_gp, float), rho_floor, rho_cap)
    rr = chain.rates(seq, rho, T_K)
    dep_prop = np.asarray(rr["dot_ep"], float) * max(float(dt_s), 0.0)
    sqrt23 = np.sqrt(2.0 / 3.0)
    dep_relax = 0.999 * seq / np.maximum(3.0 * mat.G * sqrt23, 1e-30)
    dep = np.minimum(dep_prop, dep_relax)
    if np.isfinite(max_dep_phase) and max_dep_phase > 0.0:
        dep = np.minimum(dep, max_dep_phase)
    dep = np.maximum(dep, 0.0)

    dgamma = sqrt23 * dep
    ep_new = ep_gp.copy()
    ep_new[0] += 1.5 * dgamma * nxx
    ep_new[1] += 1.5 * dgamma * nyy
    ep_new[2] += 1.5 * dgamma * nxy

    drho = k_store * np.sqrt(np.maximum(rho, 1e-30)) / max(mat.b, 1e-30) * dep - k_dyn * rho * dep
    if np.isfinite(max_rho_rel_phase) and max_rho_rel_phase > 0.0:
        drho = np.clip(drho, -max_rho_rel_phase * rho, max_rho_rel_phase * rho)
    rho_new = np.clip(rho + drho, rho_floor, rho_cap)

    seq_after = np.maximum(seq - 3.0 * mat.G * dgamma, 0.0)
    dWp = 0.5 * (seq + seq_after) * dep
    bd = chain.barrier_diagnostics(seq, rho, T_K)
    diag = {
        **rr,
        **bd,
        "seq_Pa": seq,
        "dep_eq": dep,
        "drho": rho_new - rho,
        "dWp_J_m3": dWp,
    }
    return ep_new, rho_new, dep, dWp, diag


def representative_plastic_cycle(
    mesh,
    bnd,
    mat,
    Dmat,
    ep_gp,
    rho_gp,
    Umax,
    Umin,
    T_K,
    frequency_Hz,
    n_phase,
    plast_chain,
    u_start,
    k_store,
    k_dyn,
    rho_floor,
    rho_cap,
    max_dep_phase,
    max_rho_rel_phase,
):
    phase = np.linspace(0.0, 2.0 * np.pi, int(n_phase), endpoint=False)
    frac = 0.5 * (1.0 + np.cos(phase))
    Uhist = Umin + (Umax - Umin) * frac
    dt_phase = 1.0 / max(frequency_Hz * len(Uhist), 1e-30)

    ep0, rho0 = ep_gp.copy(), rho_gp.copy()
    ep, rho, u = ep_gp.copy(), rho_gp.copy(), u_start.copy()
    dep_acc = np.zeros(mesh.ne)
    Wp_acc = np.zeros(mesh.ne)
    max_seq = np.zeros(mesh.ne)
    accum = {k: np.zeros(mesh.ne) for k in ("lambda_emit", "lambda_peierls", "lambda_taylor", "lambda_escape", "lambda_flow")}
    phi_sum = np.zeros(mesh.ne)
    Gsum = {k: np.zeros(mesh.ne) for k in ("G_emit_eV", "G_peierls_eV", "G_taylor_eV")}

    for U in Uhist:
        K, Rint, *_ = assemble_intact_mechanics(mesh, u, ep, Dmat, mat)
        u, _ = solve_symmetric_tension(K, Rint, u, bnd, mesh, U, -U)
        sig, seq, _, _ = stress_state_intact(mesh, u, ep, Dmat, mat)
        max_seq = np.maximum(max_seq, seq)
        ep, rho, dep_phase, dWp_phase, diag = arrhenius_chain_plastic_update(
            ep,
            rho,
            sig,
            mat,
            T_K,
            dt_phase,
            plast_chain,
            k_store,
            k_dyn,
            rho_floor,
            rho_cap,
            max_dep_phase,
            max_rho_rel_phase,
        )
        dep_acc += dep_phase
        Wp_acc += dWp_phase
        for k in accum:
            accum[k] += np.asarray(diag[k], float) * dt_phase
        phi_sum += np.asarray(diag["phi_taylor"], float)
        for k in Gsum:
            Gsum[k] += np.asarray(diag[k], float)
        K, Rint, *_ = assemble_intact_mechanics(mesh, u, ep, Dmat, mat)
        u, _ = solve_symmetric_tension(K, Rint, u, bnd, mesh, U, -U)

    nph = max(len(Uhist), 1)
    return {
        "dep_tensor_cycle": ep - ep0,
        "dep_eq_cycle": dep_acc,
        "drho_cycle": rho - rho0,
        "Wp_cycle_gp": Wp_acc,
        "u_end": u,
        "max_seq_gp": max_seq,
        "mu_emit_cycle_gp": accum["lambda_emit"],
        "mu_peierls_cycle_gp": accum["lambda_peierls"],
        "mu_taylor_cycle_gp": accum["lambda_taylor"],
        "mu_escape_cycle_gp": accum["lambda_escape"],
        "mu_flow_cycle_gp": accum["lambda_flow"],
        "phi_taylor_mean_gp": phi_sum / nph,
        "G_emit_mean_eV_gp": Gsum["G_emit_eV"] / nph,
        "G_peierls_mean_eV_gp": Gsum["G_peierls_eV"] / nph,
        "G_taylor_mean_eV_gp": Gsum["G_taylor_eV"] / nph,
    }


def cycle_stress_histories(mesh, bnd, mat, Dmat, ep_gp, Umax, Umin, n_phase, u_start):
    phase = np.linspace(0.0, 2.0 * np.pi, int(n_phase), endpoint=False)
    frac = 0.5 * (1.0 + np.cos(phase))
    Uhist = Umin + (Umax - Umin) * frac
    u = u_start.copy()
    sig_nodes, seq_nodes, s1_nodes, psi_nodes, Fhist, uhist = [], [], [], [], [], []
    for U in Uhist:
        K, Rint, *_ = assemble_intact_mechanics(mesh, u, ep_gp, Dmat, mat)
        u, Ft = solve_symmetric_tension(K, Rint, u, bnd, mesh, U, -U)
        sig, seq, s1, psi = stress_state_intact(mesh, u, ep_gp, Dmat, mat)
        sig_nodes.append(project_gp_to_nodes(mesh, sig))
        seq_nodes.append(project_gp_to_nodes(mesh, seq))
        s1_nodes.append(project_gp_to_nodes(mesh, s1))
        psi_nodes.append(project_gp_to_nodes(mesh, psi))
        Fhist.append(Ft)
        uhist.append(u.copy())
    return {
        "sigma_node": np.asarray(sig_nodes),
        "seq_node": np.asarray(seq_nodes),
        "s1_node": np.asarray(s1_nodes),
        "psi_node": np.asarray(psi_nodes),
        "Ftop": np.asarray(Fhist),
        "u_hist": np.asarray(uhist),
        "u_max": np.asarray(uhist[0]),
        "u_end": u,
    }


def project_plastic_state(mesh, epsp_acc_gp, rho_gp, epsp_shield_scale, epsp_damage_scale):
    epsp_node = np.maximum(project_gp_to_nodes(mesh, epsp_acc_gp), 0.0)
    rho_node = np.maximum(project_gp_to_nodes(mesh, rho_gp), 0.0)
    P = 1.0 - np.exp(-epsp_node / max(epsp_shield_scale, 1e-30))
    Dloc = 1.0 - np.exp(-epsp_node / max(epsp_damage_scale, 1e-30))
    return epsp_node, rho_node, P, Dloc


def surface_morphology_proposal(
    mesh,
    feature_nodes,
    dep_tensor_block,
    morph_band_length,
    normal_weight,
    shear_weight,
):
    idx = np.asarray(feature_nodes, dtype=int)
    t, n = feature_tangent_normal(mesh, idx)
    exx_n = project_gp_to_nodes(mesh, dep_tensor_block[0])[idx]
    eyy_n = project_gp_to_nodes(mesh, dep_tensor_block[1])[idx]
    gxy_n = project_gp_to_nodes(mesh, dep_tensor_block[2])[idx]
    E = np.zeros((len(idx), 2, 2), dtype=float)
    E[:, 0, 0] = exx_n
    E[:, 1, 1] = eyy_n
    E[:, 0, 1] = E[:, 1, 0] = 0.5 * gxy_n
    Enn = np.einsum("ni,nij,nj->n", n, E, n)
    Ent = np.einsum("ni,nij,nj->n", n, E, t)
    gamma_nt = 2.0 * Ent
    dep_eq_gp = equivalent_dep_from_tensor(dep_tensor_block)
    dep_eq_node = project_gp_to_nodes(mesh, dep_eq_gp)[idx]
    sign_reg = np.tanh(gamma_nt / np.maximum(0.05 * dep_eq_node, 1e-16))
    signed_slip = dep_eq_node * sign_reg
    dh = morph_band_length * (normal_weight * Enn + shear_weight * signed_slip)
    return dh, Enn, gamma_nt
