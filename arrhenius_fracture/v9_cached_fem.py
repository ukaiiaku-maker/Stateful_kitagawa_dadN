"""Cached linear-intact FEM operations for the v9 physical trajectory."""

from __future__ import annotations

import numpy as np
from scipy.sparse.linalg import factorized

from .sn_intact_fem import (
    _element_dofs,
    arrhenius_chain_plastic_update,
    assemble_intact_mechanics,
    project_gp_to_nodes,
    stress_state_intact,
)


class CachedIntactFEM:
    """Reuse the immutable stiffness/factorization; physics is unchanged."""

    def __init__(self, mesh, boundaries, material, Dmat):
        self.mesh = mesh
        self.boundaries = boundaries
        self.material = material
        self.Dmat = Dmat
        zeros_u = np.zeros(mesh.ndof)
        zeros_ep = np.zeros((3, mesh.ne))
        self.K = assemble_intact_mechanics(mesh, zeros_u, zeros_ep, Dmat, material)[0].tocsr()
        prescribed = np.zeros(mesh.ndof, dtype=bool)
        prescribed[2 * boundaries.top_nodes + 1] = True
        prescribed[2 * boundaries.bot_nodes + 1] = True
        x, y = mesh.nodes[:, 0], mesh.nodes[:, 1]
        self.anchor = int(np.argmin((x - x.max()) ** 2 + y**2))
        prescribed[2 * self.anchor] = True
        self.prescribed = prescribed
        self.free = ~prescribed
        self.Kff = self.K[np.ix_(self.free, self.free)].tocsc()
        self.Kfp = self.K[np.ix_(self.free, self.prescribed)].tocsr()
        self.solve_free = factorized(self.Kff)
        self.edofs = _element_dofs(mesh)

    def internal_force(self, u, ep_gp):
        ue = np.asarray(u, float)[self.edofs]
        eps_total = np.einsum("eij,ej->ei", self.mesh.B_e, ue)
        stress = (eps_total - np.asarray(ep_gp, float).T) @ self.Dmat.T
        element = np.einsum("eji,ej->ei", self.mesh.B_e, stress) * self.mesh.area_e[:, None]
        force = np.zeros(self.mesh.ndof)
        np.add.at(force, self.edofs.ravel(), element.ravel())
        return force

    def solve_tension(self, Rint, u, Uy_top, Uy_bot):
        up = np.zeros(self.mesh.ndof)
        up[2 * self.boundaries.top_nodes + 1] = Uy_top
        up[2 * self.boundaries.bot_nodes + 1] = Uy_bot
        old = np.asarray(u, float)
        delta_prescribed = up[self.prescribed] - old[self.prescribed]
        rhs = -np.asarray(Rint, float)[self.free] - self.Kfp @ delta_prescribed
        new = old.copy()
        new[self.free] = old[self.free] + self.solve_free(rhs)
        new[self.prescribed] = up[self.prescribed]
        reaction = np.asarray(Rint, float) + self.K @ (new - old)
        force_top = float(np.sum(reaction[2 * self.boundaries.top_nodes + 1]))
        return new, force_top

    def affine(self, ep_gp, sigma_max_Pa, sigma_min_Pa, u_guess):
        Rint = self.internal_force(u_guess, ep_gp)
        zero, F0 = self.solve_tension(Rint, u_guess, 0.0, 0.0)
        R0 = self.internal_force(zero, ep_gp)
        probe = max(self.mesh.hbar_tip, 1e-8) * 1e-3
        _, Fp = self.solve_tension(R0, zero, probe, -probe)
        slope = (Fp - F0) / max(probe, 1e-30)
        if not np.isfinite(slope) or abs(slope) < 1e-30:
            raise RuntimeError("cached stress-control calibration failed")
        area2d = max(float(np.ptp(self.mesh.nodes[:, 0])), 1e-30)
        Umax = (sigma_max_Pa * area2d - F0) / slope
        Umin = (sigma_min_Pa * area2d - F0) / slope
        return float(Umax), float(Umin), zero, float(F0), float(slope)

    def representative_cycle(
        self, ep_gp, rho_gp, Umax, Umin, T_K, frequency_Hz, n_phase,
        chain, u_start, k_store, k_dyn, rho_floor, rho_cap,
        max_dep_phase, max_rho_rel_phase,
    ):
        phase = np.linspace(0.0, 2.0 * np.pi, int(n_phase), endpoint=False)
        Uhist = Umin + (Umax - Umin) * 0.5 * (1.0 + np.cos(phase))
        dt_phase = 1.0 / max(frequency_Hz * len(Uhist), 1e-30)
        ep0, rho0 = ep_gp.copy(), rho_gp.copy()
        ep, rho, u = ep_gp.copy(), rho_gp.copy(), u_start.copy()
        dep_acc = np.zeros(self.mesh.ne)
        work_acc = np.zeros(self.mesh.ne)
        max_seq = np.zeros(self.mesh.ne)
        accum = {key: np.zeros(self.mesh.ne) for key in ("lambda_emit", "lambda_peierls", "lambda_taylor", "lambda_escape", "lambda_flow")}
        phi_sum = np.zeros(self.mesh.ne)
        energy_sum = {key: np.zeros(self.mesh.ne) for key in ("G_emit_eV", "G_peierls_eV", "G_taylor_eV")}
        for U in Uhist:
            Rint = self.internal_force(u, ep)
            u, _ = self.solve_tension(Rint, u, U, -U)
            stress, seq, _, _ = stress_state_intact(self.mesh, u, ep, self.Dmat, self.material)
            max_seq = np.maximum(max_seq, seq)
            ep, rho, dep, work, diag = arrhenius_chain_plastic_update(
                ep, rho, stress, self.material, T_K, dt_phase, chain,
                k_store, k_dyn, rho_floor, rho_cap, max_dep_phase, max_rho_rel_phase,
            )
            dep_acc += dep
            work_acc += work
            for key in accum:
                accum[key] += np.asarray(diag[key], float) * dt_phase
            phi_sum += np.asarray(diag["phi_taylor"], float)
            for key in energy_sum:
                energy_sum[key] += np.asarray(diag[key], float)
            Rint = self.internal_force(u, ep)
            u, _ = self.solve_tension(Rint, u, U, -U)
        nph = max(len(Uhist), 1)
        return {
            "dep_tensor_cycle": ep - ep0, "dep_eq_cycle": dep_acc,
            "drho_cycle": rho - rho0, "Wp_cycle_gp": work_acc,
            "u_end": u, "max_seq_gp": max_seq,
            **{f"mu_{key.replace('lambda_', '')}_cycle_gp": value for key, value in accum.items()},
            "phi_taylor_mean_gp": phi_sum / nph,
            **{f"{key.replace('_eV', '')}_mean_eV_gp": value / nph for key, value in energy_sum.items()},
        }

    def stress_histories(self, ep_gp, Umax, Umin, n_phase, u_start):
        phase = np.linspace(0.0, 2.0 * np.pi, int(n_phase), endpoint=False)
        Uhist = Umin + (Umax - Umin) * 0.5 * (1.0 + np.cos(phase))
        u = u_start.copy()
        sig_nodes, seq_nodes, s1_nodes, psi_nodes, forces, displacements = [], [], [], [], [], []
        for U in Uhist:
            Rint = self.internal_force(u, ep_gp)
            u, force = self.solve_tension(Rint, u, U, -U)
            sig, seq, s1, psi = stress_state_intact(
                self.mesh, u, ep_gp, self.Dmat, self.material
            )
            sig_nodes.append(project_gp_to_nodes(self.mesh, sig))
            seq_nodes.append(project_gp_to_nodes(self.mesh, seq))
            s1_nodes.append(project_gp_to_nodes(self.mesh, s1))
            psi_nodes.append(project_gp_to_nodes(self.mesh, psi))
            forces.append(force)
            displacements.append(u.copy())
        return {
            "sigma_node": np.asarray(sig_nodes), "seq_node": np.asarray(seq_nodes),
            "s1_node": np.asarray(s1_nodes), "psi_node": np.asarray(psi_nodes),
            "Ftop": np.asarray(forces), "u_hist": np.asarray(displacements),
            "u_max": displacements[0], "u_min": displacements[len(displacements)//2],
            "u_end": displacements[-1],
        }
