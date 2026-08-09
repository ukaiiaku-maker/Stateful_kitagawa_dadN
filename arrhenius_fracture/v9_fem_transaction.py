"""Side-effect-free embedded FEM/plastic/rho transaction proposals for v9."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .sn_intact_fem import (
    affine_stress_control_displacements,
    representative_plastic_cycle,
)


@dataclass(frozen=True)
class FEMPhysicalState:
    ep_gp: np.ndarray
    rho_gp: np.ndarray
    epsp_acc_gp: np.ndarray
    u: np.ndarray
    plastic_work_J_per_m: float = 0.0

    def copy_readonly(self) -> "FEMPhysicalState":
        arrays = []
        for value in (self.ep_gp, self.rho_gp, self.epsp_acc_gp, self.u):
            copied = np.asarray(value, dtype=float).copy()
            copied.setflags(write=False)
            arrays.append(copied)
        return FEMPhysicalState(*arrays, float(self.plastic_work_J_per_m))


@dataclass(frozen=True)
class FEMEmbeddedProposal:
    dN: float
    state: FEMPhysicalState
    normalized_error: float
    component_errors: dict[str, float]
    activations: tuple[str, ...]


@dataclass(frozen=True)
class FEMAdvanceResult:
    cycle: float
    state: FEMPhysicalState
    accepted_blocks: int
    rejected_blocks: int
    largest_accepted_dN: float
    next_block_dN: float


class EmbeddedFEMTransaction:
    """Heun/Euler embedded macro-step over the real representative cycle."""

    def __init__(
        self, *, mesh, boundaries, material, Dmat, plastic_chain, args,
        sigma_max_Pa: float, sigma_min_Pa: float,
        relative_tolerance: float = 2.0e-4,
        # 1e-8 strain corresponds to roughly 4.1 kPa at tungsten's modulus,
        # about four parts per million of the local GPa stress scale.
        # A near-machine-zero tolerance makes inactive tensor components, not
        # the physical plastic/rho state, control the long-life step size.
        ep_absolute_tolerance: float = 1.0e-8,
        rho_absolute_tolerance: float = 1.0,
        cached_fem=None,
    ):
        self.mesh = mesh
        self.boundaries = boundaries
        self.material = material
        self.Dmat = Dmat
        self.chain = plastic_chain
        self.args = args
        self.sigma_max_Pa = float(sigma_max_Pa)
        self.sigma_min_Pa = float(sigma_min_Pa)
        self.rtol = float(relative_tolerance)
        self.ep_atol = float(ep_absolute_tolerance)
        self.rho_atol = float(rho_absolute_tolerance)
        self.cached_fem = cached_fem

    def _cycle_derivative(self, state: FEMPhysicalState):
        if self.cached_fem is None:
            Umax, Umin, u_zero, _, _ = affine_stress_control_displacements(
                self.mesh, self.boundaries, self.material, self.Dmat,
                state.ep_gp, self.sigma_max_Pa, self.sigma_min_Pa, state.u,
            )
            cycle = representative_plastic_cycle(
                self.mesh, self.boundaries, self.material, self.Dmat,
                state.ep_gp, state.rho_gp, Umax, Umin,
                self.args.T, self.args.frequency_Hz, self.args.plastic_n_phase,
                self.chain, u_zero, self.args.k_store, self.args.k_dyn,
                self.args.rho_floor, self.args.rho_cap,
                self.args.max_dep_phase, self.args.max_rho_rel_phase,
            )
        else:
            Umax, Umin, u_zero, _, _ = self.cached_fem.affine(
                state.ep_gp, self.sigma_max_Pa, self.sigma_min_Pa, state.u
            )
            cycle = self.cached_fem.representative_cycle(
                state.ep_gp, state.rho_gp, Umax, Umin,
                self.args.T, self.args.frequency_Hz, self.args.plastic_n_phase,
                self.chain, u_zero, self.args.k_store, self.args.k_dyn,
                self.args.rho_floor, self.args.rho_cap,
                self.args.max_dep_phase, self.args.max_rho_rel_phase,
            )
        return cycle, np.asarray(cycle["u_end"], dtype=float)

    @staticmethod
    def _scaled_error(a, b, atol, rtol):
        av = np.asarray(a, dtype=float)
        bv = np.asarray(b, dtype=float)
        scale = float(atol) + float(rtol) * np.maximum(np.abs(av), np.abs(bv))
        return float(np.max(np.abs(av - bv) / scale)) if av.size else 0.0

    def propose(self, initial: FEMPhysicalState, dN: float, *, first_cycle=None) -> FEMEmbeddedProposal:
        if not math.isfinite(dN) or dN <= 0.0:
            raise ValueError("dN must be positive and finite")
        state = initial.copy_readonly()
        if first_cycle is None:
            first, _ = self._cycle_derivative(state)
        else:
            first = first_cycle
        k_ep0 = np.asarray(first["dep_tensor_cycle"], dtype=float)
        k_eq0 = np.maximum(np.asarray(first["dep_eq_cycle"], dtype=float), 0.0)
        k_rho0 = np.asarray(first["drho_cycle"], dtype=float)
        predictor = FEMPhysicalState(
            ep_gp=state.ep_gp + dN * k_ep0,
            rho_gp=state.rho_gp + dN * k_rho0,
            epsp_acc_gp=state.epsp_acc_gp + dN * k_eq0,
            u=state.u,
            plastic_work_J_per_m=state.plastic_work_J_per_m,
        )
        activations = []
        if np.any(predictor.rho_gp < self.args.rho_floor) or np.any(predictor.rho_gp > self.args.rho_cap):
            # A cap is an integration boundary, not an accepted-state clip.
            return FEMEmbeddedProposal(
                dN, state, math.inf, {"rho_cap_boundary": math.inf},
                ("rho_state_boundary_requires_localization",),
            )
        second, u_end = self._cycle_derivative(predictor)
        k_ep1 = np.asarray(second["dep_tensor_cycle"], dtype=float)
        k_eq1 = np.maximum(np.asarray(second["dep_eq_cycle"], dtype=float), 0.0)
        k_rho1 = np.asarray(second["drho_cycle"], dtype=float)
        heun_ep = state.ep_gp + 0.5 * dN * (k_ep0 + k_ep1)
        heun_rho = state.rho_gp + 0.5 * dN * (k_rho0 + k_rho1)
        heun_eq = state.epsp_acc_gp + 0.5 * dN * (k_eq0 + k_eq1)
        if np.any(heun_rho < self.args.rho_floor) or np.any(heun_rho > self.args.rho_cap):
            return FEMEmbeddedProposal(
                dN, state, math.inf, {"rho_cap_boundary": math.inf},
                ("rho_state_boundary_requires_localization",),
            )
        area = np.asarray(self.mesh.area_e, dtype=float)
        wp0 = float(np.sum(np.asarray(first["Wp_cycle_gp"], dtype=float) * area))
        wp1 = float(np.sum(np.asarray(second["Wp_cycle_gp"], dtype=float) * area))
        accepted = FEMPhysicalState(
            heun_ep, heun_rho, heun_eq, u_end,
            state.plastic_work_J_per_m + 0.5 * dN * (wp0 + wp1),
        ).copy_readonly()
        errors = {
            "ep_gp": self._scaled_error(heun_ep, predictor.ep_gp, self.ep_atol, self.rtol),
            "rho_gp": self._scaled_error(heun_rho, predictor.rho_gp, self.rho_atol, self.rtol),
            "epsp_acc_gp": self._scaled_error(heun_eq, predictor.epsp_acc_gp, self.ep_atol, self.rtol),
        }
        return FEMEmbeddedProposal(dN, accepted, max(errors.values()), errors, tuple(activations))

    def advance(
        self, initial: FEMPhysicalState, *, cycle_start: float, cycle_end: float,
        initial_block_dN: float, growth_factor: float = 2.0,
        minimum_block_dN: float = 1e-10,
    ) -> FEMAdvanceResult:
        if cycle_end < cycle_start:
            raise ValueError("cycle_end precedes cycle_start")
        state = initial.copy_readonly()
        cycle = float(cycle_start)
        candidate = float(initial_block_dN)
        accepted = rejected = 0
        largest = 0.0
        while cycle < cycle_end:
            dN = min(candidate, cycle_end - cycle)
            proposal = self.propose(state, dN)
            if proposal.normalized_error > 1.0:
                candidate = 0.5 * dN
                rejected += 1
                if candidate < minimum_block_dN:
                    raise RuntimeError("FEM transaction subdivision exhausted")
                continue
            state = proposal.state
            cycle += dN
            accepted += 1
            largest = max(largest, dN)
            if proposal.normalized_error <= 1e-16:
                factor = growth_factor
            else:
                factor = min(growth_factor, max(1.05, 0.9 / math.sqrt(proposal.normalized_error)))
            candidate = max(minimum_block_dN, dN * factor)
        return FEMAdvanceResult(cycle, state, accepted, rejected, largest, candidate)
