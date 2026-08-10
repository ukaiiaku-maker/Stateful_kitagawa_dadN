"""Higher-order deterministic quiet-tail integration with direct overlap validation."""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy.integrate import solve_ivp

from .v9_fem_transaction import FEMPhysicalState
from .v9_quiet_tail_kernel import kernel_convergence, phase_resolved_cycle


MODEL_ID = "v9_deterministic_quiet_tail_dop853_v1"


@dataclass
class DOP853Result:
    accepted: bool
    condition: object
    diagnostics: dict


def _pack(fem):
    return np.concatenate((fem.ep_gp.ravel(), fem.rho_gp.ravel(), fem.epsp_acc_gp.ravel()))


def _unpack(condition, vector, template):
    ne=condition.mesh.ne; a=3*ne; b=a+ne
    return FEMPhysicalState(np.asarray(vector[:a]).reshape(3,ne),
        np.asarray(vector[a:b]), np.asarray(vector[b:b+ne]),
        template.u.copy(), float(template.plastic_work_J_per_m))


def propagate_dop853(condition, cycles, *, rtol=2e-4, maximum_step=math.inf,
                     log_cycle_coordinate=False):
    work=condition.copy(); initial=work.fem; y0=_pack(initial); ne=work.mesh.ne
    atol=np.concatenate((np.full(3*ne,1e-8),np.full(ne,1.0),np.full(ne,1e-8)))
    evaluations=0
    N_start=max(float(work.cycles),1.0)
    N_target=N_start+float(cycles)
    def derivative(coordinate,y):
        nonlocal evaluations
        fem=_unpack(work,y,initial)
        cycle,_=work.fem_transaction._cycle_derivative(fem); evaluations+=1
        value=np.concatenate((np.asarray(cycle["dep_tensor_cycle"]).ravel(),
            np.asarray(cycle["drho_cycle"]),np.asarray(cycle["dep_eq_cycle"])))
        return value*(math.exp(coordinate) if log_cycle_coordinate else 1.0)
    span=((math.log(N_start),math.log(N_target)) if log_cycle_coordinate
          else (0.0,float(cycles)))
    rho_slice=slice(3*ne,4*ne)
    # The production constitutive derivative clips its diagnostic increment at
    # rho_cap.  An event placed at the cap can therefore be approached without
    # ever changing sign, forcing an adaptive solver into vanishing steps.
    # Stop just inside the boundary, at a margin far below the qualified state
    # tolerance, and never accept or extrapolate the capped state.
    rho_cap_guard = float(work.args.rho_cap) * (1.0 - 1.0e-10)
    def rho_cap_event(_coordinate,y):
        return rho_cap_guard-float(np.max(y[rho_slice]))
    rho_cap_event.terminal=True
    rho_cap_event.direction=-1
    def rho_floor_event(_coordinate,y):
        return float(np.min(y[rho_slice]))-float(work.args.rho_floor)
    rho_floor_event.terminal=True
    rho_floor_event.direction=-1
    sol=solve_ivp(derivative,span,y0,method="DOP853",rtol=rtol,atol=atol,
                  max_step=float(maximum_step),events=(rho_cap_event,rho_floor_event))
    if not sol.success:
        return DOP853Result(False,condition.copy(),{"model_id":MODEL_ID,"reason":sol.message})
    if sol.status == 1:
        reached=(math.exp(sol.t[-1]) if log_cycle_coordinate else N_start+sol.t[-1])
        which="rho_cap_boundary" if len(sol.t_events[0]) else "rho_floor_boundary"
        return DOP853Result(False,condition.copy(),{"model_id":MODEL_ID,"reason":which,
            "cycle_boundary":reached,"derivative_evaluations":evaluations})
    # Advance the one-way signed MPZ/hazard state over every accepted high-order
    # FEM interval using the already-qualified endpoint-average root history.
    hazard_increments=[]
    accepted_N=(np.exp(sol.t) if log_cycle_coordinate else N_start+sol.t)
    for i in range(1,len(sol.t)):
        start=_unpack(work,sol.y[:,i-1],initial); end=_unpack(work,sol.y[:,i],initial)
        root=0.5*(work._root_history(start)+work._root_history(end))
        endpoint=work.birth.advance_fem_phase_block(
            accepted_N[i]-accepted_N[i-1],work.args.frequency_Hz,work.args.T,root)
        hazard_increments.append(max(float(endpoint.get("hazard_increment",0.0)),0.0))
    work.fem=_unpack(work,sol.y[:,-1],initial)
    # Reconstruct a consistent endpoint displacement from the accepted cycle.
    _cycle,u_end=work.fem_transaction._cycle_derivative(work.fem)
    work.fem=FEMPhysicalState(work.fem.ep_gp,work.fem.rho_gp,work.fem.epsp_acc_gp,
                              u_end,work.fem.plastic_work_J_per_m)
    work.cycles += float(cycles); work.birth.time_s=work.cycles/work.args.frequency_Hz
    diagnostics={"model_id":MODEL_ID,"cycles":float(cycles),"accepted_steps":len(sol.t)-1,
                 "derivative_evaluations":evaluations,"minimum_step":float(np.min(np.diff(sol.t))),
                 "maximum_step":float(np.max(np.diff(sol.t))),"rtol":rtol,
                 "ep_atol":1e-8,"rho_atol":1.0,"epsp_atol":1e-8,
                 "hazard_increment_fsum":math.fsum(hazard_increments),
                 "integration_coordinate":"log_cycles" if log_cycle_coordinate else "cycles"}
    return DOP853Result(True,work,diagnostics)


def validate_dop853(start, cycles, *, direct_cycles=None, rtol=2e-4,
                    log_cycle_coordinate=False):
    accelerated=propagate_dop853(start,cycles,rtol=rtol,
                                 log_cycle_coordinate=log_cycle_coordinate)
    if not accelerated.accepted:
        return accelerated
    direct=start.copy(); direct.advance(float(cycles if direct_cycles is None else direct_cycles))
    if direct_cycles is not None and direct_cycles != cycles:
        comparison=propagate_dop853(start,direct_cycles,rtol=rtol,
                                    log_cycle_coordinate=log_cycle_coordinate).condition
    else:
        comparison=accelerated.condition
    kernel_direct=phase_resolved_cycle(direct); kernel_acc=phase_resolved_cycle(comparison)
    errors=kernel_convergence(kernel_direct,kernel_acc)
    fem_error=float(np.linalg.norm(_pack(direct.fem)-_pack(comparison.fem)))/max(
        float(np.linalg.norm(_pack(direct.fem))),float(np.linalg.norm(_pack(comparison.fem))),1.0)
    accepted=(fem_error<=1e-3 and errors["signed_state_relative_error"]<=1e-3 and
              errors["root_tensor_relative_error"]<=1e-3 and
              errors["backstress_relative_error"]<=1e-3 and
              errors["phase_log_rate_absolute_error"]<=1e-3 and
              errors["cycle_hazard_relative_error"]<=1e-3)
    accelerated.diagnostics.update({"direct_overlap_cycles":float(cycles if direct_cycles is None else direct_cycles),
                                    "fem_active_state_relative_error":fem_error,**errors,
                                    "qualification_passed":accepted})
    accelerated.accepted=accepted
    return accelerated
