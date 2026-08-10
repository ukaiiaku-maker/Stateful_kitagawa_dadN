"""Validated affine macro-cycle propagation for deterministic v9 quiet tails."""
from __future__ import annotations

import copy
from dataclasses import dataclass
import math

import numpy as np

from .v9_fem_transaction import FEMPhysicalState
from .v9_quiet_tail_kernel import kernel_convergence, phase_resolved_cycle


MODEL_ID = "v9_validated_affine_quiet_tail_macro_map_v1"


@dataclass
class DMDResult:
    accepted: bool
    condition: object
    cycles_advanced: float
    diagnostics: dict


def _capture(condition):
    pieces = [np.asarray(condition.fem.ep_gp).ravel(), np.asarray(condition.fem.rho_gp).ravel()]
    schema = [("fem_ep_gp", condition.fem.ep_gp.shape, pieces[0].size),
              ("fem_rho_gp", condition.fem.rho_gp.shape, pieces[1].size)]
    for name in condition.birth.mpz._CAPSULE_ARRAYS:
        value = np.asarray(getattr(condition.birth.mpz.state, name), float)
        pieces.append(value.ravel()); schema.append((name, value.shape, value.size))
    return np.concatenate(pieces), schema


def _apply(condition, vector, schema, *, projection=True):
    value = np.asarray(vector, float); offset = 0; projected = value.copy()
    decoded = {}
    for name, shape, size in schema:
        part = projected[offset:offset+size].reshape(shape)
        if projection and name == "fem_rho_gp":
            part[:] = np.clip(part, condition.args.rho_floor, condition.args.rho_cap)
        elif projection and name != "fem_ep_gp":
            part[:] = np.maximum(part, 0.0)
        decoded[name] = part.copy(); offset += size
    condition.fem = FEMPhysicalState(
        decoded.pop("fem_ep_gp"), decoded.pop("fem_rho_gp"),
        condition.fem.epsp_acc_gp.copy(), condition.fem.u.copy(),
        float(condition.fem.plastic_work_J_per_m))
    for name, array in decoded.items():
        setattr(condition.birth.mpz.state, name, array)
    modules = condition.birth.mpz.modules["signed_burgers_shared_v1025"]
    modules._sync_active(condition.birth.mpz.state); modules._sync_wake(condition.birth.mpz.state)
    return projected


def _relative(a, b):
    return float(np.linalg.norm(np.asarray(a)-np.asarray(b))) / max(
        float(np.linalg.norm(a)), float(np.linalg.norm(b)), 1.0)


def _fit(states, outputs, maximum_rank=12):
    X = np.column_stack(states); reference = X[:, 0]
    scale = np.maximum(np.max(np.abs(X), axis=1), 1.0)
    Y = (X-reference[:, None])/scale[:, None]
    sample = Y[:, 1:]
    U, singular, _ = np.linalg.svd(sample, full_matrices=False)
    if singular.size and singular[0] > 0:
        rank = max(1, min(maximum_rank, int(np.count_nonzero(singular >= singular[0]*1e-10)), len(states)-1))
        U = U[:, :rank]
    else:
        rank=1; U=np.zeros((len(reference),1))
    Z=U.T@Y; design=np.vstack([Z[:,:-1],np.ones(Z.shape[1]-1)])
    coef=np.linalg.lstsq(design.T,Z[:,1:].T,rcond=None)[0].T
    outcoef=np.linalg.lstsq(design.T,np.asarray(outputs).T,rcond=None)[0].T
    A,c=coef[:,:rank],coef[:,rank]; C,d=outcoef[:,:rank],outcoef[:,rank]
    fitted=A@Z[:,:-1]+c[:,None]
    return {"reference":reference,"scale":scale,"basis":U,"coordinates":Z,
            "A":A,"c":c,"C":C,"d":d,"rank":rank,
            "training_error":_relative(fitted,Z[:,1:]),
            "spectral_radius":float(max(abs(np.linalg.eigvals(A))))}


def _propagate(model, steps):
    A,c,C,d=(model[k] for k in ("A","c","C","d")); r=A.shape[0]; o=C.shape[0]
    M=np.zeros((r+o+1,r+o+1)); M[:r,:r]=A; M[:r,-1]=c
    M[r:r+o,:r]=C; M[r:r+o,r:r+o]=np.eye(o); M[r:r+o,-1]=d; M[-1,-1]=1
    initial=np.zeros(r+o+1); initial[:r]=model["coordinates"][:,-1]; initial[-1]=1
    final=np.linalg.matrix_power(M,int(steps))@initial
    z=final[:r]; vector=model["reference"]+model["scale"]*(model["basis"]@z)
    return vector,final[r:r+o]


def _predict_one(model, vector):
    z=model["basis"].T@((vector-model["reference"])/model["scale"])
    zn=model["A"]@z+model["c"]
    return model["reference"]+model["scale"]*(model["basis"]@zn), model["C"]@z+model["d"]


def propagate_validated_macro_dmd(condition, cycles_requested, *, macro_cycles,
                                  burst_steps=16, maximum_rank=12):
    requested=float(cycles_requested); macro=float(macro_cycles)
    total_steps=int(math.floor(requested/macro))
    if total_steps <= burst_steps+1:
        return DMDResult(False,condition.copy(),0.0,{"reason":"insufficient_projection_steps"})
    work=condition.copy(); states=[]; outputs=[]; vector,schema=_capture(work); states.append(vector)
    for _ in range(burst_steps):
        h0=phase_resolved_cycle(work)["cycle_hazard"]
        work.advance(macro)
        vector,_=_capture(work); states.append(vector)
        h1=phase_resolved_cycle(work)["cycle_hazard"]
        outputs.append([0.5*(h0+h1)*macro])
    model=_fit(states,np.asarray(outputs).T,maximum_rank)
    diagnostics={"model_id":MODEL_ID,"macro_cycles":macro,"burst_steps":burst_steps,
                 "rank":model["rank"],"training_error":model["training_error"],
                 "spectral_radius":model["spectral_radius"]}
    if model["training_error"]>5e-4:
        return DMDResult(False,work,burst_steps*macro,diagnostics|{"reason":"training_error"})
    projected_steps=total_steps-burst_steps
    raw,cumulative=_propagate(model,projected_steps)
    projected=work.copy(); physical=_apply(projected,raw,schema,projection=True)
    projection_error=_relative(raw,physical)
    projected.cycles=work.cycles+projected_steps*macro
    projected.birth.time_s=projected.cycles/projected.args.frequency_Hz
    projected.birth.mpz.state.time_s=projected.birth.time_s
    projected.birth.cumulative_cleavage_hazard += max(float(cumulative[0]),0.0)
    projected.birth.log_cumulative_cleavage_hazard=math.log(max(
        projected.birth.cumulative_cleavage_hazard,1e-300))
    predicted_next,predicted_output=_predict_one(model,physical)
    exact=projected.copy(); h0=phase_resolved_cycle(exact)["cycle_hazard"]
    # Fail closed before adaptive subdivision: a projected endpoint that cannot
    # accept the fitted macro interval is outside the qualified local map.
    trial = exact.fem_transaction.propose(exact.fem, macro)
    if trial.normalized_error > 1.0:
        diagnostics.update({"projection_error": projection_error,
                            "endpoint_macro_normalized_error": trial.normalized_error,
                            "projected_steps": projected_steps,
                            "reason": "projected_endpoint_macro_not_accepted"})
        return DMDResult(False, work, burst_steps*macro, diagnostics)
    exact.advance(macro)
    exact_vector,_=_capture(exact)
    state_error=_relative(exact_vector,predicted_next)
    h1=phase_resolved_cycle(exact)["cycle_hazard"]
    exact_dH=0.5*(h0+h1)*macro
    hazard_error=abs(exact_dH-max(float(predicted_output[0]),0.0))/max(
        exact_dH,abs(float(predicted_output[0])),1e-300)
    kernel_projected=phase_resolved_cycle(projected); kernel_exact=phase_resolved_cycle(exact)
    phase_errors=kernel_convergence(kernel_projected,kernel_exact)
    accepted=(projection_error<=1e-6 and state_error<=1e-3 and hazard_error<=1e-3 and
              phase_errors["root_tensor_relative_error"]<=1e-3 and
              phase_errors["phase_log_rate_absolute_error"]<=1e-3)
    diagnostics.update({"projection_error":projection_error,"state_validation_error":state_error,
                        "hazard_validation_error":hazard_error,**phase_errors,
                        "projected_steps":projected_steps,
                        "reason":None if accepted else "endpoint_validation_failed"})
    return DMDResult(accepted,projected if accepted else work,
                     total_steps*macro if accepted else burst_steps*macro,diagnostics)
