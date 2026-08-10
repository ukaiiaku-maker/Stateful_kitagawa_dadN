"""Phase-resolved deterministic quiet-tail and marked-renewal utilities."""
from __future__ import annotations

import copy
import math
from pathlib import Path

import numpy as np

from .config import EV_TO_J
from .v9_array_codec import AtomicArrayGenerationStore
from .v9_canonical_four_class_fem import _tensor_history
from .v9_stable_birth_energy_gate import evaluate_stable_birth_energy_gate


KERNEL_SCHEMA = "V9_CANONICAL_QUIET_TAIL_KERNEL_1"
CHECKPOINT_SCHEMA = "V9_CANONICAL_QUIET_TAIL_CHECKPOINT_1"


def phase_resolved_cycle(condition):
    """Return the accepted FEM/root/MPZ cycle map at one pre-birth boundary."""
    if hasattr(condition, "_elastic_history"):
        history = condition._elastic_history
    else:
        Umax, Umin, u_zero, _F0, _ = condition.cached_fem.affine(
            condition.fem.ep_gp, condition.sigma_max_Pa,
            condition.sigma_min_Pa, condition.fem.u)
        history = condition.cached_fem.stress_histories(
            condition.fem.ep_gp, Umax, Umin,
            condition.args.hazard_n_phase, u_zero)
    tensors = _tensor_history(history["sigma_node"], condition.root_node)
    drives = [condition.birth.mpz.resolve_root_tensor(t) for t in tensors]
    rates = [condition.birth.cleavage_rates(d["opening_stress_Pa"], condition.args.T)
             for d in drives]
    effective = np.asarray([r["cleavage_rate_effective_s"] for r in rates], float)
    dt_phase = 1.0 / (float(condition.args.frequency_Hz) * len(tensors))
    action = np.cumsum(effective * dt_phase)
    mpz = condition.birth.mpz.summary()
    return {
        "schema": KERNEL_SCHEMA,
        "cycles": float(condition.cycles),
        "phase": np.arange(len(tensors), dtype=float) / len(tensors),
        "root_tensors_Pa": tensors,
        "root_opening_Pa": np.asarray([d["opening_stress_Pa"] for d in drives]),
        "signed_channel_stress_Pa": np.stack([d["tau_signed_Pa"] for d in drives]),
        "cleavage_log_rate_raw_s": np.asarray([r["cleavage_log_rate_raw_s"] for r in rates]),
        "cleavage_log_rate_effective_s": np.asarray([r["cleavage_log_rate_effective_s"] for r in rates]),
        "cleavage_action_within_cycle": action,
        "cycle_hazard": float(action[-1]) if len(action) else 0.0,
        "u_phase": np.asarray(history["u_hist"]),
        "mobile_positive": np.asarray(mpz["mobile_positive"]),
        "mobile_negative": np.asarray(mpz["mobile_negative"]),
        "retained_positive": np.asarray(mpz["retained_positive"]),
        "retained_negative": np.asarray(mpz["retained_negative"]),
        "accumulated_slip_positive": np.asarray(mpz["accumulated_slip_positive"]),
        "accumulated_slip_negative": np.asarray(mpz["accumulated_slip_negative"]),
        "rho_back_by_system_m2": np.asarray(mpz["rho_back_by_system_m2"]),
        "tau_back_by_system_Pa": np.asarray(mpz["tau_back_by_system_Pa"]),
        "sigma_back_by_system_Pa": np.asarray(mpz["sigma_back_by_system_Pa"]),
        "signed_active_K_shield_Pa_sqrt_m": float(mpz["signed_active_K_shield_Pa_sqrt_m"]),
        "tip_radius_m": float(mpz["tip_radius_m"]),
        "wake_mobile_positive": np.asarray(condition.birth.mpz.state.wake_mobile_positive),
        "wake_mobile_negative": np.asarray(condition.birth.mpz.state.wake_mobile_negative),
        "wake_retained_positive": np.asarray(condition.birth.mpz.state.wake_retained_positive),
        "wake_retained_negative": np.asarray(condition.birth.mpz.state.wake_retained_negative),
    }


def _relative_error(a, b, floor=1e-30):
    av, bv = np.asarray(a, float), np.asarray(b, float)
    return float(np.max(np.abs(av-bv) / np.maximum(np.maximum(np.abs(av), np.abs(bv)), floor)))


def kernel_convergence(previous, current):
    signed_keys = (
        "mobile_positive", "mobile_negative", "retained_positive", "retained_negative",
        "wake_mobile_positive", "wake_mobile_negative", "wake_retained_positive",
        "wake_retained_negative", "accumulated_slip_positive", "accumulated_slip_negative",
    )
    signed = max(_relative_error(previous[k], current[k]) for k in signed_keys)
    return {
        "signed_state_relative_error": signed,
        "root_tensor_relative_error": _relative_error(previous["root_tensors_Pa"], current["root_tensors_Pa"], 1.0),
        "backstress_relative_error": _relative_error(previous["sigma_back_by_system_Pa"], current["sigma_back_by_system_Pa"], 1.0),
        "shielding_relative_error": _relative_error(previous["signed_active_K_shield_Pa_sqrt_m"], current["signed_active_K_shield_Pa_sqrt_m"], 1e-12),
        "blunting_relative_error": _relative_error(previous["tip_radius_m"], current["tip_radius_m"], 1e-15),
        "phase_log_rate_absolute_error": float(np.max(np.abs(previous["cleavage_log_rate_effective_s"]-current["cleavage_log_rate_effective_s"]))),
        "cycle_hazard_relative_error": _relative_error(previous["cycle_hazard"], current["cycle_hazard"], 1e-300),
    }


def save_condition_checkpoint(condition, root: Path, kernel, summary):
    capsule = condition.capsule()
    arrays = {}
    if "fem" in capsule:
        arrays.update({
            "fem_ep_gp": capsule["fem"].pop("ep_gp"),
            "fem_rho_gp": capsule["fem"].pop("rho_gp"),
            "fem_epsp_acc_gp": capsule["fem"].pop("epsp_acc_gp"),
            "fem_u": capsule["fem"].pop("u"),
        })
    mpz_arrays = capsule["birth"]["mpz"].pop("arrays")
    arrays.update({f"mpz_{key}": value for key, value in mpz_arrays.items()})
    for key, value in kernel.items():
        if isinstance(value, np.ndarray):
            arrays[f"kernel_{key}"] = value
    metadata = {"checkpoint_schema": CHECKPOINT_SCHEMA, "condition_capsule": capsule,
                "kernel_scalars": {k: v for k, v in kernel.items() if not isinstance(v, np.ndarray)}}
    return AtomicArrayGenerationStore(Path(root)).write(arrays, metadata, summary)


def restore_condition_checkpoint(condition, root: Path, generation: str | None = None):
    """Restore one hash-verified generation, or the atomically active one."""
    arrays, metadata, summary, manifest = AtomicArrayGenerationStore(Path(root)).load(generation)
    if metadata.get("checkpoint_schema") != CHECKPOINT_SCHEMA:
        raise RuntimeError("quiet-tail checkpoint schema mismatch")
    capsule = copy.deepcopy(metadata["condition_capsule"])
    if "fem" in capsule:
        capsule["fem"].update({
            "ep_gp": arrays["fem_ep_gp"], "rho_gp": arrays["fem_rho_gp"],
            "epsp_acc_gp": arrays["fem_epsp_acc_gp"], "u": arrays["fem_u"],
        })
    capsule["birth"]["mpz"]["arrays"] = {
        key[len("mpz_"):]: value for key, value in arrays.items() if key.startswith("mpz_")}
    condition.restore_capsule(capsule)
    kernel = dict(metadata["kernel_scalars"])
    kernel.update({key[len("kernel_"):]: value for key, value in arrays.items()
                   if key.startswith("kernel_")})
    return kernel, summary, manifest


def energy_gate_envelope(condition, kernel, xi_values):
    rows = []
    for phase_index, (phi, tensor, displacement) in enumerate(zip(
            kernel["phase"], kernel["root_tensors_Pa"], kernel["u_phase"])):
        phase_cache = {}
        drive = condition.birth.mpz.resolve_root_tensor(tensor)
        sigma_eff = condition.birth.effective_opening_stress_Pa(drive["opening_stress_Pa"])
        event_K = sigma_eff * math.sqrt(2.0*math.pi*kernel["tip_radius_m"])
        barrier_eV = float(np.asarray(
            condition.birth.mpz.state.manifest.cleavage.values_eV(sigma_eff, condition.args.T)))
        for xi in xi_values:
            gate = evaluate_stable_birth_energy_gate(
                mesh=condition.mesh, boundaries=condition.boundaries,
                displacement=displacement, ep_gp=condition.fem.ep_gp,
                Dmat=condition.cached_fem.Dmat, root_xy=condition.root_xy,
                event_direction=np.array([1.0, 0.0]), event_K_Pa_sqrt_m=event_K,
                cleavage_barrier_J=barrier_eV*EV_TO_J,
                cooperative_hits=condition.birth.m_hits,
                burgers_m=condition.birth.mpz.burgers_m,
                threshold_action=float(xi),
                plane_strain_modulus_Pa=condition.cached_fem.material.Eprime,
                evaluation_cache=phase_cache)
            trial_rows = gate["trial_rows"]
            released = max((float(row["elastic_release_J_per_m"]) for row in trial_rows), default=0.0)
            rows.append({
                "phase_index": phase_index, "phase": float(phi), "Xi": float(xi),
                "proposed_length_m": gate["stochastic_proposed_event_length_m"],
                "maximum_released_energy_J_per_m": released,
                "resistance_J_per_m2": gate["hazard_resistance_J_per_m2"],
                "admitted_length_m": gate["committed_event_length_m"],
                "reason": gate["arrest_reason"],
                "mesh_resolved": any(bool(r["topology_changed"]) for r in trial_rows),
            })
    return rows


def stationary_marked_renewal(horizon_action, cycle_hazard, phase_action, envelope_rows,
                              xi_weights=None):
    """Deterministic phase/Xi quadrature for the stationary marked-renewal law.

    The iid exponential renewal intervals are integrated with Gauss-Laguerre
    nodes.  Event phase is advanced in cumulative-hazard coordinates, retaining
    the same Xi for the threshold-scaled proposal.  The returned absorption
    probability uses a substochastic phase operator and uniformization in the
    expected attempt count.
    """
    from numpy.polynomial.laguerre import laggauss
    from scipy.linalg import expm

    h = max(float(cycle_hazard), 1e-300)
    phase_action = np.asarray(phase_action, float)
    edges = np.concatenate(([0.0], phase_action))
    nphase = len(phase_action)
    nodes, weights = laggauss(48) if xi_weights is None else xi_weights
    lookup = {(int(r["phase_index"]), round(float(r["Xi"]), 12)):
              float(r["admitted_length_m"]) > 0.0 for r in envelope_rows}
    xi_grid = np.asarray(sorted({float(r["Xi"]) for r in envelope_rows}))

    # In the VHCF tail h is commonly far below machine-useful modulo scales.
    # The wrapped exponential density in action phase differs from uniform by
    # O(h); evaluate its h->0 limit deterministically instead of taking xi % h
    # with a catastrophic quotient.  Xi remains coupled to its own proposal.
    if h < 1.0e-6:
        phase_width = np.diff(np.concatenate(([0.0], phase_action))) / h
        bounds = np.empty(len(xi_grid)+1)
        bounds[0] = 0.0
        bounds[1:-1] = 0.5*(xi_grid[:-1]+xi_grid[1:])
        bounds[-1] = math.inf
        xi_mass = np.exp(-bounds[:-1]) - np.exp(-bounds[1:])
        p_admit = 0.0
        for j, pw in enumerate(phase_width):
            flags = np.asarray([lookup[(j, round(float(x), 12))] for x in xi_grid], float)
            p_admit += float(pw) * float(xi_mass @ flags)
        p_admit = min(max(p_admit, 0.0), 1.0)
        survival = math.exp(-float(horizon_action)*p_admit)
        return {"stable_birth_survival": survival,
                "stable_birth_probability": -math.expm1(-float(horizon_action)*p_admit),
                "expected_attempt_count": float(horizon_action),
                "stationary_attempt_admission_probability": p_admit,
                "operator_spectral_radius": 1.0-p_admit,
                "wrapped_phase_approximation": "h_to_zero_uniform_action_phase",
                "wrapped_phase_total_variation_bound": 0.5*h}

    def admitted(j, xi):
        nearest = float(xi_grid[np.argmin(np.abs(xi_grid-xi))])
        return lookup[(j, round(nearest, 12))]

    Q = np.zeros((nphase, nphase))
    absorb = np.zeros(nphase)
    starts = np.concatenate(([0.0], phase_action[:-1]))
    for i, a0 in enumerate(starts):
        for xi, w in zip(nodes, weights):
            amod = (a0 + float(xi)) % h
            j = min(int(np.searchsorted(edges, amod, side="right")-1), nphase-1)
            if admitted(j, float(xi)): absorb[i] += float(w)
            else: Q[i, j] += float(w)
    # Uniformization: attempt count in action is Poisson(H); Q is the rejected
    # transition matrix conditional on an attempt.
    initial = np.zeros(nphase); initial[0] = 1.0
    survival = float(initial @ expm(float(horizon_action)*(Q-np.eye(nphase))) @ np.ones(nphase))
    expected_attempts = float(horizon_action)
    admission_probability = min(max(float(np.mean(absorb)), 0.0), 1.0)
    return {"stable_birth_survival": min(max(survival, 0.0), 1.0),
            "stable_birth_probability": min(max(1.0-survival, 0.0), 1.0),
            "expected_attempt_count": expected_attempts,
            "stationary_attempt_admission_probability": admission_probability,
            "operator_spectral_radius": float(max(abs(np.linalg.eigvals(Q))))}
