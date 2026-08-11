#!/usr/bin/env python3
"""Side-effect-free terminal exact-cycle and mechanics audit for atomic v9 PD."""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path

import numpy as np

from scripts.audit_v9_pd_persistent_sites import build_context
from arrhenius_fracture import sn_pd2d_stateful_v9_transactional as driver


def digest(value) -> str:
    return hashlib.sha256(repr(value).encode()).hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--case-dir", type=Path, required=True)
    p.add_argument("--generation", required=True)
    p.add_argument("--sigma-a-MPa", type=float, required=True)
    p.add_argument("--out", type=Path, required=True)
    ns = p.parse_args()
    args, mesh, patch, chain, crack, cached, transaction, sigma_max, sigma_min = build_context(
        ns.case_dir / "run_args.json", ns.sigma_a_MPa
    )
    stored_args = json.loads((ns.case_dir / "run_args.json").read_text())
    if "delivery_source_multiplicity" not in stored_args:
        # Reconstruct the exact predecessor signature. Rate evaluation uses
        # getattr(..., 1.0), preserving the historical per-source mapping.
        delattr(args, "delivery_source_multiplicity")
    restored = driver._load_case_checkpoint(
        ns.case_dir / "checkpoint_latest.npz", args=args, case_name="shielded",
        sigma_a_MPa=ns.sigma_a_MPa, mesh=mesh, patch=patch,
        generation=ns.generation,
    )
    state = restored["pd_state"]
    protected_before = {
        "candidate_rng": digest(copy.deepcopy(patch._candidate_rng.bit_generator.state)),
        "event_rng": digest(copy.deepcopy(patch._event_rng.bit_generator.state)),
        "threshold": hashlib.sha256(np.asarray(state.site_birth_threshold).tobytes()).hexdigest(),
        "action": hashlib.sha256(np.asarray(state.birth_cumulative_hazard).tobytes()).hexdigest(),
        "topology": hashlib.sha256(np.asarray(state.bond_damage).tobytes()).hexdigest(),
    }
    payload = driver.evaluate_dormant_exact_cycle(
        args=args, shield_on=True, mesh=mesh, patch=patch, pd_state=state,
        crack=crack, plast_chain=chain, cached_fem=cached,
        fem_transaction=transaction, sigma_max=sigma_max, sigma_min=sigma_min,
        ep_gp=restored["ep_gp"], rho_gp=restored["rho_gp"],
        epsp_acc_gp=restored["epsp_acc_gp"], u=restored["u"],
        plastic_work=restored["Wp_total"], cycles=restored["cycles"], dN=1.0,
    )
    protected_after = {
        "candidate_rng": digest(copy.deepcopy(patch._candidate_rng.bit_generator.state)),
        "event_rng": digest(copy.deepcopy(patch._event_rng.bit_generator.state)),
        "threshold": hashlib.sha256(np.asarray(state.site_birth_threshold).tobytes()).hexdigest(),
        "action": hashlib.sha256(np.asarray(state.birth_cumulative_hazard).tobytes()).hexdigest(),
        "topology": hashlib.sha256(np.asarray(state.bond_damage).tobytes()).hexdigest(),
    }
    rates = patch.last_rates
    E = float(args.E_GPa) * 1e9
    local = float(payload["diagnostics"]["fem_sigma1_cycle_max_Pa"])
    remote_ratio = float(sigma_max / E)
    local_ratio = float(local / E)
    # Linearized kinematics omit quadratic strain terms.  The explicit 5%
    # audit bound limits that neglected scale to eps^2 <= 0.0025.  It is a
    # conservative formulation bound, not a fitted fatigue threshold.
    formulation_bound = 0.05
    mechanics = "mechanics_valid_for_physical_interpretation" if local_ratio <= formulation_bound else "diagnostic_overstress_only"
    controller_path = ns.case_dir / "v9_pd_high_cycle_controller.json"
    controller = json.loads(controller_path.read_text()) if controller_path.is_file() else None
    result = {
        "schema": "V9_PD_TERMINAL_EXACT_MECHANICS_AUDIT_1",
        "cycles": float(restored["cycles"]),
        "stress_amplitude_Pa": ns.sigma_a_MPa * 1e6,
        "remote_sigma_min_Pa": float(sigma_min),
        "remote_sigma_max_Pa": float(sigma_max),
        "maximum_FEM_principal_stress_Pa": local,
        "nominal_remote_sigma_over_E": remote_ratio,
        "maximum_local_sigma_over_E": local_ratio,
        "small_strain_formulation_bound": formulation_bound,
        "small_strain_bound_basis": "neglected quadratic strain scale eps^2 <= 0.0025",
        "mechanics_classification": mechanics,
        "cleavage_critical_stress_Pa": float(crack._G0_sigc(args.T)[1]),
        "emission_critical_stress_Pa": float(
            chain.emit.sigc_Pa(args.T) if hasattr(chain.emit, "sigc_Pa") else
            max(chain.emit.sigc0_Pa + chain.emit.sT_Pa_per_K * (args.T - chain.emit.Tref_K), 1.0)
        ),
        "delivery_rate_peak_s": float(np.max(rates["delivery_rate_s"])),
        "delivery_memory_max": float(np.max(state.delivery_memory)),
        "K2_completion_max": float(np.max(state.completion)),
        "raw_cleavage_rate_peak_s": float(np.max(rates["nucleation_rate_s"])),
        "actual_site_birth_rate_peak_per_cycle": float(np.exp(np.max(payload["log_birth_action"]))),
        "exact_cycle_diagnostics": payload["diagnostics"],
        "stationary_controller": controller,
        "side_effect_free": protected_before == protected_after,
        "protected_signatures_before": protected_before,
        "protected_signatures_after": protected_after,
    }
    if not result["side_effect_free"]:
        raise RuntimeError("terminal exact diagnostic mutated protected physical state")
    ns.out.parent.mkdir(parents=True, exist_ok=True)
    serial = driver._json_safe(result)
    ns.out.write_text(json.dumps(serial, indent=2, sort_keys=True) + "\n")
    print(json.dumps(serial, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
