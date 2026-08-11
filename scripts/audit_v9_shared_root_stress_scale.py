#!/usr/bin/env python3
"""Read-only stress-scale audit of completed shared-root Peak checkpoints."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from scripts.audit_v9_pd_persistent_sites import build_context
from arrhenius_fracture import sn_pd2d_stateful_v9_transactional as driver
from arrhenius_fracture.v9_canonical_four_class_birth import (
    canonical_effective_cleavage_rate,
)
from arrhenius_fracture.v9_pd_shared_root_marked_cleavage import (
    SharedRootMarkedCleavageState,
)


DEFAULTS = {
    1000.0: Path("runs/sn_v9_shared_root_marked_pd_qualification/uninterrupted_fix2/Peak/shielded/sigmaA_1000MPa"),
    2000.0: Path("runs/sn_v9_shared_root_marked_pd_long_tail_N1e14/Peak/shielded/sigmaA_2000MPa"),
    2150.0: Path("runs/sn_v9_shared_root_marked_pd_long_tail_N1e12/Peak/shielded/sigmaA_2150MPa"),
    2500.0: Path("runs/sn_v9_shared_root_marked_pd_event_v6/Peak/shielded/sigmaA_2500MPa"),
}
QUANTILE_ACTION = {"N10": -math.log(0.9), "N50": -math.log(0.5), "N90": -math.log(0.1)}


def _clock(args, capsule, root, radius):
    source_root = getattr(args, "four_class_source_root", None)
    if source_root is None:
        source_root = args.four_class_registry_audit["source_repository"]
    state = SharedRootMarkedCleavageState(
        args.four_class_option_id, source_root,
        shear_modulus_Pa=args.E_GPa * 1e9 / (2 * (1 + args.nu)),
        poisson=args.nu, burgers_m=args.b_m, initial_tip_radius_m=radius,
        hazard_seed=args.global_cleavage_seed, mark_seed=args.spatial_mark_seed,
    )
    state.restore_capsule(capsule)
    return state


def _phase_tensors(args, mesh, cached, restored, sigma_max, sigma_min, root):
    Umax, Umin, uz, _, _ = cached.affine(restored["ep_gp"], sigma_max, sigma_min, restored["u"])
    hist = cached.stress_histories(restored["ep_gp"], Umax, Umin, args.hazard_n_phase, uz)
    node = int(np.argmin(np.linalg.norm(mesh.nodes - root[None, :], axis=1)))
    voigt = np.asarray(hist["sigma_node"][:, :, node], float)
    tensors = np.empty((len(voigt), 2, 2))
    tensors[:, 0, 0], tensors[:, 1, 1] = voigt[:, 0], voigt[:, 1]
    tensors[:, 0, 1] = tensors[:, 1, 0] = voigt[:, 2]
    return tensors, hist, node


def _ordered_one_cycle(clock, tensors, frequency, temperature):
    trial = clock.copy()
    dt = 1.0 / frequency / len(tensors)
    action = 0.0
    for tensor in tensors:
        drive = trial.mpz.resolve_root_tensor(tensor)
        trial.mpz.advance(0.5 * dt, temperature, drive["opening_stress_Pa"], drive["tau_signed_Pa"])
        raw = math.exp(trial.mpz.cleavage_log_rate_s(
            trial.effective_opening_stress_Pa(drive["opening_stress_Pa"]), temperature
        ))
        action += canonical_effective_cleavage_rate(raw, 3.0, 1e-6) * dt
        trial.mpz.advance(0.5 * dt, temperature, drive["opening_stress_Pa"], drive["tau_signed_Pa"])
    return trial, action


def audit_condition(stress, case_dir):
    args, mesh, patch, chain, crack, cached, transaction, sigma_max, sigma_min = build_context(
        case_dir / "run_args.json", stress
    )
    # Read the immutable accepted legacy capsule directly. The production
    # loader intentionally rejects predecessor source hashes after code
    # changes; this audit performs no restart or mutation.
    data = np.load(case_dir / "checkpoint_latest.npz", allow_pickle=True)
    metadata = driver._restore_nonfinite_tags(json.loads(str(data["metadata_json"].item())))
    mesh.nodes[:] = np.asarray(data["mesh_nodes"], float)
    restored = {
        "ep_gp": np.asarray(data["ep_gp"], float),
        "rho_gp": np.asarray(data["rho_gp"], float),
        "epsp_acc_gp": np.asarray(data["epsp_acc_gp"], float),
        "u": np.asarray(data["u"], float),
        "root_xy": np.asarray(data["root_xy"], float),
        "cycles": float(metadata["cycles"]),
        "shared_root_marked_clock_capsule": metadata.get("shared_root_marked_clock_capsule"),
    }
    capsule = restored.get("shared_root_marked_clock_capsule")
    if capsule is None:
        raise RuntimeError(f"{case_dir}: no shared-root capsule")
    summary = json.loads((case_dir / "summary.json").read_text())
    root = np.asarray(restored["root_xy"], float)
    radius = float(summary["root_radius_initial_m"])
    clock = _clock(args, capsule, root, radius)
    tensors, hist, root_node = _phase_tensors(
        args, mesh, cached, restored, sigma_max, sigma_min, root
    )
    drives = [clock.mpz.resolve_root_tensor(t) for t in tensors]
    nominal_open = np.asarray([d["opening_stress_Pa"] for d in drives])
    effective = np.asarray([clock.effective_opening_stress_Pa(s) for s in nominal_open])
    raw = np.exp(np.asarray([clock.mpz.cleavage_log_rate_s(s, args.T) for s in effective]))
    fem_cycle = cached.representative_cycle(
        restored["ep_gp"], restored["rho_gp"],
        *cached.affine(restored["ep_gp"], sigma_max, sigma_min, restored["u"])[0:2],
        args.T, args.frequency_Hz, args.plastic_n_phase, chain,
        cached.affine(restored["ep_gp"], sigma_max, sigma_min, restored["u"])[2],
        args.k_store, args.k_dyn, args.rho_floor, args.rho_cap,
        args.max_dep_phase, args.max_rho_rel_phase,
    )
    rates = {}
    for m in (1, 2, 3):
        eff = np.asarray([canonical_effective_cleavage_rate(x, m, 1e-6) for x in raw])
        action = float(np.mean(eff) / args.frequency_Hz)
        rates[m] = {"phase_rate_max_s": float(np.max(eff)), "per_cycle_action": action}
        rates[m].update({q: a / action if action > 0 else math.inf for q, a in QUANTILE_ACTION.items()})
    averaged = clock.copy()
    before = averaged.mpz.summary()
    avg_detail = averaged._advance_phase_block_exact(1.0, args.frequency_Hz, args.T, tensors)
    ordered, ordered_action = _ordered_one_cycle(clock, tensors, args.frequency_Hz, args.T)
    sa, so = averaged.mpz.summary(), ordered.mpz.summary()
    return {
        "sigma_a_MPa": stress, "R": args.R,
        "sigma_min_MPa": sigma_min / 1e6, "sigma_max_MPa": sigma_max / 1e6,
        "cycles_checkpoint": restored["cycles"],
        "FEM_Kt": float(summary["fem_sigma1_cycle_final_max_Pa"] / sigma_max),
        "root_sigma1_max_MPa": float(summary["fem_sigma1_cycle_final_max_Pa"] / 1e6),
        "exact_root_node_sigma1_max_MPa": float(np.max(hist["s1_node"][:, root_node]) / 1e6),
        "root_opening_max_MPa": float(np.max(nominal_open) / 1e6),
        "effective_cleavage_max_MPa": float(np.max(effective) / 1e6),
        "raw_cleavage_rate_max_s": float(np.max(raw)),
        "rates": rates,
        "epsp_acc_max": float(np.max(restored["epsp_acc_gp"])),
        "fem_dep_eq_cycle_max": float(np.max(fem_cycle["dep_eq_cycle"])),
        "fem_drho_cycle_max_m2": float(np.max(np.abs(fem_cycle["drho_cycle"]))),
        "fem_scalar_rate_maxima": {
            key: float(np.max(value)) for key, value in fem_cycle.items()
            if key.startswith("mu_")
        },
        "mpz": {
            key: before[key] for key in (
                "emitted_total", "tip_radius_m",
                "signed_active_K_shield_Pa_sqrt_m",
                "sigma_back_by_system_Pa", "rho_back_by_system_m2",
            )
        },
        "phase_order_comparison": {
            "averaged_action_m3": float(avg_detail["action_increment"]),
            "ordered_action_m3": float(ordered_action),
            "relative_action_difference": float((ordered_action - avg_detail["action_increment"]) / max(avg_detail["action_increment"], 1e-300)),
            "averaged_emitted_increment": float(sa["emitted_total"] - before["emitted_total"]),
            "ordered_emitted_increment": float(so["emitted_total"] - before["emitted_total"]),
            "averaged_shielding": float(sa["signed_active_K_shield_Pa_sqrt_m"]),
            "ordered_shielding": float(so["signed_active_K_shield_Pa_sqrt_m"]),
            "averaged_backstress_max_Pa": float(max(sa["sigma_back_by_system_Pa"])),
            "ordered_backstress_max_Pa": float(max(so["sigma_back_by_system_Pa"])),
        },
        "direct_summary_H_attempt": summary.get("H_attempt_final"),
        "mesh_root_radius_um": radius * 1e6,
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, required=True)
    ns = p.parse_args()
    rows = [audit_condition(s, pth) for s, pth in DEFAULTS.items()]
    ns.out.mkdir(parents=True, exist_ok=True)
    serial = driver._json_safe(rows)
    (ns.out / "stress_scale_audit.json").write_text(json.dumps(serial, indent=2, sort_keys=True) + "\n")
    flat = []
    for row in rows:
        for m in (1, 2, 3):
            flat.append({k: v for k, v in row.items() if not isinstance(v, dict)} | {"m_hits": m} | row["rates"][m])
    with (ns.out / "stress_scale_rates.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(flat[0])); w.writeheader(); w.writerows(flat)
    print(json.dumps(serial, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
