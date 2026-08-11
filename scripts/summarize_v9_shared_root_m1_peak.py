#!/usr/bin/env python3
"""Build compact, provenance-preserving summaries of the Peak m=1 runs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np


CASES = {
    "m1_2000_seed42_coarse_1e5_baseline": "Peak_2000/Peak/shielded/sigmaA_2000MPa",
    "m1_2000_seed42_rate_separated_v2": "rate_separated_restart_from_fine_stable_v4/Peak/shielded/sigmaA_2000MPa",
    "m1_1500_seed42": "Peak_1500_VHCF/Peak/shielded/sigmaA_1500MPa",
    "m1_1500_seed43": "Peak_1500_seed43/Peak/shielded/sigmaA_1500MPa",
    "m1_1500_seed44": "Peak_1500_seed44/Peak/shielded/sigmaA_1500MPa",
    "m1_1500_seed45": "Peak_1500_seed45/Peak/shielded/sigmaA_1500MPa",
    "m1_1500_seed46_rate_separated": "Peak_1500_seed46_rate_separated/Peak/shielded/sigmaA_1500MPa",
}

CANONICAL_BY_STRESS_SEED = {
    (2000.0, 42017, 42018): "m1_2000_seed42_rate_separated_v2",
}


def validate_condition_registry(entries):
    seen = {}
    for row in entries:
        if not row["production_analyzer_eligible"]:
            continue
        key = (float(row["sigma_a_MPa"]), int(row["hazard_seed"]), int(row["mark_seed"]))
        if key in seen:
            raise ValueError(f"ambiguous canonical stress/seed rows: {key}: {seen[key]}, {row['condition']}")
        seen[key] = row["condition"]
    return seen


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def total(value) -> float:
    return float(np.asarray(value, dtype=float).sum())


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("runs/sn_v9_shared_root_m1_peak"))
    parser.add_argument("--out", type=Path)
    args = parser.parse_args(argv)
    out = args.out or args.root / "qualification"
    out.mkdir(parents=True, exist_ok=True)

    rows, registry_rows = [], []
    manifest = {"schema": "V9_SHARED_ROOT_M1_PEAK_SUMMARY_2", "cases": {}}
    for label, rel in CASES.items():
        case = args.root / rel
        summary_path = case / "summary.json"
        summary = json.loads(summary_path.read_text())
        run_args = json.loads((case / "run_args.json").read_text())
        hazard_seed = int(run_args["global_cleavage_seed"])
        mark_seed = int(run_args["spatial_mark_seed"])
        mpz = summary["signed_mpz_state_final"]
        slip_pos = total(mpz["accumulated_slip_positive"])
        slip_neg = total(mpz["accumulated_slip_negative"])
        unsigned_slip = slip_pos + slip_neg
        # A deliberately conservative dimensional projection, not a new
        # constitutive coupling: b * line-content / the configured 2 um source zone.
        equivalent_strain_bound = 2.74e-10 * unsigned_slip / 2.0e-6
        rows.append({
            "condition": label,
            "sigma_a_MPa": summary["sigma_a_MPa"],
            "cycles_total": summary["cycles_total"],
            "N_attempt": summary["cycles_first_embryo"],
            "N_stable_seed": summary["cycles_first_stable"],
            "N_first_softening": summary["cycles_first_softening"],
            "N_root_connection": summary["cycles_root_connected"],
            "N_front_capture": summary["cycles_front_capture"],
            "physical_handoff": summary["pd_handoff_pass_final"],
            "H_attempt": summary["H_attempt_final"],
            "S_no_attempt": summary["S_no_attempt_final"],
            "attempt_count": summary["global_attempt_count_final"],
            "selected_site_id": (summary.get("last_marked_attempt") or {}).get("site_id"),
            "emitted_line_content": mpz["emitted_total"],
            "accumulated_slip_positive": slip_pos,
            "accumulated_slip_negative": slip_neg,
            "equivalent_root_plastic_strain_upper_bound": equivalent_strain_bound,
            "relaxation_upper_bound_Pa_at_E410GPa": 410e9 * equivalent_strain_bound,
            "max_MPZ_backstress_Pa": max(mpz["sigma_back_by_system_Pa"]),
            "signed_shielding_Pa_sqrt_m": mpz["signed_active_K_shield_Pa_sqrt_m"],
            "continuum_epsp_max": summary["epsp_acc_final_max"],
            "status": summary["status"],
        })
        manifest["cases"][label] = {
            "path": str(case),
            "summary_sha256": sha256(summary_path),
            "checkpoint_sha256": sha256(case / "checkpoint_latest.npz"),
            "model_id": summary["cleavage_clock_model"],
            "restart_provenance": "hash-verified atomic v9 generation store",
        }
        canonical_label = CANONICAL_BY_STRESS_SEED.get((float(summary["sigma_a_MPa"]), hazard_seed, mark_seed))
        eligible = canonical_label in (None, label)
        registry_rows.append({
            "condition": label,
            "sigma_a_MPa": float(summary["sigma_a_MPa"]),
            "hazard_seed": hazard_seed,
            "mark_seed": mark_seed,
            "numerical_protocol": ("coarse_1e5_baseline" if label.endswith("coarse_1e5_baseline")
                                   else "rate_separated_v2" if label.endswith("rate_separated_v2")
                                   else "historical_mixed_protocol_diagnostic"),
            "production_analyzer_eligible": eligible,
            "canonical_condition": canonical_label,
            "summary_path": str(summary_path),
            "summary_sha256": sha256(summary_path),
            "checkpoint_sha256": sha256(case / "checkpoint_latest.npz"),
        })

    fields = list(rows[0])
    with (out / "endpoint_results.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    seed_rows = [row for row in rows if row["sigma_a_MPa"] == 1500.0]
    n = len(seed_rows)
    stabilized = sum(row["N_stable_seed"] is not None for row in seed_rows)
    captured = sum(row["N_front_capture"] is not None for row in seed_rows)
    ensemble = {
        "schema": "V9_SHARED_ROOT_M1_PEAK_MARK_TRANSITION_ENSEMBLE_1",
        "stress_MPa": 1500.0,
        "observation_horizon_cycles_by_seed": {row["condition"]:row["cycles_total"] for row in seed_rows},
        "post_stable_exposure_cycles_by_seed": {row["condition"]:row["cycles_total"]-row["N_stable_seed"] for row in seed_rows},
        "n": n,
        "stabilized": stabilized,
        "front_captured": captured,
        "P_stabilize_given_realized_attempt_empirical": stabilized / n,
        "P_front_capture_given_stabilized_seed_empirical": captured / stabilized,
        "stable_front_survival_empirical": 1.0 - captured / n,
        "qualification": "small diagnostic ensemble with unequal censor horizons; not a converged survival estimate",
        "attempt_survival_semantics": "S_no_attempt=exp(-H_attempt), distinct from stable-front survival",
    }
    (out / "seed_ensemble.json").write_text(json.dumps(ensemble, indent=2) + "\n")

    validate_condition_registry(registry_rows)
    registry = {
        "schema": "V9_SHARED_ROOT_M1_PEAK_CONDITION_REGISTRY_1",
        "uniqueness_key": ["sigma_a_MPa", "hazard_seed", "mark_seed"],
        "selection_rule": "only production_analyzer_eligible rows may enter production S-N analysis",
        "conditions": registry_rows,
    }
    (out / "condition_registry.json").write_text(json.dumps(registry, indent=2) + "\n")

    control = args.root / CASES["m1_2000_seed42_coarse_1e5_baseline"]
    comparisons = {}
    for label, path in {
        "caller_partition_10x": args.root / "partition_2000_block1e4_v2/Peak/shielded/sigmaA_2000MPa",
        "persisted_restart_before_mark": args.root / "restart_2000_across_marked_crossing_v2/Peak/shielded/sigmaA_2000MPa",
    }.items():
        a = json.loads((control / "summary.json").read_text())
        b = json.loads((path / "summary.json").read_text())
        endpoint_keys = ["cycles_first_embryo", "cycles_first_stable",
                         "cycles_first_softening", "cycles_root_connected",
                         "cycles_front_capture"]
        deltas = {key: float(b[key] - a[key]) for key in endpoint_keys}
        with np.load(control / "pd_state_final.npz") as za, np.load(path / "pd_state_final.npz") as zb:
            arrays_exact = all(np.array_equal(za[key], zb[key], equal_nan=True)
                               for key in za.files)
        comparisons[label] = {
            "comparison_path": str(path),
            "endpoint_cycle_deltas": deltas,
            "selected_site_identity_equal": (
                a["last_marked_attempt"]["site_id"] == b["last_marked_attempt"]["site_id"]
            ),
            "attempt_action_equal": a["H_attempt_final"] == b["H_attempt_final"],
            "threshold_equal": a["global_threshold_action_final"] == b["global_threshold_action_final"],
            "final_PD_arrays_byte_equal": arrays_exact,
        }
    (out / "trajectory_equivalence.json").write_text(
        json.dumps({"schema": "V9_SHARED_ROOT_M1_TRAJECTORY_EQUIVALENCE_1",
                    "comparisons": comparisons}, indent=2) + "\n"
    )
    manifest["generated"] = {
        "endpoint_results.csv": sha256(out / "endpoint_results.csv"),
        "seed_ensemble.json": sha256(out / "seed_ensemble.json"),
        "trajectory_equivalence.json": sha256(out / "trajectory_equivalence.json"),
        "condition_registry.json": sha256(out / "condition_registry.json"),
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
