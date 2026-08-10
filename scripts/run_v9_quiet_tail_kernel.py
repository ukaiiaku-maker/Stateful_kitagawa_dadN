#!/usr/bin/env python3
"""Sequential one-class v9 quiet-tail kernel production runner."""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import tempfile

import numpy as np

from arrhenius_fracture.v9_canonical_four_class_fem import CanonicalFourClassFEMCondition
from arrhenius_fracture.v9_canonical_four_class_elastic_fem import (
    CanonicalElasticFourClassFEMCondition,
)
from arrhenius_fracture.v9_quiet_tail_kernel import (
    energy_gate_envelope, kernel_convergence, phase_resolved_cycle,
    restore_condition_checkpoint, save_condition_checkpoint, stationary_marked_renewal,
)


CLASSES = {
    "Peak": ("v913_paper_peak01_0242980_persistent_sites", 735.9214951373579),
    "DBTT": ("v913_paper_dbtt01_0202500_persistent_sites", 700.0),
    "weak-T": ("v913_paper_weakT01_0129902_persistent_sites", 350.0),
    "ceramic": ("v913_paper_ceramic01_0077080_persistent_sites", 400.0),
}
STATE_FIELDS = [
    "material_class", "option_id", "sigma_a_MPa", "N", "H_cleave",
    "cycle_hazard", "signed_state_relative_error", "root_tensor_relative_error",
    "backstress_relative_error", "shielding_relative_error", "blunting_relative_error",
    "phase_log_rate_absolute_error", "cycle_hazard_relative_error", "stationary_periodic_orbit",
    "accepted_blocks", "rejected_blocks", "checkpoint_generation",
]
GATE_FIELDS = ["material_class", "option_id", "sigma_a_MPa", "N", "phase_index", "phase", "Xi",
               "proposed_length_m", "maximum_released_energy_J_per_m", "resistance_J_per_m2",
               "minimum_energy_gate_margin_J_per_m", "admitted_length_m", "reason", "mesh_resolved"]
SURVIVAL_FIELDS = ["material_class", "option_id", "sigma_a_MPa", "N", "H_cleave",
                   "S_no_attempt", "P_at_least_one_attempt", "S_stable_birth",
                   "P_stable_birth", "expected_nonpropagating_attempt_count",
                   "conditional_attempt_admission_probability", "estimator"]
DESCENT_FIELDS = [
    "material_class", "option_id", "sigma_a_MPa", "N", "H_cleave",
    "S_no_attempt", "S_stable_birth", "conditional_attempt_admission_probability",
    "rejected_phase_xi_fraction", "minimum_admitted_length_m",
    "maximum_admitted_length_m", "minimum_energy_gate_margin_J_per_m",
]


def atomic_merge_csv(path, fields, new_rows, key_fields):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    def key_value(value):
        try:
            return ("number", format(float(value), ".17g"))
        except (TypeError, ValueError):
            return ("text", str(value))
    def row_key(row):
        return tuple(key_value(row[k]) for k in key_fields)
    rows = []
    if path.is_file():
        with path.open(newline="") as stream:
            rows = list(csv.DictReader(stream))
    keys = {row_key(row) for row in new_rows}
    rows = [row for row in rows if row_key(row) not in keys]
    rows.extend(new_rows)
    fd, name = tempfile.mkstemp(prefix="."+path.name+".", dir=path.parent)
    try:
        with os.fdopen(fd, "w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
            writer.writeheader(); writer.writerows(rows); stream.flush(); os.fsync(stream.fileno())
        os.replace(name, path)
    finally:
        if os.path.exists(name): os.unlink(name)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("material_class", choices=list(CLASSES))
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--run-args", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("runs/sn_v9_canonical_four_class"))
    parser.add_argument("--targets", default="1e6,1e7,1e8")
    parser.add_argument("--hazard-seed", type=int, default=1720)
    parser.add_argument("--stress-MPa", type=float,
                        help="override the class anchor for an adaptive descent condition")
    parser.add_argument("--resume-checkpoint", type=Path)
    parser.add_argument("--canonical-elastic-v3", action="store_true",
                        help="use elastic FEM bulk plus signed tip MPZ canonical v3")
    args = parser.parse_args()
    option, anchor_stress = CLASSES[args.material_class]
    stress = anchor_stress if args.stress_MPa is None else float(args.stress_MPa)
    condition_type = (CanonicalElasticFourClassFEMCondition
                      if args.canonical_elastic_v3 else CanonicalFourClassFEMCondition)
    condition = condition_type.from_run_args(
        args.run_args, option, args.source_root, stress, args.hazard_seed)
    condition.birth.hazard_threshold_action = 1.0e100
    targets = [float(x) for x in args.targets.split(",")]
    xi = np.asarray([.1, .3, .49, .5, .6, .8, 1., 1.25, 1.5, 2., 2.5, 3., 3.5, 3.99, 4., 8.])
    state_rows, gate_rows = [], []
    kernels = []; hazards = []
    previous = None
    if args.resume_checkpoint is not None:
        previous, _resume_summary, _manifest = restore_condition_checkpoint(
            condition, args.resume_checkpoint)
        condition.birth.hazard_threshold_action = 1.0e100
    checkpoint_class = args.material_class
    if args.stress_MPa is not None:
        checkpoint_class = f"{args.material_class}/descent/sigmaA_{stress:.12g}MPa"
    for target in targets:
        condition.advance(target-condition.cycles)
        kernel = phase_resolved_cycle(condition)
        convergence = {} if previous is None else kernel_convergence(previous, kernel)
        # Fixed, predeclared qualification tolerances.  These are reporting
        # bounds, not adjustable fitting parameters.
        stationary = bool(convergence and
            convergence["signed_state_relative_error"] <= 5e-4 and
            convergence["root_tensor_relative_error"] <= 5e-4 and
            convergence["backstress_relative_error"] <= 5e-4 and
            convergence["shielding_relative_error"] <= 5e-4 and
            convergence["blunting_relative_error"] <= 5e-4 and
            convergence["phase_log_rate_absolute_error"] <= 5e-3 and
            convergence["cycle_hazard_relative_error"] <= 5e-3)
        summary = {"material_class": args.material_class, "option_id": option,
                   "sigma_a_MPa": stress, "N": target,
                   "H_cleave": condition.birth.cumulative_cleavage_hazard,
                   "cycle_hazard": kernel["cycle_hazard"], **convergence,
                   "stationary_periodic_orbit": stationary}
        generation = save_condition_checkpoint(
            condition, args.output/"quiet_tail_checkpoints"/checkpoint_class/f"N_{target:.17g}",
            kernel, summary)
        state_rows.append(summary | {"accepted_blocks": condition.accepted_blocks,
                          "rejected_blocks": condition.rejected_blocks,
                          "checkpoint_generation": generation})
        envelope = energy_gate_envelope(condition, kernel, xi)
        for row in envelope:
            gate_rows.append({"material_class": args.material_class, "option_id": option,
                              "sigma_a_MPa": stress, "N": target, **row})
        kernels.append(kernel); hazards.append(condition.birth.cumulative_cleavage_hazard)
        previous = kernel
        atomic_merge_csv(args.output/"300K_quiet_tail_state.csv", STATE_FIELDS,
                         state_rows, ("material_class", "sigma_a_MPa", "N"))
        atomic_merge_csv(args.output/"300K_energy_gate_envelope.csv", GATE_FIELDS,
                         gate_rows, ("material_class", "sigma_a_MPa", "N", "phase_index", "Xi"))

    stable_exponent = 0.0; survival_rows = []
    last_h = 0.0
    for target, H, kernel in zip(targets, hazards, kernels):
        subset = [row for row in gate_rows if float(row["N"]) == target]
        renewal = stationary_marked_renewal(H-last_h, kernel["cycle_hazard"],
                                            kernel["cleavage_action_within_cycle"], subset)
        p = renewal["stationary_attempt_admission_probability"]
        stable_exponent += (H-last_h)*p
        stable_survival = math.exp(-stable_exponent)
        survival_rows.append({"material_class": args.material_class, "option_id": option,
            "sigma_a_MPa": stress, "N": target, "H_cleave": H,
            "S_no_attempt": math.exp(-H), "P_at_least_one_attempt": -math.expm1(-H),
            "S_stable_birth": stable_survival, "P_stable_birth": -math.expm1(-stable_exponent),
            "expected_nonpropagating_attempt_count": max(H-stable_exponent, 0.0),
            "conditional_attempt_admission_probability": p,
            "estimator": "deterministic_piecewise_stationary_marked_renewal_quadrature"})
        last_h = H

    final_stationary = bool(state_rows[-1]["stationary_periodic_orbit"])
    if final_stationary:
        final_kernel = kernels[-1]
        subset = [row for row in gate_rows if float(row["N"]) == targets[-1]]
        p = stationary_marked_renewal(1.0, final_kernel["cycle_hazard"],
                                     final_kernel["cleavage_action_within_cycle"], subset)[
                                         "stationary_attempt_admission_probability"]
        for horizon in (1e10, 1e12, 1e14):
            H = hazards[-1] + max(horizon-targets[-1], 0.0)*final_kernel["cycle_hazard"]
            exponent = stable_exponent + max(H-hazards[-1], 0.0)*p
            survival_rows.append({"material_class": args.material_class, "option_id": option,
                "sigma_a_MPa": stress, "N": horizon, "H_cleave": H,
                "S_no_attempt": math.exp(-H), "P_at_least_one_attempt": -math.expm1(-H),
                "S_stable_birth": math.exp(-exponent), "P_stable_birth": -math.expm1(-exponent),
                "expected_nonpropagating_attempt_count": max(H-exponent, 0.0),
                "conditional_attempt_admission_probability": p,
                "estimator": "qualified_stationary_marked_renewal_quadrature"})
    atomic_merge_csv(args.output/"300K_stable_birth_survival.csv", SURVIVAL_FIELDS,
                     survival_rows, ("material_class", "sigma_a_MPa", "N"))
    if args.canonical_elastic_v3 and args.stress_MPa is not None:
        descent_rows = []
        survival_by_N = {float(row["N"]): row for row in survival_rows}
        for target in targets:
            subset = [row for row in gate_rows if float(row["N"]) == target]
            admitted = np.asarray([float(row["admitted_length_m"]) for row in subset])
            survival = survival_by_N[target]
            descent_rows.append({
                "material_class": args.material_class, "option_id": option,
                "sigma_a_MPa": stress, "N": target,
                "H_cleave": survival["H_cleave"],
                "S_no_attempt": survival["S_no_attempt"],
                "S_stable_birth": survival["S_stable_birth"],
                "conditional_attempt_admission_probability":
                    survival["conditional_attempt_admission_probability"],
                "rejected_phase_xi_fraction": float(np.mean(admitted <= 0.0)),
                "minimum_admitted_length_m": float(np.min(admitted)),
                "maximum_admitted_length_m": float(np.max(admitted)),
                "minimum_energy_gate_margin_J_per_m": min(
                    float(row["minimum_energy_gate_margin_J_per_m"])
                    for row in subset),
            })
        atomic_merge_csv(
            args.output/"300K_endurance_stress_descent.csv", DESCENT_FIELDS,
            descent_rows, ("material_class", "sigma_a_MPa", "N"))
    result = {"material_class": args.material_class, "stationary": final_stationary,
              "last_cycle_hazard": kernels[-1]["cycle_hazard"],
              "last_H": hazards[-1], "survival_rows": survival_rows}
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
