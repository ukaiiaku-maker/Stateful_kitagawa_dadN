#!/usr/bin/env python3
"""Rebuild descent summaries atomically from completed physical result tables."""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np

from scripts.run_v9_quiet_tail_kernel import (
    CLASSES, DESCENT_FIELDS, GATE_FIELDS, STATE_FIELDS, SURVIVAL_FIELDS,
    atomic_merge_csv,
)


def reconstruct_direct_checkpoint_survival(states, gates):
    """Make direct checkpoint H authoritative and fail closed on partial gates."""
    gate_groups = {}
    for row in gates:
        key = (row["material_class"], float(row["sigma_a_MPa"]), float(row["N"]))
        gate_groups.setdefault(key, []).append(row)
    condition_groups = {}
    for state in states:
        condition_groups.setdefault(
            (state["material_class"], float(state["sigma_a_MPa"])), []).append(state)
    reconstructed = []
    for (_material, _stress), condition_states in condition_groups.items():
        stable_action = 0.0
        previous_H = 0.0
        for state in sorted(condition_states, key=lambda row: float(row["N"])):
            key = (state["material_class"], float(state["sigma_a_MPa"]), float(state["N"]))
            subset = gate_groups.get(key, [])
            if not subset:
                raise RuntimeError(f"missing gate envelope for direct checkpoint {key}")
            admitted = [float(row["admitted_length_m"]) > 0.0 for row in subset]
            if all(admitted):
                probability = 1.0
            elif not any(admitted):
                probability = 0.0
            else:
                raise RuntimeError(
                    f"partial gate at {key} requires checkpoint-kernel renewal reconstruction")
            H = float(state["H_cleave"])
            if H < previous_H:
                raise RuntimeError(f"nonmonotone direct checkpoint hazard for {key}")
            stable_action += (H - previous_H) * probability
            previous_H = H
            reconstructed.append({
                "material_class": state["material_class"],
                "option_id": state["option_id"],
                "sigma_a_MPa": float(state["sigma_a_MPa"]),
                "N": float(state["N"]), "H_cleave": H,
                "S_no_attempt": math.exp(-H),
                "P_at_least_one_attempt": -math.expm1(-H),
                "S_stable_birth": math.exp(-stable_action),
                "P_stable_birth": -math.expm1(-stable_action),
                "expected_nonpropagating_attempt_count": max(H-stable_action, 0.0),
                "conditional_attempt_admission_probability": probability,
                "estimator": "direct_checkpoint_H_plus_evaluated_gate_admission",
            })
    return reconstructed


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with (args.output / "300K_energy_gate_envelope.csv").open(newline="") as stream:
        gates = list(csv.DictReader(stream))
    with (args.output / "300K_quiet_tail_state.csv").open(newline="") as stream:
        states = list(csv.DictReader(stream))
    gates = list({
        (row["material_class"], float(row["sigma_a_MPa"]), float(row["N"]),
         int(row["phase_index"]), float(row["Xi"])): row
        for row in gates
    }.values())
    states = list({
        (row["material_class"], float(row["sigma_a_MPa"]), float(row["N"])): row
        for row in states
    }.values())
    survival = reconstruct_direct_checkpoint_survival(states, gates)
    atomic_merge_csv(
        args.output / "300K_stable_birth_survival.csv", SURVIVAL_FIELDS,
        survival, ("material_class", "sigma_a_MPa", "N"))
    atomic_merge_csv(
        args.output / "300K_energy_gate_envelope.csv", GATE_FIELDS, gates,
        ("material_class", "sigma_a_MPa", "N", "phase_index", "Xi"))
    atomic_merge_csv(
        args.output / "300K_quiet_tail_state.csv", STATE_FIELDS, states,
        ("material_class", "sigma_a_MPa", "N"))
    gate_groups = {}
    for row in gates:
        key = (row["material_class"], float(row["sigma_a_MPa"]), float(row["N"]))
        gate_groups.setdefault(key, []).append(row)
    rows_by_key = {}
    anchors = {name: stress for name, (_option, stress) in CLASSES.items()}
    for result in survival:
        material = result["material_class"]
        stress = float(result["sigma_a_MPa"])
        if stress == anchors[material]:
            continue
        subset = gate_groups[(material, stress, float(result["N"]))]
        admitted = np.asarray([float(row["admitted_length_m"]) for row in subset])
        row = {
            "material_class": material, "option_id": result["option_id"],
            "sigma_a_MPa": stress, "N": float(result["N"]),
            "H_cleave": float(result["H_cleave"]),
            "S_no_attempt": float(result["S_no_attempt"]),
            "S_stable_birth": float(result["S_stable_birth"]),
            "conditional_attempt_admission_probability": float(
                result["conditional_attempt_admission_probability"]),
            "rejected_phase_xi_fraction": float(np.mean(admitted <= 0.0)),
            "minimum_admitted_length_m": float(np.min(admitted)),
            "maximum_admitted_length_m": float(np.max(admitted)),
            "minimum_energy_gate_margin_J_per_m": min(
                float(row["minimum_energy_gate_margin_J_per_m"])
                for row in subset),
        }
        rows_by_key[(material, stress, float(result["N"]))] = row
    atomic_merge_csv(
        args.output / "300K_endurance_stress_descent.csv", DESCENT_FIELDS,
        list(rows_by_key.values()), ("material_class", "sigma_a_MPa", "N"))


if __name__ == "__main__":
    main()
