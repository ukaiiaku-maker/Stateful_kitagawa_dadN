#!/usr/bin/env python3
"""Rebuild descent summaries atomically from completed physical result tables."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from scripts.run_v9_quiet_tail_kernel import (
    CLASSES, DESCENT_FIELDS, GATE_FIELDS, STATE_FIELDS, SURVIVAL_FIELDS,
    atomic_merge_csv,
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with (args.output / "300K_stable_birth_survival.csv").open(newline="") as stream:
        survival = list(csv.DictReader(stream))
    with (args.output / "300K_energy_gate_envelope.csv").open(newline="") as stream:
        gates = list(csv.DictReader(stream))
    with (args.output / "300K_quiet_tail_state.csv").open(newline="") as stream:
        states = list(csv.DictReader(stream))
    survival = list({
        (row["material_class"], float(row["sigma_a_MPa"]), float(row["N"])): row
        for row in survival
    }.values())
    gates = list({
        (row["material_class"], float(row["sigma_a_MPa"]), float(row["N"]),
         int(row["phase_index"]), float(row["Xi"])): row
        for row in gates
    }.values())
    states = list({
        (row["material_class"], float(row["sigma_a_MPa"]), float(row["N"])): row
        for row in states
    }.values())
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
        key = (row["material_class"], row["sigma_a_MPa"], row["N"])
        gate_groups.setdefault(key, []).append(row)
    rows_by_key = {}
    anchors = {name: stress for name, (_option, stress) in CLASSES.items()}
    for result in survival:
        material = result["material_class"]
        stress = float(result["sigma_a_MPa"])
        if stress == anchors[material]:
            continue
        subset = gate_groups[(material, result["sigma_a_MPa"], result["N"])]
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
