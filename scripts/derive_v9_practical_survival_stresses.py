#!/usr/bin/env python3
"""Bracketed log-hazard interpolation of practical 1e14 survival stresses."""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from scripts.run_v9_quiet_tail_kernel import atomic_merge_csv


FIELDS = [
    "material_class", "target_S_stable_birth_1e14", "target_H_stable_birth",
    "estimated_sigma_a_MPa", "lower_sigma_a_MPa", "lower_H_stable_birth",
    "upper_sigma_a_MPa", "upper_H_stable_birth", "interpolation",
    "estimate_semantics",
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with (args.output / "300K_stable_birth_survival.csv").open(newline="") as stream:
        source = list(csv.DictReader(stream))
    points = {}
    for row in source:
        if float(row["N"]) != 1e14:
            continue
        key = (row["material_class"], float(row["sigma_a_MPa"]))
        stable_H = -math.log(max(float(row["S_stable_birth"]), 1e-300))
        # A numerically underflowed survival supplies only H_cleave when every
        # evaluated mark is admitted, as recorded by this campaign.
        if float(row["S_stable_birth"]) == 0.0:
            stable_H = float(row["H_cleave"])
        points[key] = stable_H
    rows = []
    for material in sorted({key[0] for key in points}):
        ordered = sorted((stress, H) for (name, stress), H in points.items()
                         if name == material)
        for survival in (0.9, 0.5, 0.1):
            target = -math.log(survival)
            bracket = next(((a, b) for a, b in zip(ordered, ordered[1:])
                            if a[1] <= target <= b[1]), None)
            if bracket is None:
                continue
            (s0, h0), (s1, h1) = bracket
            fraction = ((math.log(target) - math.log(h0)) /
                        (math.log(h1) - math.log(h0)))
            estimate = s0 + fraction * (s1 - s0)
            rows.append({
                "material_class": material,
                "target_S_stable_birth_1e14": survival,
                "target_H_stable_birth": target,
                "estimated_sigma_a_MPa": estimate,
                "lower_sigma_a_MPa": s0, "lower_H_stable_birth": h0,
                "upper_sigma_a_MPa": s1, "upper_H_stable_birth": h1,
                "interpolation": "secant_linear_sigma_vs_log_stable_birth_action",
                "estimate_semantics": "practical_1e14_survival_stress_not_infinite_endurance",
            })
    atomic_merge_csv(
        args.output / "300K_practical_survival_stress_estimates.csv", FIELDS,
        rows, ("material_class", "target_S_stable_birth_1e14"))


if __name__ == "__main__":
    main()
