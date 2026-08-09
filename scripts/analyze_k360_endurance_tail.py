#!/usr/bin/env python3
"""Reproduce the committed K360 v2.8 aggregate birth-hazard tail metrics."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np


DEFAULT_INPUT = Path("reference/k360_v2_8_tail/K360_CENSOR_DERIVED_TAIL.csv")


def analyze(path: Path, fractions=(0.50, 0.65, 0.80)) -> dict[str, object]:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    cycles = np.asarray([float(row["cycle_mid"]) for row in rows])
    hazard = np.asarray([float(row["effective_total_birth_hazard_per_cycle"]) for row in rows])
    cumulative = np.asarray([float(row["cumulative_expected_births"]) for row in rows])
    fits = []
    for fraction in fractions:
        start = int(math.floor(len(cycles) * fraction))
        x = np.log(cycles[start:])
        y = np.log(hazard[start:])
        slope, intercept = np.polyfit(x, y, 1)
        predicted = slope * x + intercept
        r_squared = 1.0 - float(np.sum((y - predicted) ** 2) / np.sum((y - np.mean(y)) ** 2))
        exponent = float(-slope)
        amplitude = float(np.exp(intercept))
        remaining = float(amplitude * cycles[-1] ** (1.0 - exponent) / (exponent - 1.0))
        fits.append({
            "fraction": fraction,
            "start_cycle": float(cycles[start]),
            "end_cycle": float(cycles[-1]),
            "exponent": exponent,
            "amplitude": amplitude,
            "r_squared": r_squared,
            "empirical_extrapolated_remaining_integrated_hazard": remaining,
        })
    conservative_remaining = max(fit["empirical_extrapolated_remaining_integrated_hazard"] for fit in fits)
    return {
        "row_count": len(rows),
        "cumulative_expected_births_end": float(cumulative[-1]),
        "nested_power_law_fits": fits,
        "naive_fit_based_survival_asymptote": math.exp(-(float(cumulative[-1]) + conservative_remaining)),
        "empirical_tail_evidence": "integrable_tail_observed",
        "physics_asymptotic_classification": "undetermined",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = analyze(args.input)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
