#!/usr/bin/env python3
"""Compare two hash-verified v9 atomic generations at a common boundary."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from arrhenius_fracture.v9_array_codec import AtomicArrayGenerationStore


EXACT_PREFIXES = (
    "pd__site_", "pd_v9__site_", "pd__bond_damage", "pd__active_front",
    "pd__front_", "pd__primary_seed", "mesh_nodes", "root_xy",
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--left", type=Path, required=True)
    parser.add_argument("--right", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    left, lm, ls, _ = AtomicArrayGenerationStore(args.left).load()
    right, rm, rs, _ = AtomicArrayGenerationStore(args.right).load()
    if set(left) != set(right):
        raise SystemExit("array schemas differ")
    rows = {}
    exact_ok = True
    max_abs = 0.0
    max_rel = 0.0
    for name in sorted(left):
        a, b = np.asarray(left[name]), np.asarray(right[name])
        exact = bool(np.array_equal(a, b, equal_nan=True))
        if np.issubdtype(a.dtype, np.number):
            finite = np.isfinite(a) & np.isfinite(b)
            absolute = float(np.max(np.abs(a[finite] - b[finite]))) if np.any(finite) else 0.0
            scale = np.maximum.reduce([np.abs(a[finite]), np.abs(b[finite]), np.ones(np.count_nonzero(finite))])
            relative = float(np.max(np.abs(a[finite] - b[finite]) / scale)) if np.any(finite) else 0.0
        else:
            absolute = relative = 0.0
        required_exact = name.startswith(EXACT_PREFIXES)
        exact_ok = exact_ok and (exact or not required_exact)
        max_abs, max_rel = max(max_abs, absolute), max(max_rel, relative)
        rows[name] = {"exact": exact, "required_exact": required_exact,
                      "max_abs": absolute, "max_relative": relative}
    report = {
        "schema": "V9_PD_ATOMIC_PARTITION_RESTART_COMPARISON_1",
        "left": str(args.left), "right": str(args.right),
        "left_cycles": lm["cycles"], "right_cycles": rm["cycles"],
        "summary_cycles_equal": float(ls["cycles"]) == float(rs["cycles"]),
        "exact_protected_state": exact_ok,
        "maximum_absolute_difference": max_abs,
        "maximum_relative_difference": max_rel,
        "arrays": rows,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: report[k] for k in (
        "left_cycles", "right_cycles", "exact_protected_state",
        "maximum_absolute_difference", "maximum_relative_difference")}, indent=2))


if __name__ == "__main__":
    main()
