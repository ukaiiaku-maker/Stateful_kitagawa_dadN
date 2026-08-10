#!/usr/bin/env python3
"""Package one hash-verified Stateful-PD atomic tail generation."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-dir", type=Path, required=True)
    parser.add_argument("--generation", required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--expected-cycles", type=float, required=True)
    parser.add_argument("--sigma-a-MPa", type=float, required=True)
    args = parser.parse_args()

    case = args.case_dir.resolve()
    source_generation = case / "v9_generations" / args.generation
    if not source_generation.is_dir():
        raise SystemExit(f"missing generation: {source_generation}")
    summary = json.loads((source_generation / "summary.json").read_text())
    if abs(float(summary["cycles"]) - args.expected_cycles) > max(1e-9, 1e-12 * abs(args.expected_cycles)):
        raise SystemExit("generation cycle does not match --expected-cycles")

    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    (out / "v9_generations").mkdir(exist_ok=True)
    target_generation = out / "v9_generations" / args.generation
    if target_generation.exists():
        raise SystemExit(f"refusing to overwrite existing package generation: {target_generation}")
    shutil.copytree(source_generation, target_generation)

    context_names = (
        "run_args.json", "checkpoint_latest.npz",
        "summary.json", "sn_stateful_pd_history.csv",
        "crack_handoff_audit_final.json", "pd_state_final.npz",
    )
    for name in context_names:
        source = case / name
        if not source.is_file():
            raise SystemExit(f"missing required context file: {source}")
        shutil.copy2(source, out / name)
    optional_context_names = (
        "v9_pd_high_cycle_controller.json", "v9_pd_high_cycle_mode_history.json",
    )
    missing_optional_context = []
    for name in optional_context_names:
        source = case / name
        if source.is_file():
            shutil.copy2(source, out / name)
        else:
            missing_optional_context.append(name)
    (out / "ACTIVE.json").write_text(
        json.dumps({"generation": args.generation}, separators=(",", ":")) + "\n"
    )

    arrays = np.load(source_generation / "state_arrays.npz", allow_pickle=False)
    site_node = arrays["pd__site_node_index"].astype(int)
    threshold = arrays["pd__site_birth_threshold"]
    status = arrays["pd__site_status"].astype(int)
    action = arrays["pd__birth_cumulative_hazard"]
    available = np.flatnonzero(status == 0)
    rows = []
    for site_id in available:
        pd_node_id = int(site_node[site_id])
        cumulative = float(action[pd_node_id])
        rows.append({
            "site_id": int(site_id), "pd_node_id": pd_node_id,
            "site_threshold": float(threshold[site_id]),
            "node_cumulative_action": cumulative,
            "remaining_threshold_action": float(threshold[site_id] - cumulative),
            "site_status": int(status[site_id]),
        })
    rows.sort(key=lambda row: (row["remaining_threshold_action"], row["site_id"]))
    table = out / "available_persistent_site_clock_table.csv"
    with table.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=tuple(rows[0]))
        writer.writeheader(); writer.writerows(rows)

    files = sorted(path for path in out.rglob("*") if path.is_file())
    manifest = {
        "schema": "V9_PD_ATOMIC_TAIL_PACKAGE_1",
        "source_case_dir": str(case),
        "generation": args.generation,
        "cycles": float(summary["cycles"]),
        "sigma_a_MPa": args.sigma_a_MPa,
        "available_site_rows": len(rows),
        "missing_optional_context": missing_optional_context,
        "files": {str(path.relative_to(out)): {"bytes": path.stat().st_size, "sha256": sha256(path)}
                  for path in files},
    }
    (out / "PACKAGE_MANIFEST.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
