#!/usr/bin/env python3
"""Build the restartable, one-condition-at-a-time v9 S-N registry."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import tempfile
from pathlib import Path
from typing import Any


CASE_GLOB = "*/solver_output/*/sigmaA_*MPa"


def canonical_hash(value: Any) -> str:
    raw = json.dumps(value, allow_nan=True, separators=(",", ":"), sort_keys=True)
    return hashlib.sha256(raw.encode()).hexdigest()


def atomic_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as stream:
            stream.write(value)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(name, path)
    finally:
        if os.path.exists(name):
            os.unlink(name)


def finite(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def read_last_history(case_dir: Path) -> tuple[dict[str, str], int, float]:
    paths = [case_dir / "sn_stateful_pd_history.csv", case_dir / "sn_stateful_pd_history_partial.csv"]
    path = next((item for item in paths if item.is_file()), None)
    if path is None:
        return {}, 0, 0.0
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        return {}, 0, 0.0
    return rows[-1], len(rows), max(float(row["dN"]) for row in rows)


def active_state(case_dir: Path) -> dict[str, Any] | None:
    pointer = case_dir / "v9_generations" / "ACTIVE.json"
    if not pointer.is_file():
        return None
    try:
        generation = json.loads(pointer.read_text())["generation"]
        root = pointer.parent / generation
        summary = json.loads((root / "summary.json").read_text())
        metadata = json.loads((root / "state_metadata.json").read_text())
        return {"generation": generation, "summary": summary, "metadata": metadata}
    except (KeyError, OSError, json.JSONDecodeError):
        return None


def classify(summary: dict[str, Any] | None, active: dict[str, Any] | None) -> str:
    if summary is None:
        return "restartable" if active is not None else "not_started"
    status = summary.get("status")
    geometry_ok = not summary.get("geometry_saturated", False) and not summary.get("geometry_invalid_reason")
    if status == "physical_handoff":
        return "physical_handoff" if geometry_ok and summary.get("pd_handoff_pass_final") else "invalid"
    if status == "right_censored":
        return "right_censored" if geometry_ok else "invalid"
    return "invalid"


def inspect_condition(case_dir: Path) -> dict[str, Any]:
    run_args_path = case_dir / "run_args.json"
    run_args = json.loads(run_args_path.read_text()) if run_args_path.is_file() else {}
    summary_path = case_dir / "summary.json"
    summary = json.loads(summary_path.read_text()) if summary_path.is_file() else None
    active = active_state(case_dir)
    last, accepted, largest = read_last_history(case_dir)
    classification = classify(summary, active)
    cycles = finite((summary or {}).get("cycles_total")) or finite(last.get("cycles_total")) or 0.0
    signature = (summary or {}).get("run_signature", run_args)
    return {
        "condition_id": case_dir.parents[2].name,
        "case": run_args.get("case", case_dir.parent.name),
        "sigma_a_MPa": finite((summary or {}).get("sigma_a_MPa")) or finite(run_args.get("sigma_a_MPa")),
        "seed": run_args.get("seed"),
        "pd_seed": run_args.get("pd_seed"),
        "request_sha256": canonical_hash(signature),
        "classification": classification,
        "cycles": cycles,
        "accepted_blocks": accepted,
        "rejected_substeps": sum(int(float(row.get("v9_fem_rejections", 0) or 0)) for row in _rows(case_dir)),
        "largest_accepted_dN": largest,
        "active_generation": None if active is None else active["generation"],
        "geometry_valid": classification not in {"invalid", "not_started"},
        "handoff_valid": classification == "physical_handoff",
        "cumulative_expected_births": finite(last.get("pd_expected_births_cumulative")),
        "instantaneous_birth_hazard_per_cycle": finite(last.get("pd_birth_rate_max_per_cycle")),
        "rho_max_m2": finite(last.get("rho_max_m2")),
        "summary_path": str(summary_path) if summary_path.is_file() else None,
        "case_dir": str(case_dir),
    }


def _rows(case_dir: Path) -> list[dict[str, str]]:
    for name in ("sn_stateful_pd_history.csv", "sn_stateful_pd_history_partial.csv"):
        path = case_dir / name
        if path.is_file():
            with path.open(newline="") as stream:
                return list(csv.DictReader(stream))
    return []


def select_next_stress(conditions: list[dict[str, Any]]) -> dict[str, Any] | None:
    finite_lives = sorted(
        (float(c["sigma_a_MPa"]), math.log10(float(c["cycles"])))
        for c in conditions
        if c["classification"] == "physical_handoff" and c["cycles"] > 0
    )
    if len(finite_lives) < 2:
        return None
    targets = [4.0, 5.0, 6.0, 7.0, 8.0]
    covered = [life for _, life in finite_lives]
    target = max(targets, key=lambda x: min(abs(x - y) for y in covered))
    by_life = sorted((life, stress) for stress, life in finite_lives)
    lo, hi = min(by_life), max(by_life)
    if not lo[0] <= target <= hi[0]:
        return None
    for (n0, s0), (n1, s1) in zip(by_life, by_life[1:]):
        if n0 <= target <= n1:
            stress = s0 + (s1 - s0) * (target - n0) / (n1 - n0)
            return {"sigma_a_MPa": stress, "target_log10_cycles": target, "method": "largest_unfilled_log_life_gap"}
    return None


def update(root: Path) -> list[dict[str, Any]]:
    conditions = [inspect_condition(path) for path in sorted(root.glob(CASE_GLOB))]
    manifest = {
        "schema": "V9_SINGLE_SEED_SN_CAMPAIGN_1",
        "policy": "one_condition_at_a_time_adaptive_log_life_gap",
        "active_conditions": [c["condition_id"] for c in conditions if c["classification"] == "restartable"],
        "conditions": conditions,
        "next_condition": select_next_stress(conditions),
    }
    atomic_text(root / "campaign_manifest.json", json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    fields = ["condition_id", "case", "sigma_a_MPa", "seed", "classification", "cycles", "geometry_valid", "handoff_valid", "request_sha256"]
    lines: list[str] = []
    from io import StringIO
    stream = StringIO(); writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader()
    for condition in conditions:
        writer.writerow({key: condition.get(key) for key in fields})
    atomic_text(root / "sn_results.csv", stream.getvalue())
    survival = StringIO(); sw = csv.writer(survival); sw.writerow(["condition_id", "cycles", "cumulative_expected_births", "survival_probability"])
    for c in conditions:
        hazard = c["cumulative_expected_births"]
        sw.writerow([c["condition_id"], c["cycles"], hazard, None if hazard is None else math.exp(-hazard)])
    atomic_text(root / "survival_results.csv", survival.getvalue())
    diagnostics = {c["condition_id"]: {
        "observation": c["classification"],
        "physics_asymptotic_classification": "undetermined",
        "cycles": c["cycles"],
        "instantaneous_birth_hazard_per_cycle": c["instantaneous_birth_hazard_per_cycle"],
        "cumulative_expected_births": c["cumulative_expected_births"],
        "note": "Finite-horizon censoring is not endurance evidence.",
    } for c in conditions}
    atomic_text(root / "endurance_diagnostics.json", json.dumps(diagnostics, indent=2, sort_keys=True) + "\n")
    return conditions


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("runs/sn_v9_single_seed_skeleton"))
    args = parser.parse_args()
    conditions = update(args.root)
    for c in conditions:
        print(f"{c['condition_id']}: {c['classification']} N={c['cycles']:.9g} blocks={c['accepted_blocks']}")


if __name__ == "__main__":
    main()
