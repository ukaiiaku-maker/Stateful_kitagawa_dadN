#!/usr/bin/env python3
"""Summarize independent Stateful-PD v8.3 case/stress jobs without stale states."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

VALID_FAILURE = {"physical_handoff"}
VALID_CENSOR = {"right_censored"}
PROVISIONAL = {"physical_handoff_geometry_saturated", "right_censored_geometry_saturated"}
INVALID_PREFIXES = ("morphology_invalid", "geometry_invalid", "invalid_")


def _read_json(path: Path):
    try:
        return json.loads(path.read_text())
    except Exception:
        return None


def _checkpoint_status(path: Path):
    try:
        with np.load(path, allow_pickle=False) as data:
            meta = json.loads(str(data["metadata_json"].item()))
        rows = meta.get("rows", [])
        last = rows[-1] if rows else {}
        return {
            "cycles_total": float(meta.get("cycles", 0.0)),
            "block": int(meta.get("next_block", 0)),
            "root_radius_over_spacing_final": last.get("root_radius_over_spacing"),
            "pd_crack_centerline_length_final_m": last.get("pd_crack_centerline_length_m"),
            "pd_crack_orientation_coherence_final": last.get("pd_crack_orientation_coherence"),
            "pd_crack_width_ratio_final": last.get("pd_crack_width_ratio"),
            "pd_handoff_offfront_broken_fraction_final": last.get("pd_off_front_broken_fraction"),
            "pd_primary_seed_reselections_final": last.get("pd_primary_seed_reselections"),
            "geometry_saturated": bool(last.get("geometry_saturated", False)),
        }
    except Exception:
        return {}


def _job_output(root: Path, job: dict) -> Path:
    stress = float(job["stress"])
    case = str(job["case"])
    seed = int(job["seed"])
    stem = root / f"seed_{seed}" / f"stress_{stress:g}MPa" / f"job_{case}"
    return stem / case / (f"sigmaA_{stress:g}MPa".replace(".", "p"))


def _category(status: str) -> str:
    if status in VALID_FAILURE:
        return "valid_failure"
    if status in VALID_CENSOR:
        return "valid_censor"
    if status in PROVISIONAL:
        return "provisional_geometry_limited"
    if status.startswith(INVALID_PREFIXES) or status == "max_blocks_reached":
        return "invalid"
    if status in {"running", "queued", "incomplete_checkpoint", "not_started"}:
        return "incomplete"
    return "unknown"


def _fmt(value, fmt):
    try:
        x = float(value)
        return format(x, fmt) if math.isfinite(x) else "-"
    except Exception:
        return "-"


def build_parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", required=True)
    p.add_argument("--model", required=True)
    p.add_argument("--cycles-max", required=True, type=float)
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    root = Path(args.root).resolve()
    campaign = _read_json(root / "campaign_definition.json") or {}
    process_status = _read_json(root / "job_process_status.json") or {}
    jobs = campaign.get("jobs", [])
    rows = []
    for job in jobs:
        key = job.get("key") or f"seed_{job['seed']}__{job['case']}__{float(job['stress']):g}MPa"
        out = _job_output(root, job)
        summary = _read_json(out / "summary.json")
        proc = process_status.get(key, {})
        if summary:
            row = dict(summary)
            status = str(summary.get("status", "unknown"))
            if str(summary.get("model", "")) != str(args.model):
                status = "invalid_model_mismatch"
            elif status == "right_censored" and float(summary.get("cycles_total", 0.0) or 0.0) < 0.999999 * float(args.cycles_max):
                status = "invalid_inconsistent_censor"
            elif status == "physical_handoff" and summary.get("cycles_connected") is None:
                status = "invalid_inconsistent_handoff"
        elif (out / "checkpoint_latest.npz").exists():
            row = _checkpoint_status(out / "checkpoint_latest.npz")
            pstat = proc.get("process_status")
            status = "running" if pstat == "running" else "incomplete_checkpoint"
        else:
            row = {}
            pstat = proc.get("process_status")
            status = "running" if pstat == "running" else "queued" if pstat == "queued" else "not_started"
        row.update({
            "job_key": key,
            "seed": int(job["seed"]),
            "case": str(job["case"]),
            "sigma_a_MPa": float(job["stress"]),
            "status": status,
            "category": _category(status),
            "process_status": proc.get("process_status", "unknown"),
            "return_code": proc.get("return_code"),
            "output_dir": str(out),
        })
        rows.append(row)

    rows.sort(key=lambda r: (r["seed"], -r["sigma_a_MPa"], r["case"]))
    keys = sorted(set().union(*(r.keys() for r in rows))) if rows else []
    with (root / "stress_sweep_cases.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader(); w.writerows(rows)

    lines = [
        "STATEFUL-PD v8.3 independent stress sweep",
        "",
        "Only unsaturated physical_handoff is a valid failure. Only unsaturated right_censored is valid censoring.",
        "Geometry-limited, morphology-invalid, and max-block outcomes are excluded from production S-N fits.",
        "",
        f"{'seed':>4} {'case':>10} {'stress':>7} {'status':>38} {'N_end':>11} {'N_handoff':>11} {'a_um':>8} {'C':>6} {'w/a':>6} {'R/h':>6} {'re':>4}",
        "-" * 120,
    ]
    for r in rows:
        lines.append(
            f"{r['seed']:4d} {r['case']:>10} {r['sigma_a_MPa']:7.0f} {r['status']:>38} "
            f"{_fmt(r.get('cycles_total'),' .3e'):>11} {_fmt(r.get('cycles_connected'),' .3e'):>11} "
            f"{_fmt(1e6*float(r.get('pd_crack_centerline_length_final_m',0) or 0),' .1f'):>8} "
            f"{_fmt(r.get('pd_crack_orientation_coherence_final'),' .2f'):>6} "
            f"{_fmt(r.get('pd_crack_width_ratio_final'),' .2f'):>6} "
            f"{_fmt(r.get('root_radius_over_spacing_final'),' .1f'):>6} "
            f"{int(r.get('pd_primary_seed_reselections_final',0) or 0):4d}"
        )

    counts = {}
    for r in rows:
        counts[r["category"]] = counts.get(r["category"], 0) + 1
    lines += ["", "category counts: " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items()))]

    # Matched-pair audit.
    pair_rows = []
    grouped = {}
    for r in rows:
        grouped.setdefault((r["seed"], r["sigma_a_MPa"]), {})[r["case"]] = r
    for (seed, stress), pair in sorted(grouped.items(), key=lambda x: (x[0][0], -x[0][1])):
        a, b = pair.get("no_shield"), pair.get("shielded")
        if not a or not b:
            continue
        ratio = None
        if a["status"] == "physical_handoff" and b["status"] == "physical_handoff":
            na = float(a.get("cycles_connected") or a.get("cycles_total"))
            nb = float(b.get("cycles_connected") or b.get("cycles_total"))
            ratio = nb / max(na, 1e-300)
        pair_rows.append({
            "seed": seed,
            "sigma_a_MPa": stress,
            "no_shield_status": a["status"],
            "shielded_status": b["status"],
            "shielded_over_no_shield_life_ratio": ratio,
        })
    with (root / "matched_pair_audit.csv").open("w", newline="") as f:
        fields = ["seed", "sigma_a_MPa", "no_shield_status", "shielded_status", "shielded_over_no_shield_life_ratio"]
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(pair_rows)

    # Same-seed monotonicity audit for valid handoffs only.
    inversions = []
    for seed in sorted({r["seed"] for r in rows}):
        for case in sorted({r["case"] for r in rows}):
            valid = [r for r in rows if r["seed"] == seed and r["case"] == case and r["status"] == "physical_handoff"]
            valid.sort(key=lambda r: r["sigma_a_MPa"])
            for low, high in zip(valid, valid[1:]):
                nlow = float(low.get("cycles_connected") or low.get("cycles_total"))
                nhigh = float(high.get("cycles_connected") or high.get("cycles_total"))
                if nhigh > nlow:
                    inversions.append({
                        "seed": seed, "case": case,
                        "lower_stress_MPa": low["sigma_a_MPa"], "lower_life": nlow,
                        "higher_stress_MPa": high["sigma_a_MPa"], "higher_life": nhigh,
                    })
    inv_fields = ["seed", "case", "lower_stress_MPa", "lower_life", "higher_stress_MPa", "higher_life"]
    with (root / "monotonicity_audit.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=inv_fields); w.writeheader(); w.writerows(inversions)
    lines.append(f"valid-handoff monotonicity inversions: {len(inversions)}")

    table = "\n".join(lines) + "\n"
    (root / "stress_sweep_status_table.txt").write_text(table)
    print(table, end="")
    payload = {
        "model": args.model,
        "cycles_max": args.cycles_max,
        "counts": counts,
        "matched_pairs": pair_rows,
        "monotonicity_inversions": inversions,
        "all_jobs_final": all(r["category"] not in {"incomplete", "unknown"} for r in rows),
    }
    (root / "stress_sweep_summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
