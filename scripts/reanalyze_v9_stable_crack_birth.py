#!/usr/bin/env python3
"""Reclassify completed v9 histories at the stable-crack-birth endpoint."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

try:
    from scripts.update_v9_sn_campaign import atomic_text, canonical_hash
except ModuleNotFoundError:  # Direct ``python scripts/...`` execution.
    from update_v9_sn_campaign import atomic_text, canonical_hash


ROOT = Path("runs/sn_v9_single_seed_skeleton")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def finite(value):
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def final_site_ledger(case_dir: Path) -> dict[str, object]:
    pointer = json.loads((case_dir / "v9_generations" / "ACTIVE.json").read_text())
    generation = case_dir / "v9_generations" / pointer["generation"]
    arrays_path = generation / "state_arrays.npz"
    with np.load(arrays_path) as arrays:
        birth = np.asarray(arrays["pd__site_birth_cycle"], float)
        stable = np.asarray(arrays["pd__site_stable_cycle"], float)
        nodes = np.asarray(arrays["pd__site_node_index"], int)
        status = np.asarray(arrays["pd__site_status"], np.uint8)
        node_birth_hazard = np.asarray(arrays["pd__birth_cumulative_hazard"], float)
        stable_ids = np.where(np.isfinite(stable))[0]
        embryo_ids = np.where(np.isfinite(birth))[0]
        stable_id = int(stable_ids[np.argmin(stable[stable_ids])]) if stable_ids.size else None
        embryo_id = int(embryo_ids[np.argmin(birth[embryo_ids])]) if embryo_ids.size else None
        return {
            "generation": pointer["generation"],
            "state_arrays_sha256": sha256(arrays_path),
            "first_embryo_site_id": embryo_id,
            "first_embryo_node": None if embryo_id is None else int(nodes[embryo_id]),
            "first_embryo_cycle_from_ledger": None if embryo_id is None else float(birth[embryo_id]),
            "stable_site_id": stable_id,
            "stable_site_node": None if stable_id is None else int(nodes[stable_id]),
            "stable_cycle_from_ledger": None if stable_id is None else float(stable[stable_id]),
            "stable_site_status_final": None if stable_id is None else int(status[stable_id]),
            "aggregate_site_birth_hazard_final": float(np.sum(node_birth_hazard[nodes])),
        }


def inspect(case_dir: Path) -> dict[str, object]:
    summary_path = case_dir / "summary.json"
    summary = json.loads(summary_path.read_text())
    ledger = final_site_ledger(case_dir)
    first_embryo = finite(summary.get("cycles_first_embryo"))
    stable = finite(summary.get("cycles_first_stable"))
    # A healed site is deliberately returned to the available pool and its
    # current birth-cycle ledger entry is cleared.  The immutable history
    # summary, not the final active-site ledger, therefore owns the exact
    # first-ever embryo time.
    if stable != ledger["stable_cycle_from_ledger"]:
        raise RuntimeError(f"first stable summary/ledger mismatch: {case_dir}")
    old_status = str(summary.get("status"))
    if stable is not None:
        classification = "stable_crack_birth"
    elif old_status == "right_censored" and finite(summary.get("cycles_total")) is not None:
        classification = "right_censored"
    else:
        classification = "invalid"
    signature = summary.get("run_signature", {})
    return {
        "condition_id": case_dir.parents[2].name,
        "class_label": "baseline_K360",
        "option_id": "legacy_baseline_not_four_class",
        "parameter_provenance_sha256": None,
        "T_K": finite(summary.get("T_K")),
        "sigma_a_MPa": finite(summary.get("sigma_a_MPa")),
        "R": finite(summary.get("R")),
        "seed": signature.get("seed"),
        "pd_seed": signature.get("pd_seed"),
        "N_first_embryo_birth": first_embryo,
        "N_stable_crack_birth": stable,
        "N_front_capture_diagnostic": finite(summary.get("cycles_front_capture")),
        "N_physical_handoff_diagnostic": finite(summary.get("cycles_connected")),
        "observation_horizon": finite(summary.get("cycles_total")),
        "classification": classification,
        "legacy_status": old_status,
        "request_sha256": canonical_hash(signature),
        "summary_sha256": sha256(summary_path),
        "source_sha256": summary.get("source_sha256"),
        # This realized path identifies exp(-H_birth) only while no embryo has
        # occurred.  After birth, feedback changes the hazards and a separate
        # no-birth counterfactual or ensemble law would be required.
        "no_embryo_survival_exp_minus_H": (
            math.exp(-float(ledger["aggregate_site_birth_hazard_final"]))
            if first_embryo is None else None
        ),
        "stable_crack_birth_survival_probability": None,
        "stable_survival_semantics": "stable_birth_survival_not_identified_from_single_seed_or_embryo_hazard_alone",
        **ledger,
    }


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    from io import StringIO
    stream = StringIO()
    writer = csv.DictWriter(stream, fields, extrasaction="ignore", lineterminator="\n")
    writer.writeheader(); writer.writerows(rows)
    atomic_text(path, stream.getvalue())


def update(root: Path = ROOT) -> list[dict[str, object]]:
    cases = sorted(root.glob("*/solver_output/*/sigmaA_*MPa"))
    rows = [inspect(case) for case in cases if (case / "summary.json").is_file()]
    primary = [
        "condition_id", "class_label", "option_id", "T_K", "sigma_a_MPa", "R", "seed", "pd_seed",
        "N_stable_crack_birth", "observation_horizon", "classification", "stable_site_id",
        "stable_site_node", "request_sha256", "summary_sha256", "state_arrays_sha256",
    ]
    stages = [
        "condition_id", "T_K", "sigma_a_MPa", "N_first_embryo_birth", "N_stable_crack_birth",
        "N_front_capture_diagnostic", "N_physical_handoff_diagnostic", "first_embryo_site_id",
        "stable_site_id", "classification",
    ]
    survival = [
        "condition_id", "observation_horizon", "classification", "no_embryo_survival_exp_minus_H",
        "stable_crack_birth_survival_probability", "stable_survival_semantics",
    ]
    legacy = [
        "condition_id", "sigma_a_MPa", "N_front_capture_diagnostic",
        "N_physical_handoff_diagnostic", "legacy_status",
    ]
    write_csv(root / "sn_stable_crack_birth_results.csv", rows, primary)
    write_csv(root / "stable_crack_birth_stage_results.csv", rows, stages)
    write_csv(root / "stable_crack_birth_survival_results.csv", rows, survival)
    write_csv(root / "legacy_physical_handoff_diagnostics.csv", rows, legacy)
    endurance = {row["condition_id"]: {
        "observation": row["classification"],
        "embryo_birth_endurance_diagnostic": "empirical_integrable_tail" if row["classification"] == "right_censored" else "event_realized",
        "stable_crack_birth_endurance_classification": "undetermined",
        "no_embryo_survival_exp_minus_H": row["no_embryo_survival_exp_minus_H"],
        "stable_crack_birth_survival_probability": None,
        "semantics": row["stable_survival_semantics"],
    } for row in rows}
    atomic_text(root / "stable_crack_birth_endurance_diagnostics.json", json.dumps(endurance, indent=2, sort_keys=True) + "\n")
    manifest_path = root / "campaign_manifest.json"
    prior = json.loads(manifest_path.read_text()) if manifest_path.is_file() else {}
    manifest = {
        "schema": "V9_STABLE_CRACK_BIRTH_CAMPAIGN_1",
        "primary_endpoint": "stable_crack_birth",
        "supersedes_primary_endpoint": prior.get("schema"),
        "existing_300K_reanalysis": rows,
        "legacy_physical_handoff_is_diagnostic_only": True,
        "four_class_registry_status": "not_yet_loaded",
    }
    atomic_text(manifest_path, json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return rows


if __name__ == "__main__":
    for row in update():
        print(row["condition_id"], row["classification"], row["N_stable_crack_birth"])
