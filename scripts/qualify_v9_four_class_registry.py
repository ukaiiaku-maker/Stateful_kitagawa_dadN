#!/usr/bin/env python3
"""Resolve and record the exact canonical four-class provenance."""
from __future__ import annotations

import json
from pathlib import Path

from arrhenius_fracture.v9_four_class_registry import EXPECTED, SCHEMA, select_canonical_option
from scripts.update_v9_sn_campaign import atomic_text


SOURCE = Path("/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1")
ROOT = Path("runs/sn_v9_single_seed_skeleton")
TEMPERATURES_K = [300, 600, 800, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300]


def main():
    records = []
    for option_id in EXPECTED:
        _, audit = select_canonical_option(option_id, SOURCE)
        records.append(audit)
    payload = {
        "schema": SCHEMA,
        "qualification": "exact_registry_identity_passed",
        "canonical_temperature_grid_K": TEMPERATURES_K,
        "temperature_grid_source": str(SOURCE / "scripts/run_v10_2_27_paper_four_class_30deg_long_rcurves.sh"),
        "records": records,
        "physical_adapter_status": "blocked_pending_explicit_signed_MPZ_to_v9_prebirth_mapping_decision",
    }
    atomic_text(ROOT / "four_class_registry_qualification.json", json.dumps(payload, indent=2, sort_keys=True) + "\n")
    manifest_path = ROOT / "campaign_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["four_class_registry_status"] = payload["qualification"]
    manifest["four_class_registry_qualification"] = payload
    atomic_text(manifest_path, json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"options": len(records), "temperatures_K": TEMPERATURES_K}))


if __name__ == "__main__":
    main()
