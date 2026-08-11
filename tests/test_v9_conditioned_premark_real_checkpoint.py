import json
from pathlib import Path

import numpy as np

from arrhenius_fracture.v9_conditioned_premark import load_verified_conditioned_capsule


ROOT = Path(__file__).resolve().parents[1]
CAPSULE_DIR = ROOT / "runs/sn_v9_shared_root_m1_peak/conditioned_premark_1500_N50_v1"
BRANCH_DIR = ROOT / "runs/sn_v9_shared_root_m1_peak/conditioned_pilot_loader_smoke5/Peak/shielded/sigmaA_1500MPa"


def test_real_3804_site_conditioned_branch_has_a_persisted_zero_age_generation():
    capsule, manifest, source, replay = load_verified_conditioned_capsule(CAPSULE_DIR)
    loc = capsule["localization"]
    assert len(loc["available_site_ids"]) == 3804
    assert np.isclose(sum(loc["normalized_mark_probability"]), 1.0, rtol=0, atol=2e-14)
    branch = json.loads((BRANCH_DIR / "conditioned_branch_manifest.json").read_text())
    assert branch["source_conditioned_capsule_sha256"] == manifest["capsule_sha256"]
    assert branch["conditioned_attempt_committed"] is True
    assert branch["analysis_only"] is False
    event = branch["branch_summary"]["event"]
    assert event["cycle"] == capsule["crossing_cycle"]
    assert event["phase_index"] == loc["phase_index"]
    generations = list((BRANCH_DIR / "v9_generations").glob("generation_*"))
    zero_age = []
    for generation in generations:
        metadata = json.loads((generation / "state_metadata.json").read_text())
        if metadata["cycles"] != capsule["crossing_cycle"]:
            continue
        with np.load(generation / "state_arrays.npz") as arrays:
            embryo = np.flatnonzero(arrays["pd__site_status"] == 1)
            assert embryo.tolist() == [event["site_id"]]
            site = embryo[0]
            assert arrays["pd__site_birth_cycle"][site] == capsule["crossing_cycle"]
            assert arrays["pd_v9__site_transition_cumulative_hazard"][site] == 0.0
            assert arrays["pd_v9__site_transition_threshold"][site] == branch["branch_summary"]["selected_site_transition_threshold"]
            assert arrays["pd_v9__site_transition_outcome_uniform"][site] == branch["branch_summary"]["selected_site_transition_outcome_uniform"]
            assert not np.any(arrays["pd__stable_sites"])
            assert not np.any(arrays["pd__bond_damage"])
        zero_age.append(generation)
    assert len(zero_age) == 1
    assert source.is_file() and replay.is_file()
