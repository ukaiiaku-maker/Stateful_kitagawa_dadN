import hashlib
from pathlib import Path

from scripts.compare_v9_fem_pd_local_tip import FEM_ROOT, PD_CASES, fem_rows, parity, pd_record
from scripts.qualify_v9_fem_pd_tensor_matrix import evaluate, fem_state, pd_state


def _hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def test_side_effect_free_cross_evaluation_has_constitutive_parity():
    watched = [FEM_ROOT / "campaign_manifest.json", FEM_ROOT / "condition_registry.json"]
    before = {path: _hash(path) for path in watched}
    result = parity()
    assert result["pass"] is True
    assert result["m3_action"]["absolute_difference"] == 0.0
    assert result["m1_raw_action"]["relative_difference"] < 5e-13
    assert before == {path: _hash(path) for path in watched}


def test_stress_concentration_definitions_are_distinct_and_root_opening_drives_transfer():
    k360 = fem_rows()[0]
    blunt = pd_record(2000.0, PD_CASES[2000.0])
    assert blunt["root_node"] == 11
    assert blunt["hotspot_node"] == 124
    assert 14.3e-6 < blunt["root_hotspot_distance_m"] < 14.5e-6
    assert 1.68 < blunt["Kt_root_opening"] < 1.69
    assert 1.70 < blunt["Kt_root_principal"] < 1.72
    assert 2.03 < blunt["Kt_hotspot_principal"] < 2.04
    assert 3.30 < k360["Kt_root_opening"] / blunt["Kt_root_opening"] < 3.31


def test_cross_geometry_saved_capsule_tensor_parity_uses_effective_tip_transform():
    result = evaluate(
        fem_state(735.9214951373579, 1e14),
        pd_state(1500.0, "generation_9a6c40a9d55b4149a6b79688650d7541"),
    )
    assert result["pass_parity"] is True
    assert result["m1_action_relative_error"] == 0.0
    assert result["m3_action_relative_error"] < 5e-13
    assert result["signed_state_m1_absolute_error"] == 0.0
