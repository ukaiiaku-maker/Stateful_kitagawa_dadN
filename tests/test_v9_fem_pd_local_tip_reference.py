import hashlib
from pathlib import Path

from scripts.compare_v9_fem_pd_local_tip import FEM_ROOT, parity


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
