"""Fail-closed bridge to the canonical audited v9.13 four-class rows."""
from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path


SCHEMA = "V9_CANONICAL_FOUR_CLASS_REGISTRY_BRIDGE_1"
EXPECTED = {
    "v913_paper_peak01_0242980_persistent_sites": ("peak", "v913_zeroD_sobol_0242980", "v10_2_25_v913_paper_campaign_registry.csv"),
    "v913_paper_dbtt01_0202500_persistent_sites": ("DBTT", "v913_zeroD_sobol_0202500", "v10_2_25_v913_paper_campaign_registry.csv"),
    "v913_paper_weakT01_0257068_persistent_sites": ("weakT", "v913_zeroD_sobol_0257068", "v10_2_26_v913_weakT_ceramic_registry.csv"),
    "v913_paper_ceramic01_0189364_persistent_sites": ("ceramic", "v913_zeroD_sobol_0189364", "v10_2_26_v913_weakT_ceramic_registry.csv"),
}
FORBIDDEN_ALTERNATES = {
    "v913_paper_weakT01_0129902_persistent_sites",
    "v913_paper_ceramic01_0077080_persistent_sites",
}
EXPECTED_HASHES = {
    "parameter_registry_v9111.py": "6946759c5786344a7111b9e25ce8b2513e5335d78fc590ba860fc730bb7d83a7",
    "v10_2_25_v913_paper_campaign_registry.csv": "d2a427a5d78047fef652ca87d769da8b41124bfc3f6f982472ceddb710371ffd",
    "v10_2_25_v913_paper_campaign_selection.json": "47e21de1c34ff220c5117caf937c907bc7b75374dc34dea0236ce97b7494af94",
    "v10_2_26_v913_weakT_ceramic_registry.csv": "04fd07113b2611f87727084bae9cffa1c740df35f0ca1ffec25fa89c0c50196a",
    "v10_2_26_v913_weakT_ceramic_selection.json": "b8a539a7b15f8aec780747450dae803000ac7a2a7e7a3c808d57a349d8d6b9b8",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _paths(source_root: str | Path) -> dict[str, Path]:
    root = Path(source_root).expanduser().resolve()
    materials = root / "arrhenius_fracture" / "data" / "materials"
    return {
        "parameter_registry_v9111.py": root / "arrhenius_fracture" / "parameter_registry_v9111.py",
        "v10_2_25_v913_paper_campaign_registry.csv": materials / "v10_2_25_v913_paper_campaign_registry.csv",
        "v10_2_25_v913_paper_campaign_selection.json": materials / "v10_2_25_v913_paper_campaign_selection.json",
        "v10_2_26_v913_weakT_ceramic_registry.csv": materials / "v10_2_26_v913_weakT_ceramic_registry.csv",
        "v10_2_26_v913_weakT_ceramic_selection.json": materials / "v10_2_26_v913_weakT_ceramic_selection.json",
    }


def select_canonical_option(option_id: str, source_root: str | Path):
    if option_id in FORBIDDEN_ALTERNATES or option_id not in EXPECTED:
        raise ValueError(f"not a canonical stable-crack-birth four-class option: {option_id!r}")
    paths = _paths(source_root)
    actual = {name: _sha256(path) for name, path in paths.items()}
    if actual != EXPECTED_HASHES:
        raise RuntimeError(f"audited four-class source hash mismatch: {actual!r}")
    selector_path = paths["parameter_registry_v9111.py"]
    spec = importlib.util.spec_from_file_location("_audited_parameter_registry_v9111", selector_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load audited selector: {selector_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    class_label, candidate_id, registry_name = EXPECTED[option_id]
    selected = module.select_option(option_id, paths[registry_name], canonical_stage3_only=False)
    if selected.candidate_id != candidate_id:
        raise RuntimeError(f"canonical candidate mismatch for {option_id}")
    selection_name = registry_name.replace("registry.csv", "selection.json")
    selection = json.loads(paths[selection_name].read_text())
    candidates = selection.get("primary_candidates", [])
    if sum(row.get("option_key") == option_id for row in candidates) != 1:
        raise RuntimeError(f"selection manifest does not uniquely authorize {option_id}")
    return selected, {
        "schema": SCHEMA,
        "class_label": class_label,
        "option_id": option_id,
        "candidate_id": candidate_id,
        "source_repository": str(Path(source_root).expanduser().resolve()),
        "selector_path": str(selector_path),
        "registry_path": str(paths[registry_name]),
        "selection_path": str(paths[selection_name]),
        "source_hashes": actual,
        "exact_registry_row": dict(selected.row),
    }
