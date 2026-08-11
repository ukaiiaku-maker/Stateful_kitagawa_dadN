import pytest

from scripts.summarize_v9_shared_root_m1_peak import validate_condition_registry


def row(name, eligible):
    return {"condition":name,"sigma_a_MPa":2000.0,"hazard_seed":42017,"mark_seed":42018,
            "production_analyzer_eligible":eligible}


def test_registry_allows_labeled_coarse_control_but_one_canonical_row():
    selected=validate_condition_registry([row("coarse",False),row("rate_separated",True)])
    assert selected == {(2000.0,42017,42018):"rate_separated"}


def test_registry_fails_closed_on_ambiguous_canonical_duplicate():
    with pytest.raises(ValueError,match="ambiguous canonical stress/seed rows"):
        validate_condition_registry([row("a",True),row("b",True)])


def test_generated_registry_has_only_rate_separated_2000_as_production():
    import json
    from pathlib import Path
    path=Path("runs/sn_v9_shared_root_m1_peak/qualification/condition_registry.json")
    data=json.loads(path.read_text())
    by_name={row["condition"]:row for row in data["conditions"]}
    assert not by_name["m1_2000_seed42_coarse_1e5_baseline"]["production_analyzer_eligible"]
    assert by_name["m1_2000_seed42_rate_separated_v2"]["production_analyzer_eligible"]
    historical=[row for row in data["conditions"] if row["numerical_protocol"]=="historical_mixed_protocol_diagnostic"]
    assert historical and all(not row["production_analyzer_eligible"] for row in historical)
    assert all(row["diagnostic_ensemble_eligible"] for row in historical)
