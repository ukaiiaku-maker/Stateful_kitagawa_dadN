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
