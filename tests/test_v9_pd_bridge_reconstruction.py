import math

import numpy as np
import pandas as pd


def test_unit_consistent_rate_separation_is_per_cycle():
    raw_s = 2.377e10
    birth_per_cycle = 2.897e-30
    separation = math.log10(raw_s / 1000.0) - math.log10(birth_per_cycle)
    assert 36.7 < separation < 37.0


def test_legacy_patch_survival_sums_realized_available_sites():
    table = pd.DataFrame({
        "cycles": [1.0, 1.0, 1.0, 2.0, 2.0, 2.0],
        "available": [True, True, False, True, True, False],
        "cumulative_action": [0.1, 0.2, 9.0, 0.3, 0.4, 9.0],
    })
    available = table[table.available]
    action = available.groupby("cycles").cumulative_action.sum()
    survival = np.exp(-action.to_numpy())
    np.testing.assert_allclose(action.to_numpy(), [0.3, 0.7])
    np.testing.assert_allclose(survival, np.exp(-np.array([0.3, 0.7])))
    assert not math.isclose(survival[0], math.exp(-0.2))
