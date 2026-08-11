import math

from scripts.analyze_v9_peak_attempt_survival import invert_piecewise


def test_action_inversion_uses_bracketing_curve_not_origin_h_over_n():
    rows=[{"cycles":10.,"H_attempt":1.},{"cycles":20.,"H_attempt":5.}]
    n,left,right=invert_piecewise(rows,3.)
    assert n == 15.
    assert (left["cycles"],right["cycles"]) == (10.,20.)
    assert not math.isclose(n,3./(5./20.))
