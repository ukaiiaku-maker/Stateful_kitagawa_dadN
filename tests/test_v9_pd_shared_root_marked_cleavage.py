import math

import numpy as np

from arrhenius_fracture.v9_pd_shared_root_marked_cleavage import (
    SharedRootMarkedCleavageState, normalized_available_site_marks,
)


class FakeMPZ:
    audit = {"source": "fake_for_clock_unit_test"}
    def __init__(self): self.value = 0.0
    def capsule(self): return {"schema": "fake", "value": self.value}
    def restore_capsule(self, c): self.value = float(c["value"])
    def summary(self):
        return {"signed_active_K_shield_Pa_sqrt_m": 0.0, "tip_radius_m": 1e-4}
    def cleavage_log_rate_s(self, sigma, T): return math.log(2.0)
    def copy(self):
        x = FakeMPZ(); x.value = self.value; return x
    def resolve_root_tensor(self, tensor):
        return {"opening_stress_Pa": float(tensor[0, 0]),
                "tau_signed_Pa": np.array([tensor[0, 1], -tensor[0, 1]])}
    def advance(self, dt, T, opening, signed): self.value += float(dt)


def make_clock():
    return SharedRootMarkedCleavageState(
        "test", ".", shear_modulus_Pa=1.0, poisson=.3, burgers_m=1e-10,
        initial_tip_radius_m=1e-4, hazard_seed=4, mark_seed=9, mpz=FakeMPZ(),
        m_hits=1.0,
    )


def test_identical_site_multiplicity_does_not_change_global_action():
    one = make_clock(); many = make_clock()
    for clock in (one, many): clock.add_log_action(math.log(.25))
    assert one.global_cumulative_action == many.global_cumulative_action
    _, p1, _ = normalized_available_site_marks([0], [True], [0.0], [1.0])
    _, pn, _ = normalized_available_site_marks([0]*8, [True]*8, [0.0], [1.0])
    np.testing.assert_allclose(p1, [1.0])
    np.testing.assert_allclose(pn, np.full(8, 1/8))


def test_nonuniform_marks_preserve_normalization_and_split_colocated_nodes():
    ids, p, _ = normalized_available_site_marks(
        [0, 0, 1], [True, True, True], [0.0, math.log(3.0)], [1.0, 1.0]
    )
    np.testing.assert_array_equal(ids, [0, 1, 2])
    np.testing.assert_allclose(p, [.125, .125, .75])
    assert np.sum(p) == 1.0


def test_capsule_preserves_both_rng_streams_and_selected_identity():
    a = make_clock(); a.global_threshold_action = .1; a.add_log_action(math.log(.1))
    capsule = a.capsule()
    first = a.select_mark([0, 1], [True, True], [0.0, 0.0], [1.0, 1.0], cycle=2, phase_index=3)
    b = make_clock(); b.restore_capsule(capsule)
    second = b.select_mark([0, 1], [True, True], [0.0, 0.0], [1.0, 1.0], cycle=2, phase_index=3)
    assert first == second
    assert a.capsule()["hazard_rng_state"] == b.capsule()["hazard_rng_state"]
    assert a.capsule()["mark_rng_state"] == b.capsule()["mark_rng_state"]


def test_no_mark_support_fails_closed():
    try:
        normalized_available_site_marks([0], [False], [0.0], [1.0])
    except RuntimeError as exc:
        assert "no available" in str(exc)
    else:
        raise AssertionError("missing mark support did not fail closed")


def test_phase_block_proposal_is_transactional_and_partition_equivalent():
    tensors = np.array([[[1., 0.], [0., 0.]], [[2., 0.], [0., 0.]]])
    base = make_clock(); before = base.capsule()
    whole = base.propose_phase_block(.01, 1.0, 300., tensors)["state"]
    assert base.capsule() == before
    half = base.propose_phase_block(.005, 1.0, 300., tensors)["state"]
    second = half.propose_phase_block(.005, 1.0, 300., tensors)["state"]
    np.testing.assert_allclose(whole.global_cumulative_action,
                               second.global_cumulative_action, rtol=1e-14)
    np.testing.assert_allclose(whole.mpz.value, second.mpz.value, rtol=0, atol=0)
