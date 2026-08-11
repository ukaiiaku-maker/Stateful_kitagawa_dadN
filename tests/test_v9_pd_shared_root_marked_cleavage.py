import math

import numpy as np
from pathlib import Path

from arrhenius_fracture.v9_canonical_four_class_birth import CanonicalFourClassBirthState
from arrhenius_fracture.v9_four_class_registry import EXPECTED
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


def test_strongly_localized_phase_crossing_is_analytically_exact():
    class PhaseMPZ(FakeMPZ):
        def copy(self):
            x = PhaseMPZ(); x.value = self.value; return x
        def cleavage_log_rate_s(self, sigma, T):
            return math.log(sigma) if sigma > 0 else -1000.0
    clock = SharedRootMarkedCleavageState(
        "test", ".", shear_modulus_Pa=1., poisson=.3, burgers_m=1e-10,
        initial_tip_radius_m=1e-4, hazard_seed=4, mark_seed=9,
        mpz=PhaseMPZ(), m_hits=1.0,
    )
    # Four quarter-cycle phases; only phase 2 has rate 8/s. Its full action is
    # 2, so a unit threshold crosses halfway through that phase: N=0.625.
    clock.global_threshold_action = 1.0
    tensors = np.zeros((4, 2, 2)); tensors[2, 0, 0] = 8.0
    result = clock.propose_phase_block(1.0, 1.0, 300.0, tensors)
    assert result["crossed"]
    assert result["phase_index"] == 2
    np.testing.assert_allclose(result["phase_fraction"], 0.5, atol=1e-14)
    np.testing.assert_allclose(result["cycles_consumed"], 0.625, atol=1e-14)
    np.testing.assert_allclose(result["state"].global_cumulative_action, 1.0)


def test_real_canonical_shared_root_action_and_mpz_parity():
    source = Path("/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1")
    option = next(iter(EXPECTED))
    kwargs = dict(shear_modulus_Pa=160.15625e9, poisson=.28,
                  burgers_m=2.74e-10, initial_tip_radius_m=1e-6)
    canonical = CanonicalFourClassBirthState(
        option, source, hazard_seed=1720, **kwargs
    )
    shared = SharedRootMarkedCleavageState(
        option, source, hazard_seed=1720, mark_seed=1721, **kwargs
    )
    canonical.hazard_threshold_action = shared.global_threshold_action = 1e100
    tensors = np.array([
        [[2.0e9, .3e9], [.3e9, 3.0e9]],
        [[1.0e9, -.2e9], [-.2e9, 1.5e9]],
        [[.5e9, .1e9], [.1e9, .8e9]],
    ])
    midpoint = canonical.copy()
    drives = [midpoint.mpz.resolve_root_tensor(t) for t in tensors]
    midpoint.mpz.advance(
        .5 * .75 / 1000., 300.,
        float(np.mean([d["opening_stress_Pa"] for d in drives])),
        np.mean(np.stack([d["tau_signed_Pa"] for d in drives]), axis=0),
    )
    expected_phase_logs = np.array([
        midpoint.cleavage_rates(d["opening_stress_Pa"], 300.)[
            "cleavage_log_rate_effective_s"
        ] for d in drives
    ])
    c = canonical._advance_phase_block_exact(.75, 1000., 300., tensors)
    s = shared.propose_phase_block(.75, 1000., 300., tensors)
    assert not s["crossed"]
    np.testing.assert_allclose(
        s["detail"]["phase_log_rate_s"],
        expected_phase_logs,
        rtol=1e-12, atol=1e-12,
    )
    np.testing.assert_allclose(s["detail"]["action_increment"], c["hazard_increment"], rtol=1e-14)
    np.testing.assert_allclose(s["state"].global_cumulative_action,
                               canonical.cumulative_cleavage_hazard, rtol=1e-14)
    for key in canonical.mpz._CAPSULE_ARRAYS:
        np.testing.assert_array_equal(getattr(s["state"].mpz.state, key),
                                      getattr(canonical.mpz.state, key))
