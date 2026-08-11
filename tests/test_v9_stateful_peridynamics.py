import unittest

import numpy as np

from arrhenius_fracture.stateful_peridynamics_v8_7_local_front_spacing import (
    _SITE_EMBRYO, _SITE_STABLE, StatefulPDConfig,
)
from arrhenius_fracture.v9_stateful_peridynamics import V9StatefulPDPatch


class _TransitionHarness(V9StatefulPDPatch):
    def __init__(self):
        self.xy = np.zeros((1, 2))
        self.cfg = StatefulPDConfig(
            patch_radius_m=1.0, horizon_m=1.0, boundary_shell_m=0.1,
            heal_return_fraction=0.0, random_seed=7,
        )
        self._event_rng = np.random.default_rng(11)

    def _sync_site_ledger(self, state):
        return None

    def _refresh_node_counts_from_site_ledger(self, state):
        return None

    def create_marked_embryo(self, state, site, cycle):
        state.site_status[site] = _SITE_EMBRYO
        state.site_birth_cycle[site] = cycle
        return int(state.site_node_index[site])


class _State:
    pass


def make_state():
    state = _State()
    state.site_node_index = np.array([0])
    state.site_status = np.array([_SITE_EMBRYO], dtype=np.uint8)
    state.site_birth_cycle = np.array([0.0])
    state.site_stable_cycle = np.array([np.nan])
    state.site_transition_threshold = np.array([0.73])
    state.site_transition_cumulative_hazard = np.array([0.0])
    state.site_transition_outcome_uniform = np.array([0.2])
    state.log_delivery_memory = np.array([-np.inf])
    state.log_birth_cumulative_hazard = np.array([-np.inf])
    state.healed_sites_cumulative = np.array([0])
    state.birth_cumulative_hazard = np.array([0.0])
    state.site_birth_threshold = np.array([1.0])
    state.born_sites_cumulative = np.array([0])
    state.available = np.array([1.0])
    state.embryo = np.array([0.0])
    state.born_cumulative = np.array([0.0])
    state.available_sites = np.array([1])
    state.embryo_sites = np.array([0])
    state.cycles_first_embryo = None
    return state


class V9StatefulPDTests(unittest.TestCase):
    def test_every_external_embryo_gets_next_persistent_transition_draw(self):
        patch=_TransitionHarness(); state=make_state(); state.site_status[0]=0
        rng=np.random.default_rng(991); expected=np.random.default_rng(991)
        first=(float(expected.exponential()),float(expected.random()))
        patch.create_external_marked_embryo(state,0,10.,rng=rng)
        self.assertEqual(state.site_transition_threshold[0],first[0]);self.assertEqual(state.site_transition_outcome_uniform[0],first[1])
        self.assertEqual(state.site_transition_cumulative_hazard[0],0.)
        state.site_status[0]=0
        second=(float(expected.exponential()),float(expected.random()))
        saved=rng.bit_generator.state
        patch.create_external_marked_embryo(state,0,20.,rng=rng)
        self.assertEqual(state.site_transition_threshold[0],second[0]);self.assertEqual(state.site_transition_outcome_uniform[0],second[1])
        replay=np.random.default_rng();replay.bit_generator.state=saved
        self.assertEqual(float(replay.exponential()),second[0]);self.assertEqual(float(replay.random()),second[1])
    def test_transition_crossing_is_partition_invariant(self):
        patch = _TransitionHarness()
        whole = make_state()
        patch._advance_discrete_embryo_transitions(whole, np.array([0.15]), np.array([0.05]), 0.0, 5.0)
        split = make_state()
        cycle = 0.0
        for dN in (0.1, 0.3, 0.7, 1.2, 2.7):
            patch._advance_discrete_embryo_transitions(split, np.array([0.15]), np.array([0.05]), cycle, dN)
            cycle += dN
        self.assertEqual(int(whole.site_status[0]), int(split.site_status[0]))
        self.assertAlmostEqual(float(whole.site_stable_cycle[0]), 3.65, places=13)
        self.assertAlmostEqual(float(split.site_stable_cycle[0]), 3.65, places=13)
        self.assertAlmostEqual(float(whole.site_transition_cumulative_hazard[0]), 0.73)
        self.assertAlmostEqual(float(split.site_transition_cumulative_hazard[0]), 0.73)

    def test_stable_endpoint_boundary_matches_full_trajectory_event(self):
        patch = _TransitionHarness()
        full = make_state()
        event_rng_before = patch._event_rng.bit_generator.state
        patch._advance_discrete_embryo_transitions(full, np.array([0.15]), np.array([0.05]), 0.0, 5.0)

        endpoint = make_state()
        wait = patch.next_embryo_transition_wait_cycles(
            endpoint, np.array([0.15]), np.array([0.05])
        )
        self.assertEqual(wait, 3.65)
        patch._advance_discrete_embryo_transitions(
            endpoint, np.array([0.15]), np.array([0.05]), 0.0, wait
        )
        self.assertEqual(int(endpoint.site_status[0]), _SITE_STABLE)
        self.assertEqual(float(endpoint.site_stable_cycle[0]), float(full.site_stable_cycle[0]))
        self.assertEqual(float(endpoint.site_transition_cumulative_hazard[0]), 0.73)
        self.assertEqual(endpoint.site_transition_threshold[0], full.site_transition_threshold[0])
        self.assertEqual(endpoint.site_transition_outcome_uniform[0], full.site_transition_outcome_uniform[0])
        self.assertEqual(patch._event_rng.bit_generator.state, event_rng_before)

    def test_periodic_large_n_acceleration_matches_partitioned_kernel(self):
        patch = _TransitionHarness()
        patch.cfg.delivery_memory_s = 1e-3
        delivery = np.array([[2e-4], [7e-5], [3e-4], [1e-5]])
        cleavage = np.array([[2e-3], [1e-3], [4e-3], [8e-4]])

        def run(memory_log, cycle, dN):
            patch._v9_log_memory_context = memory_log
            patch._v9_phase_cycle_context = cycle
            result = patch._phase_resolved_delivery_nucleation(
                delivery, cleavage, 1000.0, np.zeros(1), dN=dN
            )
            out = patch._v9_log_memory_result.copy()
            hazard = patch._v9_log_birth_increment.copy()
            for name in ("_v9_log_memory_context", "_v9_phase_cycle_context", "_v9_log_memory_result", "_v9_log_birth_increment", "_v9_periodic_remainder_bound"):
                if hasattr(patch, name):
                    delattr(patch, name)
            return out, hazard, result

        accelerated_memory, accelerated_hazard, _ = run(np.array([-np.inf]), 0.0, 2000.0)
        first_memory, first_hazard, _ = run(np.array([-np.inf]), 0.0, 1000.0)
        second_memory, second_hazard, _ = run(first_memory, 1000.0, 1000.0)
        np.testing.assert_allclose(accelerated_memory, second_memory, atol=2e-12)
        np.testing.assert_allclose(
            accelerated_hazard, np.logaddexp(first_hazard, second_hazard), atol=2e-10
        )

    def test_periodic_acceleration_matches_direct_subthreshold_spans(self):
        patch = _TransitionHarness()
        patch.cfg.delivery_memory_s = 1e-3
        delivery = np.array([[2e-4], [7e-5], [3e-4], [1e-5]])
        cleavage = np.array([[2e-3], [1e-3], [4e-3], [8e-4]])

        def advance(memory, cycle, dN):
            patch._v9_log_memory_context = memory
            patch._v9_phase_cycle_context = cycle
            patch._phase_resolved_delivery_nucleation(
                delivery, cleavage, 1000.0, np.zeros(1), dN=dN
            )
            result = patch._v9_log_memory_result.copy(), patch._v9_log_birth_increment.copy()
            for name in ("_v9_log_memory_context", "_v9_phase_cycle_context", "_v9_log_memory_result", "_v9_log_birth_increment", "_v9_periodic_remainder_bound"):
                if hasattr(patch, name):
                    delattr(patch, name)
            return result

        accelerated_memory, accelerated_hazard = advance(np.array([-np.inf]), 0.0, 100.0)
        memory = np.array([-np.inf]); hazard = np.array([-np.inf])
        for start in (0.0, 25.0, 50.0, 75.0):
            memory, increment = advance(memory, start, 25.0)
            hazard = np.logaddexp(hazard, increment)
        np.testing.assert_allclose(accelerated_memory, memory, atol=2e-12)
        np.testing.assert_allclose(accelerated_hazard, hazard, atol=2e-10)


if __name__ == "__main__":
    unittest.main()
