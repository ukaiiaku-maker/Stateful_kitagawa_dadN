import unittest

import numpy as np

from arrhenius_fracture.stateful_peridynamics_v8_7_local_front_spacing import (
    _SITE_EMBRYO, StatefulPDConfig,
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
    return state


class V9StatefulPDTests(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()
