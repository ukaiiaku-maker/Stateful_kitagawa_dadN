import math
import unittest

import numpy as np

from arrhenius_fracture.v9_physical_integrator import (
    PersistentCompetingClockLedger,
    advance_constant_phase_log,
    advance_constant_memory_log,
    advance_phase_history_log,
    cleavage_log_rate,
    log_completion_k2_from_log_memory,
    plastic_chain_log_rates,
)

from arrhenius_fracture.config import EV_TO_J, KB, ElasticProperties
from arrhenius_fracture.sn_arrhenius_chain import build_chain_from_namespace
from arrhenius_fracture.sn_pd2d_stateful_v8_7_generalized_features import (
    build_crack_barrier,
    build_parser,
)


class V9PhysicalIntegratorTests(unittest.TestCase):
    def test_k2_completion_remains_log_resolved_below_float_range(self):
        value = float(log_completion_k2_from_log_memory(-800.0))
        self.assertAlmostEqual(value, -1600.0 - math.log(2.0), places=12)

    def test_constant_memory_update_has_semigroup_property(self):
        initial = np.array([-1000.0, -4.0, -math.inf])
        rate = np.array([-900.0, -3.0, -700.0])
        whole = advance_constant_memory_log(initial, rate, 3.0, 0.7)
        split = initial
        for dt in (0.1, 0.4, 0.7, 1.8):
            split = advance_constant_memory_log(split, rate, dt, 0.7)
        np.testing.assert_allclose(whole, split, rtol=2e-14, atol=2e-14)

    def test_persistent_clock_crossing_is_partition_invariant(self):
        base = PersistentCompetingClockLedger(
            threshold=np.array([0.73, 2.0]),
            log_cumulative_hazard=np.full(2, -math.inf),
            active=np.ones(2, dtype=bool),
        )
        log_rate = np.log(np.array([0.2, 0.1]))
        whole, crossed_whole, fraction_whole = base.advance_constant(log_rate, 5.0)
        split = base
        event_cycle = None
        cycle = 0.0
        for dN in (0.07, 0.31, 0.9, 1.72, 2.0):
            updated, crossed, fraction = split.advance_constant(log_rate, dN)
            if crossed[0] and event_cycle is None:
                event_cycle = cycle + dN * fraction[0]
            split = updated
            cycle += dN
        self.assertTrue(crossed_whole[0])
        self.assertAlmostEqual(5.0 * fraction_whole[0], 3.65, places=14)
        self.assertAlmostEqual(event_cycle, 3.65, places=13)
        np.testing.assert_array_equal(whole.active, split.active)
        np.testing.assert_allclose(whole.log_cumulative_hazard, split.log_cumulative_hazard, atol=2e-15)

    def test_phase_hazard_is_invariant_to_more_than_tenfold_partition(self):
        memory0 = np.array([-8.0, -900.0])
        delivery = np.array([-5.0, -850.0])
        cleavage = np.array([-2.0, -30.0])
        whole_memory, whole_hazard = advance_constant_phase_log(
            memory0, delivery, cleavage, 1.0, 0.23
        )
        memory = memory0
        hazard = np.full(2, -math.inf)
        for _ in range(25):
            memory, increment = advance_constant_phase_log(
                memory, delivery, cleavage, 0.04, 0.23
            )
            hazard = np.logaddexp(hazard, increment)
        np.testing.assert_allclose(whole_memory, memory, rtol=0.0, atol=2e-14)
        np.testing.assert_allclose(whole_hazard, hazard, rtol=0.0, atol=2e-10)

    def test_ordered_phase_history_is_restart_equivalent(self):
        delivery = np.log(np.array([[1e-4, 2e-5], [3e-4, 7e-6], [2e-6, 9e-5]]))
        cleavage = np.log(np.array([[2e-2, 1e-4], [1e-2, 3e-4], [4e-2, 2e-4]]))
        initial = np.array([-math.inf, -math.inf])
        end_a, hazard_a = advance_phase_history_log(initial, delivery, cleavage, 1e-4, 2e-4)
        mid, hazard_b1 = advance_phase_history_log(initial, delivery[:1], cleavage[:1], 1e-4, 2e-4)
        end_b, hazard_b2 = advance_phase_history_log(mid, delivery[1:], cleavage[1:], 1e-4, 2e-4)
        np.testing.assert_allclose(end_a, end_b, atol=1e-14)
        np.testing.assert_allclose(hazard_a, np.logaddexp(hazard_b1, hazard_b2), atol=1e-14)

    def test_competing_choice_uses_persistent_uniform(self):
        ledger = PersistentCompetingClockLedger(
            np.ones(3), np.full(3, -math.inf), np.ones(3, dtype=bool)
        )
        choice = ledger.choose_competing_event(
            np.ones(3, dtype=bool),
            np.log(np.array([3.0, 1.0, 0.0 + 1e-200])),
            np.log(np.array([1.0, 3.0, 1.0])),
            np.array([0.70, 0.20, 0.0]),
        )
        np.testing.assert_array_equal(choice, np.array([True, True, True]))

    def test_real_v87_barriers_match_when_representable_and_do_not_clip(self):
        args = build_parser().parse_args([])
        # Match production preset application performed by the native main.
        from arrhenius_fracture.sn_pd2d_stateful_v8_7_generalized_features import apply_representative_fatigue_model
        apply_representative_fatigue_model(args)
        mat = ElasticProperties(E=args.E_GPa * 1e9, nu=args.nu, b=args.b_m, Tm=args.Tm_K)
        chain = build_chain_from_namespace(args, mat.b)
        stress = np.array([[0.0, 2.0e9], [1.0e9, 3.0e9]])
        rho = np.array([args.rho0, args.rho0 * 2.0])
        native = chain.rates(stress, rho, args.T)
        logs = plastic_chain_log_rates(chain, stress, rho, args.T)
        np.testing.assert_allclose(np.exp(logs["lambda_flow"]), native["lambda_flow"], rtol=2e-14)
        crack = build_crack_barrier(args)
        log_cleave = cleavage_log_rate(crack, stress, args.T)
        expected = math.log(crack.rate_prefactor) - crack.deltaG_eV(stress, args.T) / ((KB / EV_TO_J) * args.T)
        np.testing.assert_allclose(log_cleave, expected, rtol=0.0, atol=2e-14)
        self.assertTrue(np.all(log_cleave > -700.0))


if __name__ == "__main__":
    unittest.main()
