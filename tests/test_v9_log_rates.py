import math
import unittest

from arrhenius_fracture.v9_log_rates import (
    LOG_MIN_SUBNORMAL,
    arrhenius_log_rate,
    classify_exact_power_law_from_logs,
    completion_gated_log_rate,
    logistic_log_rate,
    series_log_rate,
    trace_completion_gated_rate_paths,
)


class V9LogRateTests(unittest.TestCase):
    def test_arrhenius_rate_below_float_range_remains_physical(self):
        rate = arrhenius_log_rate("cleavage", 1.0e11, 100.0, 300.0)
        self.assertLess(rate.log_rate, LOG_MIN_SUBNORMAL)
        self.assertIsNone(rate.physical_rate)
        self.assertFalse(rate.constitutive_zero)
        self.assertEqual(rate.numerical_representability, "underflow_in_float64_rate_space")
        self.assertFalse(rate.diagnostic()["numerical_floor_used_as_physics"])

    def test_exact_extinction_requires_constitutive_reason(self):
        rate = arrhenius_log_rate("delivery", 1.0e11, 1.0, 300.0, constitutive_enabled=False)
        self.assertTrue(rate.constitutive_zero)
        self.assertEqual(rate.physical_rate, 0.0)
        self.assertEqual(rate.extinction_reason, "mechanism_disabled")

    def test_series_rate_uses_log_residence_times(self):
        a = arrhenius_log_rate("a", 1.0e11, 40.0, 300.0)
        b = arrhenius_log_rate("b", 1.0e11, 41.0, 300.0)
        combined = series_log_rate("flow", (a, b))
        self.assertTrue(math.isfinite(combined.log_rate))
        self.assertLessEqual(combined.log_rate, min(a.log_rate, b.log_rate))

    def test_completion_zero_is_physical_gate_not_underflow(self):
        cleavage = arrhenius_log_rate("cleavage", 1.0e11, 1.0, 300.0)
        birth = completion_gated_log_rate(cleavage, 0.0)
        self.assertTrue(birth.constitutive_zero)
        self.assertEqual(birth.extinction_reason, "exact_zero_completion")

    def test_progression_logistic_tail_does_not_clip(self):
        rate = logistic_log_rate("stabilization", 500.0, 1000.0, -2000.0)
        self.assertEqual(rate.numerical_representability, "underflow_in_float64_rate_space")
        self.assertFalse(rate.constitutive_zero)

    def test_log_rate_endurance_invariant_to_presentation_offset(self):
        x = [math.log(10.0 ** k) for k in range(1, 10)]
        y = [-1500.0 - 1.5 * value for value in x]
        a = classify_exact_power_law_from_logs(x, y)
        b = classify_exact_power_law_from_logs(x, [value + 1200.0 for value in y])
        self.assertEqual(a[0], "endurance_supported")
        self.assertEqual(a[0], b[0])
        self.assertAlmostEqual(a[1], b[1], places=12)

    def test_shortened_v87_default_rate_path_trace(self):
        from arrhenius_fracture.sn_arrhenius_chain import build_chain_from_namespace
        from arrhenius_fracture.sn_pd2d_stateful_v8_7_generalized_features import (
            ScratchExpFloorBarrier, apply_representative_fatigue_model, build_parser,
        )

        args = build_parser().parse_args([])
        apply_representative_fatigue_model(args)
        chain = build_chain_from_namespace(args, args.b_m)
        stress = 1.0e9
        rho = 1.0e12
        phi = float(chain.taylor_phi(rho))
        crack = ScratchExpFloorBarrier(args)
        trace = trace_completion_gated_rate_paths(
            plastic_barriers=(
                ("plastic_emission", chain.emit.rate_prefactor, float(chain.emit.deltaG_eV(stress, args.T))),
                ("plastic_peierls", chain.peierls.rate_prefactor, float(chain.peierls.deltaG_eV(stress, args.T))),
                ("plastic_taylor", chain.taylor.rate_prefactor, float(chain.taylor.deltaG_eV(phi * stress, args.T))),
            ),
            cleavage_prefactor=crack.rate_prefactor,
            cleavage_barrier_eV=float(crack.deltaG_eV(stress, args.T)),
            temperature_K=args.T,
            frequency_Hz=args.frequency_Hz,
            completion_probability=0.1,
            stabilization_prefactor_s=args.nu_stabilize_s,
            healing_prefactor_s=args.nu_heal_s,
            stabilization_z=-2.0,
            growth_prefactor_s=args.nu_grow_s,
            growth_z=-3.0,
            stable_activity=0.2,
            linkage_prefactor_s=args.nu_link_s,
            linkage_z=-4.0,
            directional_activity=0.3,
            activations=("v87_crack_floor_constitutive",),
        )
        expected = {
            "plastic_emission", "plastic_peierls", "plastic_taylor", "plastic_escape",
            "plastic_completed_delivery", "cleavage", "birth", "stabilization", "healing",
            "stable_growth", "linkage_front",
        }
        self.assertEqual(set(trace), expected)
        for rate in trace.values():
            diagnostic = rate.diagnostic()
            self.assertIn("log_rate", diagnostic)
            self.assertIn("physical_rate", diagnostic)
            self.assertIn("numerical_representability", diagnostic)
            self.assertFalse(diagnostic["numerical_floor_used_as_physics"])


if __name__ == "__main__":
    unittest.main()
