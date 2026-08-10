import copy
import pickle
import tempfile
import unittest
from pathlib import Path

import numpy as np

from arrhenius_fracture.v9_canonical_four_class_elastic_fem import (
    CONTRACT, FEM_MODEL_ID, CanonicalElasticFourClassFEMCondition,
)
from arrhenius_fracture.v9_canonical_four_class_fem import CanonicalFourClassFEMCondition
from arrhenius_fracture.v9_four_class_registry import EXPECTED
from arrhenius_fracture.v9_fem_transaction import FEMPhysicalState
from arrhenius_fracture.v9_quiet_tail_kernel import (
    phase_resolved_cycle, restore_condition_checkpoint,
    save_condition_checkpoint,
)


SOURCE = Path("/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1")
RUN_ARGS = Path("runs/sn_v9_single_seed_skeleton/K360_FAILURE/solver_output/shielded/sigmaA_735p921MPa/run_args.json")


class CanonicalElasticFourClassFEMTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.peak = next(option for option in EXPECTED if "peak01" in option)
        cls.dbtt = next(option for option in EXPECTED if "dbtt01" in option)

    def make(self, option=None, stress=735.9214951373579, seed=1720):
        return CanonicalElasticFourClassFEMCondition.from_run_args(
            RUN_ARGS, option or self.peak, SOURCE, stress, seed)

    def test_generated_contract_is_explicit_and_bulk_placeholders_are_neutral(self):
        condition = self.make()
        self.assertEqual(condition.summary()["model_id"], FEM_MODEL_ID)
        for key, value in CONTRACT.items():
            self.assertEqual(condition.audit[key], value)
            self.assertEqual(condition.summary()[key], value)
        condition.assert_neutral_fem_state()
        self.assertFalse(condition.fem.ep_gp.flags.writeable)
        self.assertFalse(condition.fem.rho_gp.flags.writeable)
        self.assertTrue(np.all(condition.fem.ep_gp == 0.0))
        self.assertTrue(np.all(condition.fem.epsp_acc_gp == 0.0))

    def test_cached_elastic_root_history_matches_direct_initial_fem_scaffold(self):
        elastic = self.make()
        control = CanonicalFourClassFEMCondition.from_run_args(
            RUN_ARGS, self.peak, SOURCE, elastic.sigma_a_MPa, 1720)
        np.testing.assert_array_equal(
            elastic._root_history(), control._root_history(control.fem))
        alternate_placeholder = FEMPhysicalState(
            elastic.fem.ep_gp,
            np.full_like(elastic.fem.rho_gp, 9.9e99),
            elastic.fem.epsp_acc_gp, elastic.fem.u, 0.0)
        np.testing.assert_array_equal(
            elastic._root_history(), elastic._root_history(alternate_placeholder))
        rho_before = elastic.fem.rho_gp.copy()
        elastic.birth.hazard_threshold_action = 1e100
        elastic.advance(37.0)
        np.testing.assert_array_equal(elastic.fem.ep_gp, 0.0)
        np.testing.assert_array_equal(elastic.fem.rho_gp, rho_before)
        self.assertEqual(elastic.fem.plastic_work_J_per_m, 0.0)

    def test_tip_only_prescribed_history_parity_for_peak_and_dbtt(self):
        for option in (self.peak, self.dbtt):
            condition = self.make(option=option)
            direct = condition.birth.copy()
            condition.birth.hazard_threshold_action = 1e100
            direct.hazard_threshold_action = 1e100
            condition.advance(19.0)
            for cycles in (1.0, 2.0, 4.0, 8.0, 4.0):
                direct.advance_fem_phase_block(
                    cycles, condition.args.frequency_Hz, condition.args.T,
                    condition._root_history())
            self.assertEqual(condition.birth.cumulative_cleavage_hazard,
                             direct.cumulative_cleavage_hazard)
            self.assertEqual(condition.birth.log_cumulative_cleavage_hazard,
                             direct.log_cumulative_cleavage_hazard)
            for key in condition.birth.mpz._CAPSULE_ARRAYS:
                np.testing.assert_array_equal(
                    getattr(condition.birth.mpz.state, key),
                    getattr(direct.mpz.state, key))

    def test_eleven_way_partition_converges_and_restart_is_atomic(self):
        whole, split = self.make(), self.make()
        whole.birth.hazard_threshold_action = split.birth.hazard_threshold_action = 1e100
        whole.advance(11.0)
        for _ in range(11):
            split.advance(1.0)
        self.assertAlmostEqual(
            whole.birth.cumulative_cleavage_hazard,
            split.birth.cumulative_cleavage_hazard,
            delta=max(split.birth.cumulative_cleavage_hazard * 5e-4, 1e-300))
        for key in whole.birth.mpz._CAPSULE_ARRAYS:
            np.testing.assert_allclose(
                getattr(whole.birth.mpz.state, key),
                getattr(split.birth.mpz.state, key), rtol=5e-4, atol=1e-18)
        capsule = copy.deepcopy(split.capsule())
        resumed = self.make(); resumed.restore_capsule(capsule)
        self.assertEqual(pickle.dumps(resumed.capsule()), pickle.dumps(capsule))
        split.advance(3.0); resumed.advance(3.0)
        self.assertEqual(split.birth.capsule()["hazard_rng_state"],
                         resumed.birth.capsule()["hazard_rng_state"])
        self.assertEqual(split.birth.cumulative_cleavage_hazard,
                         resumed.birth.cumulative_cleavage_hazard)

    def test_initial_energy_gate_matches_zero_ep_control(self):
        elastic = self.make()
        control = CanonicalFourClassFEMCondition.from_run_args(
            RUN_ARGS, self.peak, SOURCE, elastic.sigma_a_MPa, 1720)
        threshold = elastic.birth.hazard_threshold_action
        pending = {"attempt_index": 0, "threshold_action": threshold,
                   "cumulative_cleavage_hazard": threshold, "time_s": 0.0}
        elastic.birth.pending_attempt = copy.deepcopy(pending)
        control.birth.pending_attempt = copy.deepcopy(pending)
        a, b = elastic._event_gate(), control._event_gate()
        for key in ("stochastic_proposed_event_length_m",
                    "committed_event_length_m", "arrest_reason"):
            self.assertEqual(a[key], b[key])

    def test_atomic_array_checkpoint_restores_v3_without_fem_state_arrays(self):
        condition = self.make(); condition.birth.hazard_threshold_action = 1e100
        condition.advance(7.0)
        kernel = phase_resolved_cycle(condition)
        self.assertGreaterEqual(kernel["cycle_hazard"], 0.0)
        with tempfile.TemporaryDirectory() as directory:
            generation = save_condition_checkpoint(
                condition, Path(directory), kernel, {"N": condition.cycles})
            resumed = self.make()
            restored_kernel, summary, manifest = restore_condition_checkpoint(
                resumed, Path(directory), generation)
        self.assertEqual(summary["N"], 7.0)
        self.assertEqual(manifest["generation"], generation)
        np.testing.assert_array_equal(
            restored_kernel["root_tensors_Pa"], kernel["root_tensors_Pa"])
        self.assertEqual(pickle.dumps(resumed.birth.capsule()),
                         pickle.dumps(condition.birth.capsule()))


if __name__ == "__main__":
    unittest.main()
