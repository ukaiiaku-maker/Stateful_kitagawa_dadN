import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from arrhenius_fracture.sn_pd2d_stateful_v8_7_generalized_features import MODEL_ID, SOURCE_SHA256, _CHECKPOINT_VERSION
from arrhenius_fracture.v9_v87_adapter import (
    BASE_ARRAYS, PD_BOND_ARRAYS, PD_PATH_ARRAYS, PD_POINT_ARRAYS, PD_SITE_ARRAYS,
    PD_VECTOR_ARRAYS, REQUIRED_PD_SCALARS, V87CheckpointValidationError,
    adapt_v87_checkpoint,
    write_v87_roundtrip_checkpoint,
)


class V87AdapterTests(unittest.TestCase):
    def make_checkpoint(self, root: Path, mutate=None):
        nn, ne, npoint, nsite, nbond = 4, 2, 3, 2, 4
        arrays = {
            "mesh_nodes": np.array([[0., 0.], [1., 0.], [0., 1.], [1., 1.]]),
            "root_xy": np.array([0., 0.]), "ep_gp": np.zeros((3, ne)),
            "rho_gp": np.full(ne, 1e12), "epsp_acc_gp": np.zeros(ne),
            "u": np.zeros(2 * nn), "last_residual": np.zeros(nn),
        }
        integer_points = {"candidate_sites", "available_sites", "embryo_sites", "stable_sites", "inactive_sites", "born_sites_cumulative", "healed_sites_cumulative"}
        for name in PD_POINT_ARRAYS:
            arrays[f"pd__{name}"] = np.zeros(npoint, dtype=np.int64 if name in integer_points else float)
        arrays["pd__candidate_sites"][:] = [1, 1, 0]
        arrays["pd__available_sites"][:] = [1, 1, 0]
        arrays["pd__available"][:] = 1.0
        arrays["pd__crack_normal_c2"][:] = 1.0
        arrays["pd__site_node_index"] = np.array([0, 1], dtype=np.int64)
        arrays["pd__site_status"] = np.zeros(nsite, dtype=np.uint8)
        arrays["pd__site_birth_threshold"] = np.array([0.5, 1.5])
        arrays["pd__site_birth_cycle"] = np.full(nsite, np.nan)
        arrays["pd__site_stable_cycle"] = np.full(nsite, np.nan)
        for name in PD_BOND_ARRAYS:
            arrays[f"pd__{name}"] = np.zeros(nbond, dtype=bool if name != "bond_damage" else float)
        for name in PD_VECTOR_ARRAYS:
            arrays[f"pd__{name}"] = np.zeros(2)
        for name in PD_PATH_ARRAYS:
            arrays[f"pd__{name}"] = np.empty((0, 2))
        scalars = {name: None for name in REQUIRED_PD_SCALARS}
        scalars.update({
            "primary_seed_node": -1, "primary_seed_stall_updates": 0,
            "primary_seed_reselections": 0, "primary_seed_last_progress": 0.0,
            "primary_seed_selected_cycles": 0.0, "active_front": False,
            "active_front_normal_c2": 1.0, "active_front_normal_s2": 0.0,
            "active_front_length_m": 0.0, "front_candidate_mode": 0,
            "front_eligible_preferred": 0, "front_eligible_fallback": 0,
            "front_stall_updates": 0, "front_last_advance_cycles": 0.0,
            "front_max_link_rate_per_cycle": 0.0, "diffuse_bad_updates": 0,
        })
        signature = {"T": 300.0, "frequency_Hz": 1000.0, "case": "shielded", "sigma_a_MPa": 700.0}
        rng = {"bit_generator": "PCG64", "state": {"state": 1, "inc": 3}, "has_uint32": 0, "uinteger": 0}
        metadata = {
            "checkpoint_version": _CHECKPOINT_VERSION, "model_id": MODEL_ID,
            "source_sha256": SOURCE_SHA256, "signature": signature,
            "next_block": 2, "cycles": 2.0, "Wp_total": 0.0,
            "pd_scalars": scalars, "candidate_rng_state": rng,
            "event_rng_state": dict(rng), "rows": [],
        }
        if mutate:
            mutate(arrays, metadata)
        arrays["metadata_json"] = np.asarray(json.dumps(metadata))
        path = root / "checkpoint.npz"
        np.savez_compressed(path, **arrays)
        return path, signature

    def test_exact_adapter_preserves_historical_arrays_clocks_and_rng(self):
        with tempfile.TemporaryDirectory() as td:
            path, signature = self.make_checkpoint(Path(td))
            adapted = adapt_v87_checkpoint(path, immutable_signature=signature, expected_cycle=2.0, expected_next_block=2)
        self.assertEqual(set(adapted.arrays), BASE_ARRAYS | {f"pd__{n}" for n in PD_POINT_ARRAYS | PD_SITE_ARRAYS | PD_BOND_ARRAYS | PD_VECTOR_ARRAYS | PD_PATH_ARRAYS})
        self.assertEqual(len(adapted.capsule.clocks), 2)
        self.assertEqual(adapted.capsule.clocks[0].threshold, 0.5)
        self.assertIn("candidate", json.loads(adapted.capsule.rng_states_json))
        self.assertIn("controller_next_block_cycles", adapted.new_v9_controller_state)
        self.assertTrue(all(not array.flags.writeable for array in adapted.arrays.values()))
        with tempfile.TemporaryDirectory() as td:
            roundtrip = Path(td) / "roundtrip.npz"
            write_v87_roundtrip_checkpoint(adapted, roundtrip)
            again = adapt_v87_checkpoint(roundtrip, immutable_signature=signature, expected_cycle=2.0, expected_next_block=2)
        for name in adapted.arrays:
            self.assertTrue(np.array_equal(adapted.arrays[name], again.arrays[name], equal_nan=True), name)
        self.assertEqual(adapted.metadata, again.metadata)

    def test_rejects_wrong_hash_model_version_signature_and_cycle(self):
        mutations = [
            lambda a, m: m.update(checkpoint_version=99),
            lambda a, m: m.update(model_id="wrong"),
            lambda a, m: m.update(source_sha256={"driver": "wrong"}),
        ]
        for mutate in mutations:
            with self.subTest(mutate=mutate), tempfile.TemporaryDirectory() as td:
                path, signature = self.make_checkpoint(Path(td), mutate)
                with self.assertRaises(V87CheckpointValidationError):
                    adapt_v87_checkpoint(path, immutable_signature=signature)
        with tempfile.TemporaryDirectory() as td:
            path, signature = self.make_checkpoint(Path(td))
            with self.assertRaisesRegex(V87CheckpointValidationError, "signature"):
                adapt_v87_checkpoint(path, immutable_signature=signature | {"T": 301.0})
            with self.assertRaisesRegex(V87CheckpointValidationError, "cycle"):
                adapt_v87_checkpoint(path, immutable_signature=signature, expected_cycle=3.0)

    def test_rejects_missing_threshold_rng_array_unknown_field_and_shape(self):
        def missing_threshold(arrays, metadata): arrays.pop("pd__site_birth_threshold")
        def missing_rng(arrays, metadata): metadata.pop("event_rng_state")
        def wrong_shape(arrays, metadata): arrays["rho_gp"] = np.zeros(3)
        def unknown(arrays, metadata): arrays["surprise"] = np.zeros(1)
        for mutation in (missing_threshold, missing_rng, wrong_shape, unknown):
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as td:
                path, signature = self.make_checkpoint(Path(td), mutation)
                with self.assertRaises(V87CheckpointValidationError):
                    adapt_v87_checkpoint(path, immutable_signature=signature)


if __name__ == "__main__":
    unittest.main()
