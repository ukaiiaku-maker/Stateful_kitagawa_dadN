import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from arrhenius_fracture.v9_array_codec import ArrayGenerationError, AtomicArrayGenerationStore


class V9ArrayCodecTests(unittest.TestCase):
    def test_atomic_roundtrip_preserves_dtype_shape_and_hashes(self):
        with tempfile.TemporaryDirectory() as td:
            store = AtomicArrayGenerationStore(Path(td))
            generation = store.write(
                {"ep_gp": np.arange(6, dtype=np.float64).reshape(3, 2), "site_status": np.array([0, 2], dtype=np.uint8)},
                {"cycle": 2.0, "schema_role": "historical_physical_state"},
                {"status": "checkpoint"},
            )
            arrays, metadata, summary, manifest = store.load()
        self.assertTrue(generation.startswith("generation_"))
        self.assertEqual(arrays["ep_gp"].dtype, np.dtype("float64"))
        self.assertEqual(arrays["ep_gp"].shape, (3, 2))
        self.assertEqual(metadata["cycle"], 2.0)
        self.assertEqual(summary["status"], "checkpoint")
        self.assertEqual(manifest["generation"], generation)

    def test_incomplete_generation_is_rejected_and_not_activated(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            (root / "generation_incomplete").mkdir()
            store = AtomicArrayGenerationStore(root)
            with self.assertRaises(ArrayGenerationError):
                store.load("generation_incomplete")
            with self.assertRaisesRegex(ArrayGenerationError, "activated"):
                store.load()

    def test_component_corruption_is_detected(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td); store = AtomicArrayGenerationStore(root)
            generation = store.write({"x": np.arange(3)}, {}, {})
            (root / generation / "summary.json").write_text(json.dumps({"corrupt": True}))
            with self.assertRaisesRegex(ArrayGenerationError, "hash"):
                store.load(generation)

    def test_pickle_object_arrays_are_forbidden(self):
        with tempfile.TemporaryDirectory() as td:
            with self.assertRaisesRegex(ArrayGenerationError, "pickle"):
                AtomicArrayGenerationStore(Path(td)).write({"bad": np.array([object()], dtype=object)}, {}, {})


if __name__ == "__main__":
    unittest.main()
