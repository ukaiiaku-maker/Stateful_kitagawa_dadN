import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from scripts.plot_v9_four_class_progress import verified_conditioned_capsules


class ConditionedProgressReconstructionTests(unittest.TestCase):
    def _capsule(self, root: Path, cls: str, name: str, *, valid_hash=True,
                 independent=True):
        directory = root / cls / "conditioned" / name
        directory.mkdir(parents=True)
        capsule = {
            "schema": "V9_CONDITIONED_ATTEMPT_PREMARK_CAPSULE_1",
            "source_generation": "generation_source",
            "replay_generation": (
                "generation_replay" if independent else "generation_source"
            ),
        }
        payload = (json.dumps(capsule, sort_keys=True) + "\n").encode()
        (directory / "conditioned_premark_capsule.json").write_bytes(payload)
        digest = hashlib.sha256(payload).hexdigest()
        if not valid_hash:
            digest = "0" * 64
        (directory / "manifest.json").write_text(
            json.dumps({"capsule_sha256": digest}) + "\n"
        )

    def test_counts_only_hash_verified_independent_capsules(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self._capsule(root, "ceramic", "premark_1000_N50_v1")
            self._capsule(root, "ceramic", "premark_1100_N50_v1", valid_hash=False)
            self._capsule(root, "ceramic", "premark_1200_N50_v1", independent=False)
            self._capsule(root, "ceramic", "premark_invalid_N50_v1")
            result = verified_conditioned_capsules(root)
            self.assertEqual(len(result["ceramic"]), 1)
            self.assertEqual(result["ceramic"][0].parent.name,
                             "premark_1000_N50_v1")


if __name__ == "__main__":
    unittest.main()
