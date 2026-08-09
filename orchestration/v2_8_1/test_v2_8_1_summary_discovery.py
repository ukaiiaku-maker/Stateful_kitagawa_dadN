from __future__ import annotations
import argparse
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
RUNNER = HERE / "run_stateful_pd_kitagawa_gated_pilot_v2_8.py"
spec = importlib.util.spec_from_file_location("pilot_v281", RUNNER)
mod = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules[spec.name] = mod
spec.loader.exec_module(mod)


class V281Tests(unittest.TestCase):
    def args(self):
        return argparse.Namespace(max_blocks=10, checkpoint_every_blocks=1, print_every=1)

    def payload(self, condition, status="physical_handoff"):
        return {
            "status": status,
            "sigma_a_MPa": condition.sigma_a_MPa,
            "cycles_total": 3.7e6 if status == "physical_handoff" else condition.cycles_max,
            "cycles_connected": 3.7e6 if status == "physical_handoff" else None,
            "pd_handoff_pass_final": True,
            "pd_production_preflight": {
                "pass": True,
                "point_spacing_root_m": 7.5e-6,
                "point_spacing_global_m": 37e-6,
            },
            "geometry_resolution_audit_final": {"pass": True},
            "pd_crack_axial_coverage_final": 1.0,
            "pd_crack_width_ratio_final": 0.24,
            "pd_crack_orientation_coherence_final": 0.98,
            "pd_handoff_offfront_broken_fraction_final": 0.0,
        }

    def make_existing_job(self, root: Path, condition, rounded_tag: str):
        job = root / condition.condition_id
        solver = job / "solver_output" / "shielded" / rounded_tag
        solver.mkdir(parents=True)
        (solver / "summary.json").write_text(json.dumps(self.payload(condition)))
        (solver / "sn_stateful_pd_history.csv").write_text("cycles_total,x\n0,0\n1,1\n")
        request = {
            "pilot_id": mod.PILOT_ID,
            "solver_model": mod.SOLVER_MODEL,
            "candidate": {"candidate_id": mod.BASELINE.candidate_id, "values": dict(mod.BASELINE.values)},
            "condition": mod.asdict(condition),
            "max_blocks": 10,
            "checkpoint_every_blocks": 1,
            "print_every": 1,
        }
        request["request_sha256"] = mod.canonical_sha(request)
        (job / "REQUEST_V2_8.json").write_text(json.dumps(request))
        (job / "EXIT_V2_8.json").write_text(json.dumps({
            "request_sha256": request["request_sha256"], "return_code": 0
        }))
        (job / "run.log").write_text("completed\n")
        return solver / "summary.json"

    def test_discovers_solver_six_digit_sigma_tag(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            condition = mod.CONDITIONS["failure"]
            expected = self.make_existing_job(root, condition, "sigmaA_735p921MPa")
            found = mod.discover_summary(root / condition.condition_id / "solver_output", condition)
            self.assertEqual(found, expected)

    def test_existing_completed_job_is_reused_without_solver(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            condition = mod.CONDITIONS["failure"]
            self.make_existing_job(root, condition, "sigmaA_735p921MPa")
            # Standalone deliberately has no importable solver. Reuse must occur
            # before any subprocess launch.
            standalone = root / "standalone"
            standalone.mkdir()
            result = mod.run_condition(root, standalone, condition, self.args())
            self.assertTrue(result.reused)
            self.assertEqual(result.classification, "failure")
            ok, errors, diagnostics = mod.validate_gate(result, condition)
            self.assertTrue(ok, errors)
            self.assertGreater(diagnostics["handoff_cycle_ratio_to_frozen_baseline"], 1.0)

    def test_mismatched_sigma_summary_is_not_reused(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            condition = mod.CONDITIONS["failure"]
            case = root / "solver_output" / "shielded" / "sigmaA_735p921MPa"
            case.mkdir(parents=True)
            payload = self.payload(condition)
            payload["sigma_a_MPa"] = condition.sigma_a_MPa + 1.0
            (case / "summary.json").write_text(json.dumps(payload))
            found = mod.discover_summary(root / "solver_output", condition)
            self.assertFalse(found.is_file())

    def test_censor_rounded_tag_discovery(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            condition = mod.CONDITIONS["censor"]
            case = root / "solver_output" / "shielded" / "sigmaA_690p443MPa"
            case.mkdir(parents=True)
            payload = self.payload(condition, status="right_censored")
            (case / "summary.json").write_text(json.dumps(payload))
            found = mod.discover_summary(root / "solver_output", condition)
            self.assertEqual(found, case / "summary.json")


if __name__ == "__main__":
    unittest.main(verbosity=2)
