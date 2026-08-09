from __future__ import annotations

import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path

from run_stateful_pd_v8_3_stress_sweep import Job, _job_paths, main as scheduler_main
import summarize_stateful_pd_v8_3_stress_sweep as summary


class StatefulPDV83SweepTests(unittest.TestCase):
    def test_independent_job_paths_are_case_and_stress_specific(self):
        root = Path('/tmp/root')
        a = _job_paths(root, Job(1, 900.0, 'no_shield'))
        b = _job_paths(root, Job(1, 900.0, 'shielded'))
        c = _job_paths(root, Job(1, 800.0, 'no_shield'))
        self.assertNotEqual(a[0], b[0])
        self.assertNotEqual(a[0], c[0])
        self.assertIn('sigmaA_900MPa', str(a[1]))

    def _write_campaign(self, root: Path, statuses):
        jobs = []
        proc = {}
        for i, (case, stress, status) in enumerate(statuses):
            job = Job(1, float(stress), case)
            jobs.append(job.__dict__ | {'key': job.key})
            proc[job.key] = {'process_status': 'completed', 'return_code': 0}
            out = _job_paths(root, job)[1]
            out.mkdir(parents=True, exist_ok=True)
            payload = {
                'model': 'M', 'status': status, 'cycles_total': 1e8,
                'cycles_connected': 5e7 if status.startswith('physical_handoff') else None,
                'pd_crack_centerline_length_final_m': 1e-4,
                'pd_crack_orientation_coherence_final': 0.9,
                'pd_crack_width_ratio_final': 0.2,
                'root_radius_over_spacing_final': 12.0,
                'pd_primary_seed_reselections_final': 1,
            }
            (out / 'summary.json').write_text(json.dumps(payload))
        (root / 'campaign_definition.json').write_text(json.dumps({'jobs': jobs}))
        (root / 'job_process_status.json').write_text(json.dumps(proc))

    def test_summary_distinguishes_valid_and_geometry_limited_results(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            self._write_campaign(root, [
                ('no_shield', 900, 'physical_handoff'),
                ('shielded', 900, 'physical_handoff_geometry_saturated'),
            ])
            summary.main(['--root', str(root), '--model', 'M', '--cycles-max', '1e9'])
            data = json.loads((root / 'stress_sweep_summary.json').read_text())
            self.assertEqual(data['counts']['valid_failure'], 1)
            self.assertEqual(data['counts']['provisional_geometry_limited'], 1)

    def test_summary_marks_queued_jobs_incomplete_not_runout(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            job = Job(1, 700.0, 'shielded')
            (root / 'campaign_definition.json').write_text(json.dumps({'jobs': [job.__dict__ | {'key': job.key}]}))
            (root / 'job_process_status.json').write_text(json.dumps({job.key: {'process_status': 'queued'}}))
            summary.main(['--root', str(root), '--model', 'M', '--cycles-max', '1e9'])
            data = json.loads((root / 'stress_sweep_summary.json').read_text())
            self.assertEqual(data['counts']['incomplete'], 1)
            self.assertFalse(data['all_jobs_final'])

    def test_monotonicity_audit_flags_higher_stress_longer_life(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            jobs = []
            proc = {}
            for stress, life in [(700.0, 1e8), (800.0, 2e8)]:
                job = Job(1, stress, 'no_shield')
                jobs.append(job.__dict__ | {'key': job.key})
                proc[job.key] = {'process_status': 'completed', 'return_code': 0}
                out = _job_paths(root, job)[1]
                out.mkdir(parents=True, exist_ok=True)
                (out / 'summary.json').write_text(json.dumps({
                    'status': 'physical_handoff', 'cycles_total': life,
                    'cycles_connected': life, 'model': 'M',
                }))
            (root / 'campaign_definition.json').write_text(json.dumps({'jobs': jobs}))
            (root / 'job_process_status.json').write_text(json.dumps(proc))
            summary.main(['--root', str(root), '--model', 'M', '--cycles-max', '1e9'])
            rows = (root / 'monotonicity_audit.csv').read_text().strip().splitlines()
            self.assertEqual(len(rows), 2)

    def test_summary_rejects_model_mismatch(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            self._write_campaign(root, [('no_shield', 900, 'physical_handoff')])
            summary.main(['--root', str(root), '--model', 'EXPECTED', '--cycles-max', '1e9'])
            data = json.loads((root / 'stress_sweep_summary.json').read_text())
            self.assertEqual(data['counts']['invalid'], 1)

    def test_summary_rejects_premature_right_censor(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            self._write_campaign(root, [('no_shield', 900, 'right_censored')])
            # Helper writes N=1e8, below the requested 1e9 cap.
            summary.main(['--root', str(root), '--model', 'M', '--cycles-max', '1e9'])
            data = json.loads((root / 'stress_sweep_summary.json').read_text())
            self.assertEqual(data['counts']['invalid'], 1)

    def test_scheduler_runs_independent_jobs_high_to_low_and_waits(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td) / 'campaign'
            fake = Path(td) / 'fake_pilot.sh'
            fake.write_text("""#!/usr/bin/env bash
set -euo pipefail
out="${OUT}/${CASES}/sigmaA_${STRESSES}MPa"
mkdir -p "$out"
python - "$out/summary.json" "$CASES" "$STRESSES" <<'PY2'
import json, sys
p, case, stress = sys.argv[1:]
json.dump({"model": "SN_2D_intact_FEM_stateful_local_peridynamics_v8_3_continuous_first_passage_fixed_geometry", "status": "right_censored", "cycles_total": 1000.0, "cycles_connected": None, "case": case, "sigma_a_MPa": float(stress)}, open(p, "w"))
PY2
""")
            fake.chmod(0o755)
            summary_script = Path(summary.__file__).resolve()
            with contextlib.redirect_stdout(io.StringIO()):
                rc = scheduler_main([
                    '--root', str(root),
                    '--stresses', '700', '900', '800',
                    '--seeds', '1',
                    '--cases', 'no_shield',
                    '--cycles-max', '1000',
                    '--block-cycles', '1000',
                    '--max-blocks', '1',
                    '--max-jobs', '1',
                    '--pilot', str(fake),
                    '--summary-script', str(summary_script),
                    '--poll-seconds', '0.01',
                ])
            self.assertEqual(rc, 0)
            status = json.loads((root / 'job_process_status.json').read_text())
            ordered = sorted(status.values(), key=lambda row: row['started_epoch'])
            self.assertEqual([row['stress_MPa'] for row in ordered], [900.0, 800.0, 700.0])
            self.assertTrue(all(row['process_status'] == 'completed' for row in ordered))
            final = json.loads((root / 'stress_sweep_summary.json').read_text())
            self.assertTrue(final['all_jobs_final'])
            self.assertEqual(final['counts']['valid_censor'], 3)


if __name__ == '__main__':
    unittest.main()
