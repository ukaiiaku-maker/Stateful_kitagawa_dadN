#!/usr/bin/env python3
"""Bounded independent-job scheduler for Stateful-PD v8.3 stress sweeps."""
from __future__ import annotations

import argparse
import json
import os
import signal
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

MODEL_ID = "SN_2D_intact_FEM_stateful_local_peridynamics_v8_3_continuous_first_passage_fixed_geometry"


@dataclass(frozen=True)
class Job:
    seed: int
    stress: float
    case: str

    @property
    def key(self) -> str:
        return f"seed_{self.seed}__{self.case}__{self.stress:g}MPa"


def _atomic_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(tmp, path)


def _job_paths(root: Path, job: Job) -> tuple[Path, Path, Path]:
    stem = root / f"seed_{job.seed}" / f"stress_{job.stress:g}MPa" / f"job_{job.case}"
    case_dir = stem / job.case / (f"sigmaA_{job.stress:g}MPa".replace(".", "p"))
    log = root / "logs" / f"{job.key}.log"
    return stem, case_dir, log


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", required=True)
    p.add_argument("--stresses", nargs="+", required=True, type=float)
    p.add_argument("--seeds", nargs="+", required=True, type=int)
    p.add_argument("--cases", nargs="+", required=True, choices=["no_shield", "shielded"])
    p.add_argument("--cycles-max", required=True, type=float)
    p.add_argument("--block-cycles", required=True, type=float)
    p.add_argument("--max-blocks", required=True, type=int)
    p.add_argument("--max-jobs", type=int, default=2)
    p.add_argument("--mesh-seed", type=int, default=1)
    p.add_argument("--resolution-profile", default="h15")
    p.add_argument("--checkpoint-every-blocks", type=int, default=10)
    p.add_argument("--snapshot-every", type=int, default=0)
    p.add_argument("--print-every", type=int, default=10)
    p.add_argument("--geometry-limit-mode", choices=["terminate", "freeze", "off"], default="terminate")
    p.add_argument("--pilot", default="run_sn_stateful_pd_v8_3_h15_pilot.sh")
    p.add_argument("--summary-script", default="summarize_stateful_pd_v8_3_stress_sweep.py")
    p.add_argument("--poll-seconds", type=float, default=1.0)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    if args.max_jobs < 1:
        raise SystemExit("--max-jobs must be positive")
    root = Path(args.root).resolve()
    project = Path.cwd().resolve()
    pilot = project / args.pilot
    summary_script = project / args.summary_script
    if not pilot.is_file() or not summary_script.is_file():
        raise SystemExit("v8.3 pilot or summary script is missing")
    root.mkdir(parents=True, exist_ok=True)
    (root / "logs").mkdir(exist_ok=True)

    stresses = sorted(set(float(x) for x in args.stresses), reverse=True)
    jobs = [Job(seed, stress, case) for seed in args.seeds for stress in stresses for case in args.cases]
    campaign = {
        "campaign": "stateful_pd_v8_3_fixed_geometry_first_passage_stress_jobs",
        "model": MODEL_ID,
        "project_root": str(project),
        "stresses_MPa_scheduled_high_to_low": stresses,
        "pd_seeds": list(args.seeds),
        "cases": list(args.cases),
        "cycles_max": args.cycles_max,
        "block_cycles_ceiling": args.block_cycles,
        "max_blocks": args.max_blocks,
        "max_concurrent_jobs": args.max_jobs,
        "mesh_seed": args.mesh_seed,
        "resolution_profile": args.resolution_profile,
        "geometry_limit_mode": args.geometry_limit_mode,
        "geometry_policy": "fixed_production_baseline",
        "birth_clock_mode": "persistent_site_first_passage_aligned",
        "one_process_per_case_stress_seed": True,
        "valid_failure_statuses": ["physical_handoff"],
        "valid_censoring_statuses": ["right_censored"],
        "provisional_geometry_limited_statuses": [
            "physical_handoff_geometry_saturated",
            "right_censored_geometry_saturated",
        ],
        "invalid_statuses": [
            "geometry_invalid_underresolved",
            "morphology_invalid_precapture_stalled",
            "morphology_invalid_stalled",
            "morphology_invalid_diffuse",
            "max_blocks_reached",
        ],
        "jobs": [job.__dict__ | {"key": job.key} for job in jobs],
    }
    _atomic_json(root / "campaign_definition.json", campaign)

    status_path = root / "job_process_status.json"
    job_status: dict[str, dict] = {
        job.key: {"seed": job.seed, "stress_MPa": job.stress, "case": job.case, "process_status": "queued"}
        for job in jobs
    }
    _atomic_json(status_path, job_status)

    pending = list(jobs)
    running: dict[str, tuple[subprocess.Popen, object, Job]] = {}
    failures = 0

    def summarize() -> None:
        subprocess.run(
            [sys.executable, str(summary_script), "--root", str(root), "--model", MODEL_ID, "--cycles-max", str(args.cycles_max)],
            check=False,
        )

    def terminate_all() -> None:
        active = list(running.values())
        for proc, _, _ in active:
            if proc.poll() is None:
                try:
                    proc.terminate()
                except Exception:
                    pass
        deadline = time.time() + 5.0
        for proc, stream, _ in active:
            if proc.poll() is None:
                try:
                    proc.wait(timeout=max(deadline - time.time(), 0.1))
                except subprocess.TimeoutExpired:
                    try:
                        proc.kill()
                        proc.wait(timeout=2.0)
                    except Exception:
                        pass
                except Exception:
                    pass
            try:
                stream.close()
            except Exception:
                pass

    def _interrupt_handler(_signum, _frame):
        raise KeyboardInterrupt

    old_handlers = {}
    for sig in (signal.SIGINT, signal.SIGTERM):
        old_handlers[sig] = signal.getsignal(sig)
        signal.signal(sig, _interrupt_handler)

    try:
        while pending or running:
            while pending and len(running) < args.max_jobs:
                job = pending.pop(0)
                out, case_dir, log = _job_paths(root, job)
                out.mkdir(parents=True, exist_ok=True)
                log.parent.mkdir(parents=True, exist_ok=True)
                env = os.environ.copy()
                env.update({
                    "OUT": str(out),
                    "RESOLUTION_PROFILE": args.resolution_profile,
                    "STRESSES": f"{job.stress:g}",
                    "CASES": job.case,
                    "CYCLES_MAX": str(args.cycles_max),
                    "BLOCK_CYCLES": str(args.block_cycles),
                    "MAX_BLOCKS": str(args.max_blocks),
                    "MESH_SEED": str(args.mesh_seed),
                    "PD_SEED": str(job.seed),
                    "CHECKPOINT_EVERY_BLOCKS": str(args.checkpoint_every_blocks),
                    "SNAPSHOT_EVERY": str(args.snapshot_every),
                    "PRINT_EVERY": str(args.print_every),
                    "GEOMETRY_LIMIT_MODE": args.geometry_limit_mode,
                    "ENABLE_GEOMETRY_EVOLUTION": "0",
                    "RESUME": "1",
                    "SKIP_EXISTING": "1",
                    "OMP_NUM_THREADS": env.get("OMP_NUM_THREADS_PER_JOB", "1"),
                    "OPENBLAS_NUM_THREADS": env.get("OPENBLAS_NUM_THREADS_PER_JOB", "1"),
                    "MKL_NUM_THREADS": env.get("MKL_NUM_THREADS_PER_JOB", "1"),
                    "VECLIB_MAXIMUM_THREADS": env.get("VECLIB_MAXIMUM_THREADS_PER_JOB", "1"),
                    "NUMEXPR_NUM_THREADS": env.get("NUMEXPR_NUM_THREADS_PER_JOB", "1"),
                })
                stream = log.open("ab", buffering=0)
                header = f"\n=== launch {job.key} at {time.strftime('%Y-%m-%d %H:%M:%S')} ===\n".encode()
                stream.write(header)
                proc = subprocess.Popen(["bash", str(pilot)], cwd=project, env=env, stdout=stream, stderr=subprocess.STDOUT)
                running[job.key] = (proc, stream, job)
                job_status[job.key].update({
                    "process_status": "running",
                    "pid": proc.pid,
                    "log": str(log),
                    "output": str(case_dir),
                    "started_epoch": time.time(),
                })
                _atomic_json(status_path, job_status)
                print(f"launched {job.key} pid={proc.pid}", flush=True)

            completed = []
            for key, (proc, stream, job) in list(running.items()):
                rc = proc.poll()
                if rc is None:
                    continue
                stream.write(f"\n=== exit {job.key} rc={rc} at {time.strftime('%Y-%m-%d %H:%M:%S')} ===\n".encode())
                stream.close()
                job_status[key].update({
                    "process_status": "completed" if rc == 0 else "process_failed",
                    "return_code": rc,
                    "finished_epoch": time.time(),
                })
                if rc != 0:
                    failures += 1
                completed.append(key)
                print(f"finished {key} rc={rc}", flush=True)
            for key in completed:
                running.pop(key, None)
            if completed:
                _atomic_json(status_path, job_status)
                summarize()
            if running and not completed:
                time.sleep(max(args.poll_seconds, 0.1))
    except KeyboardInterrupt:
        print("interrupt received; terminating active v8.3 jobs", file=sys.stderr)
        terminate_all()
        for key in running:
            job_status[key]["process_status"] = "interrupted"
        _atomic_json(status_path, job_status)
        summarize()
        return 130
    except Exception:
        terminate_all()
        for key in running:
            job_status[key]["process_status"] = "scheduler_exception"
        _atomic_json(status_path, job_status)
        summarize()
        raise
    finally:
        for sig, handler in old_handlers.items():
            signal.signal(sig, handler)

    _atomic_json(status_path, job_status)
    summarize()
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
