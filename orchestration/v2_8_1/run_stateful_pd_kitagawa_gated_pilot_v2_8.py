#!/usr/bin/env python3
"""Gated real-solver pilot for the Stateful-PD Kitagawa optimizer.

This program deliberately does not launch a parameter optimization.  It first
requires one exact baseline failure endpoint and then, in a separate explicit
invocation, one exact runout/censor endpoint.  Each gate validates the solver
output hierarchy, physical classification, geometry audit, artifacts, and
optimizer-controlled reuse before writing a hash-locked validation record.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import shlex
import subprocess
import sys
import time
from typing import Any, Iterable, Mapping

PILOT_ID = "STATEFUL_PD_KITAGAWA_GATED_PILOT_V2_8"
SOLVER_MODULE = "arrhenius_fracture.sn_pd2d_stateful_v8_7_generalized_features"
SOLVER_MODEL = (
    "SN_2D_intact_FEM_stateful_local_peridynamics_v8_7_"
    "generalized_features_local_front_spacing_fixed_geometry"
)


@dataclass(frozen=True)
class Candidate:
    candidate_id: str
    values: Mapping[str, float]


@dataclass(frozen=True)
class Condition:
    condition_id: str
    feature_type: str
    depth_um: float
    root_radius_um: float
    opening_angle_deg: float
    sigma_a_MPa: float
    cycles_max: float
    expected_class: str
    target_handoff_cycles: float | None
    pd_patch_radius_um: float
    path_refine_length_um: float
    path_refine_half_height_um: float

    @property
    def half_height_um(self) -> float:
        return math.sqrt(max(self.depth_um * self.root_radius_um, 1e-30))


@dataclass
class JobResult:
    condition_id: str
    status: str
    classification: str
    valid: bool
    cycles_total: float | None
    cycles_handoff: float | None
    geometry_pass: bool
    root_spacing_um: float | None
    global_spacing_um: float | None
    coverage: float | None
    width_ratio: float | None
    orientation_coherence: float | None
    offfront_fraction: float | None
    summary_path: str
    history_paths: list[str]
    log_path: str
    reused: bool
    return_code: int
    request_sha256: str
    summary_sha256: str | None
    error: str = ""


BASELINE = Candidate("C000", {
    "shield_chi": 0.60,
    "Gshield_eV": 0.35,
    "epsp_shield_scale": 0.0050,
    "front_link_state_shift_weight": 1.00,
    "grow_stress_GPa": 1.20,
    "link_stress_GPa": 1.10,
    "kinetic_multiplier": 1.00,
})

CONDITIONS = {
    "failure": Condition(
        condition_id="K360_FAILURE",
        feature_type="ellipse",
        depth_um=360.0,
        root_radius_um=45.0,
        opening_angle_deg=60.0,
        sigma_a_MPa=735.9214951373579,
        cycles_max=1.0e8,
        expected_class="failure",
        target_handoff_cycles=3015112.9183318703,
        pd_patch_radius_um=650.0,
        path_refine_length_um=240.0,
        path_refine_half_height_um=60.0,
    ),
    "censor": Condition(
        condition_id="K360_CENSOR",
        feature_type="ellipse",
        depth_um=360.0,
        root_radius_um=45.0,
        opening_angle_deg=60.0,
        sigma_a_MPa=690.4432004940379,
        cycles_max=1.0e8,
        expected_class="censor",
        target_handoff_cycles=None,
        pd_patch_radius_um=650.0,
        path_refine_length_um=240.0,
        path_refine_half_height_um=60.0,
    ),
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def canonical_sha(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n")
    os.replace(tmp, path)


def finite(payload: Mapping[str, Any], keys: Iterable[str]) -> float | None:
    for key in keys:
        try:
            x = float(payload.get(key))
        except (TypeError, ValueError):
            continue
        if math.isfinite(x):
            return x
    return None


def sigma_tag(value: float) -> str:
    # Retained for the legacy expected path. The solver itself currently uses
    # a shorter display tag, so output discovery must not rely on this string.
    return (f"sigmaA_{value:.12g}MPa").replace(".", "p")


def discover_summary(solver_root: Path, condition: Condition) -> Path:
    """Locate the solver summary by physics metadata, not display formatting.

    The Stateful-PD driver formats the sigma directory independently of this
    orchestrator.  A valid completed job is accepted only when exactly one
    summary under the requested case has sigma_a_MPa matching the immutable
    condition to tight floating-point tolerance.  Ambiguous or mismatched
    outputs are deliberately not reused.
    """
    expected = solver_root / "shielded" / sigma_tag(condition.sigma_a_MPa) / "summary.json"
    if expected.is_file():
        return expected

    matches: list[Path] = []
    for candidate in sorted((solver_root / "shielded").glob("sigmaA_*MPa/summary.json")):
        try:
            payload = json.loads(candidate.read_text())
            actual = float(payload.get("sigma_a_MPa"))
        except (OSError, ValueError, TypeError, json.JSONDecodeError):
            continue
        tolerance = max(1.0e-9, 1.0e-10 * abs(condition.sigma_a_MPa))
        if math.isfinite(actual) and abs(actual - condition.sigma_a_MPa) <= tolerance:
            matches.append(candidate)

    if len(matches) == 1:
        return matches[0]
    # Return the legacy expected location so parse_result reports a clear
    # missing-summary error. Ambiguity is written to the run log by caller.
    return expected


def build_command(
    python_bin: str,
    solver_root: Path,
    condition: Condition,
    *,
    max_blocks: int,
    checkpoint_every: int,
    print_every: int,
) -> list[str]:
    depth_m = condition.depth_um * 1e-6
    rho_m = condition.root_radius_um * 1e-6
    half_m = condition.half_height_um * 1e-6
    values = BASELINE.values
    kinetic = values["kinetic_multiplier"]
    return [
        python_bin, "-B", "-m", SOLVER_MODULE,
        "--out", str(solver_root),
        "--resolution-profile", "custom",
        "--fatigue-model", "plastic_shielded_case64_M1",
        "--cases", "shielded",
        "--T", "300", "--R", "0.1", "--frequency-Hz", "1000",
        "--sigma-a-MPa", f"{condition.sigma_a_MPa:.16g}",
        "--cycles-max", f"{condition.cycles_max:.16g}",
        "--block-cycles", "1e8",
        "--max-blocks", str(max_blocks),
        "--feature-type", condition.feature_type,
        "--notch-depth-m", f"{depth_m:.16g}",
        "--notch-half-height-m", f"{half_m:.16g}",
        "--notch-root-radius-m", f"{rho_m:.16g}",
        "--notch-opening-angle-deg", f"{condition.opening_angle_deg:g}",
        "--path-refine-length-m", f"{condition.path_refine_length_um * 1e-6:.16g}",
        "--path-refine-half-height-m", f"{condition.path_refine_half_height_um * 1e-6:.16g}",
        "--nx", "48", "--ny", "96", "--jitter", "0.08",
        "--root-h-fine", "7.5e-6", "--pd-horizon-m", "30e-6",
        "--pd-patch-radius-m", f"{condition.pd_patch_radius_um * 1e-6:.16g}",
        "--pd-boundary-shell-m", "60e-6",
        "--root-seed-radius-m", f"{rho_m:.16g}",
        "--seed", "1", "--pd-seed", "1",
        "--handoff-min-axial-coverage", "0.65",
        "--handoff-max-axial-gap-horizons", "1.25",
        "--checkpoint-every-blocks", str(checkpoint_every),
        "--snapshot-every", "0", "--print-every", str(print_every),
        "--shield-chi", f"{values['shield_chi']:.12g}",
        "--Gshield-eV", f"{values['Gshield_eV']:.12g}",
        "--epsp-shield-scale", f"{values['epsp_shield_scale']:.12g}",
        "--front-link-state-shift-weight", f"{values['front_link_state_shift_weight']:.12g}",
        "--grow-stress-GPa", f"{values['grow_stress_GPa']:.12g}",
        "--link-stress-GPa", f"{values['link_stress_GPa']:.12g}",
        "--nu-grow-s", f"{0.02 * kinetic:.12g}",
        "--nu-link-s", f"{0.008 * kinetic:.12g}",
    ]


def classify(payload: Mapping[str, Any], runout: float) -> tuple[str, bool]:
    status = str(payload.get("status", "missing"))
    geometry_saturated = bool(payload.get("geometry_saturated", False))
    if (
        status == "physical_handoff"
        and not geometry_saturated
        and bool(payload.get("pd_handoff_pass_final", True))
    ):
        return "failure", True
    cycles = finite(payload, ("cycles_total", "cycles_final", "cycles_end"))
    if (
        status == "right_censored"
        and cycles is not None
        and cycles >= 0.999999 * runout
        and not geometry_saturated
    ):
        return "censor", True
    return "invalid", False


def history_files(case_root: Path, payload: Mapping[str, Any]) -> list[Path]:
    candidates: list[Path] = []
    for key in ("history_csv", "pd_history_csv", "sn_history_csv"):
        value = payload.get(key)
        if value:
            p = Path(str(value))
            if not p.is_absolute():
                p = case_root / p
            candidates.append(p)
    candidates.extend(sorted(case_root.glob("*history*.csv")))
    candidates.extend(sorted(case_root.rglob("*history*.csv")))
    out: list[Path] = []
    seen: set[str] = set()
    for p in candidates:
        key = str(p.resolve()) if p.exists() else str(p)
        if p.is_file() and key not in seen:
            seen.add(key)
            out.append(p)
    return out


def parse_result(
    condition: Condition,
    summary: Path,
    log: Path,
    request_sha: str,
    *,
    reused: bool,
    return_code: int,
    error: str = "",
) -> JobResult:
    if not summary.is_file():
        return JobResult(
            condition.condition_id, "missing_summary", "invalid", False,
            None, None, False, None, None, None, None, None, None,
            str(summary), [], str(log), reused, return_code, request_sha, None,
            error or "summary missing",
        )
    try:
        payload = json.loads(summary.read_text())
    except Exception as exc:
        return JobResult(
            condition.condition_id, "invalid_summary", "invalid", False,
            None, None, False, None, None, None, None, None, None,
            str(summary), [], str(log), reused, return_code, request_sha,
            sha256_file(summary), repr(exc),
        )
    classification, physical_valid = classify(payload, condition.cycles_max)
    audit = payload.get("pd_production_preflight") or {}
    geom = payload.get("geometry_resolution_audit_final") or {}
    geometry_pass = bool(audit.get("pass", False)) and bool(geom.get("pass", True))
    case_root = summary.parent
    histories = history_files(case_root, payload)
    return JobResult(
        condition_id=condition.condition_id,
        status=str(payload.get("status", "missing")),
        classification=classification,
        valid=bool(physical_valid and geometry_pass),
        cycles_total=finite(payload, ("cycles_total", "cycles_final", "cycles_end")),
        cycles_handoff=finite(payload, ("cycles_connected", "cycles_handoff")),
        geometry_pass=geometry_pass,
        root_spacing_um=(lambda x: None if x is None else 1e6 * x)(
            finite(audit, ("point_spacing_root_m", "point_spacing_m"))
        ),
        global_spacing_um=(lambda x: None if x is None else 1e6 * x)(
            finite(audit, ("point_spacing_global_m", "point_spacing_m"))
        ),
        coverage=finite(payload, ("pd_crack_axial_coverage_final",)),
        width_ratio=finite(payload, ("pd_crack_width_ratio_final",)),
        orientation_coherence=finite(payload, ("pd_crack_orientation_coherence_final",)),
        offfront_fraction=finite(payload, (
            "pd_handoff_offfront_broken_fraction_final",
            "pd_off_front_broken_fraction_final",
        )),
        summary_path=str(summary),
        history_paths=[str(p) for p in histories],
        log_path=str(log),
        reused=reused,
        return_code=return_code,
        request_sha256=request_sha,
        summary_sha256=sha256_file(summary),
        error=error,
    )


def run_condition(
    outroot: Path,
    standalone_root: Path,
    condition: Condition,
    args: argparse.Namespace,
) -> JobResult:
    job_root = outroot / condition.condition_id
    solver_root = job_root / "solver_output"
    summary = discover_summary(solver_root, condition)
    log = job_root / "run.log"
    request_path = job_root / "REQUEST_V2_8.json"
    exit_path = job_root / "EXIT_V2_8.json"
    job_root.mkdir(parents=True, exist_ok=True)
    request: dict[str, Any] = {
        "pilot_id": PILOT_ID,
        "solver_model": SOLVER_MODEL,
        "candidate": {"candidate_id": BASELINE.candidate_id, "values": dict(BASELINE.values)},
        "condition": asdict(condition),
        "max_blocks": args.max_blocks,
        "checkpoint_every_blocks": args.checkpoint_every_blocks,
        "print_every": args.print_every,
    }
    request["request_sha256"] = canonical_sha(request)
    request_sha = str(request["request_sha256"])
    if request_path.is_file():
        prior = json.loads(request_path.read_text())
        if prior.get("request_sha256") != request_sha:
            return JobResult(
                condition.condition_id, "request_mismatch", "invalid", False,
                None, None, False, None, None, None, None, None, None,
                str(summary), [], str(log), False, -901, request_sha, None,
                "existing output root has a different immutable request",
            )
    else:
        write_json(request_path, request)

    if summary.is_file() and exit_path.is_file():
        try:
            exit_payload = json.loads(exit_path.read_text())
            if (
                exit_payload.get("request_sha256") == request_sha
                and int(exit_payload.get("return_code", 1)) == 0
            ):
                return parse_result(
                    condition, summary, log, request_sha,
                    reused=True, return_code=0,
                )
        except Exception:
            pass

    command = build_command(
        sys.executable,
        solver_root,
        condition,
        max_blocks=args.max_blocks,
        checkpoint_every=args.checkpoint_every_blocks,
        print_every=args.print_every,
    )
    (job_root / "command.txt").write_text(shlex.join(command) + "\n")
    env = os.environ.copy()
    env["PYTHONPATH"] = str(standalone_root) + (
        os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else ""
    )
    with log.open("a") as f:
        f.write("\n=== command ===\n" + shlex.join(command) + "\n")
        f.flush()
        cp = subprocess.run(
            command,
            cwd=Path(__file__).resolve().parent,
            env=env,
            stdout=f,
            stderr=subprocess.STDOUT,
        )
    write_json(exit_path, {
        "request_sha256": request_sha,
        "return_code": int(cp.returncode),
        "completed_epoch_s": time.time(),
    })
    # The solver controls its human-readable sigma directory tag. Resolve the
    # completed output again after execution rather than assuming formatting.
    summary = discover_summary(solver_root, condition)
    return parse_result(
        condition, summary, log, request_sha,
        reused=False, return_code=int(cp.returncode),
        error="" if cp.returncode == 0 else f"solver return code {cp.returncode}",
    )


def nonempty_history(paths: list[str]) -> bool:
    for value in paths:
        p = Path(value)
        try:
            if p.is_file() and p.stat().st_size > 0 and len(p.read_text(errors="replace").splitlines()) >= 2:
                return True
        except OSError:
            continue
    return False


def validate_gate(result: JobResult, condition: Condition) -> tuple[bool, list[str], dict[str, Any]]:
    errors: list[str] = []
    diagnostics: dict[str, Any] = {}
    if result.return_code != 0:
        errors.append(f"solver return code {result.return_code}")
    if result.status in {"missing_summary", "invalid_summary", "request_mismatch"}:
        errors.append(result.status)
    if result.classification != condition.expected_class:
        errors.append(
            f"expected {condition.expected_class}, obtained {result.classification} ({result.status})"
        )
    if not result.geometry_pass:
        errors.append("geometry audit failed")
    if result.root_spacing_um is None or not (6.0 <= result.root_spacing_um <= 9.0):
        errors.append(f"unexpected local root spacing: {result.root_spacing_um}")
    if not nonempty_history(result.history_paths):
        errors.append("nonempty physical history CSV not found")

    if condition.expected_class == "failure":
        if result.cycles_handoff is None:
            errors.append("failure has no cycles_handoff")
        elif condition.target_handoff_cycles:
            ratio = result.cycles_handoff / condition.target_handoff_cycles
            diagnostics["handoff_cycle_ratio_to_frozen_baseline"] = ratio
            if not (0.1 <= ratio <= 10.0):
                errors.append(f"handoff timing differs from frozen baseline by >10x: ratio={ratio:.6g}")
        if result.coverage is not None and result.coverage < 0.65:
            errors.append(f"axial coverage below handoff gate: {result.coverage}")
        if result.width_ratio is not None and result.width_ratio > 0.30:
            errors.append(f"width ratio above handoff gate: {result.width_ratio}")
        if result.orientation_coherence is not None and result.orientation_coherence < 0.65:
            errors.append(f"orientation coherence below handoff gate: {result.orientation_coherence}")
        if result.offfront_fraction is not None and result.offfront_fraction > 0.25:
            errors.append(f"off-front fraction above handoff gate: {result.offfront_fraction}")
    else:
        if result.cycles_total is None or result.cycles_total < 0.999999 * condition.cycles_max:
            errors.append(f"censor did not reach runout: cycles_total={result.cycles_total}")

    diagnostics["result"] = asdict(result)
    return not errors, errors, diagnostics


def source_identity(standalone: Path) -> dict[str, str]:
    paths = [
        standalone / "arrhenius_fracture" / "sn_pd2d_stateful_v8_7_generalized_features.py",
        standalone / "arrhenius_fracture" / "sn_feature_geometry_v8_7.py",
    ]
    return {str(p.relative_to(standalone)): sha256_file(p) for p in paths if p.is_file()}


def gate_path(outroot: Path, gate: str) -> Path:
    return outroot / f"{gate.upper()}_GATE_V2_8.json"


def load_passed_gate(outroot: Path, gate: str) -> dict[str, Any]:
    path = gate_path(outroot, gate)
    if not path.is_file():
        raise SystemExit(f"required prior gate is missing: {path}")
    payload = json.loads(path.read_text())
    if not payload.get("pass", False):
        raise SystemExit(f"required prior gate did not pass: {path}")
    return payload


def run_gate(gate: str, outroot: Path, standalone: Path, args: argparse.Namespace) -> int:
    condition = CONDITIONS[gate]
    print(f"V2.8 {gate.upper()} GATE: {condition.condition_id}", flush=True)
    first = run_condition(outroot, standalone, condition, args)
    passed, errors, diagnostics = validate_gate(first, condition)

    # A second optimizer call must reuse the completed result without touching it.
    reused = None
    reuse_errors: list[str] = []
    if passed:
        before_hash = first.summary_sha256
        second = run_condition(outroot, standalone, condition, args)
        reused = asdict(second)
        if not second.reused:
            reuse_errors.append("second call did not use optimizer-controlled reuse")
        if second.summary_sha256 != before_hash:
            reuse_errors.append("summary hash changed during reuse check")
        if second.request_sha256 != first.request_sha256:
            reuse_errors.append("request hash changed during reuse check")
        passed = passed and not reuse_errors

    payload = {
        "pilot_id": PILOT_ID,
        "gate": gate,
        "condition": asdict(condition),
        "candidate": {"candidate_id": BASELINE.candidate_id, "values": dict(BASELINE.values)},
        "pass": passed,
        "errors": errors,
        "reuse_errors": reuse_errors,
        "diagnostics": diagnostics,
        "reuse_result": reused,
        "source_identity": source_identity(standalone),
        "completed_epoch_s": time.time(),
    }
    payload["record_sha256"] = canonical_sha(payload)
    write_json(gate_path(outroot, gate), payload)
    print(json.dumps({
        "gate": gate,
        "pass": passed,
        "classification": first.classification,
        "status": first.status,
        "cycles_total": first.cycles_total,
        "cycles_handoff": first.cycles_handoff,
        "reused_on_second_call": bool(reused and reused.get("reused")),
        "errors": errors + reuse_errors,
        "record": str(gate_path(outroot, gate)),
    }, indent=2), flush=True)
    return 0 if passed else 2


def write_validation_manifest(outroot: Path, standalone: Path) -> None:
    failure = load_passed_gate(outroot, "failure")
    censor = load_passed_gate(outroot, "censor")
    payload = {
        "pilot_id": PILOT_ID,
        "validation_complete": True,
        "full_optimization_authorized": False,
        "three_candidate_sensitivity_pilot_authorized": True,
        "authorization_scope": "orchestrator and baseline endpoint validation only; the 12-candidate optimization remains prohibited until a separate three-candidate sensitivity pilot passes",
        "failure_gate_record_sha256": failure["record_sha256"],
        "censor_gate_record_sha256": censor["record_sha256"],
        "source_identity": source_identity(standalone),
        "frozen_v8_3_campaign_modified": False,
        "completed_epoch_s": time.time(),
    }
    payload["manifest_sha256"] = canonical_sha(payload)
    write_json(outroot / "PILOT_VALIDATION_MANIFEST_V2_8.json", payload)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--standalone-root", required=True)
    p.add_argument("--outroot", required=True)
    p.add_argument("--gate", choices=("failure", "censor", "status"), default="failure")
    p.add_argument("--max-blocks", type=int, default=20000)
    p.add_argument("--checkpoint-every-blocks", type=int, default=10)
    p.add_argument("--print-every", type=int, default=100)
    args = p.parse_args()

    standalone = Path(args.standalone_root).resolve()
    outroot = Path(args.outroot).resolve()
    if not (standalone / "arrhenius_fracture").is_dir():
        raise SystemExit(f"invalid standalone root: {standalone}")
    outroot.mkdir(parents=True, exist_ok=True)

    if args.gate == "status":
        payload = {
            "failure_gate": gate_path(outroot, "failure").is_file(),
            "censor_gate": gate_path(outroot, "censor").is_file(),
            "validation_manifest": (outroot / "PILOT_VALIDATION_MANIFEST_V2_8.json").is_file(),
            "outroot": str(outroot),
        }
        print(json.dumps(payload, indent=2))
        return 0

    if args.gate == "censor":
        load_passed_gate(outroot, "failure")

    rc = run_gate(args.gate, outroot, standalone, args)
    if rc == 0 and args.gate == "censor":
        write_validation_manifest(outroot, standalone)
        print(f"VALIDATION MANIFEST: {outroot / 'PILOT_VALIDATION_MANIFEST_V2_8.json'}")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
