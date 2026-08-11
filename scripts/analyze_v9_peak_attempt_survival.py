#!/usr/bin/env python3
"""Invert a recorded conditional-survival action curve without H/N scaling."""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess

ROOT=Path("runs/sn_v9_shared_root_m1_peak")
CASES={
    1172.0:ROOT/"attempt_survival_1172_v1/Peak/shielded/sigmaA_1172MPa",
    1290.0:ROOT/"attempt_survival_1290_v1/Peak/shielded/sigmaA_1290MPa",
    1460.0:ROOT/"attempt_survival_1460_v1/Peak/shielded/sigmaA_1460MPa",
    1500.0:ROOT/"attempt_survival_1500_1e12_v1/Peak/shielded/sigmaA_1500MPa",
    1630.0:ROOT/"attempt_survival_1630_v1/Peak/shielded/sigmaA_1630MPa",
    1820.0:ROOT/"attempt_survival_1820_v1/Peak/shielded/sigmaA_1820MPa",
    2000.0:ROOT/"attempt_survival_2000_v2/Peak/shielded/sigmaA_2000MPa",
}
OUT=ROOT/"attempt_sn_v1"
EXACT_BOUNDARIES={(1820.0,"N50"):ROOT/"attempt_boundary_1820_N50_v3/Peak/shielded/sigmaA_1820MPa/summary.json"}


def checkpoints(sigma_a_MPa,path):
    modes=json.loads(path.read_text())
    rows=[]
    for mode in modes:
        if "H_attempt_after_segment" not in mode:
            continue
        rows.append({
            "sigma_a_MPa":sigma_a_MPa,"cycles":float(mode["cycles_total"]),
            "H_attempt":float(mode["H_attempt_after_segment"]),
            "S_no_attempt":math.exp(-float(mode["H_attempt_after_segment"])),
            "numerical_mode":mode["mode"],
            "state_error":mode["detail"].get("state_error",0.0),
            "hazard_error":mode["detail"].get("log_hazard_error",0.0),
            "action_source":("stationary_propagation_after_certificate" if mode["mode"]=="stationary"
                             else "direct_integral_from_partition_qualified_private_window"),
        })
    if any(b["cycles"]<=a["cycles"] or b["H_attempt"]<=a["H_attempt"] for a,b in zip(rows,rows[1:])):
        raise RuntimeError("conditional action checkpoints are not strictly monotone")
    return rows


def invert_piecewise(rows,target):
    for left,right in zip(rows,rows[1:]):
        if left["H_attempt"]<=target<=right["H_attempt"]:
            # Interpolate the recorded cumulative action itself.  This is not
            # target/(H/N), and it is explicitly retained as finite-horizon
            # interpolation pending an exact endpoint trajectory.
            fraction=(target-left["H_attempt"])/(right["H_attempt"]-left["H_attempt"])
            return left["cycles"]+fraction*(right["cycles"]-left["cycles"]),left,right
    raise ValueError(f"target action {target} is outside recorded action support")


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""):h.update(chunk)
    return h.hexdigest()


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    by_stress={stress:checkpoints(stress,case/"v9_pd_high_cycle_mode_history.json")
               for stress,case in CASES.items()}
    rows=[row for stress_rows in by_stress.values() for row in stress_rows]
    with (OUT/"Peak_blunt_m1_attempt_action_checkpoints.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n");w.writeheader();w.writerows(rows)
    q=[]
    for stress,stress_rows in by_stress.items():
      for name,target in (("N10",-math.log(.9)),("N50",math.log(2.)),("N90",-math.log(.1))):
        n,left,right=invert_piecewise(stress_rows,target)
        exact_path=EXACT_BOUNDARIES.get((stress,name));exact=None
        if exact_path is not None:
            exact=json.loads(exact_path.read_text())
        interpolated=n
        reported=float(exact["cycles_total"]) if exact else interpolated
        q.append({"sigma_a_MPa":stress,"quantile":name,"target_action":target,
                  "cycles":reported,"interpolated_cycles":interpolated,
                  "absolute_interpolation_error_cycles":abs(reported-interpolated) if exact else "",
                  "relative_interpolation_error":abs(reported-interpolated)/reported if exact else "",
                  "lower_checkpoint_cycles":left["cycles"],
                  "upper_checkpoint_cycles":right["cycles"],
                  "method":"piecewise_linear_inversion_of_actual_cumulative_action",
                  "trajectory_supported":True,"exact_endpoint_run":bool(exact)})
    with (OUT/"Peak_blunt_m1_attempt_SN.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(q[0]),lineterminator="\n");w.writeheader();w.writerows(q)
    provenance={}
    for stress,case in CASES.items():
        summary=json.loads((case/"summary.json").read_text())
        args=json.loads((case/"run_args.json").read_text())
        controller=json.loads((case/"v9_pd_high_cycle_controller.json").read_text())
        modes=json.loads((case/"v9_pd_high_cycle_mode_history.json").read_text())
        totals={name:sum(float(row["accepted_cycles"]) for row in modes
                         if row["accepted"] and row["mode"]==name)
                for name in ("exact_private_window","projective","stationary","exact_burst")}
        active=json.loads((case/"v9_generations/ACTIVE.json").read_text())
        provenance[str(stress)]={
            "run_path":str(case),"model_id":summary["cleavage_clock_model"],
            "material_row_sha256":summary["four_class_registry_audit"]["active_parameter_row_sha256"],
            "material_fingerprint_sha256":summary["four_class_registry_audit"]["active_parameter_fingerprint_sha256"],
            "source_hashes":summary["source_sha256"],
            "temperature_K":summary["T_K"],"R":summary["R"],"frequency_Hz":summary["frequency_Hz"],
            "nominal_root_radius_m":summary["nominal_notch_root_radius_m"],
            "mesh_resolved_initial_root_radius_m":summary["root_radius_initial_m"],
            "conditional_survival_protocol":summary.get("conditional_survival_protocol",{
                "protocol":"conditional_no_event_action","analysis_only":True,
                "physical_threshold_draw":False,"mark_rng_consumed":False,
                "topology_continuation_permitted":False,
                "threshold_override_action":args["shared_root_survival_threshold_action"],
                "reconstructed_from_run_args":True,
            }),
            "checkpoint_generation":active["generation"],
            "checkpoint_sha256":sha256(case/"checkpoint_latest.npz"),
            "controller_sha256":sha256(case/"v9_pd_high_cycle_controller.json"),
            "mode_history_sha256":sha256(case/"v9_pd_high_cycle_mode_history.json"),
            "accepted_cycle_totals_by_mode":totals,
            "campaign_exact_map_evaluations":controller["campaign_exact_map_evaluations"],
            "mechanics_validity":"valid" if summary["geometry_resolution_audit_final"]["pass"] else "invalid",
            "geometry_resolution_audit":summary["geometry_resolution_audit_final"],
            "action_checkpoint_support_cycles":[by_stress[stress][0]["cycles"],by_stress[stress][-1]["cycles"]],
        }
    base_commit=subprocess.check_output(["git","rev-parse","HEAD"],text=True).strip()
    (OUT/"manifest.json").write_text(json.dumps({
        "schema":"V9_PEAK_BLUNT_M1_ATTEMPT_SN_1","material":"Peak","temperature_K":300.0,
        "repository_base_commit_at_analysis":base_commit,
        "geometry":"ideal smooth 600 um blunt notch","correlated_front_segment_um":10.0,
        "source_mode_histories":{str(stress):str(case/"v9_pd_high_cycle_mode_history.json")
                                 for stress,case in CASES.items()},
        "quantile_semantics":"global attempt survival only; not stable-front survival",
        "stationarity_assumed":False,
        "conditions":provenance,
    },indent=2)+"\n")


if __name__=="__main__":main()
