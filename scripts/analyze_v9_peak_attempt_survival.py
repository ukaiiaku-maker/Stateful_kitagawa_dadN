#!/usr/bin/env python3
"""Invert a recorded conditional-survival action curve without H/N scaling."""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

ROOT=Path("runs/sn_v9_shared_root_m1_peak")
CASES={
    1500.0:ROOT/"attempt_survival_1500_1e12_v1/Peak/shielded/sigmaA_1500MPa",
    1820.0:ROOT/"attempt_survival_1820_v1/Peak/shielded/sigmaA_1820MPa",
}
OUT=ROOT/"attempt_sn_v1"


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
        q.append({"sigma_a_MPa":stress,"quantile":name,"target_action":target,
                  "cycles":n,"lower_checkpoint_cycles":left["cycles"],
                  "upper_checkpoint_cycles":right["cycles"],
                  "method":"piecewise_linear_inversion_of_actual_cumulative_action",
                  "trajectory_supported":True,"exact_endpoint_run":False})
    with (OUT/"Peak_blunt_m1_attempt_SN.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(q[0]),lineterminator="\n");w.writeheader();w.writerows(q)
    (OUT/"manifest.json").write_text(json.dumps({
        "schema":"V9_PEAK_BLUNT_M1_ATTEMPT_SN_1","material":"Peak","temperature_K":300.0,
        "geometry":"ideal smooth 600 um blunt notch","correlated_front_segment_um":10.0,
        "source_mode_histories":{str(stress):str(case/"v9_pd_high_cycle_mode_history.json")
                                 for stress,case in CASES.items()},
        "quantile_semantics":"global attempt survival only; not stable-front survival",
        "stationarity_assumed":False,
    },indent=2)+"\n")


if __name__=="__main__":main()
