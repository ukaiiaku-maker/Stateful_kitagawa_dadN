#!/usr/bin/env python3
"""Extract immutable shared-root m1 event-boundary capsules without rerunning PD."""
from __future__ import annotations

import argparse, csv, hashlib, json, math
from pathlib import Path
import numpy as np

CASES = {
    "PD_1500": Path("runs/sn_v9_shared_root_m1_peak/Peak_1500_VHCF/Peak/shielded/sigmaA_1500MPa"),
    "PD_2000": Path("runs/sn_v9_shared_root_m1_peak/Peak_2000/Peak/shielded/sigmaA_2000MPa"),
}

def load(p): return json.loads(Path(p).read_text(), parse_constant=lambda _x: math.nan)
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def classify(s):
    n=float(s["cycles"])
    endpoints=[("attempt_boundary",s.get("cycles_first_embryo")),
      ("stable_seed_boundary",s.get("cycles_first_stable")),
      ("front_capture_boundary",s.get("cycles_front_capture"))]
    for name,x in endpoints:
        if x is not None and math.isfinite(float(x)) and abs(n-float(x)) < 1e-7*max(1.,abs(n)):
            return name
    if s.get("cycles_first_embryo") is None: return "pre_attempt"
    if s.get("cycles_first_stable") is None: return "post_attempt_reversible_embryo"
    if s.get("cycles_front_capture") is None: return "post_stable_pre_capture"
    return "post_front_capture"

def extract(label,base):
    found=[]
    for gd in (base/"v9_generations").glob("generation_*"):
        s=load(gd/"summary.json"); m=load(gd/"state_metadata.json")
        s.update({k:v for k,v in m.get("pd_scalars",{}).items() if k.startswith("cycles_")})
        capsule=m.get("shared_root_marked_clock_capsule") or {}
        with np.load(gd/"state_arrays.npz") as z:
            realized=z["pd__embryo_sites"]
            stable=z["pd__stable_sites"]
            damage=z["pd__bond_damage"]
        mpz=capsule.get("mpz",{}); sc=mpz.get("scalars",{}); ar=mpz.get("arrays",{})
        def total(prefix): return float(sum(np.sum(v) for k,v in ar.items() if k.startswith(prefix)))
        role=classify(s)
        found.append(dict(case=label,boundary_role=role,cycles=float(s["cycles"]),generation=gd.name,
          hazard_comparison_eligible=role in ("pre_attempt","attempt_boundary"),
          state_use=("pre-attempt hazard capsule" if role=="pre_attempt" else "exact localized hazard capsule" if role=="attempt_boundary" else "post-seed MPZ/topology diagnostic only"),
          global_action=float(capsule.get("global_cumulative_action",math.nan)),
          global_threshold=float(capsule.get("global_threshold_action",math.nan)),
          remaining_action=float(capsule.get("global_threshold_action",math.nan))-float(capsule.get("global_cumulative_action",math.nan)),
          attempt_count=int(capsule.get("attempt_count",0)),selected_site=int(capsule.get("last_selected_site_id",-1) or -1),
          embryos=float(np.sum(realized)),stable_seeds=float(np.sum(stable)),max_bond_damage=float(np.max(damage)),
          emitted_total=float(sc.get("emitted_total",math.nan)),tip_radius_m=float(sc.get("tip_radius_m",math.nan)),
          shielding_Pa_sqrt_m=float(sc.get("signed_active_K_shield_Pa_sqrt_m",math.nan)),
          mobile_positive=total("mobile_positive"),mobile_negative=total("mobile_negative"),
          retained_positive=total("retained_positive"),retained_negative=total("retained_negative"),
          manifest_sha256=sha(gd/"manifest.json"),arrays_sha256=sha(gd/"state_arrays.npz"),metadata_sha256=sha(gd/"state_metadata.json")))
    # Retain the last accepted state strictly before the first attempt and the exact
    # first occurrence of each event boundary; terminal is the latest generation.
    first_attempt=min((r for r in found if r["attempt_count"]>=1),key=lambda r:r["cycles"])
    first_attempt["boundary_role"]="attempt_boundary"
    first_attempt["hazard_comparison_eligible"]=True
    first_attempt["state_use"]="exact localized hazard capsule"
    pre=max((r for r in found if r["cycles"]<first_attempt["cycles"]),key=lambda r:r["cycles"])
    pre["boundary_role"]="pre_attempt";pre["hazard_comparison_eligible"]=True;pre["state_use"]="pre-attempt hazard capsule"
    stable=min((r for r in found if r["stable_seeds"]>=1),key=lambda r:r["cycles"])
    stable["boundary_role"]="stable_seed_boundary";stable["hazard_comparison_eligible"]=False;stable["state_use"]="post-seed MPZ/topology diagnostic only"
    terminal=max(found,key=lambda r:r["cycles"]); terminal["boundary_role"]="terminal_front_capture_state";terminal["hazard_comparison_eligible"]=False;terminal["state_use"]="post-seed MPZ/topology diagnostic only"
    initial={k:math.nan for k in pre}; initial.update(case=label,boundary_role="initial_dormant_state",cycles=0.,generation="deterministic_initial_state_from_hash_verified_run_args",hazard_comparison_eligible=True,state_use="initial hazard reference; no accepted cycle",attempt_count=0,embryos=0.,stable_seeds=0.,max_bond_damage=0.,manifest_sha256=sha(base/"run_args.json"))
    return [initial,pre,first_attempt,stable,terminal]

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--out",type=Path,default=Path("runs/sn_v9_shared_root_m1_peak/fem_pd_local_tip_reference_v2"));a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    rows=sum((extract(k,v) for k,v in CASES.items()),[])
    with (a.out/"PD_EVENT_BOUNDARY_STATES.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n");w.writeheader();w.writerows(rows)
    payload={"schema":"V9_PD_M1_EVENT_BOUNDARY_STATES_1","read_only":True,"mechanics_rerun":False,"rows":rows}
    (a.out/"PD_EVENT_BOUNDARY_STATES.json").write_text(json.dumps(payload,indent=2,allow_nan=True)+"\n")
if __name__=="__main__": main()
