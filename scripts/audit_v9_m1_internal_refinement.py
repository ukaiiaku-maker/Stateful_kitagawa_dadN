#!/usr/bin/env python3
"""Audit Peak m1 internal refinement, transition persistence, restart and attempt quantiles."""
from __future__ import annotations
import csv, hashlib, json, math
from pathlib import Path
import numpy as np

ROOT=Path("runs/sn_v9_shared_root_m1_peak")
RUNS={100000.:ROOT/"Peak_2000/Peak/shielded/sigmaA_2000MPa",10000.:ROOT/"internal_refinement_1e4/Peak/shielded/sigmaA_2000MPa",2500.:ROOT/"internal_refinement_2p5e3/Peak/shielded/sigmaA_2000MPa",625.:ROOT/"internal_refinement_625/Peak/shielded/sigmaA_2000MPa",156.25:ROOT/"internal_refinement_156p25/Peak/shielded/sigmaA_2000MPa"}
RESTART=ROOT/"internal_refinement_1e4_restart_stable/Peak/shielded/sigmaA_2000MPa"
CASES={1500.:ROOT/"Peak_1500_VHCF/Peak/shielded/sigmaA_1500MPa",2000.:ROOT/"Peak_2000/Peak/shielded/sigmaA_2000MPa"}
TRANSITIONS={39.0625:ROOT/"transition_refinement_39p0625/Peak/shielded/sigmaA_2000MPa",9.765625:ROOT/"transition_refinement_9p765625/Peak/shielded/sigmaA_2000MPa",2.44140625:ROOT/"transition_refinement_2p44140625/Peak/shielded/sigmaA_2000MPa",0.6103515625:ROOT/"transition_refinement_0p6103515625/Peak/shielded/sigmaA_2000MPa"}
OUT=ROOT/"numerical_qualification_v2"
def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    OUT.mkdir(parents=True,exist_ok=True);rows=[]
    for ceiling,p in RUNS.items():
        if not (p/"summary.json").exists():continue
        x=load(p/"summary.json")
        rows.append(dict(internal_max_cycles=ceiling,N_attempt=x["cycles_first_embryo"],N_stable=x["cycles_first_stable"],N_softening=x["cycles_first_softening"],N_root_connected=x["cycles_root_connected"],N_front_capture=x["cycles_front_capture"],H_attempt=x["H_attempt_final"],attempt_site=(x.get("last_marked_attempt") or {}).get("site_id"),summary_sha256=sha(p/"summary.json"),checkpoint_sha256=sha(p/"checkpoint_latest.npz")))
    with (OUT/"internal_macro_refinement.csv").open("w",newline="") as f:w=csv.DictWriter(f,list(rows[0]),lineterminator="\n");w.writeheader();w.writerows(rows)
    transition=[]
    for ceiling,p in TRANSITIONS.items():
      if not (p/"summary.json").exists():continue
      x=load(p/"summary.json");transition.append(dict(internal_max_cycles=ceiling,N_attempt=x["cycles_first_embryo"],N_stable=x["cycles_first_stable"],transition_delay_cycles=x["cycles_first_stable"]-x["cycles_first_embryo"],H_attempt=x["H_attempt_final"],endpoint="diagnostic_stable_seed_stop",summary_sha256=sha(p/"summary.json")))
    transition.sort(key=lambda r:r["internal_max_cycles"],reverse=True)
    with (OUT/"transition_clock_refinement.csv").open("w",newline="") as f:w=csv.DictWriter(f,list(transition[0]),lineterminator="\n");w.writeheader();w.writerows(transition)
    finest=transition[-2:]
    transition_qualification={"relative_attempt_difference_finest_pair":abs(finest[-1]["N_attempt"]-finest[-2]["N_attempt"])/finest[-1]["N_attempt"],"absolute_attempt_difference_cycles_finest_pair":abs(finest[-1]["N_attempt"]-finest[-2]["N_attempt"]),"transition_delay_spread_cycles":max(r["transition_delay_cycles"] for r in transition)-min(r["transition_delay_cycles"] for r in transition),"relative_attempt_tolerance":5e-4,"persistent_competing_clock_qualified":abs(finest[-1]["N_attempt"]-finest[-2]["N_attempt"])/finest[-1]["N_attempt"]<5e-4 and max(r["transition_delay_cycles"] for r in transition)-min(r["transition_delay_cycles"] for r in transition)<1e-9,"production_internal_transition_ceiling_cycles":finest[-1]["internal_max_cycles"],"note":"attempt coordinate retains the stated sub-cycle discretization uncertainty; stable-minus-attempt delay is partition invariant"}
    rate_path=ROOT/"rate_separated_restart_from_fine_stable_v4/Peak/shielded/sigmaA_2000MPa";rate=load(rate_path/"summary.json");fine_full=load(RUNS[156.25]/"summary.json")
    rate_separated={"source_generation":"generation_fb1f3445c7344342b6c43edc1b8c1da2","pretransition_ceiling_cycles":.6103515625,"transition_ceiling_cycles":.6103515625,"poststable_ceiling_cycles":625.,"N_attempt":rate["cycles_first_embryo"],"N_stable":rate["cycles_first_stable"],"N_softening":rate["cycles_first_softening"],"N_root_connected":rate["cycles_root_connected"],"N_front_capture":rate["cycles_front_capture"],"front_capture_relative_difference_vs_full_156p25":abs(rate["cycles_front_capture"]-fine_full["cycles_front_capture"])/fine_full["cycles_front_capture"],"action_identity":rate["H_attempt_final"]==fine_full["H_attempt_final"],"summary_sha256":sha(rate_path/"summary.json"),"checkpoint_sha256":sha(rate_path/"checkpoint_latest.npz")}
    (OUT/"rate_separated_production_qualification.json").write_text(json.dumps(rate_separated,indent=2)+"\n")
    a=RUNS[10000.];b=RESTART;sa=load(a/"summary.json");sb=load(b/"summary.json")
    endpoint_keys=("cycles_first_embryo","cycles_first_stable","cycles_first_softening","cycles_root_connected","cycles_front_capture","H_attempt_final")
    with np.load(a/"checkpoint_latest.npz",allow_pickle=True) as za,np.load(b/"checkpoint_latest.npz",allow_pickle=True) as zb:
      diffs={k:float(np.nanmax(np.abs(za[k]-zb[k]))) for k in set(za.files)&set(zb.files)-{"metadata_json"} if np.issubdtype(za[k].dtype,np.number) and not np.all(np.isnan(za[k]))}
    restart={"schema":"V9_M1_PERSISTED_RESTART_AFTER_TRANSITION_1","source_generation":"generation_9f2e88cb61d74708bec307a2dc5ea432","source_boundary_cycles":568.7777746616391,"endpoint_deltas":{k:float(sb[k]-sa[k]) for k in endpoint_keys},"endpoint_identity":all(sb[k]==sa[k] for k in endpoint_keys),"selected_site_identity":sb["last_marked_attempt"]["site_id"]==sa["last_marked_attempt"]["site_id"],"global_action_identity":sb["H_attempt_final"]==sa["H_attempt_final"],"global_threshold_identity":sb["global_threshold_action_final"]==sa["global_threshold_action_final"],"maximum_array_absolute_difference":max(diffs.values()),"array_differences_nonzero":{k:v for k,v in diffs.items() if v},"qualification":"endpoint/action/threshold identity; floating reconstruction tolerance for fields"}
    (OUT/"persisted_restart_after_transition.json").write_text(json.dumps(restart,indent=2)+"\n")
    quant=[]; targets={"N10":-math.log(.9),"N50":-math.log(.5),"N90":-math.log(.1)}
    for stress,p in CASES.items():
      x=load(p/"summary.json");H=float(x["H_attempt_final"]);N=float(x["cycles_first_embryo"]);rate=H/N
      for name,h in targets.items():quant.append(dict(sigma_a_MPa=stress,quantile=name,target_action=h,N_cycles_stationary_action_projection=h/rate,basis_H=H,basis_N=N,action_per_cycle_projection=rate,quantity="global_attempt_first_passage",not_stable_front_survival=True,representation="analytical_exponential_threshold_inversion_using_trajectory_supported_mean_action"))
    with (OUT/"analytical_attempt_quantiles.csv").open("w",newline="") as f:w=csv.DictWriter(f,list(quant[0]),lineterminator="\n");w.writeheader();w.writerows(quant)
    manifest={"schema":"V9_M1_NUMERICAL_QUALIFICATION_V2","new_stress_conditions":False,"pd_image_policy":"none","refinement_rows":len(rows),"transition_refinement_rows":len(transition),"transition_qualification":transition_qualification,"rate_separated_production":rate_separated,"restart":restart,"attempt_quantile_semantics":"S_no_attempt=exp(-H_attempt); does not include stabilization, healing, or front capture","files":{p.name:sha(p) for p in OUT.iterdir() if p.is_file() and p.name!="manifest.json"}}
    (OUT/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
if __name__=="__main__":main()
