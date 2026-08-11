#!/usr/bin/env python3
"""Cross-apply saved FEM/PD tensors and signed-MPZ capsules side-effect-free."""
from __future__ import annotations
import argparse, csv, hashlib, json
from pathlib import Path
import numpy as np
from scipy.special import logsumexp

from scripts.compare_v9_fem_pd_local_tip import FEM_ROOT, MECHANICS, PD_CASES, PEAK, SOURCE, fem_path, generation, loadjson, pd_record
from arrhenius_fracture.v9_canonical_four_class_elastic_fem import CanonicalElasticFourClassFEMCondition
from arrhenius_fracture.v9_canonical_four_class_birth import CanonicalFourClassBirthState
from arrhenius_fracture.v9_pd_shared_root_marked_cleavage import M1_PRODUCTION_MODEL_ID, M3_PARITY_MODEL_ID, SharedRootMarkedCleavageState

FEM_STATES=[(631.,1e8),(631.,1e14),(650.,1e8),(650.,1e14),(735.9214951373579,1e8),(735.9214951373579,1e14)]
PD_STATES=[(1500.,"generation_9a6c40a9d55b4149a6b79688650d7541"),(2000.,"generation_633f565da1ff4eeeaecd71df72de5f94")]

def fem_state(stress,N):
    c=CanonicalElasticFourClassFEMCondition.from_run_args(MECHANICS,PEAK,SOURCE,stress,1720)
    gd,g=generation(fem_path(stress,N)); meta=loadjson(gd/"state_metadata.json")
    with np.load(gd/"state_arrays.npz") as z: a={k:z[k].copy() for k in z.files}
    cap=c.birth.mpz.capsule();cap["arrays"].update({k[4:]:v for k,v in a.items() if k.startswith("mpz_")});cap["scalars"].update(meta["condition_capsule"]["birth"]["mpz"]["scalars"]);c.birth.mpz.restore_capsule(cap)
    return dict(label=f"K360_{stress:g}_N{N:.0e}",mpz=c.birth.mpz.copy(),tensors=a["kernel_root_tensors_Pa"],frequency=c.args.frequency_Hz,T=c.args.T,root_radius=c.root_radius_initial_m,source_generation=g)

def pd_state(stress,g):
    r=pd_record(stress,PD_CASES[stress],g)
    return dict(label=f"PD_{stress:g}_preattempt",mpz=r["_clock"].mpz.copy(),tensors=r["_tensors"],frequency=20.,T=300.,root_radius=r["_clock"].initial_tip_radius_m,source_generation=g)

def common(s):
    return dict(option_id=PEAK,source_root=SOURCE,shear_modulus_Pa=410e9/(2*(1+0.23)),poisson=s["mpz"].poisson,burgers_m=s["mpz"].burgers_m,initial_tip_radius_m=s["root_radius"],hazard_seed=1720)

def evaluate(state,tensor_source):
    tensors=tensor_source["tensors"]; kw=common(state)
    canonical=CanonicalFourClassBirthState(**kw)
    canonical.mpz=state["mpz"].copy();canonical.hazard_threshold_action=1e100
    m3=SharedRootMarkedCleavageState(**kw,mark_seed=1721,mpz=state["mpz"].copy(),m_hits=3.,model_id=M3_PARITY_MODEL_ID,endpoint_semantics="completed_cooperative_front_increment_terminate");m3.global_threshold_action=1e100
    m1=SharedRootMarkedCleavageState(**kw,mark_seed=1721,mpz=state["mpz"].copy(),m_hits=1.,model_id=M1_PRODUCTION_MODEL_ID,endpoint_semantics="elementary_attempt_to_reversible_PD_embryo");m1.global_threshold_action=1e100
    cr=canonical._advance_phase_block_exact(1.,tensor_source["frequency"],tensor_source["T"],tensors)
    sr=m3.propose_phase_block(1.,tensor_source["frequency"],tensor_source["T"],tensors)
    mr=m1.propose_phase_block(1.,tensor_source["frequency"],tensor_source["T"],tensors)
    direct=state["mpz"].copy(); drives=[direct.resolve_root_tensor(t) for t in tensors]
    opening=float(np.mean([d["opening_stress_Pa"] for d in drives]));signed=np.mean(np.stack([d["tau_signed_Pa"] for d in drives]),axis=0);dt=1./tensor_source["frequency"]
    direct.advance(.5*dt,tensor_source["T"],opening,signed)
    ds=direct.summary()
    effective=lambda x:max(max(float(x),0.)*np.sqrt(2*np.pi*state["root_radius"])-ds["signed_active_K_shield_Pa_sqrt_m"],0.)/np.sqrt(2*np.pi*max(ds["tip_radius_m"],1e-30))
    raw_logs=np.asarray([direct.cleavage_log_rate_s(effective(d["opening_stress_Pa"]),tensor_source["T"]) for d in drives])
    raw_log=float(logsumexp(raw_logs)-np.log(len(raw_logs))-np.log(tensor_source["frequency"]))
    raw=float(np.exp(raw_log)) if raw_log>-745 else 0.0
    direct.advance(.5*dt,tensor_source["T"],opening,signed)
    ca=canonical.mpz.summary();sa=sr["state"].mpz.summary();ma=mr["state"].mpz.summary();da=direct.summary()
    array_keys=("mobile_positive","mobile_negative","retained_positive","retained_negative","accumulated_slip_positive","accumulated_slip_negative")
    scalar_keys=("emitted_total","tip_radius_m","signed_active_K_shield_Pa_sqrt_m")
    state_abs=max([float(np.max(np.abs(np.asarray(ca[k])-np.asarray(sa[k])))) for k in array_keys]+[abs(float(ca[k])-float(sa[k])) for k in scalar_keys])
    m1_state_abs=max([float(np.max(np.abs(np.asarray(da[k])-np.asarray(ma[k])))) for k in array_keys]+[abs(float(da[k])-float(ma[k])) for k in scalar_keys])
    m1_rel=abs(raw-mr["detail"]["action_increment"])/max(raw,mr["detail"]["action_increment"],1e-300);m3_rel=abs(cr["hazard_increment"]-sr["detail"]["action_increment"])/max(cr["hazard_increment"],sr["detail"]["action_increment"],1e-300)
    return dict(capsule_source=state["label"],tensor_source=tensor_source["label"],root_tensor_projection_max_abs_Pa=0.,m1_raw_action_direct=raw,m1_raw_action_shared=mr["detail"]["action_increment"],m1_action_relative_error=m1_rel,m3_action_canonical=cr["hazard_increment"],m3_action_shared=sr["detail"]["action_increment"],m3_action_absolute_error=abs(cr["hazard_increment"]-sr["detail"]["action_increment"]),m3_action_relative_error=m3_rel,signed_state_m3_absolute_error=state_abs,signed_state_m1_absolute_error=m1_state_abs,aggregate_emission=float(ma["emitted_total"]),backstress_max_Pa=float(np.max(np.abs(ma["sigma_back_by_system_Pa"]))),shielding_Pa_sqrt_m=float(ma["signed_active_K_shield_Pa_sqrt_m"]),tip_radius_m=float(ma["tip_radius_m"]),pass_parity=state_abs==0. and m1_state_abs==0. and m3_rel<5e-13 and m1_rel<5e-13)

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--out",type=Path,default=Path("runs/sn_v9_shared_root_m1_peak/fem_pd_local_tip_reference_v2"));a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    states=[fem_state(*x) for x in FEM_STATES]+[pd_state(*x) for x in PD_STATES]
    rows=[evaluate(s,t) for s in states for t in states]
    with (a.out/"FEM_PD_CONSTITUTIVE_PARITY_MATRIX.csv").open("w",newline="") as f:w=csv.DictWriter(f,list(rows[0]),lineterminator="\n");w.writeheader();w.writerows(rows)
    payload={"schema":"V9_FEM_PD_CONSTITUTIVE_PARITY_MATRIX_1","read_only":True,"mechanics_trajectory_run":False,"matrix_shape":[len(states),len(states)],"all_pass":all(r["pass_parity"] for r in rows),"state_sources":[{k:v for k,v in s.items() if k not in ("mpz","tensors")}|{"tensor_sha256":hashlib.sha256(s["tensors"].tobytes()).hexdigest()} for s in states],"max_errors":{k:max(r[k] for r in rows) for k in ("m1_action_relative_error","m3_action_relative_error","m3_action_absolute_error","signed_state_m3_absolute_error","signed_state_m1_absolute_error")},"rows":rows}
    (a.out/"FEM_PD_CONSTITUTIVE_PARITY_MATRIX.json").write_text(json.dumps(payload,indent=2)+"\n")
if __name__=="__main__":main()
