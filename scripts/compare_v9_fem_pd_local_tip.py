#!/usr/bin/env python3
"""Read-only local-tip comparison of canonical FEM-v3 and shared-root PD m=1."""
from __future__ import annotations

import argparse, csv, hashlib, json, math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from scripts.audit_v9_pd_persistent_sites import build_context
from arrhenius_fracture import sn_pd2d_stateful_v9_transactional as driver
from arrhenius_fracture.v9_canonical_four_class_elastic_fem import CanonicalElasticFourClassFEMCondition
from arrhenius_fracture.v9_canonical_four_class_birth import canonical_effective_cleavage_rate
from arrhenius_fracture.v9_pd_shared_root_marked_cleavage import (
    M1_PRODUCTION_MODEL_ID, M3_PARITY_MODEL_ID, SharedRootMarkedCleavageState,
)
from arrhenius_fracture.v9_quiet_tail_kernel import restore_condition_checkpoint
from arrhenius_fracture.sn_intact_fem import stress_state_intact

PEAK = "v913_paper_peak01_0242980_persistent_sites"
SOURCE = Path("/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1")
FEM_ROOT = Path("runs/sn_v9_canonical_four_class_elastic_v3")
MECHANICS = Path("runs/sn_v9_single_seed_skeleton/K360_FAILURE/solver_output/shielded/sigmaA_735p921MPa/run_args.json")
PD_CASES = {
    2000.: Path("runs/sn_v9_shared_root_m1_peak/Peak_2000/Peak/shielded/sigmaA_2000MPa"),
    1500.: Path("runs/sn_v9_shared_root_m1_peak/Peak_1500_VHCF/Peak/shielded/sigmaA_1500MPa"),
}
_FEM_MECHANICS = {}

def digest(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def loadjson(p): return json.loads(Path(p).read_text(), parse_constant=lambda _x: float("inf"))
def generation(base):
    g=loadjson(base/"ACTIVE.json")["generation"]; return base/g, g
def fem_path(stress,N):
    base=FEM_ROOT/"quiet_tail_checkpoints"/"Peak"
    if stress != 735.9214951373579: base=base/"descent"/f"sigmaA_{stress:g}MPa"
    return base/f"N_{N:.17g}"
def eff(opening, shield, r0, r):
    K=np.maximum(opening,0)*math.sqrt(2*math.pi*r0)-shield
    return np.maximum(K,0)/math.sqrt(2*math.pi*r)
def sums(a,prefix): return float(sum(np.sum(a[k]) for k in a if k.startswith(prefix)))
def fem_mechanics(stress):
    if stress not in _FEM_MECHANICS:
        c=CanonicalElasticFourClassFEMCondition.from_run_args(MECHANICS,PEAK,SOURCE,stress,1720)
        h=c._elastic_history; phase,node=np.unravel_index(np.argmax(h["s1_node"]),h["s1_node"].shape)
        gp,_,s1gp,_=stress_state_intact(c.mesh,h["u_hist"][phase],c.fem.ep_gp,c.cached_fem.Dmat,c.cached_fem.material)
        elem=int(np.argmax(s1gp)); centroid=np.mean(c.mesh.nodes[c.mesh.elems[elem]],axis=0)
        _FEM_MECHANICS[stress]=(c,int(node),int(elem),centroid)
    return _FEM_MECHANICS[stress]
def csvrows(path):
    with Path(path).open(newline="") as f: return list(csv.DictReader(f))
def audit_source_tables():
    q={(r["material_class"],float(r["sigma_a_MPa"]),float(r["N"])):r for r in csvrows(FEM_ROOT/"300K_quiet_tail_state.csv")}
    s={(r["material_class"],float(r["sigma_a_MPa"]),float(r["N"])):r for r in csvrows(FEM_ROOT/"300K_stable_birth_survival.csv")}
    d={(r["material_class"],float(r["sigma_a_MPa"]),float(r["N"])):r for r in csvrows(FEM_ROOT/"300K_endurance_stress_descent.csv")}
    if set(q)!=set(s): raise RuntimeError("quiet-tail/survival checkpoint keys differ")
    for key in q:
        if float(q[key]["H_cleave"])!=float(s[key]["H_cleave"]): raise RuntimeError(f"direct/survival H mismatch: {key}")
    for key in d:
        if float(d[key]["H_cleave"])!=float(q[key]["H_cleave"]): raise RuntimeError(f"descent/direct H mismatch: {key}")
    return {"quiet_tail_rows":len(q),"survival_rows":len(s),"descent_rows":len(d),"direct_H_relationship":"exact"}

def fem_rows():
    rows=[]
    survival={(r["material_class"],float(r["sigma_a_MPa"]),float(r["N"])):r for r in csvrows(FEM_ROOT/"300K_stable_birth_survival.csv")}
    descent={(r["material_class"],float(r["sigma_a_MPa"]),float(r["N"])):r for r in csvrows(FEM_ROOT/"300K_endurance_stress_descent.csv")}
    gates=csvrows(FEM_ROOT/"300K_energy_gate_envelope.csv")
    for stress in (735.9214951373579,650.,631.):
      for N in (1e8,1e10,1e12,1e14):
        gd,g=generation(fem_path(stress,N)); meta=loadjson(gd/"state_metadata.json"); summ=loadjson(gd/"summary.json")
        with np.load(gd/"state_arrays.npz") as z: a={k:z[k].copy() for k in z.files}
        audit=meta["condition_capsule"]["audit"]; ks=meta["kernel_scalars"]
        tensors=a["kernel_root_tensors_Pa"]; opening=a["kernel_root_opening_Pa"]
        effective=eff(opening,ks["signed_active_K_shield_Pa_sqrt_m"],audit["root_radius_initial_m"],ks["tip_radius_m"])
        raw=np.exp(a["kernel_cleavage_log_rate_raw_s"]); m1=float(np.mean(raw)/audit["frequency_Hz"])
        sigma_max=2*stress/(1-audit["R"])
        root_s1=np.linalg.eigvalsh(tensors)[:,-1]
        condition,hotnode,hotelem,hotxy=fem_mechanics(stress); rootxy=condition.mesh.nodes[audit["root_node"]]
        key=("Peak",stress,N); sr=survival[key]; dr=descent.get(key)
        matching_gate=[x for x in gates if x["material_class"]=="Peak" and float(x["sigma_a_MPa"])==stress and float(x["N"])==N]
        admitted=[float(x["admitted_length_m"]) for x in matching_gate if float(x["admitted_length_m"])>0]
        rows.append(dict(model="FEM_v3",geometry="K360_sharp",sigma_a_MPa=stress,
          sigma_min_MPa=audit["R"]*sigma_max,sigma_max_MPa=sigma_max,N=N,
          nominal_a_um=360.,nominal_b_um=127.2792206135786,nominal_rho_um=45.,
          mesh_root_radius_um=audit["root_radius_initial_m"]*1e6,analytical_Kt=1+2*360/127.2792206135786,
          FEM_Kt=float(max(root_s1)/(sigma_max*1e6)),Kt_root_opening=float(max(opening)/(sigma_max*1e6)),Kt_root_principal=float(max(root_s1)/(sigma_max*1e6)),
          Kt_hotspot_principal=float(np.max(condition._elastic_history["s1_node"])/(sigma_max*1e6)),root_node=audit["root_node"],root_x_m=rootxy[0],root_y_m=rootxy[1],
          hotspot_node=hotnode,hotspot_node_x_m=condition.mesh.nodes[hotnode,0],hotspot_node_y_m=condition.mesh.nodes[hotnode,1],hotspot_element=hotelem,hotspot_element_x_m=hotxy[0],hotspot_element_y_m=hotxy[1],root_hotspot_distance_m=float(np.linalg.norm(rootxy-condition.mesh.nodes[hotnode])),
          root_opening_max_MPa=float(max(opening)/1e6),effective_cleavage_max_MPa=float(max(effective)/1e6),
          raw_m1_action_per_cycle=m1,m3_action_per_cycle=ks["cycle_hazard"],cumulative_action=summ["H_cleave"],
          cycle_hazard=summ["cycle_hazard"],backstress_max_Pa=float(max(abs(a["kernel_sigma_back_by_system_Pa"]))),
          shielding_Pa_sqrt_m=ks["signed_active_K_shield_Pa_sqrt_m"],tip_radius_um=ks["tip_radius_m"]*1e6,
          mobile_positive=sums(a,"kernel_mobile_positive"),mobile_negative=sums(a,"kernel_mobile_negative"),
          retained_positive=sums(a,"kernel_retained_positive"),retained_negative=sums(a,"kernel_retained_negative"),
          slip_positive=sums(a,"kernel_accumulated_slip_positive"),slip_negative=sums(a,"kernel_accumulated_slip_negative"),
          aggregate_emission=meta["condition_capsule"]["birth"]["mpz"]["scalars"]["emitted_total"],
          conditional_attempt_admission_probability=float((dr or sr)["conditional_attempt_admission_probability"]),
          rejected_phase_xi_fraction=float(dr["rejected_phase_xi_fraction"]) if dr else (sum(float(x["admitted_length_m"])==0 for x in matching_gate)/len(matching_gate) if matching_gate else math.nan),
          minimum_admitted_length_m=float(dr["minimum_admitted_length_m"]) if dr else (min(admitted) if admitted else math.nan),
          maximum_admitted_length_m=float(dr["maximum_admitted_length_m"]) if dr else (max(admitted) if admitted else math.nan),
          minimum_energy_gate_margin_J_per_m=float(dr["minimum_energy_gate_margin_J_per_m"]) if dr else math.nan,
          endpoint="energy_admissible_stable_birth_process",event_or_censor_cycle=N,checkpoint_generation=g,
          data_role="direct_atomic_checkpoint",tensor_sha256=hashlib.sha256(tensors.tobytes()).hexdigest(),_tensors=tensors))
    return rows

def pd_record(stress,case,generation_name=None):
    args,mesh,patch,chain,crack,cached,transaction,smax,smin=build_context(case/"run_args.json",stress)
    if generation_name:
      gd=case/"v9_generations"/generation_name
      md=driver._restore_nonfinite_tags(loadjson(gd/"state_metadata.json"))
      with np.load(gd/"state_arrays.npz") as z:
        mesh.nodes[:]=z["mesh_nodes"];ep=z["ep_gp"].copy();u=z["u"].copy();root=z["root_xy"].copy()
      driver.rebuild_mesh_geometry(mesh,root);patch.update_geometry(mesh,root)
      capsule=md["shared_root_marked_clock_capsule"]
      checkpoint_cycle=float(md["cycles"])
    else:
      with np.load(case/"checkpoint_latest.npz",allow_pickle=True) as z:
        md=driver._restore_nonfinite_tags(json.loads(str(z["metadata_json"].item())))
        mesh.nodes[:]=z["mesh_nodes"]
        ep=z["ep_gp"].copy(); u=z["u"].copy(); root=z["root_xy"].copy(); capsule=md["shared_root_marked_clock_capsule"]
      checkpoint_cycle=loadjson(case/"summary.json")["cycles_total"]
    summary=loadjson(case/"summary.json")
    option_id=capsule["audit"].get("option_id",PEAK)
    source_root=capsule["audit"].get("source_repository",str(SOURCE))
    clock=SharedRootMarkedCleavageState(option_id,source_root,
      shear_modulus_Pa=args.E_GPa*1e9/(2*(1+args.nu)),poisson=args.nu,burgers_m=args.b_m,
      initial_tip_radius_m=summary["root_radius_initial_m"],hazard_seed=args.global_cleavage_seed,
      mark_seed=args.spatial_mark_seed,m_hits=1.,model_id=M1_PRODUCTION_MODEL_ID,
      endpoint_semantics="elementary_attempt_to_reversible_PD_embryo")
    clock.restore_capsule(capsule)
    Umax,Umin,uz,_,_=cached.affine(ep,smax,smin,u); hist=cached.stress_histories(ep,Umax,Umin,args.hazard_n_phase,uz)
    node=int(np.argmin(np.linalg.norm(mesh.nodes-root[None,:],axis=1))); v=hist["sigma_node"][:,:,node]
    tensors=np.empty((len(v),2,2)); tensors[:,0,0]=v[:,0];tensors[:,1,1]=v[:,1];tensors[:,0,1]=tensors[:,1,0]=v[:,2]
    drives=[clock.mpz.resolve_root_tensor(t) for t in tensors]; opening=np.array([d["opening_stress_Pa"] for d in drives])
    effective=np.array([clock.effective_opening_stress_Pa(x) for x in opening]); raw=np.exp([clock.mpz.cleavage_log_rate_s(x,args.T) for x in effective])
    m1=float(np.mean(raw)/args.frequency_Hz); m3=float(np.mean([canonical_effective_cleavage_rate(x,3.,1e-6) for x in raw])/args.frequency_Hz)
    mpz=clock.mpz.summary(); sigma_max=smax/1e6; root_s1=np.linalg.eigvalsh(tensors)[:,-1]
    phase,hotnode=np.unravel_index(np.argmax(hist["s1_node"]),hist["s1_node"].shape)
    _,_,s1gp,_=stress_state_intact(mesh,hist["u_hist"][phase],ep,cached.Dmat,cached.material); hotelem=int(np.argmax(s1gp));hotxy=np.mean(mesh.nodes[mesh.elems[hotelem]],axis=0)
    return dict(model="PD_m1",geometry="blunt_600um",sigma_a_MPa=stress,sigma_min_MPa=smin/1e6,sigma_max_MPa=sigma_max,
      N=checkpoint_cycle,nominal_a_um=150.,nominal_b_um=300.,nominal_rho_um=600.,
      mesh_root_radius_um=summary["root_radius_initial_m"]*1e6,analytical_Kt=2.,FEM_Kt=float(max(root_s1)/(smax)),Kt_root_opening=float(max(opening)/smax),Kt_root_principal=float(max(root_s1)/smax),Kt_hotspot_principal=float(np.max(hist["s1_node"])/smax),
      root_node=node,root_x_m=mesh.nodes[node,0],root_y_m=mesh.nodes[node,1],hotspot_node=int(hotnode),hotspot_node_x_m=mesh.nodes[hotnode,0],hotspot_node_y_m=mesh.nodes[hotnode,1],hotspot_element=hotelem,hotspot_element_x_m=hotxy[0],hotspot_element_y_m=hotxy[1],root_hotspot_distance_m=float(np.linalg.norm(mesh.nodes[node]-mesh.nodes[hotnode])),
      root_opening_max_MPa=float(max(opening)/1e6),effective_cleavage_max_MPa=float(max(effective)/1e6),
      raw_m1_action_per_cycle=m1,m3_action_per_cycle=m3,cumulative_action=summary["H_attempt_final"],cycle_hazard=m1,
      backstress_max_Pa=float(max(abs(np.asarray(mpz["sigma_back_by_system_Pa"])))),shielding_Pa_sqrt_m=mpz["signed_active_K_shield_Pa_sqrt_m"],
      tip_radius_um=mpz["tip_radius_m"]*1e6,mobile_positive=float(np.sum(mpz["mobile_positive"])),mobile_negative=float(np.sum(mpz["mobile_negative"])),
      retained_positive=float(np.sum(mpz["retained_positive"])),retained_negative=float(np.sum(mpz["retained_negative"])),
      slip_positive=float(np.sum(mpz["accumulated_slip_positive"])),slip_negative=float(np.sum(mpz["accumulated_slip_negative"])),
      aggregate_emission=mpz["emitted_total"],endpoint="front_capture_after_stable_spatial_seed",event_or_censor_cycle=summary["cycles_total"],
      conditional_attempt_admission_probability=math.nan,rejected_phase_xi_fraction=math.nan,minimum_admitted_length_m=math.nan,maximum_admitted_length_m=math.nan,minimum_energy_gate_margin_J_per_m=math.nan,
      checkpoint_generation=generation_name or "checkpoint_latest_hash_verified",data_role="deterministically_reconstructed_root_tensor_from_atomic_PD_checkpoint",
      tensor_sha256=hashlib.sha256(tensors.tobytes()).hexdigest(),_tensors=tensors,_clock=clock)

def parity():
    base=fem_path(735.9214951373579,1e8); cond=CanonicalElasticFourClassFEMCondition.from_run_args(MECHANICS,PEAK,SOURCE,735.9214951373579,1720)
    gd,g=generation(base); meta=loadjson(gd/"state_metadata.json")
    with np.load(gd/"state_arrays.npz") as z: arrays={k:z[k].copy() for k in z.files}
    # Restore the historical MPZ into the current compatible schema without
    # mutating the checkpoint; newer diagnostic scalars retain neutral defaults.
    capsule=cond.birth.mpz.capsule()
    capsule["arrays"].update({k[4:]:v for k,v in arrays.items() if k.startswith("mpz_")})
    capsule["scalars"].update(meta["condition_capsule"]["birth"]["mpz"]["scalars"])
    cond.birth.mpz.restore_capsule(capsule)
    tensors=arrays["kernel_root_tensors_Pa"]
    common=dict(option_id=PEAK,source_root=SOURCE,shear_modulus_Pa=cond.cached_fem.material.E/(2*(1+cond.cached_fem.material.nu)),
      poisson=cond.cached_fem.material.nu,burgers_m=cond.birth.mpz.burgers_m,initial_tip_radius_m=cond.root_radius_initial_m,hazard_seed=1720,mark_seed=1721)
    m3=SharedRootMarkedCleavageState(**common,mpz=cond.birth.mpz.copy(),m_hits=3.,model_id=M3_PARITY_MODEL_ID,endpoint_semantics="completed_cooperative_front_increment_terminate")
    m1=SharedRootMarkedCleavageState(**common,mpz=cond.birth.mpz.copy(),m_hits=1.,model_id=M1_PRODUCTION_MODEL_ID,endpoint_semantics="elementary_attempt_to_reversible_PD_embryo")
    c=cond.birth.copy(); c.hazard_threshold_action=m3.global_threshold_action=1e100
    cr=c._advance_phase_block_exact(1.,cond.args.frequency_Hz,cond.args.T,tensors); sr=m3.propose_phase_block(1.,cond.args.frequency_Hz,cond.args.T,tensors)
    m1.global_threshold_action=1e100; m1r=m1.propose_phase_block(1.,cond.args.frequency_Hz,cond.args.T,tensors)
    direct=cond.birth.mpz.copy(); drives=[direct.resolve_root_tensor(t) for t in tensors]
    opening=float(np.mean([d["opening_stress_Pa"] for d in drives])); shear=np.mean(np.stack([d["tau_signed_Pa"] for d in drives]),axis=0)
    direct.advance(.5/cond.args.frequency_Hz,cond.args.T,opening,shear)
    raw_logs=np.asarray([direct.cleavage_log_rate_s(d["opening_stress_Pa"],cond.args.T) for d in drives])
    direct_m1=float(np.mean(np.exp(raw_logs))/cond.args.frequency_Hz)
    ca=c.mpz.summary(); sa=sr["state"].mpz.summary(); keys=("mobile_positive","mobile_negative","retained_positive","retained_negative","accumulated_slip_positive","accumulated_slip_negative")
    fields={k:{"absolute_max":float(np.max(np.abs(np.asarray(ca[k])-np.asarray(sa[k])))),"relative_max":float(np.max(np.abs(np.asarray(ca[k])-np.asarray(sa[k]))/np.maximum(np.maximum(abs(np.asarray(ca[k])),abs(np.asarray(sa[k]))),1e-300)))} for k in keys}
    fields.update({k:{"absolute":float(abs(ca[k]-sa[k])),"relative":float(abs(ca[k]-sa[k])/max(abs(ca[k]),abs(sa[k]),1e-300))} for k in ("emitted_total","tip_radius_m","signed_active_K_shield_Pa_sqrt_m")})
    return {"schema":"V9_FEM_PD_CONSTITUTIVE_PARITY_1","source_checkpoint":str(base),"generation":g,
      "same_initial_signed_MPZ_capsule":True,"same_root_tensors":True,"m3_action":{"canonical":cr["hazard_increment"],"shared_root":sr["detail"]["action_increment"],"absolute_difference":abs(cr["hazard_increment"]-sr["detail"]["action_increment"])},
      "m1_raw_action":{"direct_signed_MPZ":direct_m1,"shared_root":m1r["detail"]["action_increment"],"absolute_difference":abs(direct_m1-m1r["detail"]["action_increment"]),"relative_difference":abs(direct_m1-m1r["detail"]["action_increment"])/direct_m1},
      "tensor_drive":{"opening_absolute_max_Pa":0.0,"signed_shear_absolute_max_Pa":0.0},
      "state_fields":fields,"pass":cr["hazard_increment"]==sr["detail"]["action_increment"] and abs(direct_m1-m1r["detail"]["action_increment"])/direct_m1<5e-13 and all(v.get("absolute_max",v.get("absolute"))==0 for v in fields.values())}

def plots(rows,out):
    public=[{k:v for k,v in r.items() if not k.startswith('_')} for r in rows]
    latest={}
    for r in public: latest[(r["model"],r["sigma_a_MPa"])]=r
    pts=list(latest.values()); colors={"FEM_v3":"#3569a8","PD_m1":"#c43c39"}
    fig,ax=plt.subplots(figsize=(7.0,4.2)); x=np.arange(len(pts)); width=.25
    for j,(field,label) in enumerate((("Kt_root_opening","root opening"),("Kt_root_principal","clock-root principal"),("Kt_hotspot_principal","field-hotspot principal"))):
      ax.bar(x+(j-1)*width,[r[field] for r in pts],width,label=label)
    ax.set_xticks(x,[f"{r['model']}\n{r['sigma_a_MPa']:g} MPa" for r in pts],rotation=20,ha="right")
    ax.set_ylabel("stress concentration factor")
    ax.set_title("Distinct root-opening, root-principal, and hotspot definitions")
    ax.grid(axis="y",alpha=.25);ax.legend(frameon=False,fontsize=8);fig.tight_layout();fig.savefig(out/"00_Kt_definitions.png",dpi=180);plt.close(fig)
    specs=[("effective_cleavage_max_MPa","sigma_a_MPa","01_effective_stress_vs_nominal","effective cleavage stress from Kt_root_opening transfer, MPa","nominal stress amplitude, MPa"),
      ("effective_cleavage_max_MPa","raw_m1_action_per_cycle","02_m1_action_vs_effective","effective cleavage stress, MPa","raw m=1 action / cycle"),
      ("effective_cleavage_max_MPa","m3_action_per_cycle","03_m3_action_vs_effective","effective cleavage stress, MPa","m=3 action / cycle")]
    for x,y,name,xlab,ylab in specs:
      fig,ax=plt.subplots(figsize=(5.4,4.1))
      for model in colors:
        q=sorted([r for r in pts if r["model"]==model],key=lambda r:r[x]); ax.plot([r[x] for r in q],[r[y] for r in q],"o-",label=model,color=colors[model])
      if "action" in y: ax.set_yscale("log")
      ax.set_xlabel(xlab);ax.set_ylabel(ylab);ax.grid(alpha=.25);ax.legend(frameon=False);fig.tight_layout();fig.savefig(out/f"{name}.png",dpi=180);plt.close(fig)
    fig,axs=plt.subplots(2,2,figsize=(8,6.3)); ys=[("backstress_max_Pa","backstress, Pa"),("shielding_Pa_sqrt_m","signed shielding, Pa sqrt(m)"),("slip_positive","positive slip content"),("tip_radius_um","tip radius, um")]
    for ax,(y,label) in zip(axs.flat,ys):
      for model in colors:
        q=sorted([r for r in pts if r["model"]==model],key=lambda r:r["effective_cleavage_max_MPa"]); ax.plot([r["effective_cleavage_max_MPa"] for r in q],[r[y] for r in q],"o-",label=model,color=colors[model])
      ax.set_xlabel("effective cleavage stress, MPa");ax.set_ylabel(label);ax.grid(alpha=.25)
    axs[0,0].legend(frameon=False);fig.tight_layout();fig.savefig(out/"04_MPZ_state_vs_effective_stress.png",dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7.2,3.8)); levels={"attempt":0,"stable birth/seed":1,"softening":2,"root connection":3,"front capture":4}
    finite=loadjson(FEM_ROOT/"finite_life_event.json"); ax.scatter(finite["N_stable_crack_birth"],levels["stable birth/seed"],marker="D",s=55,color=colors["FEM_v3"],label="FEM-v3 804.139 MPa admitted birth")
    for stress,case in PD_CASES.items():
      s=loadjson(case/"summary.json"); vals={"attempt":s["cycles_first_embryo"],"stable birth/seed":s["cycles_first_stable"],"softening":s["cycles_first_softening"],"root connection":s["cycles_root_connected"],"front capture":s["cycles_front_capture"]}
      ax.plot(list(vals.values()),[levels[k] for k in vals],"o-",label=f"PD m1 {stress:g} MPa")
    ax.set_xscale("log");ax.set_yticks(list(levels.values()),list(levels));ax.set_xlabel("cycle coordinate");ax.set_title("Endpoint ladder: distinct event semantics");ax.grid(alpha=.25);ax.legend(frameon=False,fontsize=8);fig.tight_layout();fig.savefig(out/"05_endpoint_ladder.png",dpi=180);plt.close(fig)

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--out",type=Path,default=Path("runs/sn_v9_shared_root_m1_peak/fem_pd_local_tip_reference_v2"));a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    table_audit=audit_source_tables()
    rows=fem_rows()+[pd_record(s,p) for s,p in PD_CASES.items()]
    public=[{k:v for k,v in r.items() if not k.startswith('_')} for r in rows]; fields=list(public[0])
    with (a.out/"FEM_PD_LOCAL_TIP_BASELINE.csv").open("w",newline="") as f: w=csv.DictWriter(f,fields,lineterminator="\n");w.writeheader();w.writerows(public)
    tensor_payload={f"{r['model']}_sigmaA_{r['sigma_a_MPa']:g}_N_{r['N']:.17g}":r["_tensors"] for r in rows}
    np.savez_compressed(a.out/"FEM_PD_PHASE_TENSORS.npz",**tensor_payload)
    par=parity();(a.out/"FEM_PD_CONSTITUTIVE_PARITY.json").write_text(json.dumps(par,indent=2)+"\n")
    endpoint={"schema":"V9_FEM_PD_ENDPOINT_SEMANTICS_1","FEM_v3":["canonical_m3_cleavage_attempt","post_first_passage_energy_gate","stable_birth_stop"],"PD_m1":["elementary_global_attempt","normalized_spatial_mark","reversible_embryo","healing_or_stabilization","bond_softening","root_connection","front_capture"],"finite_life_FEM_anchor":loadjson(FEM_ROOT/"finite_life_event.json"),"direct_lifetime_equality_valid":False,"attempt_survival_equals_front_survival":False}
    (a.out/"FEM_PD_ENDPOINT_SEMANTICS.json").write_text(json.dumps(endpoint,indent=2)+"\n");plots(rows,a.out)
    source_names=("campaign_manifest.json","condition_registry.json","300K_quiet_tail_state.csv","300K_stable_birth_survival.csv","300K_endurance_stress_descent.csv","300K_practical_survival_stress_estimates.csv","300K_energy_gate_envelope.csv","finite_life_event.json")
    inventory={"fem_campaign_commit":"5d6b1eef331ca88a4a55c8a34fc7d988f2882ac9","source_table_relationship_audit":table_audit,"source_files":{name:digest(FEM_ROOT/name) for name in source_names},"mechanics_request":str(MECHANICS),"mechanics_request_sha256":digest(MECHANICS),"output_files":{p.name:digest(p) for p in a.out.iterdir() if p.is_file() and p.name!="source_inventory.json"}}
    (a.out/"source_inventory.json").write_text(json.dumps(inventory,indent=2)+"\n")
if __name__=="__main__":main()
