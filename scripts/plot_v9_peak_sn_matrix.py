#!/usr/bin/env python3
"""Read-only provisional Peak S-N matrix tables, plots, and manifest."""
from __future__ import annotations

import argparse, csv, hashlib, json, math, subprocess
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path("runs/sn_v9_shared_root_m1_peak")
STRESSES = [2500,2250,2000,1820,1630,1500,1460,1290,1172]
COLORS = {"attempt":"#0072B2","front":"#D55E00","topology":"#009E73","censor":"#CC79A7"}

def sha(path):
    h=hashlib.sha256();
    with Path(path).open("rb") as f:
        for c in iter(lambda:f.read(1<<20),b""): h.update(c)
    return h.hexdigest()

def read_csv(path):
    with path.open(newline="") as f:return list(csv.DictReader(f))

def write_csv(path, rows, fields):
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)

def find_case(branch, stress):
    m=list(Path(branch).rglob(f"sigmaA_{stress}MPa/summary.json"))
    if len(m)!=1: raise RuntimeError(f"expected one {stress} MPa summary under {branch}; found {len(m)}")
    return m[0].parent

def quantiles_from_action(stress):
    if stress in (2250,2500):
        p=ROOT/f"sn_matrix_v1/attempt_{stress}/Peak/shielded/sigmaA_{stress}MPa/sn_stateful_pd_history.csv"
        rows=read_csv(p); n=np.array([0.]+[float(r["cycles_total"]) for r in rows]); h=np.array([0.]+[float(r["H_attempt"]) for r in rows])
        vals={q:float(np.interp(t,h,n)) for q,t in (("N10",-.0+(-math.log(.9))),("N50",math.log(2)),("N90",-math.log(.1)))}
        return vals,p,"interpolated_actual_action"
    p=ROOT/"attempt_sn_v1/Peak_blunt_m1_attempt_SN_v2.csv"; rows=read_csv(p)
    own=[r for r in rows if float(r["sigma_a_MPa"])==stress]
    if len(own)!=3: raise RuntimeError(f"missing accepted attempt quantiles for {stress}")
    vals={r["quantile"]:float(r["cycles"]) for r in own}
    methods={r["method"] for r in own}
    kind="stationary_tail_supported" if any("stationary" in x for x in methods) else ("exact_action_boundary" if any(r["exact_endpoint_run"]=="True" for r in own) else "interpolated_actual_action")
    return vals,p,kind

def summary_fields(summary):
    attempt=summary.get("last_marked_attempt") or {}; fields=attempt.get("local_fields") or {}
    root=fields.get("effective_opening_stress_Pa",summary.get("root_opening_stress_max_Pa",np.nan))
    principal=fields.get("root_principal_stress_Pa",summary.get("fem_sigma1_cycle_final_max_Pa",np.nan))
    hot=summary.get("fem_sigma1_cycle_final_max_Pa",np.nan)
    mech=summary.get("mechanics_validity_classification",summary.get("mechanics_validity_class","not_recorded"))
    def mpa(value):
        try: value=float(value)
        except (TypeError,ValueError): return np.nan
        return value/1e6 if np.isfinite(value) else np.nan
    return mpa(root),mpa(principal),mpa(hot),mech

def branch_specs():
    out=[]
    for s in (2500,2250):
        for b in ("B000","B001","B002","B003"): out.append((s,b,ROOT/f"sn_matrix_v1/front_{s}/{b}"))
    for b in ("B000","B001","B002","B003"):
        out.append((1500,b,ROOT/f"conditioned_pilot_v1/{b}"))
    for b in ("B004","B005","B006"):
        out.append((1500,b,ROOT/f"conditioned_ensemble_1500_v1/{b}"))
    # B007 is intentionally excluded from the montage but is valid matrix data.
    out.append((1500,"B007",ROOT/"conditioned_ensemble_1500_v1/B007"))
    return out

def km_quantile(delays, events, q):
    order=np.argsort(delays); t=np.asarray(delays)[order]; e=np.asarray(events)[order]
    surv=1.; n=len(t)
    for x in np.unique(t):
        at=np.sum(t>=x); d=np.sum((t==x)&e)
        if d: surv*=1-d/at
        if 1-surv>=q:return "identified",float(x)
    return "lower_bound",float(np.max(t))

def style(ax):
    ax.grid(True,which="both",alpha=.22); ax.spines[["top","right"]].set_visible(False)
    ax.set_xlabel("cycles to event, N"); ax.set_ylabel(r"stress amplitude $\sigma_a$ [MPa]")

def save(fig,out,name):
    fig.savefig(out/"current"/(name+".png"),dpi=190,bbox_inches="tight",pil_kwargs={"compress_level":9})
    fig.savefig(out/"current"/(name+".pdf"),bbox_inches="tight");plt.close(fig)

def main(argv=None):
    ap=argparse.ArgumentParser();ap.add_argument("--out",type=Path,required=True);a=ap.parse_args(argv)
    out=a.out; current=out/"current"; data=out/"figure_data"; final=out/"final"; end=out/"end_states"
    for d in (current,data,final,end):d.mkdir(parents=True,exist_ok=True)
    commit=subprocess.check_output(["git","rev-parse","HEAD"],text=True).strip()
    attempt=[]; aq={}
    for s in STRESSES:
        vals,source,kind=quantiles_from_action(s); aq[s]=vals
        # Prefer a physical summary, otherwise the accepted action summary.
        candidates=list((ROOT/f"sn_matrix_v1/front_{s}").rglob(f"sigmaA_{s}MPa/summary.json"))
        if not candidates:candidates=list((ROOT/f"attempt_survival_{s}_v1").rglob(f"sigmaA_{s}MPa/summary.json"))
        if not candidates and s in (2250,2500):candidates=list((ROOT/f"sn_matrix_v1/attempt_{s}").rglob(f"sigmaA_{s}MPa/summary.json"))
        summ=json.loads(candidates[0].read_text()) if candidates else {}; root,_,hot,mech=summary_fields(summ)
        attempt.append(dict(sigma_a_MPa=s,sigma_min_MPa=2*s*.1/.9,sigma_max_MPa=2*s/.9,N10_attempt=vals["N10"],N50_attempt=vals["N50"],N90_attempt=vals["N90"],root_opening_max_MPa=root,FEM_sigma1_hotspot_max_MPa=hot,mechanics_validity_class=mech,action_source_path=source,exact_or_interpolated=kind,source_commit=commit))
    af=list(attempt[0]);write_csv(data/"attempt_SN_current.csv",attempt,af)
    branches=[]
    for s,b,path in branch_specs():
        try:case=find_case(path,s)
        except RuntimeError as exc: print(f"SKIP optional branch: {exc}");continue
        p=case/"summary.json";d=json.loads(p.read_text()); attempt_info=d.get("last_marked_attempt") or {}; stable=d.get("cycles_first_stable"); front=d.get("cycles_front_capture")
        status="captured" if front is not None else ("administrative_right_censor" if d.get("status") in ("right_censored","captured") or stable is not None else "partial_restartable")
        local=attempt_info.get("local_fields") or {}; node=attempt_info.get("selected_node",d.get("stable_crack_birth_node")); site=attempt_info.get("selected_site",d.get("stable_crack_birth_site_id"))
        prob=attempt_info.get("selected_probability",attempt_info.get("mark_probability",np.nan)); seed=d.get("pd_primary_seed_node_final",-1)
        morphology="captured_root_connected_front" if front is not None else ("root_connected_uncaptured" if d.get("cycles_root_connected") is not None else "damaged_uncaptured")
        _,_,_,mech=summary_fields(d)
        branches.append(dict(sigma_a_MPa=s,branch_id=b,source_path=case,selected_site=site,selected_node=node,mark_probability=prob,attempt_cycle=attempt_info.get("cycle",d.get("cycles_first_embryo")),stable_cycle=stable,first_softening_cycle=d.get("cycles_first_softening"),root_connection_cycle=d.get("cycles_root_connected"),front_capture_cycle=front,post_stable_exposure_cycles=(float(d["cycles_total"])-float(stable)) if stable is not None else "",observation_status=status,max_bond_damage=d.get("pd_bond_damage_final_max"),broken_bond_count=d.get("pd_broken_bonds_final"),connected_extent_m=d.get("pd_connected_extent_final_m"),active_front_length_m=d.get("pd_active_front_length_final_m"),primary_seed=seed,morphology_class=morphology,mechanics_validity_class=mech))
    bf=list(branches[0]);write_csv(data/"front_branches_current.csv",branches,bf)
    summaries=[]
    for s in STRESSES:
        bb=[x for x in branches if x["sigma_a_MPa"]==s]; cap=[x for x in bb if x["observation_status"]=="captured"]; cens=[x for x in bb if x["observation_status"]=="administrative_right_censor"]; partial=[x for x in bb if x["observation_status"]=="partial_restartable"]
        vals={}
        if bb:
            observations=[float(x["front_capture_cycle"]) if x["observation_status"]=="captured"
                          else float(x["stable_cycle"])+float(x["post_stable_exposure_cycles"])
                          for x in bb]
            events=[x["observation_status"]=="captured" for x in bb]
            for label,q in (("N10",.1),("N50",.5),("N90",.9)):vals[label]=km_quantile(observations,events,q)
        else:
            for label in ("N10","N50","N90"):vals[label]=("unavailable","")
        common=[x for x in bb if x["observation_status"]=="captured" or float(x["post_stable_exposure_cycles"] or 0)>=1e6-1]
        frac=(sum(x["observation_status"]=="captured" and float(x["post_stable_exposure_cycles"])<=1e6+1 for x in common)/len(common)) if common else ""
        summaries.append(dict(sigma_a_MPa=s,physical_branch_count=len(bb),capture_count=len(cap),administrative_censor_count=len(cens),partial_count=len(partial),capture_fraction_at_1e6=frac,N10_front_status=vals["N10"][0],N10_front_or_bound=vals["N10"][1],N50_front_status=vals["N50"][0],N50_front_or_bound=vals["N50"][1],N90_front_status=vals["N90"][0],N90_front_or_bound=vals["N90"][1],median_identifiable=vals["N50"][0]=="identified"))
    write_csv(data/"front_summary_current.csv",summaries,list(summaries[0]))
    # 1 attempt S-N
    fig,ax=plt.subplots(figsize=(7.2,5.2)); ss=np.array(STRESSES); n10=np.array([aq[s]["N10"] for s in STRESSES]);n50=np.array([aq[s]["N50"] for s in STRESSES]);n90=np.array([aq[s]["N90"] for s in STRESSES])
    ax.fill_betweenx(ss,n10,n90,color=COLORS["attempt"],alpha=.15,label="N10–N90");ax.plot(n50,ss,"-o",color=COLORS["attempt"],label="N50")
    for row in attempt:
        marker={"exact_action_boundary":"s","stationary_tail_supported":"D"}.get(row["exact_or_interpolated"],"o");ax.scatter(row["N50_attempt"],row["sigma_a_MPa"],marker=marker,s=42,color=COLORS["attempt"])
    ax.set_xscale("log");style(ax);ax.legend();ax.set_title("Peak Elementary-Attempt S–N Response\n300 K, R=0.1, 1000 Hz; 600 µm blunt notch; one 10 µm segment");save(fig,out,"01_attempt_SN_current")
    # 2 physical branches
    fig,ax=plt.subplots(figsize=(7.2,5.2));
    for x in branches:
        y=x["sigma_a_MPa"]+(int(hashlib.sha1(x["branch_id"].encode()).hexdigest()[:4],16)%17-8)*.65
        if x["observation_status"]=="captured":ax.scatter(x["front_capture_cycle"],y,c=COLORS["front"],s=45,marker="o")
        elif x["observation_status"]=="administrative_right_censor":ax.scatter(float(x["stable_cycle"])+float(x["post_stable_exposure_cycles"]),y,facecolors="none",edgecolors=COLORS["censor"],s=55,marker=">")
        else:ax.scatter(float(x["attempt_cycle"]),y,c="gray",s=18)
    ax.set_xscale("log");style(ax);ax.set_title("Peak Physical Front-Capture Branches\nfilled: capture; open >: administrative right-censor");save(fig,out,"02_front_branches_current")
    # 3 provisional conditional front response/bounds plus attempt reference
    fig,ax=plt.subplots(figsize=(7.2,5.2));ax.plot(n50,ss,"--",lw=1.1,color=COLORS["attempt"],label="elementary-attempt N50")
    for row in summaries:
        s=row["sigma_a_MPa"]
        for label,marker in (("N10","v"),("N50","o"),("N90","^")):
            status=row[f"{label}_front_status"];v=row[f"{label}_front_or_bound"]
            if v=="":continue
            x=float(v); fill=COLORS["front"] if status=="identified" else "none"
            ax.scatter(x,s,marker=marker,s=52,facecolors=fill,edgecolors=COLORS["front"])
            if status!="identified":ax.annotate("bound →",(x,s),xytext=(4,4),textcoords="offset points",fontsize=7)
    ax.set_xscale("log");style(ax);ax.legend();ax.set_title("PROVISIONAL — INCOMPLETE FRONT MATRIX\ncensor-aware absolute front-time quantiles/bounds; dashed = attempt N50")
    save(fig,out,"03_combined_front_SN_provisional")
    # 4 decomposition
    fig,ax=plt.subplots(figsize=(7.2,5.2));ax.plot(n50,ss,"--o",color=COLORS["attempt"],label="median attempt life")
    for x in branches:
        if x["observation_status"]=="captured":ax.scatter(float(x["post_stable_exposure_cycles"]),x["sigma_a_MPa"],color=COLORS["topology"],alpha=.8,s=32)
    ax.set_xscale("log");style(ax);ax.legend();ax.set_title("Attempt versus Captured Topology Delay\n(no sum-of-medians total is shown)");save(fig,out,"04_attempt_topology_decomposition")
    # 5 capture fractions with Wilson intervals
    fig,ax=plt.subplots(figsize=(6.8,4.7));
    for row in summaries:
        if row["capture_fraction_at_1e6"]=="":continue
        s=row["sigma_a_MPa"]; n=row["capture_count"]+row["administrative_censor_count"];k=row["capture_count"];p=k/n;z=1.96;den=1+z*z/n;cen=(p+z*z/(2*n))/den;half=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
        ax.errorbar(s,p,yerr=[[p-max(0,cen-half)],[min(1,cen+half)-p]],fmt="o",color=COLORS["front"],capsize=4)
    ax.set_ylim(-.04,1.04);ax.set_xlabel(r"stress amplitude $\sigma_a$ [MPa]");ax.set_ylabel(r"P(front by $10^6$ post-stable cycles | stable seed)");ax.grid(alpha=.25);ax.set_title("Peak Front Capture at Common Observation Horizon\n95% Wilson intervals")
    save(fig,out,"05_capture_by_1e6")
    # 6 stress mapping
    fig,ax=plt.subplots(figsize=(7.2,5.0));x=np.array([r["sigma_a_MPa"] for r in attempt]);root=np.array([float(r["root_opening_max_MPa"]) for r in attempt]);hot=np.array([float(r["FEM_sigma1_hotspot_max_MPa"]) for r in attempt]);
    ax.plot(x,root,"o-",label="root opening maximum");ax.plot(x,hot,"s-",label="field principal-stress hotspot");ax.plot(x,hot,"^--",alpha=.55,label="clock-node principal (fallback: recorded root hotspot)")
    ax.set_xlabel(r"stress amplitude $\sigma_a$ [MPa]");ax.set_ylabel("local stress [MPa]");ax.grid(alpha=.25);ax.legend(fontsize=8);ax.set_title("Peak Local Root Stress Mapping\nexplicit stress definitions; missing values remain unfilled")
    save(fig,out,"06_root_stress_vs_sigma_a")
    readme=("# Peak S–N provisional visualization package\n\n"
      "Endpoint: first captured front after a stable spatial seed. The elementary-attempt curve describes the conserved global cleavage clock; it is not the front-capture curve. Filled branch markers are captures and open right-pointing markers are finite administrative right-censors, never arrests.\n\n"
      "Representative end states: captured/censored B001/B000 at 2500 and 2250 MPa, and captured B004/long-censored B005 at 1500 MPa. All figures are post-processed from existing summaries and terminal arrays. No simulation or image sequence is produced.\n\n"
      "Fixed context: Peak, 300 K, R=0.1, 1000 Hz, ideal smooth 600 µm-radius blunt notch, one 10 µm correlated front segment, stress axis = sigma_a, no specimen-width scaling, no microflaw calibration.\n")
    (out/"README.md").write_text(readme)
    # Manifest after renderer has run; include every figure in current and end_states.
    records=[]
    sources=sorted({Path(x["source_path"]) for x in branches})
    source_hashes={str(p):sha(p/"summary.json") for p in sources if (p/"summary.json").is_file()}
    now=datetime.now(timezone.utc).isoformat()
    for p in sorted(list(current.glob("*"))+list(end.glob("*.png"))):
        records.append(dict(filename=str(p.relative_to(out)),figure_type=("representative_pd_end_state" if p.parent==end else "summary_plot"),status="provisional",source_branch_paths=[str(x) for x in sources],source_summary_hashes=source_hashes,source_commit=commit,creation_script=("scripts/render_v9_peak_pd_end_states.py" if p.parent==end or "07_" in p.name else "scripts/plot_v9_peak_sn_matrix.py"),creation_timestamp=now,file_size=p.stat().st_size))
    (out/"FIGURE_MANIFEST.json").write_text(json.dumps({"schema":"V9_PEAK_SN_FIGURE_MANIFEST_1","figures":records},indent=2)+"\n")

if __name__=="__main__":main()
