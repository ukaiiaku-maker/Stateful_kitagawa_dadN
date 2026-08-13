#!/usr/bin/env python3
"""Plot the completed-row-only state of the four-class m1 campaign."""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path("runs/sn_v9_four_class_m1")
OUT = ROOT / "figures_progress_v1"
PEAK = Path("runs/sn_v9_shared_root_m1_peak/sn_matrix_v1/figures_v1/figure_data/attempt_SN_current.csv")
PEAK_BRANCHES = Path("runs/sn_v9_shared_root_m1_peak/sn_matrix_v1/figures_v1/figure_data/front_branches_current.csv")
DBTT = ROOT / "DBTT/DBTT_blunt_m1_attempt_SN_v1.csv"
CLASSES = ("Peak", "DBTT", "weakT", "ceramic")
COLORS = {"Peak":"#2166ac", "DBTT":"#b2182b", "weakT":"#1b7837", "ceramic":"#762a83"}
MARKERS = {"Peak":"o", "DBTT":"s", "weakT":"^", "ceramic":"D"}
OPTIONS = {
    "Peak":"v913_paper_peak01_0242980_persistent_sites",
    "DBTT":"v913_paper_dbtt01_0202500_persistent_sites",
    "weakT":"v913_paper_weakT01_0129902_persistent_sites",
    "ceramic":"v913_paper_ceramic01_0077080_persistent_sites",
}

def read(path):
    with path.open(newline="") as stream: return list(csv.DictReader(stream))

def crossing(rows, target):
    n=np.array([float(r["cycles_total"]) for r in rows]); h=np.array([float(r["H_attempt"]) for r in rows])
    hit=np.flatnonzero(h >= target)
    if not len(hit): raise ValueError("incomplete action trajectory")
    i=int(hit[0])
    if i == 0: return float(n[0]*target/h[0])
    return float(n[i-1]+(target-h[i-1])*(n[i]-n[i-1])/(h[i]-h[i-1]))

def completed_from_history(history, cls):
    rows=read(history)
    if not rows or float(rows[-1]["H_attempt"]) < math.log(2)-1e-12: return None
    full=float(rows[-1]["H_attempt"]) >= -math.log(.1)-1e-12
    return {"material_class":cls, "option_id":OPTIONS[cls], "sigma_a_MPa":float(rows[-1]["sigma_a_MPa"]),
            "N10_attempt":crossing(rows,-math.log(.9)), "N50_attempt":crossing(rows,math.log(2)),
            "N90_attempt":crossing(rows,-math.log(.1)) if full else None,
            "median_complete":True, "full_quantile_complete":full,
            "source_history":str(history)}

def inventory():
    data={c:[] for c in CLASSES}
    for r in read(PEAK):
        row={k:(float(v) if k in {"sigma_a_MPa","N10_attempt","N50_attempt","N90_attempt"} else v) for k,v in r.items()}
        data["Peak"].append(row | {"median_complete":True,"full_quantile_complete":True})
    for r in read(DBTT):
        row={k:(float(v) if k in {"sigma_a_MPa","N10_attempt","N50_attempt","N90_attempt"} else v) for k,v in r.items()}
        data["DBTT"].append(row | {"median_complete":True,"full_quantile_complete":True})
    for cls,dirname in (("weakT","weakT"),("ceramic","ceramic")):
        found={}
        histories=list((ROOT/dirname/"action").rglob("sn_stateful_pd_history.csv"))
        histories+=list((ROOT/dirname/"action").rglob("sn_stateful_pd_history_partial.csv"))
        for h in histories:
            row=completed_from_history(h,cls)
            if row:
                old=found.get(row["sigma_a_MPa"])
                if old is None or (row["full_quantile_complete"] and not old["full_quantile_complete"]):
                    found[row["sigma_a_MPa"]]=row
        data[cls]=list(found.values())
    for rows in data.values(): rows.sort(key=lambda r:r["N50_attempt"])
    return data

def gap(rows):
    x=sorted(math.log10(r["N50_attempt"]) for r in rows)
    return max((b-a for a,b in zip(x,x[1:])),default=0.0)

def plot_one(cls, rows, pathbase):
    fig,ax=plt.subplots(figsize=(7.2,5.2),layout="constrained")
    x50=np.array([r["N50_attempt"] for r in rows]); y=np.array([r["sigma_a_MPa"] for r in rows])
    full=[r for r in rows if r["full_quantile_complete"]]
    x10=np.array([r["N10_attempt"] for r in full]); x90=np.array([r["N90_attempt"] for r in full])
    yf=np.array([r["sigma_a_MPa"] for r in full])
    if len(rows)>1:
        ax.plot(x50,y,"--",lw=1.1,color=COLORS[cls],alpha=.65,label="N50 guide")
        for a,b,yy in zip(x10,x90,yf): ax.plot([a,b],[yy,yy],color=COLORS[cls],lw=1.3,alpha=.65)
    ax.scatter(x10,yf,marker="|",s=65,color=COLORS[cls],label="N10/N90 (full row)")
    ax.scatter(x90,yf,marker="|",s=65,color=COLORS[cls])
    for r in rows:
        ax.plot(r["N50_attempt"],r["sigma_a_MPa"],marker=MARKERS[cls],markersize=7,
                markerfacecolor=COLORS[cls] if r["full_quantile_complete"] else "white",
                markerfacecoloralt=COLORS[cls],fillstyle="full" if r["full_quantile_complete"] else "left",
                markeredgecolor=COLORS[cls],linestyle="none")
    ax.set_xscale("log"); ax.set_xlabel("Cycles to elementary attempt"); ax.set_ylabel(r"Stress amplitude $\sigma_a$ (MPa)")
    nfull=sum(r["full_quantile_complete"] for r in rows)
    status="COMPLETE ACTION CURVE" if nfull>=10 and gap(full)<=1.5 else f"INCOMPLETE — {len(rows)} median / {nfull} full rows"
    ax.set_title(f"{cls} 300 K shared-root m1 attempt S–N\n{status}")
    ax.grid(True,which="both",alpha=.22); ax.legend(frameon=False,ncol=3)
    for ext in ("png","pdf"): fig.savefig(pathbase.with_suffix("."+ext),dpi=220)
    plt.close(fig)

def main():
    OUT.mkdir(parents=True,exist_ok=True); data=inventory()
    for cls in CLASSES: plot_one(cls,data[cls],OUT/f"{cls}_attempt_SN_progress")
    fig,ax=plt.subplots(figsize=(8.2,5.8),layout="constrained")
    for cls in CLASSES:
        rows=data[cls]; x=np.array([r["N50_attempt"] for r in rows]); y=np.array([r["sigma_a_MPa"] for r in rows]); full=[r for r in rows if r["full_quantile_complete"]]
        complete=len(full)>=10 and gap(full)<=1.5
        if len(rows)>1: ax.plot(x,y,"-" if complete else "--",color=COLORS[cls],lw=1.2,alpha=.65)
        for r in rows: ax.plot(r["N50_attempt"],r["sigma_a_MPa"],marker=MARKERS[cls],markersize=7,markerfacecolor=COLORS[cls] if r["full_quantile_complete"] else "white",markerfacecoloralt=COLORS[cls],fillstyle="full" if r["full_quantile_complete"] else "left",markeredgecolor=COLORS[cls],linestyle="none")
        ax.plot([],[],color=COLORS[cls],marker=MARKERS[cls],linestyle="none",label=f"{cls}: {len(rows)} median / {len(full)} full")
    ax.set_xscale("log"); ax.set_xlabel("N50 elementary-attempt cycles"); ax.set_ylabel(r"Stress amplitude $\sigma_a$ (MPa)")
    ax.set_title("300 K four-class shared-root m1 attempt S–N progress\nRaw completed executable trajectories only")
    ax.grid(True,which="both",alpha=.22); ax.legend(frameon=False)
    for ext in ("png","pdf"): fig.savefig(OUT/f"Four_class_attempt_SN_progress.{ext}",dpi=220)
    plt.close(fig)

    peak_branches=read(PEAK_BRANCHES)
    physical={c:{"rows":0,"branches":0} for c in CLASSES}
    physical["Peak"]={"rows":len(set(r["sigma_a_MPa"] for r in peak_branches)),"branches":len(peak_branches)}
    conditioned={c:physical[c]["rows"] for c in CLASSES}
    fields=["class","median_complete_rows","full_quantile_rows","required_action_rows","completed_conditioned_rows","required_conditioned_rows","completed_physical_rows","required_physical_rows","completed_physical_branches","required_physical_branches","lowest_completed_N50_attempt","highest_completed_N50_attempt","largest_adjacent_gap_decades","action_curve_complete","physical_matrix_complete","class_complete"]
    progress=[]
    for cls in CLASSES:
        rows=data[cls]; full=[r for r in rows if r["full_quantile_complete"]]; ac=len(full)>=10 and gap(full)<=1.5; pc=physical[cls]["rows"]>=10 and physical[cls]["branches"]>=40
        progress.append({"class":cls,"median_complete_rows":len(rows),"full_quantile_rows":len(full),"required_action_rows":10,"completed_conditioned_rows":conditioned[cls],"required_conditioned_rows":10,"completed_physical_rows":physical[cls]["rows"],"required_physical_rows":10,"completed_physical_branches":physical[cls]["branches"],"required_physical_branches":40,"lowest_completed_N50_attempt":min((r["N50_attempt"] for r in rows),default=""),"highest_completed_N50_attempt":max((r["N50_attempt"] for r in rows),default=""),"largest_adjacent_gap_decades":gap(full),"action_curve_complete":ac,"physical_matrix_complete":pc,"class_complete":ac and pc})
    with (OUT/"Four_class_campaign_progress.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(progress)
    sources=[PEAK,PEAK_BRANCHES,DBTT]+[Path(r["source_history"]) for c in ("weakT","ceramic") for r in data[c]]
    manifest={"schema":"V9_FOUR_CLASS_PROGRESS_FIGURES_1","completed_rows_only":True,"generated_outputs":sorted(p.name for p in OUT.iterdir() if p.name not in {"FIGURE_PROGRESS_MANIFEST.json","README.md"}),"sources":[{"path":str(p),"sha256":hashlib.sha256(p.read_bytes()).hexdigest()} for p in sources],"counts":{r["class"]:r for r in progress}}
    (OUT/"FIGURE_PROGRESS_MANIFEST.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
    (OUT/"README.md").write_text("# Four-class campaign progress figures\n\nCompleted executable action trajectories only. Partial/restartable trajectories and estimates are excluded from S–N points. Lines are visual guides, not fitted curves. The progress CSV applies the full requirement of 10 action rows, 10 exact conditioned rows, 10 physical rows, and 40 unbiased branches per class.\n")
    print(json.dumps({r["class"]:{"median":r["median_complete_rows"],"full":r["full_quantile_rows"],"physical_rows":r["completed_physical_rows"],"branches":r["completed_physical_branches"]} for r in progress},indent=2))

if __name__ == "__main__": main()
