#!/usr/bin/env python3
"""Audit the standardized Peak 1500 MPa conditioned four-branch pilot."""
import csv, hashlib, json, math
from pathlib import Path
import numpy as np

ROOT=Path("runs/sn_v9_shared_root_m1_peak")

def sha(p):
 h=hashlib.sha256();h.update(Path(p).read_bytes());return h.hexdigest()
def case(root,b):return root/b/"Peak/shielded/sigmaA_1500MPa"
def active(p):
 a=json.loads((p/"v9_generations/ACTIVE.json").read_text())
 return p/"v9_generations"/a["generation"]
def load_arrays(g):
 with np.load(g/"state_arrays.npz") as z:return {k:z[k].copy() for k in z.files}

def main():
 cap=json.loads((ROOT/"conditioned_premark_1500_N50_v1/conditioned_premark_capsule.json").read_text())
 probs=np.asarray(cap["localization"]["normalized_mark_probability"],float);order=np.sort(probs)[::-1]
 pilot=ROOT/"conditioned_pilot_v1"; rows=[]
 for b in ("B000","B001","B002","B003"):
  p=case(pilot,b);m=json.loads((p/"conditioned_branch_manifest.json").read_text());e=m["branch_summary"]["event"]
  g=active(p);md=json.loads((g/"state_metadata.json").read_text());z=load_arrays(g)
  rows.append({"branch_id":b,"selected_site":e["site_id"],"selected_node":e["pd_node_id"],
   "mark_probability":e["mark_probability"],"transition_threshold":m["branch_summary"]["selected_site_transition_threshold"],
   "transition_outcome_uniform":m["branch_summary"]["selected_site_transition_outcome_uniform"],
   "cycles_checkpoint":md["cycles"],"cycles_first_stable":md["pd_scalars"].get("cycles_first_stable"),
   "stable_sites":int(z["pd__stable_sites"].sum()),"front_capture":md["pd_scalars"].get("cycles_front_capture"),
   "right_censored":md["pd_scalars"].get("cycles_front_capture") is None})
 comparisons=[]
 for b in ("B001","B002"):
  a=load_arrays(active(case(pilot,b)));r=load_arrays(active(case(ROOT/"conditioned_restart_v2",b+"_restart")))
  keys=[k for k in a if k in r and np.issubdtype(a[k].dtype,np.number) and a[k].shape==r[k].shape]
  exact={k:bool(np.array_equal(a[k],r[k],equal_nan=True)) for k in keys}
  differences={}
  for k in keys:
   delta=np.abs(a[k].astype(float)-r[k].astype(float));finite=delta[np.isfinite(delta)]
   differences[k]=float(np.max(finite)) if finite.size else 0.0
  comparisons.append({"branch_id":b,"all_numeric_arrays_exact":all(exact.values()),
                      "maximum_absolute_error":max(differences.values(),default=0.),"fields_compared":len(differences),
                      "nonidentical_fields":[k for k,v in exact.items() if not v]})
 out=ROOT/"conditioned_pilot_audit_v1";out.mkdir(exist_ok=True)
 with (out/"branch_results.csv").open("w",newline="") as f:w=csv.DictWriter(f,rows[0],lineterminator="\n");w.writeheader();w.writerows(rows)
 result={"schema":"V9_CONDITIONED_PEAK_1500_PILOT_AUDIT_1","branches":len(rows),
  "mark_support":{"entropy_nats":cap["localization"]["mark_entropy_nats"],
   "effective_support":math.exp(cap["localization"]["mark_entropy_nats"]),"maximum_site_probability":float(order[0]),
   "top5_cumulative_probability":float(order[:5].sum()),"top10_cumulative_probability":float(order[:10].sum()),
   "unique_selected_sites":len(set(r["selected_site"] for r in rows)),"duplicate_selections":len(rows)-len(set(r["selected_site"] for r in rows))},
  "restart_comparisons":comparisons,"B000_long_continuation_is_restartable_partial":True,
  "B000_not_claimed_1e6_complete":True,"branch_results_sha256":sha(out/"branch_results.csv")}
 (out/"audit.json").write_text(json.dumps(result,indent=2)+"\n")

if __name__=="__main__":main()
