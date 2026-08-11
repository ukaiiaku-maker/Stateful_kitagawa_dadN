#!/usr/bin/env python3
"""Fail-closed conversion of an analysis boundary into a conditioned premark capsule."""
from __future__ import annotations
import argparse,hashlib,json,math
from pathlib import Path
import numpy as np

TARGET=math.log(2.0)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):h.update(chunk)
    return h.hexdigest()

def main(argv=None):
    p=argparse.ArgumentParser();p.add_argument("--source",type=Path,required=True)
    p.add_argument("--replay",type=Path,required=True);p.add_argument("--out",type=Path,required=True)
    p.add_argument("--branch-namespace",required=True);a=p.parse_args(argv)
    source=json.loads((a.source/"summary.json").read_text())
    replay=json.loads((a.replay/"summary.json").read_text());loc=replay.get("analysis_boundary_localization")
    if loc is None:raise RuntimeError("exact phase localization is absent from replay")
    for label,data in (("source",source),("replay",replay)):
        if not math.isclose(float(data["H_attempt_final"]),TARGET,rel_tol=0,abs_tol=2e-15):raise RuntimeError(f"{label} H != ln2")
        if int(data["global_attempt_count_final"])!=0:raise RuntimeError(f"{label} attempt ledger is nonzero")
        if data.get("cycles_first_embryo") is not None or data.get("cycles_first_stable") is not None:raise RuntimeError(f"{label} contains embryo/stable state")
    if source["cycles_total"]!=replay["cycles_total"]:raise RuntimeError("replay crossing cycle differs from approved boundary")
    protocol=replay["conditional_survival_protocol"]
    if protocol.get("mark_rng_consumed") is not False:raise RuntimeError("mark RNG provenance is not untouched")
    with np.load(a.replay/"pd_state_final.npz") as z:
        if np.any(z["embryo_sites"]) or np.any(z["stable_sites"]):raise RuntimeError("nonzero realized embryo/stable sites")
        if np.any(z["bond_damage"]):raise RuntimeError("conditioned boundary has bond damage")
        if int(z["primary_seed_node"])!=-1 or bool(z["active_front"]):raise RuntimeError("conditioned boundary has seed/front topology")
        candidate_hash=hashlib.sha256(np.ascontiguousarray(z["candidate_sites"]).view(np.uint8)).hexdigest()
        geometry_hash=hashlib.sha256(np.ascontiguousarray(z["xy"]).view(np.uint8)).hexdigest()
    modes=json.loads((a.replay/"v9_pd_high_cycle_mode_history.json").read_text())
    bracket=[{"cycles":row["cycles_total"],"H":row["H_attempt_after_segment"]}
             for row in modes if "H_attempt_after_segment" in row][-2:]
    if not bracket or bracket[-1]["H"]>=TARGET:raise RuntimeError("replay lacks a strict pre-crossing action bracket")
    probs=np.asarray(loc["normalized_mark_probability"],float)
    if np.any(probs<0) or not math.isclose(float(probs.sum()),1.,rel_tol=0,abs_tol=2e-14):raise RuntimeError("invalid mark normalization")
    capsule={
      "schema":"V9_CONDITIONED_ATTEMPT_PREMARK_CAPSULE_1",
      "protocol":"conditioned_attempt_premark","source_analysis_checkpoint":sha(a.source/"checkpoint_latest.npz"),
      "source_replay_checkpoint":sha(a.replay/"checkpoint_latest.npz"),"conditioned_action":TARGET,
      "physical_threshold_draw":False,"conditioned_first_attempt_pending":True,
      "mark_rng_consumed":False,"topology_continuation_permitted":True,
      "exact_phase_localization_present":True,"branch_seed_namespace":a.branch_namespace,
      "crossing_cycle":replay["cycles_total"],"action_bracket":bracket,
      "candidate_population_sha256":candidate_hash,"geometry_sha256":geometry_hash,
      "source_summary_sha256":sha(a.source/"summary.json"),"replay_summary_sha256":sha(a.replay/"summary.json"),
      "source_generation":json.loads((a.source/"v9_generations/ACTIVE.json").read_text())["generation"],
      "replay_generation":json.loads((a.replay/"v9_generations/ACTIVE.json").read_text())["generation"],
      "localization":loc,
    }
    a.out.mkdir(parents=True,exist_ok=False)
    (a.out/"conditioned_premark_capsule.json").write_text(json.dumps(capsule,indent=2)+"\n")
    (a.out/"manifest.json").write_text(json.dumps({"capsule_sha256":sha(a.out/"conditioned_premark_capsule.json"),
      "source_files":{str(x):sha(x) for x in (a.source/"checkpoint_latest.npz",a.source/"summary.json",a.replay/"checkpoint_latest.npz",a.replay/"summary.json")}},indent=2)+"\n")

if __name__=="__main__":main()
