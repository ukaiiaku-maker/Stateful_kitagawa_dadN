#!/usr/bin/env python3
"""Reproducible synthetic qualification of the v9 PD high-cycle state machine."""
from __future__ import annotations
import argparse, csv, json, math
from pathlib import Path
import numpy as np
from arrhenius_fracture.v9_pd_high_cycle import (
    ActiveState, CycleEvaluation, DormantPDHighCycleEngine, HighCycleConfig,
    ProtectedSignatures, _digest,
)

class LinearAdapter:
    def __init__(self, rate, threshold, drift=0.0):
        self.x=np.zeros(1); self.rate=rate; self.threshold=np.array([threshold]); self.action=np.zeros(1)
        self.cycles=0.0; self.ledger=0.0; self.log_action=np.full(1,-math.inf); self.status=np.zeros(1,np.uint8); self.rng=np.random.default_rng(91); self.drift=drift
    def dormant_eligibility(self): return (not np.any(self.status), "dormant")
    def active_state(self): return ActiveState(self.x, (("x",(1,),"float64"),))
    def restore_active_state(self,s,v): self.x=np.asarray(v,float).copy()
    def protected_signatures(self): return ProtectedSignatures(_digest(self.ledger),_digest((self.threshold,self.action,self.rng.bit_generator.state)),_digest(self.status))
    def exact_private_cycle(self):
        s=self.active_state(); e=ActiveState(s.vector+self.drift,s.specification); lr=np.array([math.log(self.rate)])
        return CycleEvaluation(s,e,lr,{"ledger":1.0},np.array([0.,1.]),lr[None,:],{},"dormant",self.protected_signatures().topology)
    def commit_private_cycle(self,e): self.x=e.state_end.vector.copy(); self.commit_ledger_increments(e.ledger_increments)
    def commit_ledger_increments(self,x): self.ledger+=float(x.get("ledger",0.0))
    def remaining_birth_actions(self): return self.threshold-self.action
    def commit_birth_action(self,q,n): self.action+=q
    def commit_log_birth_action(self,q,n): self.log_action=np.logaddexp(self.log_action,q); self.action=np.exp(self.log_action)
    def physical_cycles(self): return self.cycles
    def set_physical_cycles(self,n): self.cycles=float(n)

def main():
    p=argparse.ArgumentParser(); p.add_argument("--out",default="runs/v9_pd_high_cycle_qualification"); a=p.parse_args()
    out=Path(a.out); out.mkdir(parents=True,exist_ok=True)
    cases=[("stationary_1e12",1e12,1e-20,10.,0.),("stationary_1e14",1e14,1e-20,10.,0.),
           ("projective_1e8",1e8,1e-30,10.,1e-9),
           ("event_guard_near_1e12",1e14,1e-12,(1e12+.25)*1e-12,0.)]
    rows=[]; modes={}
    for name,horizon,rate,threshold,drift in cases:
        model=LinearAdapter(rate,threshold,drift)
        cfg=HighCycleConfig(periodic_max_iterations=2,periodic_relative_tolerance=1e-12,
            periodic_admission_distance=1e-12,projective_state_tolerance=1e-12,
            projective_log_hazard_tolerance=1e-12,projective_initial_cycles=10**8)
        engine=DormantPDHighCycleEngine(model,cfg); result=engine.advance(horizon)
        rows.append({"case":name,"requested_cycles":horizon,"consumed_cycles":result.cycles_consumed,
            "event_guard":result.event_guard_reached,"cumulative_action":model.action[0],
            "exact_cycle_maps":result.exact_map_evaluations,"accepted_projected_cycles":result.accepted_projected_cycles})
        modes[name]=[{"mode":m.mode,"cycles":m.cycles,"accepted":m.accepted,"detail":m.detail} for m in result.modes]
    with (out/"comparison.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=rows[0]); w.writeheader(); w.writerows(rows)
    (out/"mode_history.json").write_text(json.dumps(modes,indent=2)+"\n")
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig,ax=plt.subplots(figsize=(6.2,3.7)); ax.bar([r["case"] for r in rows],[max(r["exact_cycle_maps"],1) for r in rows])
        ax.set_yscale("log"); ax.set_ylabel("exact one-cycle map evaluations"); ax.tick_params(axis="x",rotation=25)
        ax.set_title("v9 PD high-cycle kernel qualification"); fig.tight_layout(); fig.savefig(out/"map_count.png",dpi=180); plt.close(fig)
    except ImportError: pass
    print(json.dumps(rows,indent=2))
if __name__=="__main__": main()
