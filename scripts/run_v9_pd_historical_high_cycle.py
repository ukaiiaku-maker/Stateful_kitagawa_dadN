#!/usr/bin/env python3
"""Replay an exact historical Stateful-PD request with v9 high-cycle orchestration."""
from __future__ import annotations
import argparse, json
from pathlib import Path
from arrhenius_fracture import sn_pd2d_stateful_v9_transactional as driver

def main():
    p=argparse.ArgumentParser(); p.add_argument("--source-run-args",required=True)
    p.add_argument("--out",required=True); p.add_argument("--cycles-max",type=float,required=True)
    p.add_argument("--max-blocks",type=int,default=20000); p.add_argument("--resume",action="store_true")
    p.add_argument("--high-cycle-start",type=float,default=1e4)
    a=p.parse_args(); source=json.loads(Path(a.source_run_args).read_text())
    target=driver.build_parser().parse_args([])
    valid=set(vars(target))
    for key,value in source.items():
        if key in valid: setattr(target,key,value)
    for key in ("fatigue_model_preset_applied",):
        if key in source: setattr(target,key,source[key])
    target.out=a.out; target.cycles_max=a.cycles_max; target.max_blocks=a.max_blocks
    # Never allow an archived absolute/relative checkpoint target to point the
    # replay back into its immutable source campaign.
    target.checkpoint_path=""
    target.resume=a.resume; target.skip_existing=False; target.pd_image_policy="none"
    target.pd_high_cycle=True; target.pd_high_cycle_max_segment=min(a.cycles_max,1e9)
    target.pd_high_cycle_start_cycles=min(a.high_cycle_start,a.cycles_max)
    target.checkpoint_every_blocks=10; target.print_every=10
    case=str(source.get("case","shielded")); sigma=float(source["sigma_a_MPa"])
    print(json.dumps({"source":a.source_run_args,"case":case,"sigma_a_MPa":sigma,
                      "cycles_max":a.cycles_max,"pd_high_cycle":True},indent=2))
    driver.run_case_stress(target,case,sigma)
if __name__=="__main__": main()
