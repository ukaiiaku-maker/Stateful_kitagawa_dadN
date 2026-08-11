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
    p.add_argument("--disable-high-cycle", action="store_true")
    p.add_argument("--sigma-a-MPa", type=float, default=None)
    p.add_argument("--stress-step-source-checkpoint", default="")
    p.add_argument("--stress-step-source-generation", default="")
    p.add_argument("--stress-step-source-sigma-a-MPa", type=float, default=None)
    p.add_argument("--protocol-label", default="")
    a=p.parse_args(); source=json.loads(Path(a.source_run_args).read_text())
    target=driver.build_parser().parse_args([])
    valid=set(vars(target))
    for key,value in source.items():
        if key in valid: setattr(target,key,value)
    metadata_path = None
    if a.stress_step_source_checkpoint and a.stress_step_source_generation:
        metadata_path = (
            Path(a.stress_step_source_checkpoint).parent / "v9_generations"
            / a.stress_step_source_generation / "state_metadata.json"
        )
    elif a.resume:
        generation_root = Path(a.source_run_args).resolve().parent / "v9_generations"
        active = json.loads((generation_root / "ACTIVE.json").read_text())
        metadata_path = generation_root / active["generation"] / "state_metadata.json"
    if metadata_path is not None:
        metadata = driver._restore_nonfinite_tags(json.loads(metadata_path.read_text()))
        # The signed atomic signature, not a later presentation run_args file,
        # is authoritative for exact historical/four-class reconstruction.
        for key, value in metadata["signature"].items():
            if key not in {"case", "sigma_a_MPa"}:
                setattr(target, key, value)
        if "delivery_source_multiplicity" not in metadata["signature"]:
            # Preserve predecessor signatures and their historical unit
            # per-source delivery mapping. Corrected production is launched
            # separately and never silently substituted during resume.
            delattr(target, "delivery_source_multiplicity")
    for key in ("fatigue_model_preset_applied",):
        if key in source: setattr(target,key,source[key])
    target.out=a.out; target.cycles_max=a.cycles_max; target.max_blocks=a.max_blocks
    # Never allow an archived absolute/relative checkpoint target to point the
    # replay back into its immutable source campaign.
    target.checkpoint_path=""
    target.stress_step_source_checkpoint=a.stress_step_source_checkpoint
    target.stress_step_source_generation=a.stress_step_source_generation
    target.stress_step_source_sigma_a_MPa=a.stress_step_source_sigma_a_MPa
    target.protocol_label=a.protocol_label
    target.resume=a.resume; target.skip_existing=False; target.pd_image_policy="none"
    # A certified stationary map is horizon-independent and must not be split
    # into an arbitrary 1e9-cycle orchestration loop. Event guarding remains
    # inside the engine and bounds the actually accepted segment.
    target.pd_high_cycle=not a.disable_high_cycle; target.pd_high_cycle_max_segment=a.cycles_max
    target.pd_high_cycle_start_cycles=min(a.high_cycle_start,a.cycles_max)
    target.checkpoint_every_blocks=10; target.print_every=10
    case=str(source.get("case","shielded")); sigma=(float(source["sigma_a_MPa"])
        if a.sigma_a_MPa is None else float(a.sigma_a_MPa))
    print(json.dumps({"source":a.source_run_args,"case":case,"sigma_a_MPa":sigma,
                      "cycles_max":a.cycles_max,"pd_high_cycle":target.pd_high_cycle},indent=2))
    driver.run_case_stress(target,case,sigma)
if __name__=="__main__": main()
