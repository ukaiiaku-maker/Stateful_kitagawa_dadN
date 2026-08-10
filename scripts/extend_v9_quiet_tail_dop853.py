#!/usr/bin/env python3
"""Qualify and extend one 300 K quiet-tail condition without brute-force FEM."""
from __future__ import annotations

import argparse, csv, json, math, os, tempfile
from pathlib import Path
import numpy as np

from arrhenius_fracture.v9_canonical_four_class_fem import CanonicalFourClassFEMCondition
from arrhenius_fracture.v9_quiet_tail_dop853 import propagate_dop853, validate_dop853
from arrhenius_fracture.v9_quiet_tail_kernel import (
    energy_gate_envelope, kernel_convergence, phase_resolved_cycle,
    restore_condition_checkpoint, save_condition_checkpoint, stationary_marked_renewal)
from scripts.run_v9_quiet_tail_kernel import CLASSES, GATE_FIELDS, SURVIVAL_FIELDS, atomic_merge_csv


def main():
    p=argparse.ArgumentParser(); p.add_argument("material_class",choices=list(CLASSES))
    p.add_argument("--source-root",type=Path,required=True); p.add_argument("--run-args",type=Path,required=True)
    p.add_argument("--checkpoint",type=Path,required=True); p.add_argument("--output",type=Path,required=True)
    p.add_argument("--horizons",default="1e10,1e12,1e14"); a=p.parse_args()
    option,stress=CLASSES[a.material_class]
    c=CanonicalFourClassFEMCondition.from_run_args(a.run_args,option,a.source_root,stress,1720)
    base_kernel,base_summary,_=restore_condition_checkpoint(c,a.checkpoint)
    c.birth.hazard_threshold_action=1e100
    qualifications=[]
    for window in (1e6,1e7):
        q=validate_dop853(c,window,log_cycle_coordinate=True); qualifications.append(q.diagnostics)
        if not q.accepted:
            raise RuntimeError(f"DOP853 direct-overlap qualification failed at {window:g}: {q.diagnostics}")
    xi=np.asarray([.1,.3,.49,.5,.6,.8,1.,1.25,1.5,2.,2.5,3.,3.5,3.99,4.,8.])
    existing=[]; survival_path=a.output/"300K_stable_birth_survival.csv"
    if survival_path.is_file(): existing=list(csv.DictReader(survival_path.open()))
    base_candidates=[r for r in existing if r["material_class"]==a.material_class and float(r["N"])<=c.cycles]
    if base_candidates:
        base=max(base_candidates,key=lambda r:float(r["N"])); total_H=float(base["H_cleave"])
        stable_exponent=-math.log(max(float(base["S_stable_birth"]),1e-300))
    else:
        total_H=float(c.birth.cumulative_cleavage_hazard); stable_exponent=0.0
    prior_kernel=base_kernel; envelope_rows=[]; survival_rows=[]; horizon_kernels=[]
    for horizon in [float(x) for x in a.horizons.split(",")]:
        result=propagate_dop853(c,horizon-c.cycles,log_cycle_coordinate=True)
        if not result.accepted: raise RuntimeError(result.diagnostics)
        c=result.condition; dH=float(result.diagnostics["hazard_increment_fsum"]); total_H+=dH
        kernel=phase_resolved_cycle(c); envelope=energy_gate_envelope(c,kernel,xi)
        decorated=[{"material_class":a.material_class,"option_id":option,"sigma_a_MPa":stress,
                    "N":horizon,**row} for row in envelope]; envelope_rows.extend(decorated)
        renewal=stationary_marked_renewal(dH,kernel["cycle_hazard"],
                                         kernel["cleavage_action_within_cycle"],decorated)
        p_admit=renewal["stationary_attempt_admission_probability"]
        stable_exponent+=dH*p_admit
        survival_rows.append({"material_class":a.material_class,"option_id":option,
            "sigma_a_MPa":stress,"N":horizon,"H_cleave":total_H,
            "S_no_attempt":math.exp(-total_H),"P_at_least_one_attempt":-math.expm1(-total_H),
            "S_stable_birth":math.exp(-stable_exponent),"P_stable_birth":-math.expm1(-stable_exponent),
            "expected_nonpropagating_attempt_count":max(total_H-stable_exponent,0.0),
            "conditional_attempt_admission_probability":p_admit,
            "estimator":"direct-qualified_DOP853_plus_deterministic_marked_renewal"})
        conv=kernel_convergence(prior_kernel,kernel); horizon_kernels.append((horizon,kernel,conv,result.diagnostics))
        save_condition_checkpoint(c,a.output/"quiet_tail_checkpoints"/a.material_class/f"N_{horizon:.17g}",
                                  kernel,{"material_class":a.material_class,"N":horizon,"H_cleave":total_H,
                                          "acceleration":result.diagnostics,"kernel_convergence":conv})
        prior_kernel=kernel
    atomic_merge_csv(a.output/"300K_energy_gate_envelope.csv",GATE_FIELDS,envelope_rows,
                     ("material_class","N","phase_index","Xi"))
    atomic_merge_csv(survival_path,SURVIVAL_FIELDS,survival_rows,("material_class","N"))
    final_conv=horizon_kernels[-1][2]
    stationary=(final_conv["signed_state_relative_error"]<=5e-4 and
                final_conv["root_tensor_relative_error"]<=5e-4 and
                final_conv["backstress_relative_error"]<=5e-4 and
                final_conv["phase_log_rate_absolute_error"]<=5e-3 and
                final_conv["cycle_hazard_relative_error"]<=5e-3)
    final_rows=[r for r in envelope_rows if float(r["N"])==horizon_kernels[-1][0]]
    admissible=sum(float(r["admitted_length_m"])>0 for r in final_rows)
    classification=("no_endurance" if stationary and admissible>0 else
                    "endurance_supported" if stationary and admissible==0 else "undetermined")
    print(json.dumps({"material_class":a.material_class,"qualification":qualifications,
        "stationary_at_final":stationary,"final_convergence":final_conv,
        "final_admissible_points":admissible,"final_gate_points":len(final_rows),
        "stable_crack_birth_asymptotic_classification":classification,
        "survival":survival_rows},indent=2))

if __name__=="__main__": main()
