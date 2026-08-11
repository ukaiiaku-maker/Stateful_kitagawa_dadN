#!/usr/bin/env python3
"""Executable Peak emission -> delivery -> K=2 bridge ledger at one checkpoint."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from scripts.audit_v9_pd_persistent_sites import build_context
from arrhenius_fracture import sn_pd2d_stateful_v9_transactional as driver
from arrhenius_fracture.v9_four_class_signed_mpz import SignedMPZPreBirthState


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--case-dir", type=Path, required=True)
    p.add_argument("--generation", required=True)
    p.add_argument("--sigma-a-MPa", type=float, required=True)
    p.add_argument("--source-root", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    ns = p.parse_args()
    stored = driver._restore_nonfinite_tags(json.loads((ns.case_dir / "run_args.json").read_text()))
    args, mesh, patch, chain, crack, cached, transaction, sigma_max, sigma_min = build_context(
        ns.case_dir / "run_args.json", ns.sigma_a_MPa
    )
    if "delivery_source_multiplicity" not in stored:
        delattr(args, "delivery_source_multiplicity")
    restored = driver._load_case_checkpoint(
        ns.case_dir / "checkpoint_latest.npz", args=args, case_name="shielded",
        sigma_a_MPa=ns.sigma_a_MPa, mesh=mesh, patch=patch, generation=ns.generation,
    )
    state = restored["pd_state"]
    payload = driver.evaluate_dormant_exact_cycle(
        args=args, shield_on=True, mesh=mesh, patch=patch, pd_state=state,
        crack=crack, plast_chain=chain, cached_fem=cached,
        fem_transaction=transaction, sigma_max=sigma_max, sigma_min=sigma_min,
        ep_gp=restored["ep_gp"], rho_gp=restored["rho_gp"],
        epsp_acc_gp=restored["epsp_acc_gp"], u=restored["u"],
        plastic_work=restored["Wp_total"], cycles=restored["cycles"], dN=1.0,
    )
    log_action = np.asarray(payload["log_birth_action"], float)
    node = int(np.argmax(log_action))
    rates = patch.last_rates
    global_node = int(patch.global_nodes[node])
    phase = int(np.argmax(
        np.asarray(payload["diagnostic_fields"]["equivalent_stress_mid_global_Pa"])[
            :, global_node
        ]
    ))
    voigt = np.asarray(payload["diagnostic_fields"]["sigma_mid_global_voigt_Pa"])[phase, :, global_node]
    tensor = np.array([[voigt[0], voigt[2]], [voigt[2], voigt[1]]], float)

    option = stored["four_class_option_id"]
    signed = SignedMPZPreBirthState(
        option, ns.source_root, shear_modulus_Pa=float(args.E_GPa) * 1e9 / (2 * (1 + args.nu)),
        poisson=args.nu, burgers_m=args.b_m, initial_tip_radius_m=600e-6,
    )
    resolved = signed.resolve_root_tensor(tensor)
    s = signed.state
    geom = dict(s.persistent_site_last_geometry)
    back = np.asarray(signed.summary()["sigma_back_by_system_Pa"], float)
    opening = float(resolved["opening_stress_Pa"])
    drives = np.abs(np.asarray(resolved["tau_signed_Pa"], float))
    effective = np.maximum(drives - back, 0.0)
    per_source = np.array([math.exp(signed.emission_log_rate_per_site_s(x, 300.0)) for x in effective])
    M = float(geom["multiplicity_per_system"])
    authoritative_aggregate = float(M * np.sum(per_source))
    scalar_per_source = float(rates["delivery_rate_s"][node])
    scalar_aggregate_corrected = M * scalar_per_source
    raw = float(rates["nucleation_rate_s"][node])
    Q = float(state.completion[node])
    actual = float(math.exp(log_action[node]))
    pd_area = float(patch.area[node])
    candidate_density = float(args.site_density_m2)
    candidate_area = 1.0 / candidate_density
    row = stored["four_class_registry_audit"]["exact_registry_row"]
    result = {
        "schema": "V9_PD_EXECUTABLE_EMISSION_DELIVERY_BRIDGE_LEDGER_1",
        "scientific_label": "Peak_blunt_4000MPa_K2_suppressed_runout_diagnostic",
        "cycles": restored["cycles"], "pd_node_id": node, "global_fem_node_id": global_node,
        "phase_index": phase,
        "global_FEM_stress_tensor_Pa": tensor.tolist(),
        "opening_stress_Pa": opening,
        "equivalent_stress_Pa": float(payload["diagnostic_fields"]["equivalent_stress_mid_global_Pa"][phase, global_node]),
        "crystallographic_signed_emission_drives_Pa": drives.tolist(),
        "emission_backstress_by_system_Pa": back.tolist(),
        "effective_emission_drive_by_system_Pa": effective.tolist(),
        "emission_barrier_eV_by_system": [float(s.manifest.emission.values_eV(x, 300.0)) for x in effective],
        "per_source_lambda_emit_s-1_by_system": per_source.tolist(),
        "source_density_m-2": float(row["rho_source0_m2"]),
        "reference_source_area_m2": float(row["reference_source_area_um2"]) * 1e-12,
        "source_multiplicity_per_system_M": M,
        "M_times_lambda_sum_authoritative_s-1": authoritative_aggregate,
        "legacy_scalar_per_source_delivery_s-1": scalar_per_source,
        "corrected_scalar_M_times_lambda_s-1": scalar_aggregate_corrected,
        "PD_lumped_point_area_m2": pd_area,
        "PD_candidate_site_density_m-2": candidate_density,
        "PD_candidate_site_area_m2": candidate_area,
        "area_conversion_statement": "M already contains rho_source*reference_source_area; PD point area and candidate density are reported but not multiplied into delivery",
        "delivery_deposited_per_cycle_legacy": float(rates["delivery_events_per_cycle"][node]),
        "delivery_memory_Lambda_legacy": float(state.delivery_memory[node]),
        "K2_completion_Q_legacy": Q,
        "raw_cleavage_barrier_eV": float(crack.deltaG_eV(rates["effective_opening_stress_Pa"][node], 300.0)),
        "raw_cleavage_rate_s-1": raw,
        "state_shift_eV": float(rates["state_shift_eV"][node]),
        "effective_gated_birth_rate_per_cycle": actual,
        "persistent_site_birth_action_per_cycle": actual,
        "suppression_decades_raw_s-1_to_birth_cycle-1": math.log10(raw) - math.log10(actual),
        "omitted_multiplicity_decades": math.log10(M),
        "small_Lambda_K2_birth_boost_decades_if_only_M_corrected": 2 * math.log10(M),
        "authoritative_projection_differs_from_scalar_equivalent_projection": True,
        "unresolved_spatial_state_mapping": "authoritative signed backstress/MPZ state is root-local; repository does not specify broadcast versus node-local signed state for the spatial PD patch",
    }
    ns.out.parent.mkdir(parents=True, exist_ok=True)
    ns.out.write_text(json.dumps(driver._json_safe(result), indent=2, sort_keys=True) + "\n")
    print(json.dumps(driver._json_safe(result), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
