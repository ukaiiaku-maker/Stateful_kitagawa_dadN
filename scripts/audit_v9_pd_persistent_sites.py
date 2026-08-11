#!/usr/bin/env python3
"""Reconstruct matched persistent-site diagnostics at atomic v9 checkpoints."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from arrhenius_fracture.config import ElasticProperties
from arrhenius_fracture.sn_arrhenius_chain import build_chain_from_namespace
from arrhenius_fracture.sn_feature_geometry_v8_7 import (
    BluntNotchGeometry, identify_feature_surface_nodes, local_root_xy,
    make_blunt_edge_notch_mesh,
)
from arrhenius_fracture.sn_intact_fem import plane_strain_D
from arrhenius_fracture import sn_pd2d_stateful_v9_transactional as driver
from arrhenius_fracture.v9_cached_fem import CachedIntactFEM
from arrhenius_fracture.v9_fem_transaction import EmbeddedFEMTransaction
from arrhenius_fracture.v9_stateful_peridynamics import V9StatefulPDPatch


def build_context(run_args: Path, sigma_a_MPa: float):
    source = driver._restore_nonfinite_tags(json.loads(run_args.read_text()))
    args = driver.build_parser().parse_args([])
    valid = set(vars(args)) | {
        "four_class_option_id", "four_class_transfer_contract",
        "four_class_registry_audit", "audited_four_class_barriers",
        "fatigue_model_preset_applied",
    }
    for key, value in source.items():
        if key in valid:
            setattr(args, key, value)
    if source.get("fatigue_model_preset_applied"):
        args.fatigue_model_preset_applied = True
    mat = ElasticProperties(E=args.E_GPa * 1e9, nu=args.nu, b=args.b_m, Tm=args.Tm_K)
    Dmat = plane_strain_D(mat)
    geom = BluntNotchGeometry(
        args.Lx, args.Ly, args.notch_depth_m, args.notch_half_height_m,
        feature_type=args.feature_type, root_radius_m=args.notch_root_radius_m,
        opening_angle_deg=args.notch_opening_angle_deg,
        path_refine_length_m=args.path_refine_length_m,
        path_refine_half_height_m=args.path_refine_half_height_m,
    )
    mesh, bnd, _ = make_blunt_edge_notch_mesh(
        geom, nx=args.nx, ny=args.ny, jitter=args.jitter,
        root_h_fine=args.root_h_fine, seed=args.seed,
    )
    feature = identify_feature_surface_nodes(mesh, geom)
    root = local_root_xy(mesh, feature)
    patch = V9StatefulPDPatch(
        mesh, geom, root, mat, driver._pd_config_from_args(args),
        feature_surface_global_nodes=feature,
    )
    chain = build_chain_from_namespace(args, mat.b)
    crack = driver.build_crack_barrier(args)
    cached = CachedIntactFEM(mesh, bnd, mat, Dmat)
    sigma_max = 2.0 * sigma_a_MPa * 1e6 / (1.0 - args.R)
    sigma_min = args.R * sigma_max
    transaction = EmbeddedFEMTransaction(
        mesh=mesh, boundaries=bnd, material=mat, Dmat=Dmat,
        plastic_chain=chain, args=args, sigma_max_Pa=sigma_max,
        sigma_min_Pa=sigma_min, cached_fem=cached,
    )
    return args, mesh, patch, chain, crack, cached, transaction, sigma_max, sigma_min


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-args", type=Path, required=True)
    parser.add_argument("--sigma-a-MPa", type=float, required=True)
    parser.add_argument("--case", default="shielded")
    parser.add_argument("--checkpoint", action="append", nargs=3,
                        metavar=("LABEL", "STORE", "GENERATION"), required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--table-prefix", default="690MPa_persistent_site")
    args_cli = parser.parse_args()
    rows = []
    for label, store_text, generation in args_cli.checkpoint:
        args, mesh, patch, chain, crack, cached, transaction, sigma_max, sigma_min = build_context(
            args_cli.run_args, args_cli.sigma_a_MPa
        )
        if "delivery_source_multiplicity" not in driver._restore_nonfinite_tags(
            json.loads(args_cli.run_args.read_text())
        ):
            delattr(args, "delivery_source_multiplicity")
        restored = driver._load_case_checkpoint(
            Path(store_text) / "checkpoint_latest.npz", args=args,
            case_name=args_cli.case, sigma_a_MPa=args_cli.sigma_a_MPa,
            mesh=mesh, patch=patch, generation=generation,
        )
        state = restored["pd_state"]
        payload = driver.evaluate_dormant_exact_cycle(
            args=args, shield_on=args_cli.case == "shielded", mesh=mesh, patch=patch,
            pd_state=state, crack=crack, plast_chain=chain, cached_fem=cached,
            fem_transaction=transaction, sigma_max=sigma_max, sigma_min=sigma_min,
            ep_gp=restored["ep_gp"], rho_gp=restored["rho_gp"],
            epsp_acc_gp=restored["epsp_acc_gp"], u=restored["u"],
            plastic_work=restored["Wp_total"], cycles=restored["cycles"], dN=1.0,
        )
        rates = patch.last_rates
        nodes = np.asarray(state.site_node_index, int)
        available = np.asarray(state.site_status, int) == 0
        log_rate_node = np.asarray(payload["log_birth_action"], float)
        cumulative_node = np.asarray(state.birth_cumulative_hazard, float)
        threshold = np.asarray(state.site_birth_threshold, float)
        for site_id, node in enumerate(nodes):
            log_rate = float(log_rate_node[node])
            rate = math.exp(log_rate) if log_rate >= math.log(np.nextafter(0.0, 1.0)) else 0.0
            remaining = max(float(threshold[site_id] - cumulative_node[node]), 0.0)
            rows.append({
                "checkpoint": label, "cycles": restored["cycles"],
                "site_id": site_id, "pd_node_id": int(node),
                "x_m": float(patch.xy[node, 0]), "y_m": float(patch.xy[node, 1]),
                "available": bool(available[site_id]),
                "cumulative_action": float(cumulative_node[node]),
                "remaining_threshold_action": remaining,
                "log_birth_rate_per_cycle": log_rate,
                "birth_rate_per_cycle": rate,
                "instantaneous_wait_cycles": remaining / rate if available[site_id] and rate > 0 else math.inf,
                "delivery_rate_s": float(rates["delivery_rate_s"][node]),
                "delivery_memory": float(state.delivery_memory[node]),
                "K2_completion": float(state.completion[node]),
                "raw_cleavage_rate_s": float(rates["nucleation_rate_s"][node]),
                "effective_opening_stress_Pa": float(rates["effective_opening_stress_Pa"][node]),
                "emission_drive_equivalent_stress_Pa": float(
                    payload["diagnostics"]["emission_drive_equivalent_stress_Pa"][node]
                ),
                "backstress_Pa": float(rates["backstress_Pa"][node]),
                "state_shift_eV": float(rates["state_shift_eV"][node]),
            })
    args_cli.out.mkdir(parents=True, exist_ok=True)
    table = args_cli.out / f"{args_cli.table_prefix}_matched_checkpoints.csv"
    with table.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=tuple(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    extrema = {}
    for label in dict.fromkeys(row["checkpoint"] for row in rows):
        subset = [row for row in rows if row["checkpoint"] == label and row["available"]]
        extrema[label] = {
            "minimum_remaining_action_site": min(subset, key=lambda r: (r["remaining_threshold_action"], r["site_id"])),
            "maximum_rate_site": max(subset, key=lambda r: (r["birth_rate_per_cycle"], -r["site_id"])),
            "minimum_wait_site": min(subset, key=lambda r: (r["instantaneous_wait_cycles"], r["site_id"])),
            "maximum_cumulative_action_site": max(subset, key=lambda r: (r["cumulative_action"], -r["site_id"])),
        }
    (args_cli.out / f"{args_cli.table_prefix}_extrema.json").write_text(
        json.dumps(extrema, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(extrema, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
