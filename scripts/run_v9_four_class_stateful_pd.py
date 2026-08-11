#!/usr/bin/env python3
"""Sequential 300 K four-class Stateful-PD stable-spatial-birth driver."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from arrhenius_fracture.sn_pd2d_stateful_v9_transactional import (
    apply_resolution_profile,
    build_parser as build_legacy_parser,
    run_sweep,
)
from arrhenius_fracture.v9_four_class_registry import select_canonical_option
from arrhenius_fracture.v9_four_class_signed_mpz import load_audited_modules


MODEL_ID = "v9_four_class_stateful_PD_shared_root_MPZ_marked_cleavage_v1"
OPTIONS = {
    "Peak": "v913_paper_peak01_0242980_persistent_sites",
    "DBTT": "v913_paper_dbtt01_0202500_persistent_sites",
    "weak-T": "v913_paper_weakT01_0129902_persistent_sites",
    "ceramic": "v913_paper_ceramic01_0077080_persistent_sites",
}


def _surface(barrier):
    return {
        "G00_eV": float(barrier.G00_eV),
        "gT_eV_per_K": float(barrier.gT_eV_per_K),
        "sigc0_Pa": float(barrier.sigc0_Pa),
        "sT_Pa_per_K": float(barrier.sT_Pa_per_K),
        "alpha": float(barrier.alpha),
        "exponent": float(barrier.exponent),
        "floor_fraction": float(barrier.floor_fraction),
        "floor_min_eV": float(barrier.floor_min_eV),
        "floor_max_fraction": float(barrier.floor_max_fraction),
        "Tref_K": float(barrier.Tref_K),
        "rate_prefactor": float(barrier.attempt_frequency_s),
    }


def configure_four_class(args, material_class, source_root):
    option_id = OPTIONS[material_class]
    selected, registry_audit = select_canonical_option(option_id, source_root)
    modules, source_hashes, _paths = load_audited_modules(source_root)
    manifest = modules["reduced_campaign_v1024"].manifest_from_row(selected.row)

    args.T = 300.0
    args.cases = ["shielded"]
    args.fatigue_model = "custom"
    args.fatigue_model_preset_applied = False
    args.fatigue_endpoint = "stable_spatial_crack_birth"
    args.pd_image_policy = args.pd_image_policy
    args.cleavage_clock_model = MODEL_ID
    args.four_class_source_root = str(Path(source_root).resolve())
    args.align_blocks_to_birth_clock = False
    args.birth_scale = 0.0

    # Historical, intentionally blunt S-N geometry. These are fixed contract
    # values, not stress/life tuning parameters.
    args.Lx = 2.0e-3
    args.Ly = 4.0e-3
    args.feature_type = "ellipse"
    args.notch_depth_m = 150.0e-6
    args.notch_half_height_m = 300.0e-6
    args.notch_root_radius_m = None

    cleavage = manifest.cleavage
    args.crack_G00_eV = float(cleavage.G00_eV)
    args.crack_gT_eV_per_K = float(cleavage.gT_eV_per_K)
    args.crack_sigc0_GPa = float(cleavage.sigc0_Pa) * 1.0e-9
    args.crack_sT_GPa_per_K = float(cleavage.sT_Pa_per_K) * 1.0e-9
    args.crack_exp_a = float(cleavage.alpha)
    args.crack_exp_n = float(cleavage.exponent)
    args.crack_floor_frac = float(cleavage.floor_fraction)
    args.crack_T_mode = "audited_linear"
    args.crack_Tref_K = float(cleavage.Tref_K)
    args.nu0_crack = float(cleavage.attempt_frequency_s)
    args.S_crack_kB = 0.0

    # Local fields below are retained only for relative spatial marks and the
    # existing stabilization/healing/topology response.  They do not multiply
    # the conserved root-local canonical attempt intensity.
    args.audited_four_class_barriers = {
        "emission": _surface(manifest.emission),
        "peierls": _surface(manifest.peierls.as_surface(manifest.emission)),
        "taylor": _surface(manifest.taylor.as_surface(manifest.emission)),
    }
    # Legacy scalar delivery remains diagnostic-only in this driver mode.
    args.delivery_source = "emission"
    source_density_m2 = float(selected.row["rho_source0_m2"])
    reference_source_area_m2 = float(selected.row["reference_source_area_um2"]) * 1e-12
    args.delivery_source_multiplicity = source_density_m2 * reference_source_area_m2
    args.four_class_option_id = option_id
    args.four_class_registry_audit = registry_audit | {
        "constitutive_source_hashes": source_hashes,
    }
    args.four_class_transfer_contract = {
        "model_id": MODEL_ID,
        "exactly_transferred": [
            "cleavage_EXP_floor_surface", "emission_EXP_floor_surface_and_delivery",
            "Peierls_transport_surface", "Taylor_transport_surface",
            "attempt_frequencies",
        ],
        "parallel_spatial_PD_physics": [
            "candidate_site_population", "normalized_spatial_mark",
            "stabilization_healing",
            "directional_bond_damage", "local_PD_redistribution",
            "stable_seed_competition", "front_capture", "root_topology",
        ],
        "delivery_bridge": {
            "model": "audited_v10221_aggregate_persistent_emission",
            "per_source_rate": "lambda_emit_s-1",
            "source_density_m2": source_density_m2,
            "reference_source_area_m2": reference_source_area_m2,
            "multiplicity_per_system": args.delivery_source_multiplicity,
            "aggregate_rate": "multiplicity_per_system*lambda_emit_s-1",
            "candidate_site_density_applied_to_delivery": False,
            "qualification": "relative_spatial_mark_only_not_total_hazard",
        },
        "not_claimed_exact_signed_MPZ_parity": [
            "scalar_FEM_plastic_state", "spatial_mark_field",
            "PD_stabilization_healing", "PD_bond_growth_linkage",
            "signed_mobile_retained_populations",
        ],
        "primary_endpoint": "first_front_capture_after_stable_spatial_seed",
        "transient_embryo_is_failure": False,
        "post_endpoint_propagation": False,
        "mechanics_treatment": "stateful_PD",
    }
    return args


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--material-class", choices=tuple(OPTIONS), required=True)
    parser.add_argument("--sigma-a-MPa", type=float, required=True)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=Path("runs/sn_v9_four_class_stateful_pd"))
    parser.add_argument("--cycles-max", type=float, default=1.0e8)
    parser.add_argument("--block-cycles", type=float, default=1.0e5)
    parser.add_argument("--min-block-cycles", type=float, default=1.0e-6)
    parser.add_argument("--max-blocks", type=int, default=3000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--pd-seed", type=int)
    parser.add_argument("--global-cleavage-seed", type=int, default=42017)
    parser.add_argument("--spatial-mark-seed", type=int, default=42018)
    parser.add_argument("--site-density-m2", type=float, default=5e10)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--checkpoint-every-blocks", type=int, default=25)
    parser.add_argument("--print-every", type=int, default=10)
    parser.add_argument("--pd-image-policy", choices=("none", "event_only", "selected"), default="none")
    parser.add_argument("--resolution-profile", choices=("h15", "h10", "custom"), default="h15")
    parser.add_argument("--pd-high-cycle", action="store_true")
    parser.add_argument("--pd-high-cycle-start-cycles", type=float, default=1.0e4)
    parser.add_argument("--pd-high-cycle-max-segment", type=float, default=1.0e8)
    return parser


def main(argv=None):
    cli = build_parser().parse_args(argv)
    args = build_legacy_parser().parse_args([])
    for name in (
        "cycles_max", "block_cycles", "min_block_cycles", "max_blocks", "seed",
        "pd_seed", "global_cleavage_seed", "spatial_mark_seed",
        "site_density_m2", "resume", "checkpoint_every_blocks", "print_every", "pd_image_policy",
        "resolution_profile",
        "pd_high_cycle", "pd_high_cycle_start_cycles", "pd_high_cycle_max_segment",
    ):
        setattr(args, name, getattr(cli, name))
    args.sigma_a_MPa = [float(cli.sigma_a_MPa)]
    args.out = str(cli.out / cli.material_class)
    configure_four_class(args, cli.material_class, cli.source_root)
    apply_resolution_profile(args)
    cli.out.mkdir(parents=True, exist_ok=True)
    contract = {
        "schema": "V9_FOUR_CLASS_STATEFUL_PD_CAMPAIGN_CONTRACT_1",
        "model_id": MODEL_ID,
        "temperature_K": 300.0,
        "material_class": cli.material_class,
        "option_id": args.four_class_option_id,
        "analytic_notch_root_radius_m": 600.0e-6,
        "pd_image_policy": args.pd_image_policy,
        "transfer": args.four_class_transfer_contract,
        "registry": args.four_class_registry_audit,
    }
    (cli.out / "campaign_contract.json").write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n")
    run_sweep(args)


if __name__ == "__main__":
    main()
