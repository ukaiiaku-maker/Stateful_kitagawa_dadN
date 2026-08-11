import math

import numpy as np

from arrhenius_fracture.config import ElasticProperties
from arrhenius_fracture.sn_arrhenius_chain import build_chain_from_namespace
from arrhenius_fracture import sn_pd2d_stateful_v9_transactional as driver
from scripts.run_v9_four_class_stateful_pd import configure_four_class


SOURCE = "/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1"


def test_peak_delivery_is_audited_aggregate_m_times_per_source_lambda():
    args = driver.build_parser().parse_args([])
    configure_four_class(args, "Peak", SOURCE)
    material = ElasticProperties(E=args.E_GPa * 1e9, nu=args.nu, b=args.b_m, Tm=args.Tm_K)
    chain = build_chain_from_namespace(args, material.b)
    stress = np.array([[4.0e9, 5.0e9]])
    rho = np.array([5.0e12, 5.0e12])
    per_source = chain.rates(stress, rho, 300.0)["lambda_emit"]
    aggregate = driver.phase_resolved_delivery_rate(args, chain, stress, rho, 300.0)
    np.testing.assert_allclose(
        aggregate, args.delivery_source_multiplicity * per_source, rtol=2e-15
    )
    log_aggregate = driver.phase_resolved_delivery_log_rate(
        args, chain, stress, rho, 300.0
    )
    np.testing.assert_allclose(np.exp(log_aggregate), aggregate, rtol=2e-12)
    row = args.four_class_registry_audit["exact_registry_row"]
    expected = float(row["rho_source0_m2"]) * float(row["reference_source_area_um2"]) * 1e-12
    assert math.isclose(args.delivery_source_multiplicity, expected, rel_tol=1e-15)
    assert args.four_class_transfer_contract["delivery_bridge"][
        "candidate_site_density_applied_to_delivery"
    ] is False


def test_non_four_class_delivery_retains_unit_multiplicity():
    args = driver.build_parser().parse_args([])
    args.delivery_source = "emission"
    material = ElasticProperties(E=args.E_GPa * 1e9, nu=args.nu, b=args.b_m, Tm=args.Tm_K)
    chain = build_chain_from_namespace(args, material.b)
    stress = np.array([[1.0e9]])
    rho = np.array([5.0e12])
    expected = chain.rates(stress, rho, args.T)["lambda_emit"]
    actual = driver.phase_resolved_delivery_rate(args, chain, stress, rho, args.T)
    np.testing.assert_allclose(actual, expected, rtol=2e-15)
