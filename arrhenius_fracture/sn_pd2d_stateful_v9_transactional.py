"""Transactional v9 blunt-scratch S-N initiation driver.

Active architecture
-------------------
1. Intact 2-D FEM resolves cyclic stress, Arrhenius plastic eigenstrain,
   dislocation-density evolution, residual stress, and optional ALE scratch
   reshaping.
2. A local nonlocal spring/peridynamic patch receives the FEM shell motion and
   plastic eigenstrain and re-equilibrates as cohesive bonds soften.
3. Candidate sites undergo finite-memory plastic-event delivery, independent cleavage nucleation, reversible
   embryo formation, stabilization/healing, stable-defect growth, and spatial
   bond linkage.
4. Crack formation is a graph-connectivity event in the softened/broken bond
   network.  No phase-field or AT1/AT2 module is imported or called.

This v8.3 production candidate is one-way coupled: PD bond degradation redistributes the local
patch deformation but does not yet modify the global FEM stiffness.
"""
from __future__ import annotations

import argparse
import csv
import json
import hashlib
import math
import os
import sys
from copy import deepcopy
from dataclasses import fields
from pathlib import Path

import numpy as np

from .config import EV_TO_J, KB, ElasticProperties
from .sn_arrhenius_chain import build_chain_from_namespace
from .sn_feature_geometry_v8_7 import (
    BluntNotchGeometry,
    apply_local_ale_surface_update,
    identify_feature_surface_nodes,
    local_root_radius,
    local_root_xy,
    make_blunt_edge_notch_mesh,
    rebuild_mesh_geometry,
)
from .sn_intact_fem import (
    affine_stress_control_displacements,
    cycle_stress_histories,
    plane_strain_D,
    project_gp_to_nodes,
    project_plastic_state,
    representative_plastic_cycle,
    stress_state_intact,
    surface_morphology_proposal,
)
from .stateful_peridynamics_v8_7_local_front_spacing import StatefulPDConfig
from .v9_stateful_peridynamics import V9StatefulPDPatch as StatefulPDPatch
from .v9_fem_transaction import EmbeddedFEMTransaction, FEMPhysicalState
from .v9_physical_integrator import plastic_chain_log_rates
from .v9_array_codec import AtomicArrayGenerationStore
from .v9_cached_fem import CachedIntactFEM
from .v9_pd_high_cycle import DormantPDHighCycleEngine, HighCycleConfig
from .v9_pd_high_cycle_adapter import SpatialPDDormantAdapter


MODEL_ID = "SN_2D_intact_FEM_stateful_local_peridynamics_v9_transactional_log_domain"


def _logdiffexp(log_total, log_old):
    total = np.asarray(log_total, float); old = np.asarray(log_old, float)
    out = np.full_like(total, -math.inf)
    only = np.isfinite(total) & ~np.isfinite(old); out[only] = total[only]
    valid = np.isfinite(total) & np.isfinite(old) & (total > old)
    out[valid] = total[valid] + np.log1p(-np.exp(old[valid] - total[valid]))
    return out


def evaluate_dormant_exact_cycle(*, args, shield_on, mesh, patch, pd_state, crack,
        plast_chain, cached_fem, fem_transaction, sigma_max, sigma_min,
        ep_gp, rho_gp, epsp_acc_gp, u, plastic_work, cycles, dN=1.0):
    """Authoritative private exact/macro-cycle operation extracted from the block loop.

    The same cached FEM, pre/post midpoint construction, PD equilibrium, and
    v9 finite-memory update are used. Persistent thresholds are made unreachable
    on the private state, so the operation computes action without comparing or
    consuming a physical first-passage threshold.
    """
    trial = deepcopy(pd_state)
    trial.site_birth_threshold = np.full_like(trial.site_birth_threshold, math.inf)
    stochastic_before = (
        deepcopy(patch._candidate_rng.bit_generator.state),
        deepcopy(patch._event_rng.bit_generator.state),
    )
    fem0 = FEMPhysicalState(ep_gp, rho_gp, epsp_acc_gp, u, plastic_work)
    Umax, Umin, u_zero, _, _ = cached_fem.affine(ep_gp, sigma_max, sigma_min, u)
    first = cached_fem.representative_cycle(
        ep_gp, rho_gp, Umax, Umin, args.T, args.frequency_Hz,
        args.plastic_n_phase, plast_chain, u_zero, args.k_store, args.k_dyn,
        args.rho_floor, args.rho_cap, args.max_dep_phase, args.max_rho_rel_phase,
    )
    dN = float(dN)
    if not math.isfinite(dN) or dN <= 0.0:
        raise ValueError("private dormant window dN must be positive and finite")
    proposal = fem_transaction.propose(fem0, dN, first_cycle=first)
    if proposal.normalized_error > 1.0:
        raise RuntimeError("private dormant window fails the accepted FEM transaction tolerance")

    def coupled_fields(ep, rho, epsacc, displacement):
        Uhi, Ulo, uzero, _, _ = cached_fem.affine(ep, sigma_max, sigma_min, displacement)
        eps_node, rho_node, P, Dloc = project_plastic_state(
            mesh, epsacc, rho, args.epsp_shield_scale, args.epsp_damage_scale)
        chi = args.shield_chi if shield_on else 0.0
        Gsh = args.Gshield_eV if shield_on else 0.0
        sigma_back = args.sigma_back_max_GPa * 1e9 * P
        state_shift = Gsh * P - args.Gstored_eV * Dloc
        hist = cached_fem.stress_histories(ep, Uhi, Ulo, args.hazard_n_phase, uzero)
        delivery = phase_resolved_delivery_rate(args, plast_chain, hist["seq_node"], rho_node, args.T)
        log_delivery = phase_resolved_delivery_log_rate(args, plast_chain, hist["seq_node"], rho_node, args.T)
        ep_node = project_gp_to_nodes(mesh, ep)
        _, _, _, bond_amp, point_amp = patch.solve_local_mechanics(trial, hist["u_max"], ep_node)
        return hist, delivery, log_delivery, P, state_shift, sigma_back, chi, point_amp, bond_amp

    pre = coupled_fields(ep_gp, rho_gp, epsp_acc_gp, u)
    state1 = proposal.state
    post = coupled_fields(state1.ep_gp, state1.rho_gp, state1.epsp_acc_gp, state1.u)
    sigma_mid = 0.5 * (pre[0]["sigma_node"] + post[0]["sigma_node"])
    delivery_mid = 0.5 * (pre[1] + post[1])
    log_delivery_mid = np.logaddexp(pre[2], post[2]) - math.log(2.0)
    point_amp = 0.5 * (pre[7] + post[7]); bond_amp = 0.5 * (pre[8] + post[8])
    log_birth0 = np.asarray(trial.log_birth_cumulative_hazard, float).copy()
    born0 = np.asarray(trial.born_cumulative, float).copy()
    healed0 = np.asarray(trial.healed_cumulative, float).copy()
    diagnostics = patch.update(
        trial, crack, sigma_mid, delivery_mid, args.T, args.frequency_Hz, dN,
        cycles, post[4], post[5], post[6], post[3], point_amp, bond_amp,
        log_delivery_rate_phase_global=log_delivery_mid,
    )
    if (patch._candidate_rng.bit_generator.state != stochastic_before[0]
            or patch._event_rng.bit_generator.state != stochastic_before[1]):
        raise RuntimeError("private dormant cycle consumed RNG state")
    if np.any(trial.site_status != pd_state.site_status) or np.any(trial.bond_damage != pd_state.bond_damage):
        raise RuntimeError("private dormant cycle changed discrete PD state")
    log_node_action = _logdiffexp(trial.log_birth_cumulative_hazard, log_birth0)
    ledger = {
        "plastic_work": float(state1.plastic_work_J_per_m - plastic_work),
        "born_cumulative": np.maximum(trial.born_cumulative - born0, 0.0),
        "healed_cumulative": np.maximum(trial.healed_cumulative - healed0, 0.0),
    }
    s1_mid = 0.5 * (pre[0]["s1_node"] + post[0]["s1_node"])
    sigma_tensor_mid = 0.5 * (pre[0]["sigma_node"] + post[0]["sigma_node"])
    seq_mid = 0.5 * (pre[0]["seq_node"] + post[0]["seq_node"])
    iph, inode = np.unravel_index(int(np.argmax(s1_mid)), s1_mid.shape)
    return {
        "ep_gp": np.asarray(state1.ep_gp).copy(), "rho_gp": np.asarray(state1.rho_gp).copy(),
        "epsp_acc_gp": np.asarray(state1.epsp_acc_gp).copy(), "u": np.asarray(post[0]["u_end"]).copy(),
        "log_delivery_memory": np.asarray(trial.log_delivery_memory).copy(),
        "available": np.asarray(trial.available).copy(), "embryo": np.asarray(trial.embryo).copy(),
        "stable": np.asarray(trial.stable).copy(), "inactive": np.asarray(trial.inactive).copy(),
        "completion": np.asarray(trial.completion).copy(),
        "log_birth_action": log_node_action, "ledger_increments": ledger,
        "phase": np.linspace(0.0, 1.0, args.hazard_n_phase, endpoint=False),
        "phase_log_birth_rate": np.empty((0, len(log_node_action))),
        "diagnostic_fields": {
            "sigma_mid_global_voigt_Pa": sigma_tensor_mid,
            "equivalent_stress_mid_global_Pa": seq_mid,
        },
        "diagnostics": {"max_effective_stress_Pa": diagnostics.max_effective_stress_Pa,
                        "max_delivery_memory": diagnostics.max_delivery_memory,
                        "max_completion": diagnostics.max_completion,
                        "fem_embedded_error": proposal.normalized_error,
                        "private_window_cycles": dN,
                        "fem_sigma1_cycle_max_Pa": float(s1_mid[iph, inode]),
                        "fem_sigma1_hotspot_phase_index": int(iph),
                        "fem_sigma1_hotspot_global_node": int(inode),
                        "fem_sigma1_hotspot_x_m": float(mesh.nodes[inode, 0]),
                        "fem_sigma1_hotspot_y_m": float(mesh.nodes[inode, 1]),
                        "fem_hotspot_stress_voigt_xx_yy_xy_Pa": np.asarray(
                            sigma_tensor_mid[iph, :, inode], float
                        ),
                        "fem_hotspot_equivalent_stress_Pa": float(seq_mid[iph, inode]),
                        "local_sigma1_over_remote_sigma_max": float(
                            s1_mid[iph, inode] / max(abs(sigma_max), 1e-300)
                        ),
                        "max_raw_cleavage_rate_s": float(
                            np.max(diagnostics.max_raw_cleavage_rate_s)
                            if np.ndim(diagnostics.max_raw_cleavage_rate_s) else
                            diagnostics.max_raw_cleavage_rate_s
                        ) if hasattr(diagnostics, "max_raw_cleavage_rate_s") else float(
                            np.max(patch.last_rates["nucleation_rate_s"])
                        ),
                        "max_actual_site_birth_rate_per_cycle": float(
                            np.exp(np.max(log_node_action))
                        ),
                        "emission_drive_equivalent_stress_Pa": np.max(
                            0.5 * (pre[0]["seq_node"] + post[0]["seq_node"]), axis=0
                        )[patch.global_nodes]},
        "transition_signature": "dormant_fixed_topology",
    }


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _active_source_sha256() -> dict[str, str]:
    driver_file = Path(__file__)
    modules = {
        "pd_module": StatefulPDPatch.__module__,
        "pd_base_module": StatefulPDPatch.__mro__[1].__module__,
        "fem_transaction": EmbeddedFEMTransaction.__module__,
        "physical_integrator": plastic_chain_log_rates.__module__,
        "cached_fem": CachedIntactFEM.__module__,
        "array_codec": AtomicArrayGenerationStore.__module__,
        "pd_high_cycle_engine": DormantPDHighCycleEngine.__module__,
        "pd_high_cycle_adapter": SpatialPDDormantAdapter.__module__,
    }
    paths = {name: Path(getattr(sys.modules.get(module), "__file__", "")) for name, module in modules.items()}
    if not driver_file.is_file() or any(not path.is_file() for path in paths.values()):
        raise RuntimeError("cannot resolve active v9 source files for provenance audit")
    return {"driver": _sha256_file(driver_file)} | {
        name: _sha256_file(path) for name, path in paths.items()
    }


SOURCE_SHA256 = _active_source_sha256()
VERIFIED_COMPATIBLE_PREDECESSOR_SOURCES = (
    {
        # Peak 4 GPa legacy diagnostic stationary extension to N=1e14. This
        # exact generation predates only read-only audit-script improvements.
        "array_codec": "e99fc4a65d0b1343c7e945124ecd3dd69703345dcddc8b607e77a386372a8f3a",
        "cached_fem": "e5679ac0a613b0bcaefe7013c874671edc5829ac405f86738e3fac404d6a8490",
        "driver": "e69a5b087684cee66efa6113c23e1685d9609787ff239b8dc67a5c172a3bb8ae",
        "fem_transaction": "5c8c5467bf7043c4d8ccaae59ab1ad2ea4f2e043459b9cf3aa4b7023d9be9d7e",
        "pd_base_module": "ab38cbd7a9db8f5a5052d7c3b8490cc5d0470720ba2b8176bcd3f23e45c1c21b",
        "pd_high_cycle_adapter": "bf2a05c31244e50bc52ebab5b3e78b9b31bccd3bb3310b05c058afac94775837",
        "pd_high_cycle_engine": "7af93b75ae0df2594c1dd9455ed989cf9600cd31ba27c0de0c2209a15e9af2e5",
        "pd_module": "1d164d367994b8119cfc48552e89221adec26aaf5259820b32a731ecc67185e7",
        "physical_integrator": "a087d2dacdcf52497de5964f0ed9170f44f7a5a77daa15a90cc9774f3bc97fe3",
    },
    {
        # Peak 4 GPa legacy diagnostic stationary extension to N=1e12. The
        # subsequent change only made its absent multiplicity explicit as the
        # historical default 1.0 in presentation output.
        "array_codec": "e99fc4a65d0b1343c7e945124ecd3dd69703345dcddc8b607e77a386372a8f3a",
        "cached_fem": "e5679ac0a613b0bcaefe7013c874671edc5829ac405f86738e3fac404d6a8490",
        "driver": "dfa8db7cd3bbca116e5ca6ee05d4f266cc3dc4e3aeb41cecd8e08ec201ee62ee",
        "fem_transaction": "5c8c5467bf7043c4d8ccaae59ab1ad2ea4f2e043459b9cf3aa4b7023d9be9d7e",
        "pd_base_module": "ab38cbd7a9db8f5a5052d7c3b8490cc5d0470720ba2b8176bcd3f23e45c1c21b",
        "pd_high_cycle_adapter": "bf2a05c31244e50bc52ebab5b3e78b9b31bccd3bb3310b05c058afac94775837",
        "pd_high_cycle_engine": "7af93b75ae0df2594c1dd9455ed989cf9600cd31ba27c0de0c2209a15e9af2e5",
        "pd_module": "1d164d367994b8119cfc48552e89221adec26aaf5259820b32a731ecc67185e7",
        "physical_integrator": "a087d2dacdcf52497de5964f0ed9170f44f7a5a77daa15a90cc9774f3bc97fe3",
    },
    {
        # Repaired Peak 12 GPa direct topology diagnostic at N=3e6.  Its
        # physical modules match the accepted slow-seed repair; only the later
        # HC-020 engine fingerprint differs from the 4 GPa tail package.
        "array_codec": "e99fc4a65d0b1343c7e945124ecd3dd69703345dcddc8b607e77a386372a8f3a",
        "cached_fem": "e5679ac0a613b0bcaefe7013c874671edc5829ac405f86738e3fac404d6a8490",
        "driver": "147ffbac6ef75976a018580e4e040272e59aee57baa0e8cd58f46c898a283e4a",
        "fem_transaction": "5c8c5467bf7043c4d8ccaae59ab1ad2ea4f2e043459b9cf3aa4b7023d9be9d7e",
        "pd_base_module": "ab38cbd7a9db8f5a5052d7c3b8490cc5d0470720ba2b8176bcd3f23e45c1c21b",
        "pd_high_cycle_adapter": "bf2a05c31244e50bc52ebab5b3e78b9b31bccd3bb3310b05c058afac94775837",
        "pd_high_cycle_engine": "8c4537283a9302e4806f4742507b0a8bf45a4c4ac209f17ba38552749eea6fae",
        "pd_module": "1d164d367994b8119cfc48552e89221adec26aaf5259820b32a731ecc67185e7",
        "physical_integrator": "a087d2dacdcf52497de5964f0ed9170f44f7a5a77daa15a90cc9774f3bc97fe3",
    },
    {
        # Accepted 29786fc Peak diagnostic generations.  The subsequent
        # aggregate-emission bridge and terminal-audit additions do not alter
        # these preserved arrays; they are loaded only for explicit diagnostic
        # reconstruction and never silently resumed as corrected production.
        "array_codec": "e99fc4a65d0b1343c7e945124ecd3dd69703345dcddc8b607e77a386372a8f3a",
        "cached_fem": "e5679ac0a613b0bcaefe7013c874671edc5829ac405f86738e3fac404d6a8490",
        "driver": "147ffbac6ef75976a018580e4e040272e59aee57baa0e8cd58f46c898a283e4a",
        "fem_transaction": "5c8c5467bf7043c4d8ccaae59ab1ad2ea4f2e043459b9cf3aa4b7023d9be9d7e",
        "pd_base_module": "ab38cbd7a9db8f5a5052d7c3b8490cc5d0470720ba2b8176bcd3f23e45c1c21b",
        "pd_high_cycle_adapter": "bf2a05c31244e50bc52ebab5b3e78b9b31bccd3bb3310b05c058afac94775837",
        "pd_high_cycle_engine": "7af93b75ae0df2594c1dd9455ed989cf9600cd31ba27c0de0c2209a15e9af2e5",
        "pd_module": "1d164d367994b8119cfc48552e89221adec26aaf5259820b32a731ecc67185e7",
        "physical_integrator": "a087d2dacdcf52497de5964f0ed9170f44f7a5a77daa15a90cc9774f3bc97fe3",
    },
    {
        # Peak 12 GPa accepted boundary immediately before repairing the
        # controller-partition-dependent primary-seed stall counter. The
        # newly fingerprinted base PD module was already the executable base.
        "array_codec": "e99fc4a65d0b1343c7e945124ecd3dd69703345dcddc8b607e77a386372a8f3a",
        "cached_fem": "e5679ac0a613b0bcaefe7013c874671edc5829ac405f86738e3fac404d6a8490",
        "driver": "c4ac205a5d2492b24e1a8fdda3f59ca280e96ddb115a2f3d390745f8823fa243",
        "fem_transaction": "5c8c5467bf7043c4d8ccaae59ab1ad2ea4f2e043459b9cf3aa4b7023d9be9d7e",
        "pd_high_cycle_adapter": "bf2a05c31244e50bc52ebab5b3e78b9b31bccd3bb3310b05c058afac94775837",
        "pd_high_cycle_engine": "8c4537283a9302e4806f4742507b0a8bf45a4c4ac209f17ba38552749eea6fae",
        "pd_module": "1d164d367994b8119cfc48552e89221adec26aaf5259820b32a731ecc67185e7",
        "physical_integrator": "a087d2dacdcf52497de5964f0ed9170f44f7a5a77daa15a90cc9774f3bc97fe3",
    },
    {
        # Hash-verified terminal 690.443 MPa generation at
        # N=168634947.0289538, accepted for the explicit 690->840 MPa
        # state-preserving stress-step protocol.
        "array_codec": "e99fc4a65d0b1343c7e945124ecd3dd69703345dcddc8b607e77a386372a8f3a",
        "cached_fem": "e5679ac0a613b0bcaefe7013c874671edc5829ac405f86738e3fac404d6a8490",
        "driver": "858610e76be701fc0a2db5c79b3d5f9e72181f26ac5ee272115da8c7d439ea5e",
        "fem_transaction": "5c8c5467bf7043c4d8ccaae59ab1ad2ea4f2e043459b9cf3aa4b7023d9be9d7e",
        "pd_high_cycle_adapter": "c314f55d18bfa929c8345ec0c893a0c3554230f6e672188d125b416f93cc3a6e",
        "pd_high_cycle_engine": "d5d92d616a9ba8dc3726bb544070c3f1fdd7a2bdef59e49edacfe5eb2f8085cc",
        "pd_module": "1d164d367994b8119cfc48552e89221adec26aaf5259820b32a731ecc67185e7",
        "physical_integrator": "a087d2dacdcf52497de5964f0ed9170f44f7a5a77daa15a90cc9774f3bc97fe3",
    },
    {
        # Completed 690.443 MPa N=1e8 boundary before making the diagnostic
        # high-cycle mode history append-only across atomic restarts.
        "array_codec": "e99fc4a65d0b1343c7e945124ecd3dd69703345dcddc8b607e77a386372a8f3a",
        "cached_fem": "e5679ac0a613b0bcaefe7013c874671edc5829ac405f86738e3fac404d6a8490",
        "driver": "454c00588a077a3be3e9d9fe3261f643156ff7eca4b92831f5697a8e46ef7fc8",
        "fem_transaction": "5c8c5467bf7043c4d8ccaae59ab1ad2ea4f2e043459b9cf3aa4b7023d9be9d7e",
        "pd_high_cycle_adapter": "c314f55d18bfa929c8345ec0c893a0c3554230f6e672188d125b416f93cc3a6e",
        "pd_high_cycle_engine": "d5d92d616a9ba8dc3726bb544070c3f1fdd7a2bdef59e49edacfe5eb2f8085cc",
        "pd_module": "1d164d367994b8119cfc48552e89221adec26aaf5259820b32a731ecc67185e7",
        "physical_integrator": "a087d2dacdcf52497de5964f0ed9170f44f7a5a77daa15a90cc9774f3bc97fe3",
    },
    {
        # Accepted 690.443 MPa N=1.6023057e7 boundary before adding an
        # orchestration-only projective efficiency budget and fingerprinting
        # the already active high-cycle engine/adapter sources.
        "array_codec": "e99fc4a65d0b1343c7e945124ecd3dd69703345dcddc8b607e77a386372a8f3a",
        "cached_fem": "e5679ac0a613b0bcaefe7013c874671edc5829ac405f86738e3fac404d6a8490",
        "driver": "e05cc9d486d76c3520dc3963dfc5f89f4300b22be1dda62b525fd9f7f6b10de0",
        "fem_transaction": "5c8c5467bf7043c4d8ccaae59ab1ad2ea4f2e043459b9cf3aa4b7023d9be9d7e",
        "pd_module": "1d164d367994b8119cfc48552e89221adec26aaf5259820b32a731ecc67185e7",
        "physical_integrator": "a087d2dacdcf52497de5964f0ed9170f44f7a5a77daa15a90cc9774f3bc97fe3",
    },
    {
        # Historical K360 accepted boundary immediately before sub-ULP event
        # localization repair; physical arrays and stochastic state unchanged.
        "array_codec": "e99fc4a65d0b1343c7e945124ecd3dd69703345dcddc8b607e77a386372a8f3a",
        "cached_fem": "e5679ac0a613b0bcaefe7013c874671edc5829ac405f86738e3fac404d6a8490",
        "driver": "8ca9eef649b781e501d9edf439ef47f5e00683926d6ef24e0f6e77651a83d2fb",
        "fem_transaction": "5c8c5467bf7043c4d8ccaae59ab1ad2ea4f2e043459b9cf3aa4b7023d9be9d7e",
        "pd_module": "1d164d367994b8119cfc48552e89221adec26aaf5259820b32a731ecc67185e7",
        "physical_integrator": "a087d2dacdcf52497de5964f0ed9170f44f7a5a77daa15a90cc9774f3bc97fe3",
    },
    {
        # Isolated historical K360 N=1e5 accepted boundary before recording
        # rejected high-cycle qualifications in the mode history.
        "array_codec": "e99fc4a65d0b1343c7e945124ecd3dd69703345dcddc8b607e77a386372a8f3a",
        "cached_fem": "e5679ac0a613b0bcaefe7013c874671edc5829ac405f86738e3fac404d6a8490",
        "driver": "2b9cb4687aa5d2635cee3ba742391ef5d40b3486c2032b738a29694c48ae4c37",
        "fem_transaction": "5c8c5467bf7043c4d8ccaae59ab1ad2ea4f2e043459b9cf3aa4b7023d9be9d7e",
        "pd_module": "1d164d367994b8119cfc48552e89221adec26aaf5259820b32a731ecc67185e7",
        "physical_integrator": "a087d2dacdcf52497de5964f0ed9170f44f7a5a77daa15a90cc9774f3bc97fe3",
    },
    {
        # Accepted f7c8fb1 synthetic high-cycle milestone immediately before
        # extracting the side-effect-free real one-cycle callback. The frozen
        # physical modules are identical; only driver orchestration changed.
        "array_codec": "e99fc4a65d0b1343c7e945124ecd3dd69703345dcddc8b607e77a386372a8f3a",
        "cached_fem": "e5679ac0a613b0bcaefe7013c874671edc5829ac405f86738e3fac404d6a8490",
        "driver": "37352c88a763492cdafa47d51d587c52285baf86e52840fa7c2c0485a6a80ff5",
        "fem_transaction": "5c8c5467bf7043c4d8ccaae59ab1ad2ea4f2e043459b9cf3aa4b7023d9be9d7e",
        "pd_module": "1d164d367994b8119cfc48552e89221adec26aaf5259820b32a731ecc67185e7",
        "physical_integrator": "a087d2dacdcf52497de5964f0ed9170f44f7a5a77daa15a90cc9774f3bc97fe3",
    },
    {
        # S002 accepted boundary before lowering the exact periodic-map
        # activation threshold from 1024 to 32 cycles.
        "array_codec": "e99fc4a65d0b1343c7e945124ecd3dd69703345dcddc8b607e77a386372a8f3a",
        "cached_fem": "e5679ac0a613b0bcaefe7013c874671edc5829ac405f86738e3fac404d6a8490",
        "driver": "89ccd8288c3a31bd1421cef465b9ee7e6e5829b6f728e35a987cdcebde2594b1",
        "fem_transaction": "5c8c5467bf7043c4d8ccaae59ab1ad2ea4f2e043459b9cf3aa4b7023d9be9d7e",
        "pd_module": "3eeb5707625062d16d4325beab5392638ccff8f96f6c3465039229ad16ce14b0",
        "physical_integrator": "a087d2dacdcf52497de5964f0ed9170f44f7a5a77daa15a90cc9774f3bc97fe3",
    },
    {
        # Real K360 accepted boundary with 1e-10 absolute strain control;
        # 1e-8 successor was separately checked against 1e-9.
        "array_codec": "e99fc4a65d0b1343c7e945124ecd3dd69703345dcddc8b607e77a386372a8f3a",
        "cached_fem": "e5679ac0a613b0bcaefe7013c874671edc5829ac405f86738e3fac404d6a8490",
        "driver": "e8b95f80249c6ff4deba9ab63fc4b701aec657677cc9bf94486fcf6945a781cf",
        "fem_transaction": "be5c8f52a4147e140ce5f803bf4bf6745468d299baf40085a6dfb88e2b02e4e2",
        "pd_module": "3eeb5707625062d16d4325beab5392638ccff8f96f6c3465039229ad16ce14b0",
        "physical_integrator": "a087d2dacdcf52497de5964f0ed9170f44f7a5a77daa15a90cc9774f3bc97fe3",
    },
    {
        # Real K360 accepted boundary with 1e-12 absolute strain control;
        # 1e-10 successor was separately checked against 1e-11.
        "array_codec": "e99fc4a65d0b1343c7e945124ecd3dd69703345dcddc8b607e77a386372a8f3a",
        "cached_fem": "e5679ac0a613b0bcaefe7013c874671edc5829ac405f86738e3fac404d6a8490",
        "driver": "f38bf43c8370822d91f41946f5b9ceed4120f014a6bf12ff37449a42fb4a5199",
        "fem_transaction": "f342bbd3d20f77be6d40901aca02ac0c695ffdb74c99fdfda44cf2809f38a98e",
        "pd_module": "3eeb5707625062d16d4325beab5392638ccff8f96f6c3465039229ad16ce14b0",
        "physical_integrator": "a087d2dacdcf52497de5964f0ed9170f44f7a5a77daa15a90cc9774f3bc97fe3",
    },
    {
        # Diagnostic-only component-error source preceding the validated
        # absolute strain tolerance correction in v9_fem_transaction.py.
        "driver": "7b97d38560021600c1059b03a36c424e65269d6b3f9b9cb6c99e7baa0a278296",
        "pd_module": "3eeb5707625062d16d4325beab5392638ccff8f96f6c3465039229ad16ce14b0",
    },
    {
        # c8ab0dc: adds persistent controller scheduling; diagnostic-only
        # successor adds componentwise embedded error reporting.
        "driver": "5f035212efdfeb2cb31b783cb0679360ea2467062a1fa18a2dd0f85a335cf018",
        "pd_module": "3eeb5707625062d16d4325beab5392638ccff8f96f6c3465039229ad16ce14b0",
    },
    {
        # ebb0bd7/11a64e2: identical physics and accepted state; predecessor
        # recomputed the embedded candidate from the broad physical bound.
        "driver": "82029215d84b9543d41583003ddade17054afed8145028434731416fd8d6bd03",
        "pd_module": "3eeb5707625062d16d4325beab5392638ccff8f96f6c3465039229ad16ce14b0",
    },
    {
        "driver": "cf73e1b71c50730ee8c026c2de8894109d23730a4cca26ba3c92f5effeabee01",
        "pd_module": "3eeb5707625062d16d4325beab5392638ccff8f96f6c3465039229ad16ce14b0",
    },
    {
        "driver": "f7678fa38a6cf18a400251effa93f8ed115688b3de6086e6d941f0ef000b685e",
        "pd_module": "3eeb5707625062d16d4325beab5392638ccff8f96f6c3465039229ad16ce14b0",
    },
)


class ScratchExpFloorBarrier:
    """Cleavage EXP-floor barrier matching the production six-case closure."""

    def __init__(self, args):
        self.G00_eV = float(args.crack_G00_eV)
        self.sigc0_Pa = float(args.crack_sigc0_GPa) * 1.0e9
        self.a = float(args.crack_exp_a)
        self.n = float(args.crack_exp_n)
        self.floor_frac = float(args.crack_floor_frac)
        self.T_mode = str(args.crack_T_mode)
        self.Tref_K = float(args.crack_Tref_K)
        self.mu_dlnmu_dT_per_K = float(args.crack_mu_dlnmu_dT_per_K)
        self.G0_mu_power = float(args.crack_G0_mu_power)
        self.sigc_mu_power = float(args.crack_sigc_mu_power)
        self.S_crack_kB = float(args.S_crack_kB)
        self.gT_eV_per_K = getattr(args, "crack_gT_eV_per_K", None)
        self.sT_Pa_per_K = (
            None if getattr(args, "crack_sT_GPa_per_K", None) is None
            else float(args.crack_sT_GPa_per_K) * 1.0e9
        )
        self.rate_prefactor = float(args.nu0_crack)

    def _G0_sigc(self, T_K):
        if self.T_mode == "mu_scale":
            ratio = max(
                0.05,
                1.0 + self.mu_dlnmu_dT_per_K * (float(T_K) - self.Tref_K),
            )
            G0 = self.G00_eV * ratio ** self.G0_mu_power
            sigc = self.sigc0_Pa * ratio ** self.sigc_mu_power
        elif self.T_mode == "audited_linear":
            if self.gT_eV_per_K is None or self.sT_Pa_per_K is None:
                raise RuntimeError("audited_linear cleavage requires exact temperature slopes")
            dT = float(T_K) - self.Tref_K
            G0 = self.G00_eV + float(self.gT_eV_per_K) * dT
            sigc = self.sigc0_Pa + float(self.sT_Pa_per_K) * dT
        else:
            G0 = self.G00_eV
            sigc = self.sigc0_Pa
        return max(G0, 1e-12), max(sigc, 1e6)

    def deltaG_eV(self, sigma_Pa, T_K):
        G0, sigc = self._G0_sigc(T_K)
        floor = min(0.95 * G0, max(1e-4, self.floor_frac * G0))
        x = np.maximum(np.asarray(sigma_Pa, float), 0.0) / sigc
        G = floor + (G0 - floor) * np.exp(-self.a * np.power(x, self.n))
        G -= (float(T_K) - 300.0) * self.S_crack_kB * (KB / EV_TO_J)
        return np.maximum(G, 1e-12)


REPRESENTATIVE_FATIGUE_MODELS = {
    # G00_eV, sigc0_GPa, a, n, floor, emit_S, peierls_S, taylor_S
    "FCC_like_case29": (1.0, 2.5, 0.70, 0.60, 0.020, 0.375, 0.001875, 0.0075),
    "shifted_ductile_case64": (1.0, 3.0, 0.70, 0.60, 0.010, 0.375, 0.001875, 0.0075),
    "steep_cleavage_case35": (1.0, 2.5, 0.70, 1.00, 0.020, 0.0, 0.0, 0.0),
    "slow_threshold_case101": (1.0, 3.5, 0.70, 0.60, 0.020, 0.0, 0.0, 0.0),
    "higher_barrier_case171": (1.1, 2.5, 0.70, 0.60, 0.005, 0.0, 0.0, 0.0),
    "plastic_shielded_case64_M1": (1.0, 3.0, 0.70, 0.60, 0.010, 0.75, 0.00375, 0.015),
}


def apply_representative_fatigue_model(args):
    """Apply one of the six production fatigue barrier parameterizations.

    The cleavage EXP-floor surface and the emission/Peierls/Taylor entropy
    scales are copied from the prior six-case production runners.  The
    scratch-specific embryo/growth/linkage closure remains a spatial
    translation of that model and is reported separately in the output.
    """
    label = str(args.fatigue_model)
    if label == "custom":
        args.fatigue_model_preset_applied = False
        return args
    vals = REPRESENTATIVE_FATIGUE_MODELS[label]
    (
        args.crack_G00_eV,
        args.crack_sigc0_GPa,
        args.crack_exp_a,
        args.crack_exp_n,
        args.crack_floor_frac,
        args.emit_entropy_scale,
        args.peierls_entropy_scale,
        args.taylor_entropy_scale,
    ) = vals
    # The production six-case family used these common energy scales.
    args.emit_energy_scale = 0.75
    args.peierls_energy_scale = 0.00375
    args.taylor_energy_scale = 0.015
    args.crack_T_mode = "mu_scale"
    args.crack_Tref_K = 481.33
    args.crack_mu_dlnmu_dT_per_K = -1.5e-4
    args.crack_G0_mu_power = 1.0
    args.crack_sigc_mu_power = 1.0
    args.fatigue_model_preset_applied = True
    return args


def build_crack_barrier(args):
    return ScratchExpFloorBarrier(args)


def phase_resolved_delivery_rate(args, plast_chain, seq_node_phase, rho_node, T_K):
    """Return the plastic-event delivery rate on every cycle phase and node.

    The default ``completed_flow`` source uses the completed Arrhenius
    emission--Peierls--Taylor chain.  Alternative sources are retained for
    mechanism sensitivity, but the cleavage hazard is never reused as the
    delivery process.
    """
    seq = np.maximum(np.asarray(seq_node_phase, float), 0.0)
    rho = np.maximum(np.asarray(rho_node, float), 1.0e6)
    rr = plast_chain.rates(seq, rho, float(T_K))
    source = str(args.delivery_source)
    if source == "completed_flow":
        base = np.asarray(rr["lambda_flow"], float)
    elif source == "emission":
        base = np.asarray(rr["lambda_emit"], float)
    elif source == "plastic_strain_rate":
        event_strain = max(float(args.delivery_event_strain), 1e-30)
        base = np.asarray(rr["dot_ep"], float) / event_strain
    elif source == "weighted_mechanisms":
        base = (
            float(args.delivery_weight_emit) * np.asarray(rr["lambda_emit"], float)
            + float(args.delivery_weight_peierls) * np.asarray(rr["lambda_peierls"], float)
            + float(args.delivery_weight_taylor) * np.asarray(rr["lambda_taylor"], float)
        )
    else:
        raise ValueError(f"unknown delivery source: {source}")
    # ``base`` is a per-source Arrhenius rate.  Canonical four-class rows use
    # the audited v10.2.21 aggregate persistent-emission observable
    # M*lambda.  The multiplicity is one local source-zone measure; PD
    # candidate density is a separate first-passage population and must not be
    # folded into delivery a second time.
    multiplicity = max(float(getattr(args, "delivery_source_multiplicity", 1.0)), 0.0)
    rate = (max(float(args.delivery_scale), 0.0) * multiplicity
            * np.maximum(base, 0.0))
    cap = float(args.delivery_rate_cap_s)
    if np.isfinite(cap) and cap > 0.0:
        rate = np.minimum(rate, cap)
    return rate


def phase_resolved_delivery_log_rate(args, plast_chain, seq_node_phase, rho_node, T_K):
    """V9 companion to the legacy presentation-rate array."""
    logs = plastic_chain_log_rates(
        plast_chain, np.maximum(np.asarray(seq_node_phase, float), 0.0),
        np.maximum(np.asarray(rho_node, float), 1.0e6), float(T_K),
    )
    source = str(args.delivery_source)
    if source == "completed_flow":
        out = logs["lambda_flow"]
    elif source == "emission":
        out = logs["lambda_emit"]
    elif source == "plastic_strain_rate":
        out = logs["dot_ep"] - math.log(max(float(args.delivery_event_strain), 1e-300))
    elif source == "weighted_mechanisms":
        terms = []
        for weight, name in (
            (args.delivery_weight_emit, "lambda_emit"),
            (args.delivery_weight_peierls, "lambda_peierls"),
            (args.delivery_weight_taylor, "lambda_taylor"),
        ):
            if float(weight) > 0.0:
                terms.append(logs[name] + math.log(float(weight)))
        out = np.logaddexp.reduce(np.stack(terms), axis=0) if terms else np.full_like(logs["lambda_flow"], -math.inf)
    else:
        raise ValueError(f"unknown delivery source: {source}")
    scale = (float(args.delivery_scale)
             * float(getattr(args, "delivery_source_multiplicity", 1.0)))
    if scale <= 0.0:
        return np.full_like(out, -math.inf)
    out = out + math.log(scale)
    cap = float(args.delivery_rate_cap_s)
    if np.isfinite(cap) and cap > 0.0:
        out = np.minimum(out, math.log(cap))
    return out


def apply_resolution_profile(args):
    """Apply versioned local-resolution defaults without hiding overrides.

    h15 targets an effective scratch-root spacing near 15 micrometres and uses
    delta=45 micrometres. h10 is the convergence profile with delta=30
    micrometres. The mesh generator's input spacing is smaller than the
    measured hbar_tip because the final Delaunay connectivity also contains the
    graded background mesh.
    """
    profile = str(args.resolution_profile)
    if profile == "custom":
        return args
    if profile == "h15":
        args.nx, args.ny = 36, 72
        args.root_h_fine = 12e-6
        args.pd_horizon_m = 45e-6
        args.pd_patch_radius_m = 0.50e-3
        args.pd_boundary_shell_m = 90e-6
        args.root_seed_radius_m = 60e-6
    elif profile == "h10":
        args.nx, args.ny = 48, 96
        args.root_h_fine = 5e-6
        args.pd_horizon_m = 30e-6
        args.pd_patch_radius_m = 0.50e-3
        args.pd_boundary_shell_m = 60e-6
        args.root_seed_radius_m = 45e-6
    else:
        raise ValueError(f"unknown resolution profile: {profile}")
    return args


def _write_csv(path: Path, rows):
    if not rows:
        return
    keys = sorted(set().union(*(r.keys() for r in rows)))
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)



_CHECKPOINT_VERSION = 9
_CHECKPOINT_EXCLUDED_ARGS = {
    "out", "resume", "skip_existing", "checkpoint_every_blocks",
    "checkpoint_path", "snapshot_every", "print_every", "max_blocks",
    "cycles_max",
    "pd_high_cycle", "pd_high_cycle_max_segment", "pd_high_cycle_checkpoint_decades",
    "pd_high_cycle_start_cycles",
    "stress_step_source_checkpoint", "stress_step_source_generation",
    "stress_step_source_sigma_a_MPa", "protocol_label",
}


def _json_safe(value):
    if isinstance(value, dict):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def _strict_json_safe(value):
    value = _json_safe(value)
    if isinstance(value, dict):
        return {key: _strict_json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_strict_json_safe(item) for item in value]
    if isinstance(value, float) and not np.isfinite(value):
        return "Infinity" if value > 0.0 else "-Infinity" if value < 0.0 else "NaN"
    return value


def _restore_nonfinite_tags(value):
    if isinstance(value, dict):
        return {key: _restore_nonfinite_tags(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_restore_nonfinite_tags(item) for item in value]
    if value == "Infinity":
        return float("inf")
    if value == "-Infinity":
        return float("-inf")
    if value == "NaN":
        return float("nan")
    return value


def _checkpoint_signature(args, case_name, sigma_a_MPa):
    payload = {
        k: _json_safe(v)
        for k, v in vars(args).items()
        if k not in _CHECKPOINT_EXCLUDED_ARGS
    }
    payload["case"] = str(case_name)
    payload["sigma_a_MPa"] = float(sigma_a_MPa)
    return payload


def _save_case_checkpoint(
    path,
    *,
    args,
    case_name,
    sigma_a_MPa,
    next_block,
    controller_next_block_cycles,
    cycles,
    Wp_total,
    mesh,
    root_xy,
    ep_gp,
    rho_gp,
    epsp_acc_gp,
    u,
    last_residual,
    pd_state,
    patch,
    rows,
):
    """Atomically save all evolving state needed for exact continuation."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    arrays = {
        "mesh_nodes": np.asarray(mesh.nodes, float),
        "root_xy": np.asarray(root_xy, float),
        "ep_gp": np.asarray(ep_gp, float),
        "rho_gp": np.asarray(rho_gp, float),
        "epsp_acc_gp": np.asarray(epsp_acc_gp, float),
        "u": np.asarray(u, float),
        "last_residual": np.asarray(last_residual, float),
    }
    pd_scalars = {}
    for field in fields(pd_state):
        value = getattr(pd_state, field.name)
        if isinstance(value, np.ndarray):
            arrays[f"pd__{field.name}"] = value
        else:
            pd_scalars[field.name] = _json_safe(value)
    for name, value in patch.v9_extra_state_arrays(pd_state).items():
        arrays[f"pd_v9__{name}"] = value
    metadata = {
        "checkpoint_version": _CHECKPOINT_VERSION,
        "model_id": MODEL_ID,
        "source_sha256": SOURCE_SHA256,
        "signature": _checkpoint_signature(args, case_name, sigma_a_MPa),
        "next_block": int(next_block),
        "controller_next_block_cycles": float(controller_next_block_cycles),
        "cycles": float(cycles),
        "Wp_total": float(Wp_total),
        "pd_scalars": pd_scalars,
        "candidate_rng_state": _json_safe(patch._candidate_rng.bit_generator.state),
        "event_rng_state": _json_safe(patch._event_rng.bit_generator.state),
        "rows": _json_safe(rows),
    }
    arrays["metadata_json"] = np.asarray(json.dumps(metadata, allow_nan=True))
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("wb") as f:
        np.savez_compressed(f, **arrays)
    os.replace(tmp, path)
    generation_arrays = {
        name: np.asarray(value) for name, value in arrays.items()
        if name != "metadata_json"
    }
    AtomicArrayGenerationStore(path.parent / "v9_generations").write(
        generation_arrays,
        _strict_json_safe(metadata) | {
            "legacy_checkpoint_sha256": _sha256_file(path),
            "checkpoint_role": "accepted_physical_boundary",
        },
        {
            "status": "restartable",
            "cycles": float(cycles),
            "next_block": int(next_block),
            "model_id": MODEL_ID,
        },
    )


def _load_case_checkpoint(
    path,
    *,
    args,
    case_name,
    sigma_a_MPa,
    mesh,
    patch,
    generation=None,
):
    """Restore the newest hash-verified atomic v9 array generation."""
    path = Path(path)
    store = AtomicArrayGenerationStore(path.parent / "v9_generations")
    data, encoded_metadata, summary, _ = store.load(generation=generation)
    metadata = _restore_nonfinite_tags(encoded_metadata)
    if metadata.get("checkpoint_version") != _CHECKPOINT_VERSION:
        raise RuntimeError("unsupported STATEFUL_PD_V9 checkpoint version")
    if metadata.get("model_id") != MODEL_ID:
        raise RuntimeError("checkpoint model identifier does not match v9")
    checkpoint_sources = metadata.get("source_sha256")
    if (
        checkpoint_sources != SOURCE_SHA256
        and checkpoint_sources not in VERIFIED_COMPATIBLE_PREDECESSOR_SOURCES
    ):
        raise RuntimeError("checkpoint source hashes do not match active or verified-compatible v9 code")
    expected = _checkpoint_signature(args, case_name, sigma_a_MPa)
    if metadata.get("signature") != expected:
        raise RuntimeError("checkpoint arguments do not match this case")
    if float(summary.get("cycles", -1.0)) != float(metadata["cycles"]):
        raise RuntimeError("generation summary/capsule cycle mismatch")
    mesh.nodes[:] = np.asarray(data["mesh_nodes"], float)
    root_xy = np.asarray(data["root_xy"], float).copy()
    rebuild_mesh_geometry(mesh, root_xy)
    patch.update_geometry(mesh, root_xy)
    state = patch.initial_state()
    for field in fields(state):
        key = f"pd__{field.name}"
        if key in data:
            setattr(state, field.name, np.asarray(data[key]).copy())
        elif field.name in metadata["pd_scalars"]:
            setattr(state, field.name, metadata["pd_scalars"][field.name])
    for name in (
        "site_transition_threshold", "site_transition_cumulative_hazard",
        "site_transition_outcome_uniform", "log_delivery_memory",
        "log_birth_cumulative_hazard",
    ):
        key = f"pd_v9__{name}"
        if key not in data:
            raise RuntimeError(f"v9 checkpoint missing persistent transition state: {name}")
        setattr(state, name, np.asarray(data[key]).copy())
    patch._candidate_rng.bit_generator.state = metadata["candidate_rng_state"]
    patch._event_rng.bit_generator.state = metadata["event_rng_state"]
    restored = {
        "next_block": int(metadata["next_block"]), "cycles": float(metadata["cycles"]),
        "controller_next_block_cycles": float(metadata.get("controller_next_block_cycles", 0.0)),
        "Wp_total": float(metadata["Wp_total"]), "root_xy": root_xy,
        "ep_gp": np.asarray(data["ep_gp"], float).copy(),
        "rho_gp": np.asarray(data["rho_gp"], float).copy(),
        "epsp_acc_gp": np.asarray(data["epsp_acc_gp"], float).copy(),
        "u": np.asarray(data["u"], float).copy(),
        "last_residual": np.asarray(data["last_residual"], float).copy(),
        "pd_state": state, "rows": list(metadata.get("rows", [])),
    }
    return restored


def _plot_fem_fields(mesh, fields, out_png, root_xy, feature_nodes=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri

    n = len(fields)
    fig, axes = plt.subplots(1, n, figsize=(4.3 * n, 4.2), constrained_layout=True)
    axes = np.atleast_1d(axes)
    tri = mtri.Triangulation(mesh.nodes[:, 0] * 1e3, mesh.nodes[:, 1] * 1e3, mesh.elems)
    for ax, (title, data) in zip(axes, fields):
        pc = ax.tripcolor(tri, np.asarray(data), shading="gouraud")
        if feature_nodes is not None:
            q = mesh.nodes[np.asarray(feature_nodes, int)] * 1e3
            ax.plot(q[:, 0], q[:, 1], "k-", lw=0.8)
        ax.plot(root_xy[0] * 1e3, root_xy[1] * 1e3, "rx", ms=7)
        ax.set_aspect("equal")
        ax.set_title(title)
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        fig.colorbar(pc, ax=ax, fraction=0.046, pad=0.04)
    fig.savefig(out_png, dpi=220)
    plt.close(fig)


def _geometry_resolution_audit(mesh, feature_nodes, point_spacing_m, initial_min_area, args):
    """Audit whether the evolved scratch remains resolved by the current mesh."""
    radius = float(local_root_radius(mesh, feature_nodes))
    spacing = max(float(point_spacing_m), 1e-30)
    ratio = radius / spacing if np.isfinite(radius) else float("nan")
    min_area = float(np.min(mesh.area_e)) if len(mesh.area_e) else float("nan")
    area_fraction = min_area / max(float(initial_min_area), 1e-300)
    finite = bool(np.isfinite(radius) and np.isfinite(ratio) and np.isfinite(area_fraction))
    radius_pass = finite and ratio >= float(args.geometry_min_valid_radius_spacing)
    area_pass = finite and area_fraction >= float(args.geometry_min_global_area_fraction)
    reasons = []
    if not finite:
        reasons.append("nonfinite_geometry")
    if finite and not radius_pass:
        reasons.append("root_radius_underresolved")
    if finite and not area_pass:
        reasons.append("cumulative_element_area_loss")
    return {
        "root_radius_m": radius,
        "root_radius_over_spacing": ratio,
        "minimum_element_area_m2": min_area,
        "minimum_element_area_over_initial": area_fraction,
        "radius_pass": bool(radius_pass),
        "area_pass": bool(area_pass),
        "pass": bool(finite and radius_pass and area_pass),
        "failure_reasons": reasons,
    }


def _geometry_taper_factor(audit, args):
    """Smoothly reduce ALE motion as the root approaches the resolution floor."""
    if args.geometry_limit_mode != "freeze":
        return 1.0
    ratio = float(audit.get("root_radius_over_spacing", float("nan")))
    if not np.isfinite(ratio):
        return 0.0
    start = float(args.geometry_saturation_start_radius_spacing)
    floor = float(args.geometry_min_valid_radius_spacing)
    if start <= floor:
        return 1.0 if ratio > floor else 0.0
    return float(np.clip((ratio - floor) / (start - floor), 0.0, 1.0))


def _pd_config_from_args(args):
    return StatefulPDConfig(
        patch_radius_m=args.pd_patch_radius_m,
        horizon_m=args.pd_horizon_m,
        boundary_shell_m=args.pd_boundary_shell_m,
        residual_bond_stiffness=args.pd_residual_stiffness,
        initiation_radius_m=args.pd_initiation_radius_m,
        initiation_taper_m=args.pd_initiation_taper_m,
        initiation_back_extent_m=args.pd_initiation_back_extent_m,
        site_density_m2=args.site_density_m2,
        delivery_hit_count=args.delivery_hit_count,
        delivery_memory_s=args.delivery_memory_s,
        birth_scale=args.birth_scale,
        nu_stabilize_s=args.nu_stabilize_s,
        nu_heal_s=args.nu_heal_s,
        stabilize_stress_Pa=args.stabilize_stress_GPa * 1e9,
        stabilize_width_Pa=args.stabilize_width_GPa * 1e9,
        stabilize_plastic_gain=args.stabilize_plastic_gain,
        heal_return_fraction=args.heal_return_fraction,
        nu_grow_s=args.nu_grow_s,
        grow_stress_Pa=args.grow_stress_GPa * 1e9,
        grow_width_Pa=args.grow_width_GPa * 1e9,
        stable_count_scale=args.stable_count_scale,
        nu_link_s=args.nu_link_s,
        link_stress_Pa=args.link_stress_GPa * 1e9,
        link_width_Pa=args.link_width_GPa * 1e9,
        link_orientation_power=args.link_orientation_power,
        link_orientation_floor=args.link_orientation_floor,
        directional_band_horizons=args.directional_band_horizons,
        seed_influence_horizons=args.seed_influence_horizons,
        front_neighbor_spacing_factor=args.front_neighbor_spacing_factor,
        front_orientation_tolerance_deg=args.front_orientation_tolerance_deg,
        neighbor_link_gain=args.neighbor_link_gain,
        front_capture_enabled=args.front_capture_enabled,
        front_capture_min_bonds=args.front_capture_min_bonds,
        front_capture_min_length_horizons=args.front_capture_min_length_horizons,
        front_capture_min_orientation_coherence=args.front_capture_min_orientation_coherence,
        front_capture_max_surface_gap_horizons=args.front_capture_max_surface_gap_horizons,
        front_process_ahead_horizons=args.front_process_ahead_horizons,
        front_process_behind_horizons=args.front_process_behind_horizons,
        front_band_horizons=args.front_band_horizons,
        front_wake_band_horizons=args.front_wake_band_horizons,
        front_tip_seed_gain=args.front_tip_seed_gain,
        front_preferred_orientation_tolerance_deg=args.front_preferred_orientation_tolerance_deg,
        front_fallback_orientation_tolerance_deg=args.front_fallback_orientation_tolerance_deg,
        front_preferred_band_horizons=args.front_preferred_band_horizons,
        front_fallback_band_horizons=args.front_fallback_band_horizons,
        front_fallback_activity_scale=args.front_fallback_activity_scale,
        front_direction_smoothing=args.front_direction_smoothing,
        front_max_turn_deg=args.front_max_turn_deg,
        front_recent_segment_horizons=args.front_recent_segment_horizons,
        front_crack_tube_horizons=args.front_crack_tube_horizons,
        front_backbone_band_horizons=args.front_backbone_band_horizons,
        front_min_advance_spacing_factor=args.front_min_advance_spacing_factor,
        front_stall_enabled=args.front_stall_enabled,
        front_stall_patience_updates=args.front_stall_patience_updates,
        primary_seed_reselection_enabled=args.primary_seed_reselection_enabled,
        primary_seed_reselection_patience_updates=args.primary_seed_reselection_patience_updates,
        primary_seed_progress_damage_increment=args.primary_seed_progress_damage_increment,
        primary_seed_max_reselections=args.primary_seed_max_reselections,
        primary_seed_surface_score_horizons=args.primary_seed_surface_score_horizons,
        primary_seed_stress_score_weight=args.primary_seed_stress_score_weight,
        primary_seed_population_score_weight=args.primary_seed_population_score_weight,
        primary_seed_damage_score_weight=args.primary_seed_damage_score_weight,
        pre_capture_birth_scale_with_primary=args.pre_capture_birth_scale_with_primary,
        post_capture_birth_scale=args.post_capture_birth_scale,
        post_capture_stabilization_scale=args.post_capture_stabilization_scale,
        off_front_growth_scale=args.off_front_growth_scale,
        front_amplification_ahead_horizons=args.front_amplification_ahead_horizons,
        front_amplification_behind_horizons=args.front_amplification_behind_horizons,
        front_amplification_band_horizons=args.front_amplification_band_horizons,
        front_amplification_cap=args.front_amplification_cap,
        front_link_state_shift_weight=args.front_link_state_shift_weight,
        front_link_state_shift_scale_eV=args.front_link_state_shift_scale_eV,
        front_link_state_shift_z_clip=args.front_link_state_shift_z_clip,
        diffuse_abort_enabled=args.diffuse_abort_enabled,
        diffuse_abort_min_bonds=args.diffuse_abort_min_bonds,
        diffuse_abort_min_offfront_fraction=args.diffuse_abort_min_offfront_fraction,
        diffuse_abort_patience_updates=args.diffuse_abort_patience_updates,
        surface_connection_horizons=args.surface_connection_horizons,
        topology_neighbor_spacing_factor=args.topology_neighbor_spacing_factor,
        topology_orientation_tolerance_deg=args.topology_orientation_tolerance_deg,
        handoff_min_orientation_coherence=args.handoff_min_orientation_coherence,
        handoff_min_axial_coverage=args.handoff_min_axial_coverage,
        handoff_max_axial_gap_horizons=args.handoff_max_axial_gap_horizons,
        handoff_min_slenderness=args.handoff_min_slenderness,
        handoff_require_active_front=args.handoff_require_active_front,
        handoff_max_offfront_broken_fraction=args.handoff_max_offfront_broken_fraction,
        max_transition_probability=args.max_transition_probability,
        broken_damage=args.broken_damage,
        root_seed_radius_m=args.root_seed_radius_m,
        established_extent_m=args.established_extent_m,
        handoff_mode=args.handoff_mode,
        handoff_min_length_m=args.handoff_min_length_m,
        handoff_min_length_horizons=args.handoff_min_length_horizons,
        handoff_max_width_ratio=args.handoff_max_width_ratio,
        handoff_max_tip_width_horizons=args.handoff_max_tip_width_horizons,
        handoff_tip_window_horizons=args.handoff_tip_window_horizons,
        handoff_width_root_exclusion_horizons=args.handoff_width_root_exclusion_horizons,
        handoff_min_boundary_clearance_horizons=args.handoff_min_boundary_clearance_horizons,
        handoff_min_connected_bonds=args.handoff_min_connected_bonds,
        handoff_edge_geometry_factor=args.handoff_edge_geometry_factor,
        handoff_min_remote_K_MPam05=args.handoff_min_remote_K_MPam05,
        pd_amplification_cap=args.pd_amplification_cap,
        amplification_damage_scale=args.pd_amplification_damage_scale,
        random_seed=args.pd_seed if args.pd_seed is not None else args.seed,
        softening_damage=args.softening_damage,
    )


def run_case_stress(args, case_name: str, sigma_a_MPa: float):
    if args.enable_geometry_evolution:
        raise RuntimeError("v9 cached transactional driver currently requires fixed geometry")
    shield_on = case_name == "shielded"
    mat = ElasticProperties(E=args.E_GPa * 1e9, nu=args.nu, b=args.b_m, Tm=args.Tm_K)
    Dmat = plane_strain_D(mat)
    geom = BluntNotchGeometry(
        args.Lx,
        args.Ly,
        args.notch_depth_m,
        args.notch_half_height_m,
        feature_type=args.feature_type,
        root_radius_m=args.notch_root_radius_m,
        opening_angle_deg=args.notch_opening_angle_deg,
        path_refine_length_m=args.path_refine_length_m,
        path_refine_half_height_m=args.path_refine_half_height_m,
    )
    mesh, bnd, _ = make_blunt_edge_notch_mesh(
        geom,
        nx=args.nx,
        ny=args.ny,
        jitter=args.jitter,
        root_h_fine=args.root_h_fine,
        seed=args.seed,
    )
    feature_nodes = identify_feature_surface_nodes(mesh, geom)
    root_xy = local_root_xy(mesh, feature_nodes)
    root_radius0 = local_root_radius(mesh, feature_nodes)
    initial_min_area = float(np.min(mesh.area_e))
    fixed_mesh_nodes = np.unique(np.r_[bnd.top_nodes, bnd.bot_nodes])

    plast_chain = build_chain_from_namespace(args, mat.b)
    crack = build_crack_barrier(args)
    pd_cfg = _pd_config_from_args(args)
    patch = StatefulPDPatch(
        mesh, geom, root_xy, mat, pd_cfg,
        feature_surface_global_nodes=feature_nodes,
    )
    preflight_audit = patch.production_preflight_audit()
    print(
        "STATEFUL_PD_V8_7 resolution audit: "
        f"h={preflight_audit['point_spacing_m']*1e6:.2f} um, "
        f"delta/h={preflight_audit['horizon_over_spacing']:.2f}, "
        f"surface_nodes={preflight_audit['surface_node_count']}, "
        f"clearance/delta={preflight_audit['active_region_clearance_m']/max(pd_cfg.horizon_m,1e-30):.2f}"
    )
    if not preflight_audit["pass"]:
        message = (
            "STATEFUL_PD_V8_7 production preflight failed: "
            + ", ".join(preflight_audit["failure_reasons"])
        )
        if args.resolution_profile in {"h15", "h10"}:
            raise RuntimeError(message)
        print("WARNING: " + message + " (custom profile allowed for smoke/debug only)")
    sigma_max = 2.0 * sigma_a_MPa * 1e6 / max(1.0 - args.R, 1e-30)
    sigma_min = args.R * sigma_max
    patch.remote_sigma_max_Pa = float(sigma_max)
    cached_fem = CachedIntactFEM(mesh, bnd, mat, Dmat)
    fem_transaction = EmbeddedFEMTransaction(
        mesh=mesh,
        boundaries=bnd,
        material=mat,
        Dmat=Dmat,
        plastic_chain=plast_chain,
        args=args,
        sigma_max_Pa=sigma_max,
        sigma_min_Pa=sigma_min,
        cached_fem=cached_fem,
    )

    outdir = Path(args.out) / case_name / (f"sigmaA_{sigma_a_MPa:g}MPa".replace(".", "p"))
    outdir.mkdir(parents=True, exist_ok=True)
    checkpoint_path = Path(args.checkpoint_path) if args.checkpoint_path else outdir / "checkpoint_latest.npz"
    stress_step_source = str(getattr(args, "stress_step_source_checkpoint", "") or "")

    pd_state = patch.initial_state()
    u = np.zeros(mesh.ndof)
    ep_gp = np.zeros((3, mesh.ne))
    rho_gp = np.full(mesh.ne, args.rho0)
    epsp_acc_gp = np.zeros(mesh.ne)
    cycles = 0.0
    Wp_total = 0.0
    rows = []
    last_residual = np.zeros(mesh.nn)
    start_block = 0
    resumed = False
    controller_next_block_cycles = 0.0
    geometry_saturated = False
    geometry_saturation_cycles = None
    geometry_invalid_reason = None
    geometry_audit = _geometry_resolution_audit(
        mesh, feature_nodes, patch.point_spacing_m, initial_min_area, args
    )
    if stress_step_source:
        source_sigma = getattr(args, "stress_step_source_sigma_a_MPa", None)
        source_generation = str(getattr(args, "stress_step_source_generation", "") or "") or None
        if source_sigma is None:
            raise RuntimeError("stress-step restart requires stress_step_source_sigma_a_MPa")
        restored = _load_case_checkpoint(
            Path(stress_step_source),
            args=args,
            case_name=case_name,
            sigma_a_MPa=float(source_sigma),
            mesh=mesh,
            patch=patch,
            generation=source_generation,
        )
        start_block = restored["next_block"]
        cycles = restored["cycles"]
        Wp_total = restored["Wp_total"]
        root_xy = restored["root_xy"]
        ep_gp = restored["ep_gp"]
        rho_gp = restored["rho_gp"]
        epsp_acc_gp = restored["epsp_acc_gp"]
        u = restored["u"]
        last_residual = restored["last_residual"]
        pd_state = restored["pd_state"]
        rows = restored["rows"]
        controller_next_block_cycles = restored["controller_next_block_cycles"]
        resumed = True
        geometry_audit = _geometry_resolution_audit(
            mesh, feature_nodes, patch.point_spacing_m, initial_min_area, args
        )
        print(
            f"STATEFUL_PD_V9 stress-step restored {case_name} "
            f"from sigma_a={float(source_sigma):g} MPa at N={cycles:.9e}; "
            f"continuing at sigma_a={sigma_a_MPa:g} MPa"
        )
    elif args.resume and checkpoint_path.exists():
        restored = _load_case_checkpoint(
            checkpoint_path,
            args=args,
            case_name=case_name,
            sigma_a_MPa=sigma_a_MPa,
            mesh=mesh,
            patch=patch,
        )
        start_block = restored["next_block"]
        cycles = restored["cycles"]
        Wp_total = restored["Wp_total"]
        root_xy = restored["root_xy"]
        ep_gp = restored["ep_gp"]
        rho_gp = restored["rho_gp"]
        epsp_acc_gp = restored["epsp_acc_gp"]
        u = restored["u"]
        last_residual = restored["last_residual"]
        pd_state = restored["pd_state"]
        rows = restored["rows"]
        controller_next_block_cycles = restored["controller_next_block_cycles"]
        resumed = True
        if rows:
            geometry_saturated = bool(rows[-1].get("geometry_saturated", False))
            raw_sat_cycles = rows[-1].get("geometry_saturation_cycles", np.nan)
            if raw_sat_cycles is not None and np.isfinite(raw_sat_cycles):
                geometry_saturation_cycles = float(raw_sat_cycles)
        geometry_audit = _geometry_resolution_audit(
            mesh, feature_nodes, patch.point_spacing_m, initial_min_area, args
        )
        if args.geometry_limit_mode != "off" and not geometry_audit["pass"]:
            geometry_invalid_reason = ",".join(geometry_audit["failure_reasons"])
        print(
            f"STATEFUL_PD_V8_7 resumed {case_name} sigma_a={sigma_a_MPa:g} MPa "
            f"from block={start_block} N={cycles:.6e}"
        )

    with (outdir / "run_args.json").open("w") as f:
        json.dump(
            vars(args)
            | {
                "case": case_name,
                "sigma_a_MPa": sigma_a_MPa,
                "root_radius_initial_m": root_radius0,
                "pd_points": len(patch.xy),
                "pd_bonds": len(patch.bonds),
                "pd_production_preflight": preflight_audit,
                "checkpoint_path": str(checkpoint_path),
                "source_sha256": SOURCE_SHA256,
                "resumed": resumed,
                "protocol_label": str(getattr(args, "protocol_label", "") or ""),
                "stress_step_source_checkpoint": stress_step_source,
                "stress_step_source_generation": str(getattr(args, "stress_step_source_generation", "") or ""),
                "stress_step_source_sigma_a_MPa": getattr(args, "stress_step_source_sigma_a_MPa", None),
            },
            f,
            indent=2,
            sort_keys=True,
        )

    last_hist = None
    last_diag = None
    next_block = start_block
    stable_endpoint = args.fatigue_endpoint == "stable_crack_birth"
    spatial_birth_endpoint = args.fatigue_endpoint == "stable_spatial_crack_birth"
    high_cycle_mode_path = outdir / "v9_pd_high_cycle_mode_history.json"
    high_cycle_mode_rows = []
    if resumed and high_cycle_mode_path.is_file():
        existing_mode_rows = json.loads(high_cycle_mode_path.read_text())
        if not isinstance(existing_mode_rows, list):
            raise RuntimeError("high-cycle mode history is not a JSON list")
        high_cycle_mode_rows = existing_mode_rows
    high_cycle_retry_after = float(args.pd_high_cycle_start_cycles)

    def build_high_cycle_adapter():
        def evaluator(adapter, dN=1.0):
            payload = evaluate_dormant_exact_cycle(
                args=args, shield_on=shield_on, mesh=adapter.mesh, patch=adapter.patch,
                pd_state=adapter.pd_state, crack=crack, plast_chain=plast_chain,
                cached_fem=cached_fem, fem_transaction=fem_transaction,
                sigma_max=sigma_max, sigma_min=sigma_min, ep_gp=adapter.ep_gp,
                rho_gp=adapter.rho_gp, epsp_acc_gp=adapter.epsp_acc_gp, u=adapter.u,
                plastic_work=adapter.plastic_work, cycles=adapter.cycles, dN=dN,
            )
            adapter.ep_gp = payload.pop("ep_gp"); adapter.rho_gp = payload.pop("rho_gp")
            adapter.epsp_acc_gp = payload.pop("epsp_acc_gp"); adapter.u = payload.pop("u")
            adapter.pd_state.log_delivery_memory = payload.pop("log_delivery_memory")
            adapter.pd_state.delivery_memory = np.where(np.isfinite(adapter.pd_state.log_delivery_memory),
                np.exp(np.minimum(adapter.pd_state.log_delivery_memory,709.0)),0.0)
            for field in ("available","embryo","stable","inactive","completion"):
                setattr(adapter.pd_state,field,payload.pop(field))
            return payload
        return SpatialPDDormantAdapter(
            patch=patch, pd_state=pd_state, mesh=mesh, ep_gp=ep_gp, rho_gp=rho_gp,
            epsp_acc_gp=epsp_acc_gp, u=u, cycles=cycles, plastic_work=Wp_total,
            cycle_evaluator=evaluator, window_evaluator=evaluator,
        )

    for ib in range(start_block, args.max_blocks):
        if (
            cycles >= args.cycles_max
            or pd_state.cycles_connected is not None
            or pd_state.cycles_diffuse_abort is not None
            or pd_state.cycles_front_stalled is not None
            or pd_state.cycles_precapture_stalled is not None
            or geometry_invalid_reason is not None
            or (stable_endpoint and pd_state.cycles_first_stable is not None)
            or (spatial_birth_endpoint and pd_state.cycles_front_capture is not None)
        ):
            break

        if args.pd_high_cycle and cycles >= high_cycle_retry_after:
            adapter = build_high_cycle_adapter()
            eligible, _ = adapter.dormant_eligibility()
            if eligible:
                decade = 10.0 ** math.ceil(math.log10(max(cycles, 1.0)) - 1e-14)
                if decade <= cycles * (1.0 + 1e-14): decade *= 10.0
                segment = min(args.cycles_max - cycles, args.pd_high_cycle_max_segment,
                              max(decade - cycles, 1.0))
                engine = DormantPDHighCycleEngine(adapter, HighCycleConfig(
                    projective_max_cycles=max(int(args.pd_high_cycle_max_segment), 2),
                    exact_retry_cycles=0,
                    minimum_projected_cycles_per_exact_map=16.0))
                hc = engine.advance(segment)
                high_cycle_mode_rows.extend({
                    "cycles_total": adapter.cycles, "requested_segment": segment,
                    "mode": mode.mode, "accepted_cycles": mode.cycles,
                    "exact_map_evaluations": mode.exact_map_evaluations,
                    "accepted": mode.accepted, "detail": mode.detail,
                } for mode in hc.modes)
                mode_path = high_cycle_mode_path
                tmp_mode = mode_path.with_name(mode_path.name + ".tmp")
                tmp_mode.write_text(json.dumps(high_cycle_mode_rows, indent=2, default=_json_safe) + "\n")
                os.replace(tmp_mode, mode_path)
                engine.write_atomic_mode_checkpoint(outdir / "v9_pd_high_cycle_controller.json")
                if hc.cycles_consumed > 0.0:
                    ep_gp = adapter.ep_gp.copy(); rho_gp = adapter.rho_gp.copy()
                    epsp_acc_gp = adapter.epsp_acc_gp.copy(); u = adapter.u.copy()
                    Wp_total = adapter.plastic_work; cycles = adapter.cycles
                    _, _, u_zero_hc, _, _ = cached_fem.affine(ep_gp, sigma_max, sigma_min, u)
                    _, _, s1_res_hc, _ = stress_state_intact(mesh, u_zero_hc, ep_gp, Dmat, mat)
                    last_residual = project_gp_to_nodes(mesh, s1_res_hc)
                    efficiency_limited = any(mode.mode == "efficiency_budget" for mode in hc.modes)
                    if efficiency_limited:
                        high_cycle_retry_after = 10.0 ** math.ceil(
                            math.log10(max(cycles + 1.0, 10.0))
                        )
                    if not hc.event_guard_reached and not efficiency_limited:
                        continue
                else:
                    # The authoritative embedded macro-stepper is the efficient
                    # exact transient path. Avoid repeating expensive private
                    # qualifications on every accepted macro block.
                    high_cycle_retry_after = max(
                        cycles + 8.0 * max(controller_next_block_cycles, 1.0),
                        10.0 ** math.ceil(math.log10(max(cycles + 1.0, 10.0))),
                    )

        geometry_audit = _geometry_resolution_audit(
            mesh, feature_nodes, patch.point_spacing_m, initial_min_area, args
        )
        if args.geometry_limit_mode != "off" and not geometry_audit["pass"]:
            geometry_invalid_reason = ",".join(geometry_audit["failure_reasons"])
            break
        geometry_taper = _geometry_taper_factor(geometry_audit, args)
        if (
            args.geometry_limit_mode == "freeze"
            and args.enable_geometry_evolution
            and geometry_taper <= 0.0
        ):
            geometry_saturated = True
            if geometry_saturation_cycles is None:
                geometry_saturation_cycles = float(cycles)

        Umax, Umin, u_zero, F0, _ = cached_fem.affine(
            ep_gp, sigma_max, sigma_min, u
        )
        cyc = cached_fem.representative_cycle(
            ep_gp, rho_gp, Umax, Umin, args.T, args.frequency_Hz,
            args.plastic_n_phase, plast_chain, u_zero, args.k_store,
            args.k_dyn, args.rho_floor, args.rho_cap, args.max_dep_phase,
            args.max_rho_rel_phase,
        )
        dep_tensor_cycle = cyc["dep_tensor_cycle"]
        dep_eq_cycle = np.maximum(cyc["dep_eq_cycle"], 0.0)
        drho_cycle = cyc["drho_cycle"]

        epsp_node_pre, rho_node_pre, P, Dloc = project_plastic_state(
            mesh, epsp_acc_gp, rho_gp, args.epsp_shield_scale, args.epsp_damage_scale
        )
        chi = args.shield_chi if shield_on else 0.0
        Gsh = args.Gshield_eV if shield_on else 0.0
        sigma_back = args.sigma_back_max_GPa * 1e9 * P
        state_shift = Gsh * P - args.Gstored_eV * Dloc
        hist_pre = cached_fem.stress_histories(
            ep_gp, Umax, Umin, args.hazard_n_phase, u_zero
        )
        delivery_pre = phase_resolved_delivery_rate(
            args, plast_chain, hist_pre["seq_node"], rho_node_pre, args.T
        )
        log_delivery_pre = phase_resolved_delivery_log_rate(
            args, plast_chain, hist_pre["seq_node"], rho_node_pre, args.T
        )
        ep_node_pre = project_gp_to_nodes(mesh, ep_gp)
        _, _, _, bond_amp_pre, point_amp_pre = patch.solve_local_mechanics(
            pd_state, hist_pre["u_max"], ep_node_pre
        )
        pre_rates = patch.preview_rates(
            pd_state,
            crack,
            hist_pre["sigma_node"],
            delivery_pre,
            args.T,
            args.frequency_Hz,
            state_shift,
            sigma_back,
            chi,
            P,
            point_amp_pre,
            bond_amp_pre,
        )

        geometry_active = (
            args.enable_geometry_evolution
            and not geometry_saturated
            and geometry_invalid_reason is None
            and not (args.freeze_ale_after_front_capture and pd_state.active_front)
        )
        remaining = args.cycles_max - cycles
        dN = min(args.block_cycles, remaining)
        max_dep_cycle = float(np.max(dep_eq_cycle)) if dep_eq_cycle.size else 0.0
        if max_dep_cycle > 0.0:
            dN = min(dN, args.target_dep_eq_block / max_dep_cycle)
        rel_rho = float(np.max(np.abs(drho_cycle) / np.maximum(rho_gp, args.rho0)))
        if rel_rho > 0.0:
            dN = min(dN, args.target_rho_rel_block / rel_rho)
        max_delivery = float(np.max(pre_rates["delivery_events_per_cycle"]))
        if (
            np.isfinite(args.target_delivery_events)
            and args.target_delivery_events > 0.0
            and max_delivery > 0.0
        ):
            dN = min(dN, args.target_delivery_events / max_delivery)
        if pre_rates["max_rate_per_cycle"] > 0.0:
            target_hazard = -math.log(max(1.0 - args.max_transition_probability, 1e-12))
            dN = min(dN, target_hazard / pre_rates["max_rate_per_cycle"])

        next_birth_wait = float("inf")
        block_limited_by_birth_clock = False
        if args.align_blocks_to_birth_clock:
            next_birth_wait = patch.next_birth_wait_cycles(
                pd_state, pre_rates["mu_birth_bound"]
            )
            if np.isfinite(next_birth_wait) and next_birth_wait <= dN:
                dN = max(
                    next_birth_wait * args.birth_clock_safety_factor,
                    min(args.min_block_cycles, 1e-6),
                )
                block_limited_by_birth_clock = True

        # Apply this exact event limiter after the birth-clock safety factor:
        # that factor intentionally overshoots a birth threshold and must
        # never overshoot an already-active embryo's stable/heal transition.
        next_embryo_transition_wait = patch.next_embryo_transition_wait_cycles(
            pd_state, pre_rates["mu_stab"], pre_rates["mu_heal"]
        )
        block_limited_by_embryo_transition = bool(
            np.isfinite(next_embryo_transition_wait)
            and next_embryo_transition_wait <= dN
        )
        if np.isfinite(next_embryo_transition_wait):
            dN = min(dN, next_embryo_transition_wait)

        if geometry_active and max_dep_cycle > 0.0:
            dh_cycle, _, _ = surface_morphology_proposal(
                mesh,
                feature_nodes,
                dep_tensor_cycle,
                args.morph_band_length_m,
                args.morph_normal_weight,
                args.morph_shear_weight,
            )
            dh_cycle = np.asarray(dh_cycle, float) * float(geometry_taper)
            max_dh_cycle = float(np.max(np.abs(dh_cycle))) if dh_cycle.size else 0.0
            if max_dh_cycle > 0.0:
                dN = min(dN, args.target_surface_move_fraction * mesh.hbar_tip / max_dh_cycle)

        # Reuse the accepted embedded controller's next proposal.  Starting
        # every macro-step from the broad physical bound caused 5--8 discarded
        # FEM solves per accepted interval in quiet long-life trajectories.
        controller_before_event_localization = float(controller_next_block_cycles)
        if controller_next_block_cycles > 0.0:
            dN = min(dN, controller_next_block_cycles)

        if block_limited_by_embryo_transition:
            # Exact first-passage boundaries take precedence over the normal
            # macro-step floor, including when the residual wait is sub-floor.
            dN = min(dN, remaining)
            representable_step = float(np.nextafter(float(cycles), math.inf) - float(cycles))
            if 0.0 < dN < representable_step:
                # A positive sub-ULP wait cannot advance the absolute cycle
                # coordinate and otherwise repeats forever. The next
                # representable boundary is the tightest localizable event.
                dN = min(representable_step, remaining)
        elif block_limited_by_birth_clock:
            dN = max(min(dN, remaining), min(args.min_block_cycles, 1e-6))
        else:
            dN = max(min(dN, remaining), args.min_block_cycles)
        dN = min(dN, remaining)
        if dN <= 0.0:
            break

        # Side-effect-free embedded proposal. Rejected attempts mutate neither
        # physical arrays nor stochastic state and therefore cannot consume RNG.
        fem_initial = FEMPhysicalState(ep_gp, rho_gp, epsp_acc_gp, u, Wp_total)
        fem_rejections = 0
        while True:
            fem_proposal = fem_transaction.propose(fem_initial, dN, first_cycle=cyc)
            if fem_proposal.normalized_error <= 1.0:
                break
            dN *= 0.5
            fem_rejections += 1
            if dN < args.min_block_cycles or fem_rejections > 128:
                raise RuntimeError(
                    "v9 embedded FEM transaction could not satisfy tolerance "
                    f"before minimum block; error={fem_proposal.normalized_error:g}"
                )
        if fem_proposal.normalized_error <= 1e-16:
            controller_factor = 2.0
        else:
            controller_factor = min(
                2.0, max(1.05, 0.9 / math.sqrt(fem_proposal.normalized_error))
            )
        controller_next_block_cycles = max(
            args.min_block_cycles, dN * controller_factor
        )
        if block_limited_by_embryo_transition and controller_before_event_localization > 0.0:
            controller_next_block_cycles = max(
                controller_next_block_cycles, controller_before_event_localization
            )
        dep_tensor_block = fem_proposal.state.ep_gp - ep_gp
        dep_eq_block = fem_proposal.state.epsp_acc_gp - epsp_acc_gp
        ep_gp = np.asarray(fem_proposal.state.ep_gp).copy()
        rho_gp = np.asarray(fem_proposal.state.rho_gp).copy()
        epsp_acc_gp = np.asarray(fem_proposal.state.epsp_acc_gp).copy()
        u = np.asarray(fem_proposal.state.u).copy()
        Wp_total = float(fem_proposal.state.plastic_work_J_per_m)

        mesh_scale = 0.0
        dh_block = np.zeros(len(feature_nodes))
        gamma_nt = np.zeros(len(feature_nodes))
        if geometry_active:
            dh_block, _, gamma_nt = surface_morphology_proposal(
                mesh,
                feature_nodes,
                dep_tensor_block,
                args.morph_band_length_m,
                args.morph_normal_weight,
                args.morph_shear_weight,
            )
            dh_block = np.asarray(dh_block, float) * float(geometry_taper)
            old_nodes = mesh.nodes.copy()
            old_root_xy = np.asarray(root_xy, float).copy()
            mesh_scale = apply_local_ale_surface_update(
                mesh,
                feature_nodes,
                dh_block,
                decay_length=args.morph_decay_length_m,
                fixed_nodes=fixed_mesh_nodes,
                max_move=args.max_surface_move_fraction * mesh.hbar_tip,
                min_area_fraction=args.min_area_fraction,
            )
            if mesh_scale > 0.0:
                root_xy = local_root_xy(mesh, feature_nodes)
                rebuild_mesh_geometry(mesh, root_xy)
                patch.update_geometry(mesh, root_xy)
                proposed_audit = _geometry_resolution_audit(
                    mesh, feature_nodes, patch.point_spacing_m, initial_min_area, args
                )
                if args.geometry_limit_mode != "off" and not proposed_audit["pass"]:
                    # Roll back the geometry; material-state evolution for this
                    # accepted fatigue block is retained and the termination or
                    # freeze decision is recorded explicitly.
                    mesh.nodes[:] = old_nodes
                    root_xy = old_root_xy
                    rebuild_mesh_geometry(mesh, root_xy)
                    patch.update_geometry(mesh, root_xy)
                    mesh_scale = 0.0
                    dh_block[:] = 0.0
                    if args.geometry_limit_mode == "freeze":
                        geometry_saturated = True
                        if geometry_saturation_cycles is None:
                            geometry_saturation_cycles = float(cycles + dN)
                    else:
                        geometry_invalid_reason = ",".join(proposed_audit["failure_reasons"])
            geometry_audit = _geometry_resolution_audit(
                mesh, feature_nodes, patch.point_spacing_m, initial_min_area, args
            )

        Umax2, Umin2, u_zero2, F02, _ = cached_fem.affine(
            ep_gp, sigma_max, sigma_min, u_zero
        )
        _, _, s1_res, _ = stress_state_intact(mesh, u_zero2, ep_gp, Dmat, mat)
        residual_node = project_gp_to_nodes(mesh, s1_res)
        last_residual = residual_node.copy()

        epsp_node, rho_node_post, P, Dloc = project_plastic_state(
            mesh, epsp_acc_gp, rho_gp, args.epsp_shield_scale, args.epsp_damage_scale
        )
        sigma_back = args.sigma_back_max_GPa * 1e9 * P
        state_shift = Gsh * P - args.Gstored_eV * Dloc
        hist_post = cached_fem.stress_histories(
            ep_gp, Umax2, Umin2, args.hazard_n_phase, u_zero2
        )
        delivery_post = phase_resolved_delivery_rate(
            args, plast_chain, hist_post["seq_node"], rho_node_post, args.T
        )
        log_delivery_post = phase_resolved_delivery_log_rate(
            args, plast_chain, hist_post["seq_node"], rho_node_post, args.T
        )
        ep_node_post = project_gp_to_nodes(mesh, ep_gp)
        _, _, _, bond_amp_post, point_amp_post = patch.solve_local_mechanics(
            pd_state, hist_post["u_max"], ep_node_post
        )
        # Mid-block representative cycle: preserve the physical phase order
        # while averaging the pre/post state, rather than concatenating two
        # cycles and compressing them into one period.
        sigma_combined = 0.5 * (hist_pre["sigma_node"] + hist_post["sigma_node"])
        delivery_combined = 0.5 * (delivery_pre + delivery_post)
        log_delivery_combined = np.logaddexp(log_delivery_pre, log_delivery_post) - math.log(2.0)
        point_amp = 0.5 * (point_amp_pre + point_amp_post)
        bond_amp = 0.5 * (bond_amp_pre + bond_amp_post)
        diag = patch.update(
            pd_state,
            crack,
            sigma_combined,
            delivery_combined,
            args.T,
            args.frequency_Hz,
            dN,
            cycles,
            state_shift,
            sigma_back,
            chi,
            P,
            point_amp,
            bond_amp,
            log_delivery_rate_phase_global=log_delivery_combined,
        )
        cycles += dN
        last_hist = hist_post
        last_diag = diag

        root_radius = local_root_radius(mesh, feature_nodes)
        row = {
            "block": ib,
            "case": case_name,
            "fatigue_model": args.fatigue_model,
            "fatigue_model_preset_applied": bool(args.fatigue_model_preset_applied),
            "delivery_source": args.delivery_source,
            "delivery_hit_count": args.delivery_hit_count,
            "delivery_memory_s": args.delivery_memory_s,
            "delivery_scale": args.delivery_scale,
            "delivery_source_multiplicity": float(
                getattr(args, "delivery_source_multiplicity", 1.0)
            ),
            "sigma_a_MPa": sigma_a_MPa,
            "cycles_total": cycles,
        "cycles_target": float(args.cycles_max),
        "max_blocks_requested": int(args.max_blocks),
            "dN": dN,
            "v9_fem_embedded_error": fem_proposal.normalized_error,
            "v9_fem_error_ep_gp": fem_proposal.component_errors.get("ep_gp", np.nan),
            "v9_fem_error_rho_gp": fem_proposal.component_errors.get("rho_gp", np.nan),
            "v9_fem_error_epsp_acc_gp": fem_proposal.component_errors.get("epsp_acc_gp", np.nan),
            "v9_fem_rejections": fem_rejections,
            "next_birth_wait_cycles_pre": next_birth_wait,
            "block_limited_by_birth_clock": bool(block_limited_by_birth_clock),
            "Umax_m": Umax2,
            "Umin_m": Umin2,
            "F0_residual_N_per_m": F02,
            "dep_eq_cycle_max": max_dep_cycle,
            "dep_eq_block_max": float(np.max(dep_eq_block)),
            "rho_max_m2": float(np.max(rho_gp)),
            "rho_mean_m2": float(np.mean(rho_gp)),
            "epsp_acc_max": float(np.max(epsp_acc_gp)),
            "residual_sigma1_max_Pa": float(np.max(residual_node)),
            "sigma1_cycle_max_Pa": float(np.max(hist_post["s1_node"])),
            "remote_sigma_max_Pa": float(sigma_max),
            "fem_scratch_Kt_max": float(np.max(hist_post["s1_node"]) / max(sigma_max, 1e-30)),
            "effective_scratch_Kt_max": float(diag.max_effective_stress_Pa / max(sigma_max, 1e-30)),
            "P_max": float(np.max(P)),
            "Dloc_max": float(np.max(Dloc)),
            "pd_delivery_memory_max": diag.max_delivery_memory,
            "pd_completion_max": diag.max_completion,
            "pd_embryo_max": diag.max_embryo,
            "pd_stable_max": diag.max_stable,
            "pd_growth_max": diag.max_growth,
            "pd_expected_embryos": diag.expected_embryos,
            "pd_expected_births_cumulative": diag.expected_births_cumulative,
            "pd_expected_stable": diag.expected_stable,
            "pd_realized_embryos": diag.realized_embryos,
            "pd_realized_births_cumulative": diag.realized_births_cumulative,
            "pd_realized_stable": diag.realized_stable,
            "pd_bond_damage_max": diag.max_bond_damage,
            "pd_broken_bonds": diag.broken_bonds,
            "pd_connected_extent_m": diag.connected_extent_m,
            "pd_connected_bonds": diag.connected_bonds,
            "pd_crack_centerline_length_m": diag.crack_centerline_length_m,
            "pd_crack_width_m": diag.crack_width_m,
            "pd_crack_tip_width_m": diag.crack_tip_width_m,
            "pd_crack_width_ratio": diag.crack_width_ratio,
            "pd_crack_tip_radius_eff_m": diag.crack_tip_radius_eff_m,
            "pd_crack_orientation_deg": diag.crack_orientation_deg,
            "pd_crack_boundary_clearance_m": diag.crack_boundary_clearance_m,
            "pd_crack_remote_K_MPam05": diag.crack_remote_K_MPam05,
            "pd_crack_local_upper_K_MPam05": diag.crack_local_upper_K_MPam05,
            "pd_crack_remote_KI_MPam05": diag.crack_remote_KI_MPam05,
            "pd_crack_remote_KII_MPam05": diag.crack_remote_KII_MPam05,
            "pd_crack_orientation_coherence": diag.crack_orientation_coherence,
            "pd_crack_axial_coverage": diag.crack_axial_coverage,
            "pd_crack_max_axial_gap_m": diag.crack_max_axial_gap_m,
            "pd_crack_surface_gap_m": diag.crack_surface_gap_m,
            "pd_crack_slenderness": diag.crack_slenderness,
            "pd_handoff_pass": diag.handoff_pass,
            "pd_handoff_length_pass": diag.handoff_length_pass,
            "pd_handoff_slenderness_pass": diag.handoff_slenderness_pass,
            "pd_handoff_tip_pass": diag.handoff_tip_pass,
            "pd_handoff_boundary_pass": diag.handoff_boundary_pass,
            "pd_handoff_K_pass": diag.handoff_K_pass,
            "pd_handoff_orientation_pass": diag.handoff_orientation_pass,
            "pd_handoff_coverage_pass": diag.handoff_coverage_pass,
            "pd_handoff_surface_pass": diag.handoff_surface_pass,
            "pd_handoff_front_pass": diag.handoff_front_pass,
            "pd_handoff_competition_pass": diag.handoff_competition_pass,
            "pd_active_front": diag.active_front,
            "pd_active_front_length_m": diag.active_front_length_m,
            "pd_active_front_bonds": diag.active_front_bonds,
            "pd_front_backbone_bonds": diag.front_backbone_bonds,
            "pd_front_wake_bonds": diag.front_wake_bonds,
            "pd_front_process_bonds": diag.front_process_bonds,
            "pd_front_candidate_mode": diag.front_candidate_mode,
            "pd_front_eligible_preferred": diag.front_eligible_preferred,
            "pd_front_eligible_fallback": diag.front_eligible_fallback,
            "pd_front_stall_updates": diag.front_stall_updates,
            "pd_front_stalled": diag.front_stalled,
            "pd_front_max_link_rate_per_cycle": diag.front_max_link_rate_per_cycle,
            "pd_primary_seed_node": diag.primary_seed_node,
            "pd_primary_seed_stall_updates": diag.primary_seed_stall_updates,
            "pd_primary_seed_reselections": diag.primary_seed_reselections,
            "pd_precapture_stalled": diag.precapture_stalled,
            "pd_off_front_broken_bonds": diag.off_front_broken_bonds,
            "pd_off_front_broken_fraction": diag.off_front_broken_fraction,
            "pd_diffuse_abort": diag.diffuse_abort,
            "pd_amplification_max": diag.max_pd_amplification,
            "pd_effective_stress_max_Pa": diag.max_effective_stress_Pa,
            "pd_max_rate_per_cycle": diag.max_rate_per_cycle,
            "pd_delivery_rate_max_s": diag.max_delivery_rate_s,
            "pd_delivery_events_cycle_max": diag.max_delivery_events_per_cycle,
            "pd_nucleation_rate_max_s": diag.max_nucleation_rate_s,
            "pd_nucleation_hazard_cycle_max": diag.max_nucleation_hazard_per_cycle,
            "pd_temporal_transient_cycles": diag.temporal_transient_cycles,
            "pd_birth_rate_max_per_cycle": diag.max_birth_rate_per_cycle,
            "pd_expected_candidate_sites": diag.expected_candidate_sites,
            "pd_realized_candidate_sites": diag.realized_candidate_sites,
            "root_x_m": root_xy[0],
            "root_y_m": root_xy[1],
            "root_radius_m": root_radius,
            "root_radius_over_initial": root_radius / max(root_radius0, 1e-30) if np.isfinite(root_radius) else np.nan,
            "root_radius_over_spacing": float(geometry_audit.get("root_radius_over_spacing", np.nan)),
            "minimum_element_area_over_initial": float(geometry_audit.get("minimum_element_area_over_initial", np.nan)),
            "geometry_limit_mode": args.geometry_limit_mode,
            "geometry_taper_factor": float(geometry_taper),
            "geometry_saturated": bool(geometry_saturated),
            "geometry_saturation_cycles": geometry_saturation_cycles if geometry_saturation_cycles is not None else np.nan,
            "geometry_invalid_reason": geometry_invalid_reason or "",
            "surface_move_max_m": float(np.max(np.abs(dh_block))) if dh_block.size else 0.0,
            "surface_shear_gamma_nt_max": float(np.max(np.abs(gamma_nt))) if gamma_nt.size else 0.0,
            "ale_accept_scale": mesh_scale,
            "ale_geometry_active": bool(geometry_active),
            "plastic_work_J_per_m": Wp_total,
            "cycles_first_embryo": pd_state.cycles_first_embryo if pd_state.cycles_first_embryo is not None else np.nan,
            "cycles_first_stable": pd_state.cycles_first_stable if pd_state.cycles_first_stable is not None else np.nan,
            "cycles_first_expected_embryo": pd_state.cycles_first_expected_embryo if pd_state.cycles_first_expected_embryo is not None else np.nan,
            "cycles_first_expected_stable": pd_state.cycles_first_expected_stable if pd_state.cycles_first_expected_stable is not None else np.nan,
            "cycles_first_softening": pd_state.cycles_first_softening if pd_state.cycles_first_softening is not None else np.nan,
            "cycles_root_connected": pd_state.cycles_root_connected if pd_state.cycles_root_connected is not None else np.nan,
            "cycles_two_horizon_crack": pd_state.cycles_two_horizon_crack if pd_state.cycles_two_horizon_crack is not None else np.nan,
            "cycles_front_capture": pd_state.cycles_front_capture if pd_state.cycles_front_capture is not None else np.nan,
            "cycles_front_stalled": pd_state.cycles_front_stalled if pd_state.cycles_front_stalled is not None else np.nan,
            "cycles_diffuse_abort": pd_state.cycles_diffuse_abort if pd_state.cycles_diffuse_abort is not None else np.nan,
            "cycles_connected": pd_state.cycles_connected if pd_state.cycles_connected is not None else np.nan,
        }
        rows.append(row)
        next_block = ib + 1

        if args.print_every and ib % args.print_every == 0:
            print(
                f"STATEFUL_PD {case_name} sigma_a={sigma_a_MPa:g}MPa block={ib} "
                f"N={cycles:.3e} dN={dN:.3g} dep={row['dep_eq_block_max']:.2e} "
                f"Dmem={diag.max_delivery_memory:.2e} Q={diag.max_completion:.2e} "
                f"Ebirth={diag.expected_births_cumulative:.2g} Estab={diag.expected_stable:.2g} "
                f"Rbirth={diag.realized_births_cumulative:d} Rstab={diag.realized_stable:d} "
                f"omega={diag.max_bond_damage:.3g} broken={diag.broken_bonds} "
                f"a={diag.crack_centerline_length_m*1e6:.1f}um "
                f"C={diag.crack_orientation_coherence:.2f} "
                f"cov={diag.crack_axial_coverage:.2f} "
                f"front={int(diag.active_front)} off={diag.off_front_broken_fraction:.2f} "
                f"mode={diag.front_candidate_mode} pref={diag.front_eligible_preferred} "
                f"fall={diag.front_eligible_fallback} seed={diag.primary_seed_node} "
                f"reselect={diag.primary_seed_reselections} handoff={int(diag.handoff_pass)} "
                f"stall={int(diag.front_stalled)} pre={int(diag.precapture_stalled)} "
                f"geomR/h={geometry_audit.get('root_radius_over_spacing', float('nan')):.1f} "
                f"gsat={int(geometry_saturated)} abort={int(diag.diffuse_abort)}"
            )

        if getattr(args, "pd_image_policy", "selected") == "selected" and (
            ib == 0 or (args.snapshot_every > 0 and ib % args.snapshot_every == 0)
        ):
            patch.plot_snapshot(
                pd_state,
                outdir / f"pd_patch_block_{ib:05d}.png",
                title=f"{case_name}, sigma_a={sigma_a_MPa:g} MPa, N={cycles:.3e}",
            )
            _plot_fem_fields(
                mesh,
                [
                    ("accumulated eps_p", epsp_node),
                    ("rho (m^-2)", rho_node_post),
                    ("residual sigma1 (MPa)", residual_node * 1e-6),
                ],
                outdir / f"fem_fields_block_{ib:05d}.png",
                root_xy,
                feature_nodes,
            )
        u = hist_post["u_end"].copy()
        checkpoint_due = (
            args.checkpoint_every_blocks > 0
            and (
                (ib + 1) % args.checkpoint_every_blocks == 0
                or pd_state.cycles_connected is not None
                or pd_state.cycles_diffuse_abort is not None
                or pd_state.cycles_front_stalled is not None
                or pd_state.cycles_precapture_stalled is not None
                or geometry_invalid_reason is not None
                or (stable_endpoint and pd_state.cycles_first_stable is not None)
                or (spatial_birth_endpoint and pd_state.cycles_front_capture is not None)
            )
        )
        if checkpoint_due:
            _write_csv(outdir / "sn_stateful_pd_history_partial.csv", rows)
            _save_case_checkpoint(
                checkpoint_path,
                args=args,
                case_name=case_name,
                sigma_a_MPa=sigma_a_MPa,
                next_block=next_block,
                controller_next_block_cycles=controller_next_block_cycles,
                cycles=cycles,
                Wp_total=Wp_total,
                mesh=mesh,
                root_xy=root_xy,
                ep_gp=ep_gp,
                rho_gp=rho_gp,
                epsp_acc_gp=epsp_acc_gp,
                u=u,
                last_residual=last_residual,
                pd_state=pd_state,
                patch=patch,
                rows=rows,
            )

    _write_csv(outdir / "sn_stateful_pd_history.csv", rows)
    if args.checkpoint_every_blocks > 0:
        _save_case_checkpoint(
            checkpoint_path,
            args=args,
            case_name=case_name,
            sigma_a_MPa=sigma_a_MPa,
            next_block=next_block,
            controller_next_block_cycles=controller_next_block_cycles,
            cycles=cycles,
            Wp_total=Wp_total,
            mesh=mesh,
            root_xy=root_xy,
            ep_gp=ep_gp,
            rho_gp=rho_gp,
            epsp_acc_gp=epsp_acc_gp,
            u=u,
            last_residual=last_residual,
            pd_state=pd_state,
            patch=patch,
            rows=rows,
        )
    image_policy = getattr(args, "pd_image_policy", "selected")
    save_terminal_image = image_policy == "selected" or (
        image_policy == "event_only"
        and (pd_state.cycles_first_stable is not None or pd_state.cycles_front_capture is not None)
    )
    if save_terminal_image:
        patch.plot_snapshot(pd_state, outdir / "pd_patch_final.png", title=f"final N={cycles:.3e}")
        patch.plot_initiation_diagnostics(
            pd_state,
            outdir / "pd_initiation_diagnostics_final.png",
            title=f"{case_name}, sigma_a={sigma_a_MPa:g} MPa, final N={cycles:.3e}",
        )
    epsp_node, rho_node, P, Dloc = project_plastic_state(
        mesh, epsp_acc_gp, rho_gp, args.epsp_shield_scale, args.epsp_damage_scale
    )
    if save_terminal_image:
        _plot_fem_fields(
            mesh,
            [
                ("accumulated eps_p", epsp_node),
                ("rho (m^-2)", rho_node),
                ("residual sigma1 (MPa)", last_residual * 1e-6),
            ],
            outdir / "fem_fields_final.png",
            root_xy,
            feature_nodes,
        )
    np.savez_compressed(
        outdir / "pd_state_final.npz",
        global_nodes=patch.global_nodes,
        xy=patch.xy,
        bonds=patch.bonds,
        boundary=patch.boundary,
        initiation_weight=patch.initiation_weight,
        mean_candidate_sites=patch.mean_candidate_sites,
        available=pd_state.available,
        embryo=pd_state.embryo,
        stable=pd_state.stable,
        inactive=pd_state.inactive,
        candidate_sites=pd_state.candidate_sites,
        available_sites=pd_state.available_sites,
        embryo_sites=pd_state.embryo_sites,
        stable_sites=pd_state.stable_sites,
        inactive_sites=pd_state.inactive_sites,
        born_sites_cumulative=pd_state.born_sites_cumulative,
        healed_sites_cumulative=pd_state.healed_sites_cumulative,
        delivery_memory=pd_state.delivery_memory,
        completion=pd_state.completion,
        growth=pd_state.growth,
        crack_normal_c2=pd_state.crack_normal_c2,
        crack_normal_s2=pd_state.crack_normal_s2,
        crack_orientation_weight=pd_state.crack_orientation_weight,
        feature_surface_xy=patch.feature_surface_xy,
        bond_midpoints=patch.bond_midpoints,
        born_cumulative=pd_state.born_cumulative,
        healed_cumulative=pd_state.healed_cumulative,
        bond_damage=pd_state.bond_damage,
        primary_seed_node=np.asarray(pd_state.primary_seed_node),
        primary_seed_rejected=pd_state.primary_seed_rejected,
        primary_seed_stall_updates=np.asarray(pd_state.primary_seed_stall_updates),
        primary_seed_reselections=np.asarray(pd_state.primary_seed_reselections),
        primary_seed_last_progress=np.asarray(pd_state.primary_seed_last_progress),
        primary_seed_selected_cycles=np.asarray(pd_state.primary_seed_selected_cycles),
        active_front=np.asarray(pd_state.active_front),
        active_front_normal_c2=np.asarray(pd_state.active_front_normal_c2),
        active_front_normal_s2=np.asarray(pd_state.active_front_normal_s2),
        active_front_contact_xy=pd_state.active_front_contact_xy,
        active_front_tip_xy=pd_state.active_front_tip_xy,
        active_front_length_m=np.asarray(pd_state.active_front_length_m),
        active_front_bonds=pd_state.active_front_bonds,
        front_backbone_bonds=pd_state.front_backbone_bonds,
        front_wake_bonds=pd_state.front_wake_bonds,
        front_process_bonds=pd_state.front_process_bonds,
        active_front_path_xy=pd_state.active_front_path_xy,
        front_candidate_mode=np.asarray(pd_state.front_candidate_mode),
        front_stall_updates=np.asarray(pd_state.front_stall_updates),
        diffuse_bad_updates=np.asarray(pd_state.diffuse_bad_updates),
        point_amplification=getattr(patch, "last_point_amp", np.ones(len(patch.xy))),
        bond_amplification=getattr(patch, "last_bond_amp", np.ones(len(patch.bonds))),
        delivery_rate_s=getattr(patch, "last_rates", {}).get("delivery_rate_s", np.zeros(len(patch.xy))),
        delivery_events_per_cycle=getattr(patch, "last_rates", {}).get("delivery_events_per_cycle", np.zeros(len(patch.xy))),
        nucleation_rate_s=getattr(patch, "last_rates", {}).get("nucleation_rate_s", np.zeros(len(patch.xy))),
        nucleation_hazard_per_cycle=getattr(patch, "last_rates", {}).get("nucleation_hazard_per_cycle", np.zeros(len(patch.xy))),
        birth_intensity_per_cycle=getattr(patch, "last_rates", {}).get("mu_birth", np.zeros(len(patch.xy))),
        effective_stress_Pa=getattr(patch, "last_rates", {}).get("smax", np.zeros(len(patch.xy))),
    )

    final_connected_extent, final_connected_bonds = patch.connected_crack(pd_state)
    final_effective = float(last_diag.max_effective_stress_Pa) if last_diag is not None else 0.0
    final_audit = patch.crack_handoff_audit(
        pd_state, remote_sigma_max_Pa=sigma_max, local_effective_stress_Pa=final_effective
    )
    geometry_audit_final = _geometry_resolution_audit(
        mesh, feature_nodes, patch.point_spacing_m, initial_min_area, args
    )
    if spatial_birth_endpoint and pd_state.cycles_front_capture is not None:
        status = "stable_spatial_crack_birth"
    elif stable_endpoint and pd_state.cycles_first_stable is not None:
        status = "stable_crack_birth"
    elif geometry_invalid_reason is not None:
        status = "geometry_invalid_underresolved"
    elif pd_state.cycles_precapture_stalled is not None:
        status = "morphology_invalid_precapture_stalled"
    elif pd_state.cycles_connected is not None:
        status = (
            "physical_handoff_geometry_saturated"
            if geometry_saturated else "physical_handoff"
        )
    elif pd_state.cycles_front_stalled is not None:
        status = "morphology_invalid_stalled"
    elif pd_state.cycles_diffuse_abort is not None:
        status = "morphology_invalid_diffuse"
    elif cycles >= args.cycles_max:
        status = (
            "right_censored_geometry_saturated"
            if geometry_saturated else "right_censored"
        )
    else:
        status = "max_blocks_reached"
    final_rates = getattr(patch, "last_rates", {})
    birth_final = np.asarray(final_rates.get("mu_birth", np.zeros(len(patch.xy))), float)
    delivery_final = np.asarray(final_rates.get("delivery_rate_s", np.zeros(len(patch.xy))), float)
    nucleation_final = np.asarray(final_rates.get("nucleation_rate_s", np.zeros(len(patch.xy))), float)
    effective_final = np.asarray(final_rates.get("smax", np.zeros(len(patch.xy))), float)
    site_nodes_final = np.asarray(pd_state.site_node_index, dtype=np.int64)
    stable_site_ids = np.where(np.isfinite(np.asarray(pd_state.site_stable_cycle, float)))[0]
    first_stable_site_id = (
        int(stable_site_ids[np.argmin(np.asarray(pd_state.site_stable_cycle, float)[stable_site_ids])])
        if stable_site_ids.size else None
    )
    site_available_final = np.asarray(pd_state.site_status, dtype=np.uint8) == 0
    if np.any(site_available_final):
        birth_residual_final = (
            np.asarray(pd_state.site_birth_threshold, float)[site_available_final]
            - np.asarray(pd_state.birth_cumulative_hazard, float)[site_nodes_final[site_available_final]]
        )
        min_birth_residual_final = float(np.min(np.maximum(birth_residual_final, 0.0)))
    else:
        min_birth_residual_final = float("inf")
    birth_hot = int(np.argmax(birth_final)) if birth_final.size else 0
    effective_hot = int(np.argmax(effective_final)) if effective_final.size else 0
    if last_hist is not None:
        iph, inode = np.unravel_index(np.argmax(last_hist["s1_node"]), last_hist["s1_node"].shape)
        fem_hot_xy = mesh.nodes[int(inode)]
    else:
        iph, inode = 0, 0
        fem_hot_xy = np.array([np.nan, np.nan])
    summary = {
        "model": MODEL_ID,
        "source_sha256": SOURCE_SHA256,
        "run_signature": _checkpoint_signature(args, case_name, sigma_a_MPa),
        "coupling": "one_way_FEM_to_PD_with_local_PD_redistribution",
        "case": case_name,
        "mechanics_treatment": "stateful_PD",
        "fatigue_endpoint": args.fatigue_endpoint,
        "pd_image_policy": getattr(args, "pd_image_policy", "selected"),
        "four_class_option_id": getattr(args, "four_class_option_id", None),
        "four_class_registry_audit": getattr(args, "four_class_registry_audit", None),
        "four_class_transfer_contract": getattr(args, "four_class_transfer_contract", None),
        "fatigue_model": args.fatigue_model,
        "fatigue_model_preset_applied": bool(args.fatigue_model_preset_applied),
        "fatigue_model_role": (
            ("representative" if case_name == "shielded" else "no_shield_ablation")
            if args.fatigue_model == "plastic_shielded_case64_M1"
            else ("representative" if case_name == "no_shield" else "added_shielding_ablation")
        ),
        "crack_G00_eV": float(args.crack_G00_eV),
        "crack_sigc0_GPa": float(args.crack_sigc0_GPa),
        "crack_exp_a": float(args.crack_exp_a),
        "crack_exp_n": float(args.crack_exp_n),
        "crack_floor_frac": float(args.crack_floor_frac),
        "crack_T_mode": args.crack_T_mode,
        "crack_Tref_K": float(args.crack_Tref_K),
        "crack_mu_dlnmu_dT_per_K": float(args.crack_mu_dlnmu_dT_per_K),
        "crack_G0_at_run_T_eV": float(crack._G0_sigc(args.T)[0]),
        "crack_sigc_at_run_T_GPa": float(crack._G0_sigc(args.T)[1] / 1e9),
        "emit_energy_scale": float(args.emit_energy_scale),
        "emit_entropy_scale": float(args.emit_entropy_scale),
        "peierls_energy_scale": float(args.peierls_energy_scale),
        "peierls_entropy_scale": float(args.peierls_entropy_scale),
        "taylor_energy_scale": float(args.taylor_energy_scale),
        "taylor_entropy_scale": float(args.taylor_entropy_scale),
        "sigma_a_MPa": sigma_a_MPa,
        "remote_sigma_max_Pa": float(sigma_max),
        "fem_sigma1_cycle_final_max_Pa": float(np.max(last_hist["s1_node"])) if last_hist is not None else None,
        "fem_scratch_Kt_final": float(np.max(last_hist["s1_node"]) / max(sigma_max, 1e-30)) if last_hist is not None else None,
        "effective_scratch_Kt_final": float(last_diag.max_effective_stress_Pa / max(sigma_max, 1e-30)) if last_diag is not None else None,
        "fem_sigma1_hotspot_x_m": float(fem_hot_xy[0]),
        "fem_sigma1_hotspot_y_m": float(fem_hot_xy[1]),
        "fem_sigma1_hotspot_phase_index": int(iph),
        "pd_effective_stress_hotspot_x_m": float(patch.xy[effective_hot, 0]) if len(patch.xy) else None,
        "pd_effective_stress_hotspot_y_m": float(patch.xy[effective_hot, 1]) if len(patch.xy) else None,
        "pd_birth_hotspot_x_m": float(patch.xy[birth_hot, 0]) if len(patch.xy) else None,
        "pd_birth_hotspot_y_m": float(patch.xy[birth_hot, 1]) if len(patch.xy) else None,
        "pd_delivery_rate_final_max_s": float(np.max(delivery_final)) if delivery_final.size else 0.0,
        "pd_nucleation_rate_final_max_s": float(np.max(nucleation_final)) if nucleation_final.size else 0.0,
        "pd_birth_intensity_final_max_per_cycle": float(np.max(birth_final)) if birth_final.size else 0.0,
        "pd_birth_cumulative_hazard_final_max": float(np.max(pd_state.birth_cumulative_hazard)),
        "pd_birth_threshold_min_remaining_final": min_birth_residual_final,
        "pd_birth_clock_available_sites_final": int(np.count_nonzero(site_available_final)),
        "pd_effective_link_traction_final_max_Pa": float(np.max(final_rates.get("tn", np.zeros(1)))),
        "pd_raw_link_traction_final_max_Pa": float(np.max(final_rates.get("tn_raw", np.zeros(1)))),
        "pd_link_state_shift_z_final_max": float(np.max(final_rates.get("link_state_shift_z", np.zeros(1)))),
        "T_K": args.T,
        "R": args.R,
        "frequency_Hz": args.frequency_Hz,
        "cycles_total": cycles,
        "cycles_first_embryo": pd_state.cycles_first_embryo,
        "cycles_first_stable": pd_state.cycles_first_stable,
        "cycles_first_expected_embryo": pd_state.cycles_first_expected_embryo,
        "cycles_first_expected_stable": pd_state.cycles_first_expected_stable,
        "cycles_first_softening": pd_state.cycles_first_softening,
        "cycles_root_connected": pd_state.cycles_root_connected,
        "cycles_two_horizon_crack": pd_state.cycles_two_horizon_crack,
        "cycles_front_capture": pd_state.cycles_front_capture,
        "cycles_primary_seed_reselected": pd_state.cycles_primary_seed_reselected,
        "cycles_precapture_stalled": pd_state.cycles_precapture_stalled,
        "cycles_front_stalled": pd_state.cycles_front_stalled,
        "cycles_diffuse_abort": pd_state.cycles_diffuse_abort,
        "cycles_connected": pd_state.cycles_connected,
        "status": status,
        "resumed_from_checkpoint": bool(resumed),
        "checkpoint_path": str(checkpoint_path),
        "pd_points": len(patch.xy),
        "pd_bonds": len(patch.bonds),
        "pd_production_preflight": preflight_audit,
        "pd_horizon_m": pd_cfg.horizon_m,
        "pd_patch_radius_m": pd_cfg.patch_radius_m,
        "pd_initiation_radius_m": pd_cfg.initiation_radius_m,
        "pd_initiation_taper_m": pd_cfg.initiation_taper_m,
        "pd_initiation_back_extent_m": pd_cfg.initiation_back_extent_m,
        "pd_established_extent_criterion_m": pd_cfg.established_extent_m,
        "pd_front_capture_enabled": bool(pd_cfg.front_capture_enabled),
        "freeze_ale_after_front_capture": bool(args.freeze_ale_after_front_capture),
        "pd_front_capture_min_bonds": int(pd_cfg.front_capture_min_bonds),
        "pd_front_capture_min_length_horizons": float(pd_cfg.front_capture_min_length_horizons),
        "pd_front_capture_min_orientation_coherence": float(pd_cfg.front_capture_min_orientation_coherence),
        "pd_front_capture_max_surface_gap_horizons": float(pd_cfg.front_capture_max_surface_gap_horizons),
        "pd_front_process_ahead_horizons": float(pd_cfg.front_process_ahead_horizons),
        "pd_front_process_behind_horizons": float(pd_cfg.front_process_behind_horizons),
        "pd_front_band_horizons": float(pd_cfg.front_band_horizons),
        "pd_front_wake_band_horizons": float(pd_cfg.front_wake_band_horizons),
        "pd_front_orientation_tolerance_deg": float(pd_cfg.front_orientation_tolerance_deg),
        "pd_front_preferred_orientation_tolerance_deg": float(pd_cfg.front_preferred_orientation_tolerance_deg),
        "pd_front_fallback_orientation_tolerance_deg": float(pd_cfg.front_fallback_orientation_tolerance_deg),
        "pd_front_preferred_band_horizons": float(pd_cfg.front_preferred_band_horizons),
        "pd_front_fallback_band_horizons": float(pd_cfg.front_fallback_band_horizons),
        "pd_front_fallback_activity_scale": float(pd_cfg.front_fallback_activity_scale),
        "pd_front_direction_smoothing": float(pd_cfg.front_direction_smoothing),
        "pd_front_max_turn_deg": float(pd_cfg.front_max_turn_deg),
        "pd_front_crack_tube_horizons": float(pd_cfg.front_crack_tube_horizons),
        "pd_front_stall_patience_updates": int(pd_cfg.front_stall_patience_updates),
        "pd_primary_seed_reselection_enabled": bool(pd_cfg.primary_seed_reselection_enabled),
        "pd_primary_seed_reselection_patience_updates": int(pd_cfg.primary_seed_reselection_patience_updates),
        "pd_primary_seed_progress_damage_increment": float(pd_cfg.primary_seed_progress_damage_increment),
        "pd_primary_seed_max_reselections": int(pd_cfg.primary_seed_max_reselections),
        "pd_pre_capture_birth_scale_with_primary": float(pd_cfg.pre_capture_birth_scale_with_primary),
        "pd_post_capture_birth_scale": float(pd_cfg.post_capture_birth_scale),
        "pd_post_capture_stabilization_scale": float(pd_cfg.post_capture_stabilization_scale),
        "pd_off_front_growth_scale": float(pd_cfg.off_front_growth_scale),
        "pd_diffuse_abort_enabled": bool(pd_cfg.diffuse_abort_enabled),
        "pd_diffuse_abort_min_bonds": int(pd_cfg.diffuse_abort_min_bonds),
        "pd_diffuse_abort_min_offfront_fraction": float(pd_cfg.diffuse_abort_min_offfront_fraction),
        "pd_diffuse_abort_patience_updates": int(pd_cfg.diffuse_abort_patience_updates),
        "pd_front_amplification_cap": float(pd_cfg.front_amplification_cap),
        "pd_handoff_mode": pd_cfg.handoff_mode,
        "pd_effective_point_spacing_m": float(patch.point_spacing_m),
        "pd_handoff_min_length_m": float(max(pd_cfg.handoff_min_length_m, pd_cfg.handoff_min_length_horizons * pd_cfg.horizon_m)),
        "pd_handoff_min_length_horizons": float(pd_cfg.handoff_min_length_horizons),
        "pd_handoff_max_width_ratio": float(pd_cfg.handoff_max_width_ratio),
        "pd_handoff_max_tip_width_horizons": float(pd_cfg.handoff_max_tip_width_horizons),
        "pd_handoff_min_boundary_clearance_horizons": float(pd_cfg.handoff_min_boundary_clearance_horizons),
        "pd_surface_connection_horizons": float(pd_cfg.surface_connection_horizons),
        "pd_handoff_min_orientation_coherence": float(pd_cfg.handoff_min_orientation_coherence),
        "pd_handoff_min_axial_coverage": float(pd_cfg.handoff_min_axial_coverage),
        "pd_handoff_max_axial_gap_horizons": float(pd_cfg.handoff_max_axial_gap_horizons),
        "pd_handoff_min_slenderness": float(pd_cfg.handoff_min_slenderness),
        "pd_handoff_require_active_front": bool(pd_cfg.handoff_require_active_front),
        "pd_handoff_max_offfront_broken_fraction": float(pd_cfg.handoff_max_offfront_broken_fraction),
        "pd_connected_extent_final_m": float(final_connected_extent),
        "pd_crack_centerline_length_final_m": float(final_audit.centerline_length_m),
        "pd_crack_width_final_m": float(final_audit.width_m),
        "pd_crack_tip_width_final_m": float(final_audit.tip_width_m),
        "pd_crack_width_ratio_final": float(final_audit.width_ratio),
        "pd_crack_tip_radius_eff_final_m": float(final_audit.tip_radius_eff_m),
        "pd_crack_orientation_final_deg": float(final_audit.orientation_deg),
        "pd_crack_boundary_clearance_final_m": float(final_audit.boundary_clearance_m),
        "pd_crack_remote_K_final_MPam05": float(final_audit.remote_K_MPam05),
        "pd_crack_local_upper_K_final_MPam05": float(final_audit.local_upper_K_MPam05),
        "pd_crack_remote_KI_final_MPam05": float(final_audit.remote_KI_MPam05),
        "pd_crack_remote_KII_final_MPam05": float(final_audit.remote_KII_MPam05),
        "pd_crack_orientation_coherence_final": float(final_audit.orientation_coherence),
        "pd_crack_axial_coverage_final": float(final_audit.axial_coverage),
        "pd_crack_max_axial_gap_final_m": float(final_audit.max_axial_gap_m),
        "pd_crack_surface_gap_final_m": float(final_audit.surface_gap_m),
        "pd_crack_surface_contact_x_final_m": float(final_audit.surface_contact_x_m),
        "pd_crack_surface_contact_y_final_m": float(final_audit.surface_contact_y_m),
        "pd_crack_slenderness_final": float(final_audit.slenderness),
        "pd_handoff_length_pass_final": bool(final_audit.length_pass),
        "pd_handoff_slenderness_pass_final": bool(final_audit.slenderness_pass),
        "pd_handoff_tip_pass_final": bool(final_audit.tip_pass),
        "pd_handoff_boundary_pass_final": bool(final_audit.boundary_pass),
        "pd_handoff_K_pass_final": bool(final_audit.K_pass),
        "pd_handoff_orientation_pass_final": bool(final_audit.orientation_pass),
        "pd_handoff_coverage_pass_final": bool(final_audit.coverage_pass),
        "pd_handoff_surface_pass_final": bool(final_audit.surface_pass),
        "pd_handoff_front_pass_final": bool(final_audit.front_pass),
        "pd_handoff_competition_pass_final": bool(final_audit.competition_pass),
        "pd_handoff_offfront_broken_fraction_final": float(final_audit.off_front_broken_fraction),
        "pd_handoff_failure_reasons_final": list(final_audit.failure_reasons),
        "pd_handoff_pass_final": bool(final_audit.handoff_pass),
        "pd_connected_bonds_final": int(final_connected_bonds),
        "pd_active_front_final": bool(pd_state.active_front),
        "pd_active_front_length_final_m": float(pd_state.active_front_length_m),
        "pd_active_front_bonds_final": int(np.count_nonzero(pd_state.active_front_bonds)),
        "pd_front_backbone_bonds_final": int(np.count_nonzero(pd_state.front_backbone_bonds)),
        "pd_front_wake_bonds_final": int(np.count_nonzero(pd_state.front_wake_bonds)),
        "pd_front_process_bonds_final": int(np.count_nonzero(pd_state.front_process_bonds)),
        "pd_front_candidate_mode_final": int(pd_state.front_candidate_mode),
        "pd_front_eligible_preferred_final": int(pd_state.front_eligible_preferred),
        "pd_front_eligible_fallback_final": int(pd_state.front_eligible_fallback),
        "pd_front_stall_updates_final": int(pd_state.front_stall_updates),
        "pd_front_stalled_final": bool(pd_state.cycles_front_stalled is not None),
        "pd_front_max_link_rate_per_cycle_final": float(pd_state.front_max_link_rate_per_cycle),
        "pd_active_front_path_xy_final": np.asarray(pd_state.active_front_path_xy, float).tolist(),
        "pd_primary_seed_node_final": int(pd_state.primary_seed_node),
        "pd_primary_seed_stall_updates_final": int(pd_state.primary_seed_stall_updates),
        "pd_primary_seed_reselections_final": int(pd_state.primary_seed_reselections),
        "pd_primary_seed_rejected_nodes_final": int(np.count_nonzero(pd_state.primary_seed_rejected)),
        "pd_precapture_stalled_final": bool(pd_state.cycles_precapture_stalled is not None),
        "pd_broken_bonds_final": int(np.count_nonzero(pd_state.bond_damage >= pd_cfg.broken_damage)),
        "pd_off_front_broken_bonds_final": int(np.count_nonzero(
            (pd_state.bond_damage >= pd_cfg.broken_damage) & ~np.asarray(pd_state.active_front_bonds, bool)
        )),
        "pd_bond_damage_final_max": float(np.max(pd_state.bond_damage)),
        "pd_expected_candidate_sites": float(np.sum(patch.mean_candidate_sites)),
        "pd_expected_embryos_final": float(np.sum(pd_state.embryo * patch.mean_candidate_sites)),
        "pd_expected_births_cumulative_final": float(np.sum(pd_state.born_cumulative * patch.mean_candidate_sites)),
        "pd_expected_stable_final": float(np.sum(pd_state.stable * patch.mean_candidate_sites)),
        "pd_candidate_sites_realized": int(np.sum(pd_state.candidate_sites)),
        "pd_realized_births_cumulative_final": int(np.sum(pd_state.born_sites_cumulative)),
        "pd_realized_embryos_final": int(np.sum(pd_state.embryo_sites)),
        "pd_realized_stable_final": int(np.sum(pd_state.stable_sites)),
        "pd_delivery_source": str(args.delivery_source),
        "pd_delivery_scale": float(args.delivery_scale),
        "pd_delivery_source_multiplicity": float(
            getattr(args, "delivery_source_multiplicity", 1.0)
        ),
        "pd_delivery_event_strain": float(args.delivery_event_strain),
        "pd_delivery_weight_emit": float(args.delivery_weight_emit),
        "pd_delivery_weight_peierls": float(args.delivery_weight_peierls),
        "pd_delivery_weight_taylor": float(args.delivery_weight_taylor),
        "pd_delivery_rate_cap_s": float(args.delivery_rate_cap_s),
        "pd_delivery_hit_count": float(args.delivery_hit_count),
        "pd_delivery_memory_s": float(args.delivery_memory_s),
        "pd_delivery_memory_cycles": float(args.frequency_Hz * args.delivery_memory_s),
        "pd_temporal_multihit_scheme": "separated_delivery_plus_persistent_site_first_passage_plus_shielded_active_front_v8_3",
        "pd_hazard_n_phase": int(args.hazard_n_phase),
        "pd_amplification_damage_scale": float(pd_cfg.amplification_damage_scale),
        "pd_softening_damage_threshold": float(pd_cfg.softening_damage),
        "pd_stochastic_transition_probability_cap": float(args.max_transition_probability),
        "pd_rng_streams": "independent_candidate_and_event_streams_from_pd_seed",
        "root_radius_initial_m": root_radius0,
        "analytic_notch_root_radius_m": float(geom.root_radius),
        "root_radius_final_m": local_root_radius(mesh, feature_nodes),
        "root_radius_over_spacing_final": float(geometry_audit_final.get("root_radius_over_spacing", np.nan)),
        "minimum_element_area_over_initial_final": float(geometry_audit_final.get("minimum_element_area_over_initial", np.nan)),
        "geometry_policy": ("ale_sensitivity" if args.enable_geometry_evolution else "fixed_production_baseline"),
        "geometry_evolution_enabled": bool(args.enable_geometry_evolution),
        "birth_clock_mode": ("persistent_site_first_passage_aligned" if args.align_blocks_to_birth_clock else "persistent_site_first_passage_unaligned"),
        "front_link_state_shift_weight": float(args.front_link_state_shift_weight),
        "front_link_state_shift_scale_eV": float(args.front_link_state_shift_scale_eV),
        "geometry_limit_mode": args.geometry_limit_mode,
        "geometry_saturation_start_radius_spacing": float(args.geometry_saturation_start_radius_spacing),
        "geometry_min_valid_radius_spacing": float(args.geometry_min_valid_radius_spacing),
        "geometry_min_global_area_fraction": float(args.geometry_min_global_area_fraction),
        "geometry_saturated": bool(geometry_saturated),
        "geometry_saturation_cycles": geometry_saturation_cycles,
        "geometry_invalid_reason": geometry_invalid_reason,
        "geometry_resolution_audit_final": geometry_audit_final,
        "epsp_acc_final_max": float(np.max(epsp_acc_gp)),
        "rho_final_max_m2": float(np.max(rho_gp)),
        "residual_sigma1_final_max_Pa": float(np.max(last_residual)),
        "plastic_work_final_J_per_m": Wp_total,
        "history_csv": str(outdir / "sn_stateful_pd_history.csv"),
        "fatigue_endpoint": args.fatigue_endpoint,
        "stable_crack_birth_site_id": first_stable_site_id,
        "stable_crack_birth_node": (
            int(site_nodes_final[first_stable_site_id]) if first_stable_site_id is not None else None
        ),
        "stable_crack_birth_embryo_cycle": (
            float(pd_state.site_birth_cycle[first_stable_site_id]) if first_stable_site_id is not None else None
        ),
        "stable_crack_birth_cycle": pd_state.cycles_first_stable,
        "stable_crack_birth_transition_threshold": (
            float(pd_state.site_transition_threshold[first_stable_site_id]) if first_stable_site_id is not None else None
        ),
        "stable_crack_birth_transition_cumulative_hazard": (
            float(pd_state.site_transition_cumulative_hazard[first_stable_site_id]) if first_stable_site_id is not None else None
        ),
        "stable_crack_birth_outcome_uniform": (
            float(pd_state.site_transition_outcome_uniform[first_stable_site_id]) if first_stable_site_id is not None else None
        ),
        "stable_crack_birth_at_accepted_boundary": bool(
            pd_state.cycles_first_stable is not None
            and abs(float(cycles) - float(pd_state.cycles_first_stable))
            <= 16.0 * np.finfo(float).eps * max(1.0, abs(float(cycles)))
        ),
    }
    audit_payload = {
        k: (bool(v) if isinstance(v, (np.bool_, bool)) else float(v) if isinstance(v, (np.floating, float)) else int(v) if isinstance(v, (np.integer, int)) else v)
        for k, v in vars(final_audit).items()
    }
    with (outdir / "crack_handoff_audit_final.json").open("w") as f:
        json.dump(audit_payload, f, indent=2, sort_keys=True)
    with (outdir / "summary.json").open("w") as f:
        json.dump(summary, f, indent=2, sort_keys=True)
    return summary


def run_sweep(args):
    Path(args.out).mkdir(parents=True, exist_ok=True)
    summaries = []
    for case in args.cases:
        for s in args.sigma_a_MPa:
            sigma = float(s)
            case_dir = Path(args.out) / case / (f"sigmaA_{sigma:g}MPa".replace(".", "p"))
            summary_path = case_dir / "summary.json"
            if args.skip_existing and summary_path.exists():
                try:
                    cached = json.loads(summary_path.read_text())
                    complete_statuses = {
                        "physical_handoff",
                        "physical_handoff_geometry_saturated",
                        "morphology_invalid_diffuse",
                        "morphology_invalid_stalled",
                        "morphology_invalid_precapture_stalled",
                        "geometry_invalid_underresolved",
                        "right_censored",
                        "right_censored_geometry_saturated",
                        "stable_crack_birth",
                        "stable_spatial_crack_birth",
                    }
                    complete = (
                        cached.get("status") in complete_statuses
                        or float(cached.get("cycles_total", 0.0)) >= float(args.cycles_max)
                    )
                    same_model = cached.get("model") == MODEL_ID
                    same_source = cached.get("source_sha256") == SOURCE_SHA256
                    same_signature = cached.get("run_signature") == _checkpoint_signature(
                        args, case, sigma
                    )
                    if complete and same_model and same_source and same_signature:
                        print(
                            f"=== STATEFUL PD V8.3 skip existing case={case} "
                            f"sigma_a={sigma:g} MPa ==="
                        )
                        summaries.append(cached)
                        continue
                except Exception:
                    pass
            print(f"=== STATEFUL PD V8.3 case={case} sigma_a={sigma:g} MPa ===")
            summaries.append(run_case_stress(args, case, sigma))
    _write_csv(Path(args.out) / "sn_stateful_pd_summary.csv", summaries)
    return summaries


def build_parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", default="runs/sn_stateful_pd_v8_7_generalized_features")
    p.add_argument("--skip-existing", action="store_true", dest="skip_existing")
    p.add_argument("--resolution-profile", choices=["h15", "h10", "custom"], default="h15", dest="resolution_profile")
    p.add_argument(
        "--fatigue-model",
        choices=list(REPRESENTATIVE_FATIGUE_MODELS) + ["custom"],
        default="plastic_shielded_case64_M1",
        dest="fatigue_model",
        help="one of the six prior representative fatigue barrier sets",
    )
    p.add_argument("--cases", nargs="+", choices=["no_shield", "shielded"], default=["no_shield", "shielded"])
    p.add_argument("--T", type=float, default=300.0)
    p.add_argument("--sigma-a-MPa", nargs="+", type=float, default=[500, 600, 700], dest="sigma_a_MPa")
    p.add_argument("--R", type=float, default=0.1)
    p.add_argument("--frequency-Hz", type=float, default=1000.0, dest="frequency_Hz")
    p.add_argument("--cycles-max", type=float, default=1e9, dest="cycles_max")
    p.add_argument("--block-cycles", type=float, default=1e7, dest="block_cycles")
    p.add_argument("--min-block-cycles", type=float, default=1e-6, dest="min_block_cycles")
    p.add_argument("--max-blocks", type=int, default=3000, dest="max_blocks")
    p.add_argument("--pd-high-cycle", action="store_true", dest="pd_high_cycle",
                   help="enable the versioned dormant fixed-topology event-to-event engine")
    p.add_argument("--pd-high-cycle-max-segment", type=float, default=1e9, dest="pd_high_cycle_max_segment")
    p.add_argument("--pd-high-cycle-start-cycles", type=float, default=1e4,
                   dest="pd_high_cycle_start_cycles",
                   help="use the efficient transactional macro-stepper through the initial transient")
    p.add_argument("--pd-high-cycle-checkpoint-decades", action="store_true", default=True,
                   dest="pd_high_cycle_checkpoint_decades")
    p.add_argument("--target-dep-eq-block", type=float, default=2e-4, dest="target_dep_eq_block")
    p.add_argument("--target-rho-rel-block", type=float, default=0.05, dest="target_rho_rel_block")
    p.add_argument("--target-delivery-events", type=float, default=float("inf"), dest="target_delivery_events", help="optional cap on expected plastic-delivery events per adaptive block")
    p.add_argument("--max-transition-probability", type=float, default=0.08, dest="max_transition_probability")
    p.add_argument("--align-blocks-to-birth-clock", action="store_true", default=True, dest="align_blocks_to_birth_clock")
    p.add_argument("--no-align-blocks-to-birth-clock", action="store_false", dest="align_blocks_to_birth_clock")
    p.add_argument("--birth-clock-safety-factor", type=float, default=1.000001, dest="birth_clock_safety_factor")
    p.add_argument("--plastic-n-phase", type=int, default=12, dest="plastic_n_phase")
    p.add_argument("--hazard-n-phase", type=int, default=32, dest="hazard_n_phase")

    p.add_argument("--Lx", type=float, default=2e-3)
    p.add_argument("--Ly", type=float, default=4e-3)
    p.add_argument("--notch-depth-m", type=float, default=0.15e-3, dest="notch_depth_m")
    p.add_argument("--notch-half-height-m", type=float, default=0.30e-3, dest="notch_half_height_m")
    p.add_argument("--feature-type", choices=["ellipse", "rounded_v"], default="ellipse", dest="feature_type")
    p.add_argument("--notch-root-radius-m", type=float, default=None, dest="notch_root_radius_m")
    p.add_argument("--notch-opening-angle-deg", type=float, default=60.0, dest="notch_opening_angle_deg")
    p.add_argument("--path-refine-length-m", type=float, default=0.24e-3, dest="path_refine_length_m")
    p.add_argument("--path-refine-half-height-m", type=float, default=60e-6, dest="path_refine_half_height_m")
    p.add_argument("--nx", type=int, default=36)
    p.add_argument("--ny", type=int, default=72)
    p.add_argument("--jitter", type=float, default=0.08)
    p.add_argument("--root-h-fine", type=float, default=30e-6, dest="root_h_fine")
    p.add_argument("--seed", type=int, default=42)

    p.add_argument("--E-GPa", type=float, default=410.0, dest="E_GPa")
    p.add_argument("--nu", type=float, default=0.28)
    p.add_argument("--b-m", type=float, default=2.74e-10, dest="b_m")
    p.add_argument("--Tm-K", type=float, default=3695.0, dest="Tm_K")
    p.add_argument("--rho0", type=float, default=1e12)
    p.add_argument("--rho-floor", type=float, default=1e8, dest="rho_floor")
    p.add_argument("--rho-cap", type=float, default=1e17, dest="rho_cap")
    p.add_argument("--k-store", type=float, default=np.sqrt(2.0), dest="k_store")
    p.add_argument("--k-dyn", type=float, default=1.0, dest="k_dyn")

    p.add_argument("--exp-system", default="W[100]", dest="exp_system")
    p.add_argument("--exp-G00-eV", type=float, default=None, dest="exp_G00_eV")
    p.add_argument("--exp-gT-eV-per-K", type=float, default=None, dest="exp_gT_eV_per_K")
    p.add_argument("--exp-sigc0-GPa", type=float, default=None, dest="exp_sigc0_GPa")
    p.add_argument("--exp-sT-MPa-per-K", type=float, default=None, dest="exp_sT_MPa_per_K")
    p.add_argument("--exp-Tref-K", type=float, default=None, dest="exp_Tref_K")
    p.add_argument("--exp-floor-frac", type=float, default=None, dest="exp_floor_frac")
    p.add_argument("--exp-a", type=float, default=None, dest="exp_a")
    p.add_argument("--exp-n", type=float, default=None, dest="exp_n")
    p.add_argument("--emit-energy-scale", type=float, default=0.75, dest="emit_energy_scale")
    p.add_argument("--emit-entropy-scale", type=float, default=0.75, dest="emit_entropy_scale")
    p.add_argument("--emit-stress-scale", type=float, default=1.0, dest="emit_stress_scale")
    p.add_argument("--peierls-energy-scale", type=float, default=0.00375, dest="peierls_energy_scale")
    p.add_argument("--peierls-entropy-scale", type=float, default=0.00375, dest="peierls_entropy_scale")
    p.add_argument("--peierls-stress-scale", type=float, default=1.0, dest="peierls_stress_scale")
    p.add_argument("--taylor-energy-scale", type=float, default=0.015, dest="taylor_energy_scale")
    p.add_argument("--taylor-entropy-scale", type=float, default=0.015, dest="taylor_entropy_scale")
    p.add_argument("--taylor-stress-scale", type=float, default=1.0, dest="taylor_stress_scale")
    p.add_argument("--nu0-emit-pz", type=float, default=1e11, dest="nu0_emit_pz")
    p.add_argument("--nu0-peierls", type=float, default=1e11, dest="nu0_peierls")
    p.add_argument("--nu0-taylor", type=float, default=1e11, dest="nu0_taylor")
    p.add_argument("--plastic-event-strain", type=float, default=1e-5, dest="plastic_event_strain")
    p.add_argument("--phi-taylor-max", type=float, default=20.0, dest="phi_taylor_max")
    p.add_argument("--max-dep-phase", type=float, default=2e-5, dest="max_dep_phase")
    p.add_argument("--max-rho-rel-phase", type=float, default=0.02, dest="max_rho_rel_phase")

    p.add_argument("--epsp-shield-scale", type=float, default=5e-3, dest="epsp_shield_scale")
    p.add_argument("--epsp-damage-scale", type=float, default=2e-2, dest="epsp_damage_scale")
    p.add_argument("--crack-G00-eV", type=float, default=1.0, dest="crack_G00_eV")
    p.add_argument("--crack-sigc0-GPa", type=float, default=3.0, dest="crack_sigc0_GPa")
    p.add_argument("--crack-exp-a", type=float, default=0.70, dest="crack_exp_a")
    p.add_argument("--crack-exp-n", type=float, default=0.60, dest="crack_exp_n")
    p.add_argument("--crack-floor-frac", type=float, default=0.010, dest="crack_floor_frac")
    p.add_argument("--crack-T-mode", choices=["linear", "mu_scale", "audited_linear"], default="mu_scale", dest="crack_T_mode")
    p.add_argument("--crack-Tref-K", type=float, default=481.33, dest="crack_Tref_K")
    p.add_argument("--crack-gT-eV-per-K", type=float, default=None, dest="crack_gT_eV_per_K")
    p.add_argument("--crack-sT-GPa-per-K", type=float, default=None, dest="crack_sT_GPa_per_K")
    p.add_argument("--crack-mu-dlnmu-dT-per-K", type=float, default=-1.5e-4, dest="crack_mu_dlnmu_dT_per_K")
    p.add_argument("--crack-G0-mu-power", type=float, default=1.0, dest="crack_G0_mu_power")
    p.add_argument("--crack-sigc-mu-power", type=float, default=1.0, dest="crack_sigc_mu_power")
    p.add_argument("--nu0-crack", type=float, default=1e11, dest="nu0_crack")
    p.add_argument("--S-crack-kB", type=float, default=0.0, dest="S_crack_kB")
    p.add_argument("--sigma-back-max-GPa", type=float, default=1.0, dest="sigma_back_max_GPa")
    p.add_argument("--shield-chi", type=float, default=0.6, dest="shield_chi")
    p.add_argument("--Gshield-eV", type=float, default=0.35, dest="Gshield_eV")
    p.add_argument("--Gstored-eV", type=float, default=0.25, dest="Gstored_eV")

    p.add_argument("--enable-geometry-evolution", action="store_true", default=False, dest="enable_geometry_evolution", help="explicit ALE sensitivity; fixed geometry is the v8.3 production baseline")
    p.add_argument("--disable-geometry-evolution", action="store_false", dest="enable_geometry_evolution")
    p.add_argument("--freeze-ale-after-front-capture", action="store_true", default=True, dest="freeze_ale_after_front_capture")
    p.add_argument("--continue-ale-after-front-capture", action="store_false", dest="freeze_ale_after_front_capture")
    p.add_argument("--morph-band-length-m", type=float, default=100e-6, dest="morph_band_length_m")
    p.add_argument("--morph-normal-weight", type=float, default=0.25, dest="morph_normal_weight")
    p.add_argument("--morph-shear-weight", type=float, default=1.0, dest="morph_shear_weight")
    p.add_argument("--morph-decay-length-m", type=float, default=150e-6, dest="morph_decay_length_m")
    p.add_argument("--target-surface-move-fraction", type=float, default=0.10, dest="target_surface_move_fraction")
    p.add_argument("--max-surface-move-fraction", type=float, default=0.20, dest="max_surface_move_fraction")
    p.add_argument("--min-area-fraction", type=float, default=0.15, dest="min_area_fraction")
    p.add_argument(
        "--geometry-limit-mode",
        choices=["terminate", "freeze", "off"],
        default="terminate",
        dest="geometry_limit_mode",
        help="handling when ALE root curvature approaches the local resolution limit",
    )
    p.add_argument("--geometry-saturation-start-radius-spacing", type=float, default=8.0, dest="geometry_saturation_start_radius_spacing")
    p.add_argument("--geometry-min-valid-radius-spacing", type=float, default=5.0, dest="geometry_min_valid_radius_spacing")
    p.add_argument("--geometry-min-global-area-fraction", type=float, default=0.05, dest="geometry_min_global_area_fraction")

    p.add_argument("--pd-patch-radius-m", type=float, default=0.45e-3, dest="pd_patch_radius_m")
    p.add_argument("--pd-horizon-m", type=float, default=90e-6, dest="pd_horizon_m")
    p.add_argument("--pd-boundary-shell-m", type=float, default=100e-6, dest="pd_boundary_shell_m")
    p.add_argument("--pd-residual-stiffness", type=float, default=1e-7, dest="pd_residual_stiffness")
    p.add_argument("--pd-initiation-radius-m", type=float, default=240e-6, dest="pd_initiation_radius_m")
    p.add_argument("--pd-initiation-taper-m", type=float, default=60e-6, dest="pd_initiation_taper_m")
    p.add_argument("--pd-initiation-back-extent-m", type=float, default=60e-6, dest="pd_initiation_back_extent_m")
    p.add_argument("--pd-amplification-cap", type=float, default=4.0, dest="pd_amplification_cap")
    p.add_argument("--pd-amplification-damage-scale", type=float, default=0.05, dest="pd_amplification_damage_scale")
    p.add_argument("--pd-seed", type=int, default=None, dest="pd_seed",
                   help="seed for the discrete candidate-site realization; defaults to --seed")
    p.add_argument("--site-density-m2", type=float, default=5e10, dest="site_density_m2")
    p.add_argument(
        "--delivery-source",
        choices=["completed_flow", "emission", "plastic_strain_rate", "weighted_mechanisms"],
        default="completed_flow",
        dest="delivery_source",
        help="plastic-event process used to populate the finite-memory delivery clock",
    )
    p.add_argument("--delivery-scale", type=float, default=1.0, dest="delivery_scale")
    p.add_argument("--delivery-source-multiplicity", type=float, default=1.0,
                   dest="delivery_source_multiplicity")
    p.add_argument(
        "--delivery-event-strain",
        type=float,
        default=1e-5,
        dest="delivery_event_strain",
        help="strain increment represented by one delivery event when delivery-source=plastic_strain_rate",
    )
    p.add_argument("--delivery-weight-emit", type=float, default=1.0, dest="delivery_weight_emit")
    p.add_argument("--delivery-weight-peierls", type=float, default=0.0, dest="delivery_weight_peierls")
    p.add_argument("--delivery-weight-taylor", type=float, default=0.0, dest="delivery_weight_taylor")
    p.add_argument("--delivery-rate-cap-s", type=float, default=float("inf"), dest="delivery_rate_cap_s")
    p.add_argument("--delivery-hit-count", type=float, default=2.0, dest="delivery_hit_count")
    p.add_argument("--delivery-memory-s", type=float, default=1e-3, dest="delivery_memory_s")
    p.add_argument("--birth-scale", type=float, default=1.0, dest="birth_scale")
    p.add_argument("--nu-stabilize-s", type=float, default=5e2, dest="nu_stabilize_s")
    p.add_argument("--nu-heal-s", type=float, default=2e2, dest="nu_heal_s")
    p.add_argument("--stabilize-stress-GPa", type=float, default=1.4, dest="stabilize_stress_GPa")
    p.add_argument("--stabilize-width-GPa", type=float, default=0.25, dest="stabilize_width_GPa")
    p.add_argument("--stabilize-plastic-gain", type=float, default=2.0, dest="stabilize_plastic_gain")
    p.add_argument("--heal-return-fraction", type=float, default=0.9, dest="heal_return_fraction")
    p.add_argument("--nu-grow-s", type=float, default=2e-2, dest="nu_grow_s")
    p.add_argument("--grow-stress-GPa", type=float, default=1.2, dest="grow_stress_GPa")
    p.add_argument("--grow-width-GPa", type=float, default=0.25, dest="grow_width_GPa")
    p.add_argument("--stable-count-scale", type=float, default=2.0, dest="stable_count_scale")
    p.add_argument("--nu-link-s", type=float, default=8e-3, dest="nu_link_s")
    p.add_argument("--link-stress-GPa", type=float, default=1.1, dest="link_stress_GPa")
    p.add_argument("--link-width-GPa", type=float, default=0.20, dest="link_width_GPa")
    p.add_argument("--link-orientation-power", type=float, default=8.0, dest="link_orientation_power")
    p.add_argument("--link-orientation-floor", type=float, default=0.0, dest="link_orientation_floor")
    p.add_argument("--directional-band-horizons", type=float, default=0.40, dest="directional_band_horizons")
    p.add_argument("--seed-influence-horizons", type=float, default=1.25, dest="seed_influence_horizons")
    p.add_argument("--front-neighbor-spacing-factor", type=float, default=1.75, dest="front_neighbor_spacing_factor")
    p.add_argument("--front-orientation-tolerance-deg", type=float, default=22.5, dest="front_orientation_tolerance_deg")
    p.add_argument("--neighbor-link-gain", type=float, default=2.0, dest="neighbor_link_gain")
    p.add_argument("--front-capture-enabled", action="store_true", default=True, dest="front_capture_enabled")
    p.add_argument("--disable-front-capture", action="store_false", dest="front_capture_enabled")
    p.add_argument("--front-capture-min-bonds", type=int, default=6, dest="front_capture_min_bonds")
    p.add_argument("--front-capture-min-length-horizons", type=float, default=0.50, dest="front_capture_min_length_horizons")
    p.add_argument("--front-capture-min-orientation-coherence", type=float, default=0.55, dest="front_capture_min_orientation_coherence")
    p.add_argument("--front-capture-max-surface-gap-horizons", type=float, default=1.0, dest="front_capture_max_surface_gap_horizons")
    p.add_argument("--front-process-ahead-horizons", type=float, default=1.50, dest="front_process_ahead_horizons")
    p.add_argument("--front-process-behind-horizons", type=float, default=0.35, dest="front_process_behind_horizons")
    p.add_argument("--front-band-horizons", type=float, default=0.30, dest="front_band_horizons")
    p.add_argument("--front-wake-band-horizons", type=float, default=0.25, dest="front_wake_band_horizons")
    p.add_argument("--front-tip-seed-gain", type=float, default=1.0, dest="front_tip_seed_gain")
    p.add_argument("--front-preferred-orientation-tolerance-deg", type=float, default=25.0, dest="front_preferred_orientation_tolerance_deg")
    p.add_argument("--front-fallback-orientation-tolerance-deg", type=float, default=42.0, dest="front_fallback_orientation_tolerance_deg")
    p.add_argument("--front-preferred-band-horizons", type=float, default=0.35, dest="front_preferred_band_horizons")
    p.add_argument("--front-fallback-band-horizons", type=float, default=0.70, dest="front_fallback_band_horizons")
    p.add_argument("--front-fallback-activity-scale", type=float, default=0.35, dest="front_fallback_activity_scale")
    p.add_argument("--front-direction-smoothing", type=float, default=0.50, dest="front_direction_smoothing")
    p.add_argument("--front-max-turn-deg", type=float, default=18.0, dest="front_max_turn_deg")
    p.add_argument("--front-recent-segment-horizons", type=float, default=1.25, dest="front_recent_segment_horizons")
    p.add_argument("--front-crack-tube-horizons", type=float, default=0.80, dest="front_crack_tube_horizons")
    p.add_argument("--front-backbone-band-horizons", type=float, default=0.40, dest="front_backbone_band_horizons")
    p.add_argument("--front-min-advance-spacing-factor", type=float, default=0.25, dest="front_min_advance_spacing_factor")
    p.add_argument("--front-stall-enabled", action="store_true", default=True, dest="front_stall_enabled")
    p.add_argument("--disable-front-stall-audit", action="store_false", dest="front_stall_enabled")
    p.add_argument("--front-stall-patience-updates", type=int, default=12, dest="front_stall_patience_updates")
    p.add_argument("--primary-seed-reselection-enabled", action="store_true", default=True, dest="primary_seed_reselection_enabled")
    p.add_argument("--disable-primary-seed-reselection", action="store_false", dest="primary_seed_reselection_enabled")
    p.add_argument("--primary-seed-reselection-patience-updates", type=int, default=24, dest="primary_seed_reselection_patience_updates")
    p.add_argument("--primary-seed-progress-damage-increment", type=float, default=0.02, dest="primary_seed_progress_damage_increment")
    p.add_argument("--primary-seed-max-reselections", type=int, default=64, dest="primary_seed_max_reselections")
    p.add_argument("--primary-seed-surface-score-horizons", type=float, default=2.0, dest="primary_seed_surface_score_horizons")
    p.add_argument("--primary-seed-stress-score-weight", type=float, default=1.0, dest="primary_seed_stress_score_weight")
    p.add_argument("--primary-seed-population-score-weight", type=float, default=0.35, dest="primary_seed_population_score_weight")
    p.add_argument("--primary-seed-damage-score-weight", type=float, default=1.5, dest="primary_seed_damage_score_weight")
    p.add_argument("--pre-capture-birth-scale-with-primary", type=float, default=0.05, dest="pre_capture_birth_scale_with_primary")
    p.add_argument("--post-capture-birth-scale", type=float, default=0.0, dest="post_capture_birth_scale")
    p.add_argument("--post-capture-stabilization-scale", type=float, default=0.0, dest="post_capture_stabilization_scale")
    p.add_argument("--off-front-growth-scale", type=float, default=0.0, dest="off_front_growth_scale")
    p.add_argument("--front-amplification-ahead-horizons", type=float, default=1.50, dest="front_amplification_ahead_horizons")
    p.add_argument("--front-amplification-behind-horizons", type=float, default=0.25, dest="front_amplification_behind_horizons")
    p.add_argument("--front-amplification-band-horizons", type=float, default=1.00, dest="front_amplification_band_horizons")
    p.add_argument("--front-amplification-cap", type=float, default=3.0, dest="front_amplification_cap")
    p.add_argument("--front-link-state-shift-weight", type=float, default=1.0, dest="front_link_state_shift_weight")
    p.add_argument("--front-link-state-shift-scale-eV", type=float, default=0.35, dest="front_link_state_shift_scale_eV")
    p.add_argument("--front-link-state-shift-z-clip", type=float, default=3.0, dest="front_link_state_shift_z_clip")
    p.add_argument("--diffuse-abort-enabled", action="store_true", default=True, dest="diffuse_abort_enabled")
    p.add_argument("--disable-diffuse-abort", action="store_false", dest="diffuse_abort_enabled")
    p.add_argument("--diffuse-abort-min-bonds", type=int, default=200, dest="diffuse_abort_min_bonds")
    p.add_argument("--diffuse-abort-min-offfront-fraction", type=float, default=0.50, dest="diffuse_abort_min_offfront_fraction")
    p.add_argument("--diffuse-abort-patience-updates", type=int, default=5, dest="diffuse_abort_patience_updates")
    p.add_argument("--surface-connection-horizons", type=float, default=1.0, dest="surface_connection_horizons")
    p.add_argument("--topology-neighbor-spacing-factor", type=float, default=1.75, dest="topology_neighbor_spacing_factor")
    p.add_argument("--topology-orientation-tolerance-deg", type=float, default=30.0, dest="topology_orientation_tolerance_deg")
    p.add_argument("--handoff-min-orientation-coherence", type=float, default=0.65, dest="handoff_min_orientation_coherence")
    p.add_argument("--handoff-min-axial-coverage", type=float, default=0.70, dest="handoff_min_axial_coverage")
    p.add_argument("--handoff-max-axial-gap-horizons", type=float, default=0.75, dest="handoff_max_axial_gap_horizons")
    p.add_argument("--handoff-min-slenderness", type=float, default=5.0, dest="handoff_min_slenderness")
    p.add_argument("--handoff-require-active-front", action="store_true", default=True, dest="handoff_require_active_front")
    p.add_argument("--handoff-allow-uncaptured-front", action="store_false", dest="handoff_require_active_front")
    p.add_argument("--handoff-max-offfront-broken-fraction", type=float, default=0.25, dest="handoff_max_offfront_broken_fraction")
    p.add_argument("--broken-damage", type=float, default=0.95, dest="broken_damage")
    p.add_argument("--softening-damage", type=float, default=1e-3, dest="softening_damage",
                   help="cohesive-damage level used to record first measurable softening")
    p.add_argument("--root-seed-radius-m", type=float, default=120e-6, dest="root_seed_radius_m")
    p.add_argument("--established-extent-m", type=float, default=240e-6, dest="established_extent_m", help="legacy length-only criterion used only with --handoff-mode legacy_extent")
    p.add_argument("--handoff-mode", choices=["physical", "legacy_extent", "disabled"], default="physical", dest="handoff_mode")
    p.add_argument("--handoff-min-length-m", type=float, default=0.0, dest="handoff_min_length_m")
    p.add_argument("--handoff-min-length-horizons", type=float, default=2.0, dest="handoff_min_length_horizons")
    p.add_argument("--handoff-max-width-ratio", type=float, default=0.30, dest="handoff_max_width_ratio")
    p.add_argument("--handoff-max-tip-width-horizons", type=float, default=1.0, dest="handoff_max_tip_width_horizons")
    p.add_argument("--handoff-tip-window-horizons", type=float, default=1.0, dest="handoff_tip_window_horizons")
    p.add_argument("--handoff-width-root-exclusion-horizons", type=float, default=0.5, dest="handoff_width_root_exclusion_horizons")
    p.add_argument("--handoff-min-boundary-clearance-horizons", type=float, default=3.0, dest="handoff_min_boundary_clearance_horizons")
    p.add_argument("--handoff-min-connected-bonds", type=int, default=3, dest="handoff_min_connected_bonds")
    p.add_argument("--handoff-edge-geometry-factor", type=float, default=1.12, dest="handoff_edge_geometry_factor")
    p.add_argument("--handoff-min-remote-K-MPam05", type=float, default=0.0, dest="handoff_min_remote_K_MPam05")

    p.add_argument("--checkpoint-every-blocks", type=int, default=25, dest="checkpoint_every_blocks")
    p.add_argument("--checkpoint-path", default="", dest="checkpoint_path")
    p.add_argument("--resume", action="store_true", dest="resume")
    p.add_argument("--stress-step-source-checkpoint", default="", dest="stress_step_source_checkpoint")
    p.add_argument("--stress-step-source-generation", default="", dest="stress_step_source_generation")
    p.add_argument("--stress-step-source-sigma-a-MPa", type=float, default=None,
                   dest="stress_step_source_sigma_a_MPa")
    p.add_argument("--protocol-label", default="", dest="protocol_label")
    p.add_argument("--snapshot-every", type=int, default=25, dest="snapshot_every")
    p.add_argument("--pd-image-policy", choices=("none", "event_only", "selected"),
                   default="selected", dest="pd_image_policy")
    p.add_argument("--print-every", type=int, default=1, dest="print_every")
    p.add_argument(
        "--fatigue-endpoint", choices=("physical_handoff", "stable_crack_birth", "stable_spatial_crack_birth"),
        default="physical_handoff", dest="fatigue_endpoint",
    )
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    apply_resolution_profile(args)
    apply_representative_fatigue_model(args)
    print(f"STATEFUL_PD_V8_7 fatigue model: {args.fatigue_model} (preset={args.fatigue_model_preset_applied})")
    run_sweep(args)


if __name__ == "__main__":
    main()
