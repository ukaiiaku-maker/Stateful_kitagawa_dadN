import copy
import numpy as np

from arrhenius_fracture.config import ElasticProperties
from arrhenius_fracture.sn_arrhenius_chain import build_chain_from_namespace
from arrhenius_fracture.sn_feature_geometry_v8_7 import (
    BluntNotchGeometry, identify_feature_surface_nodes, local_root_xy,
    make_blunt_edge_notch_mesh,
)
from arrhenius_fracture.sn_intact_fem import plane_strain_D
from arrhenius_fracture.sn_pd2d_stateful_v9_transactional import (
    _pd_config_from_args, apply_representative_fatigue_model, build_crack_barrier,
    build_parser, evaluate_dormant_exact_cycle,
)
from arrhenius_fracture.v9_cached_fem import CachedIntactFEM
from arrhenius_fracture.v9_fem_transaction import EmbeddedFEMTransaction
from arrhenius_fracture.v9_stateful_peridynamics import V9StatefulPDPatch
from arrhenius_fracture.v9_pd_high_cycle import private_cycle
from arrhenius_fracture.v9_pd_high_cycle_adapter import SpatialPDDormantAdapter


def test_real_callback_uses_fem_pd_path_without_consuming_physical_state():
    args = build_parser().parse_args([
        "--resolution-profile", "custom", "--nx", "14", "--ny", "28",
        "--jitter", "0", "--root-h-fine", "2.5e-5",
        "--pd-horizon-m", "1.05e-4", "--pd-patch-radius-m", "4.6e-4",
        "--pd-boundary-shell-m", "1e-4", "--plastic-n-phase", "8",
        "--hazard-n-phase", "12", "--fatigue-model", "plastic_shielded_case64_M1",
    ])
    apply_representative_fatigue_model(args)
    mat = ElasticProperties(E=args.E_GPa*1e9, nu=args.nu, b=args.b_m, Tm=args.Tm_K)
    Dmat = plane_strain_D(mat)
    geom = BluntNotchGeometry(args.Lx,args.Ly,args.notch_depth_m,args.notch_half_height_m,
        feature_type=args.feature_type,root_radius_m=args.notch_root_radius_m,
        opening_angle_deg=args.notch_opening_angle_deg,
        path_refine_length_m=args.path_refine_length_m,path_refine_half_height_m=args.path_refine_half_height_m)
    mesh,bnd,_=make_blunt_edge_notch_mesh(geom,nx=args.nx,ny=args.ny,jitter=args.jitter,
        root_h_fine=args.root_h_fine,seed=args.seed)
    feature=identify_feature_surface_nodes(mesh,geom); root=local_root_xy(mesh,feature)
    patch=V9StatefulPDPatch(mesh,geom,root,mat,_pd_config_from_args(args),feature_surface_global_nodes=feature)
    state=patch.initial_state(); state_before=copy.deepcopy(state)
    candidate_rng=copy.deepcopy(patch._candidate_rng.bit_generator.state); event_rng=copy.deepcopy(patch._event_rng.bit_generator.state)
    chain=build_chain_from_namespace(args,mat.b); crack=build_crack_barrier(args)
    cached=CachedIntactFEM(mesh,bnd,mat,Dmat); sigma_max=2*700e6/(1-args.R); sigma_min=args.R*sigma_max
    transaction=EmbeddedFEMTransaction(mesh=mesh,boundaries=bnd,material=mat,Dmat=Dmat,
        plastic_chain=chain,args=args,sigma_max_Pa=sigma_max,sigma_min_Pa=sigma_min,cached_fem=cached)
    ep=np.zeros((3,mesh.ne)); rho=np.full(mesh.ne,args.rho0); eps=np.zeros(mesh.ne); u=np.zeros(mesh.ndof)
    result=evaluate_dormant_exact_cycle(args=args,shield_on=True,mesh=mesh,patch=patch,pd_state=state,
        crack=crack,plast_chain=chain,cached_fem=cached,fem_transaction=transaction,
        sigma_max=sigma_max,sigma_min=sigma_min,ep_gp=ep,rho_gp=rho,epsp_acc_gp=eps,u=u,
        plastic_work=0.0,cycles=0.0)
    assert result["ep_gp"].shape==ep.shape and result["log_birth_action"].shape==state.birth_cumulative_hazard.shape
    assert result["root_phase_tensors_Pa"].shape==(args.hazard_n_phase,2,2)
    assert np.all(np.isfinite(result["rho_gp"])) and result["ledger_increments"]["plastic_work"] >= 0.0
    np.testing.assert_array_equal(state.site_status,state_before.site_status)
    np.testing.assert_array_equal(state.birth_cumulative_hazard,state_before.birth_cumulative_hazard)
    np.testing.assert_array_equal(state.bond_damage,state_before.bond_damage)
    assert patch._candidate_rng.bit_generator.state==candidate_rng
    assert patch._event_rng.bit_generator.state==event_rng

    window=evaluate_dormant_exact_cycle(args=args,shield_on=True,mesh=mesh,patch=patch,pd_state=state,
        crack=crack,plast_chain=chain,cached_fem=cached,fem_transaction=transaction,
        sigma_max=sigma_max,sigma_min=sigma_min,ep_gp=ep,rho_gp=rho,epsp_acc_gp=eps,u=u,
        plastic_work=0.0,cycles=0.0,dN=2.0)
    assert window["diagnostics"]["private_window_cycles"] == 2.0
    assert np.all(np.isfinite(window["rho_gp"]))
    np.testing.assert_array_equal(state.site_status,state_before.site_status)
    np.testing.assert_array_equal(state.birth_cumulative_hazard,state_before.birth_cumulative_hazard)
    assert patch._candidate_rng.bit_generator.state==candidate_rng
    assert patch._event_rng.bit_generator.state==event_rng

    def evaluator(adapter):
        payload=evaluate_dormant_exact_cycle(args=args,shield_on=True,mesh=adapter.mesh,patch=adapter.patch,
            pd_state=adapter.pd_state,crack=crack,plast_chain=chain,cached_fem=cached,
            fem_transaction=transaction,sigma_max=sigma_max,sigma_min=sigma_min,
            ep_gp=adapter.ep_gp,rho_gp=adapter.rho_gp,epsp_acc_gp=adapter.epsp_acc_gp,u=adapter.u,
            plastic_work=adapter.plastic_work,cycles=adapter.cycles)
        adapter.ep_gp=payload.pop("ep_gp"); adapter.rho_gp=payload.pop("rho_gp")
        adapter.epsp_acc_gp=payload.pop("epsp_acc_gp"); adapter.u=payload.pop("u")
        adapter.pd_state.log_delivery_memory=payload.pop("log_delivery_memory")
        adapter.pd_state.delivery_memory=np.where(np.isfinite(adapter.pd_state.log_delivery_memory),np.exp(adapter.pd_state.log_delivery_memory),0.0)
        for field in ("available","embryo","stable","inactive","completion"):
            setattr(adapter.pd_state,field,payload.pop(field))
        return payload
    adapter=SpatialPDDormantAdapter(patch=patch,pd_state=state,mesh=mesh,ep_gp=ep,rho_gp=rho,
        epsp_acc_gp=eps,u=u,cycles=0.0,plastic_work=0.0,cycle_evaluator=evaluator)
    signatures=adapter.protected_signatures(); evaluated=private_cycle(adapter)
    assert evaluated.state_end.vector.shape==evaluated.state_start.vector.shape
    assert adapter.protected_signatures()==signatures and adapter.cycles==0.0
