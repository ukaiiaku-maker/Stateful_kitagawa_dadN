"""v9 adapter for the immutable audited signed persistent-site MPZ closure.

The constitutive implementation is loaded unchanged under an isolated package
name.  v9 owns transaction/capsule semantics; the audited source owns every
local pre-stable rate and state update.
"""
from __future__ import annotations

import copy
import hashlib
import importlib
import json
from pathlib import Path
import sys
import types

import numpy as np

from .config import EV_TO_J, KB

from .v9_four_class_registry import select_canonical_option


AUDITED_PACKAGE = "_v9_audited_signed_mpz_v10221"
FAMILY_RELATIVE = "runs/v10_2_28_kernel_cache/1447653d199f0b43cb475951092d69444c9b785f6fdf518c723792abb3b1f5e5/family.json"
MODULES = (
    "material_manifest", "unified_mpz", "signed_kernel_family_v10214",
    "signed_burgers_shared_v1025", "persistent_site_source_v10221",
    "reduced_campaign_v1024", "campaign_calibrated_tip",
    "anisotropic_emission_v10174",
)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_audited_modules(source_root: str | Path):
    root = Path(source_root).expanduser().resolve()
    package_root = root / "arrhenius_fracture"
    package = sys.modules.get(AUDITED_PACKAGE)
    if package is None:
        package = types.ModuleType(AUDITED_PACKAGE)
        package.__path__ = [str(package_root)]
        package.__package__ = AUDITED_PACKAGE
        sys.modules[AUDITED_PACKAGE] = package
    loaded = {
        name: importlib.import_module(f"{AUDITED_PACKAGE}.{name}")
        for name in MODULES
    }
    paths = {name: Path(module.__file__).resolve() for name, module in loaded.items()}
    family = root / FAMILY_RELATIVE
    paths["signed_kernel_family_json"] = family
    return loaded, {name: _sha256(path) for name, path in paths.items()}, paths


class SignedMPZPreBirthState:
    """Transaction-copyable authoritative signed MPZ state and diagnostics."""

    def __init__(
        self, option_id: str, source_root: str | Path, *,
        shear_modulus_Pa: float, poisson: float, burgers_m: float,
        initial_tip_radius_m: float,
    ):
        selected, registry_audit = select_canonical_option(option_id, source_root)
        modules, hashes, paths = load_audited_modules(source_root)
        row = selected.row
        manifest = modules["reduced_campaign_v1024"].manifest_from_row(row)
        length_m = float(row["L_pz_um_recommended"]) * 1e-6
        n_bins = int(float(row["n_bins_recommended"]))
        cfg = modules["unified_mpz"].MPZConfig(
            length_m=length_m,
            n_bins=n_bins,
            n_systems=int(float(row["n_slip_channels"])),
            source_bin_count=max(1, int(np.ceil(2e-6 / (length_m / n_bins)))),
            mobile_shield_fraction=float(row["mobile_shield_fraction"]),
            forest_density_floor_m2=float(row["rho_forest_floor_m2"]),
            peierls_stress_fraction=float(row["peierls_stress_fraction"]),
            taylor_stress_fraction=float(row["taylor_stress_fraction"]),
            mobile_recovery_rate_s=0.0,
            wake_shielding=False,
        )
        state = modules["unified_mpz"].UnifiedMPZState(manifest, cfg)
        family = modules["signed_kernel_family_v10214"].ActiveOnlySigned2DShieldingKernelFamily.from_json(paths["signed_kernel_family_json"])
        family = family.bind_to_state_grid(state)
        modules["signed_burgers_shared_v1025"].install_signed_burgers_population(
            state, family, modules["signed_burgers_shared_v1025"].VALIDATED_SCALAR_TRANSPORT,
        )
        family.resolve(r_eff_over_r0=1.0, opening_strength_fraction=0.0, crack_extension_m=0.0)
        state._signed_kernel = family
        state._campaign_G_Pa = float(shear_modulus_Pa)
        state._campaign_b = float(burgers_m)
        state._campaign_backstress_scale = 1.0
        state._anisotropic_drive_reliable = True
        state._anisotropic_drive_factors = np.ones(state.n_systems)
        state._anisotropic_tau_signed_Pa = np.ones(state.n_systems)
        persistent_cfg = modules["persistent_site_source_v10221"].PersistentSiteConfig(
            rho_site0_m2=float(row["rho_source0_m2"]),
            reference_source_area_m2=float(row["reference_source_area_um2"]) * 1e-12,
            reference_front_width_m=float(row["reference_front_width_um"]) * 1e-6,
            reference_density_m2=float(row["rho_forest_floor_m2"]),
            source_zone_length_m=float(row["source_zone_length_um"]) * 1e-6,
            minimum_front_width_m=0.0,
            maximum_front_width_m=float(row["L_pz_um_recommended"]) * 1e-6,
            implicit_tolerance=1e-30,
            implicit_max_iterations=512,
        )
        modules["persistent_site_source_v10221"].install_persistent_site_source(
            state, config=persistent_cfg, r0_m=initial_tip_radius_m, b_m=burgers_m,
        )
        self.state = state
        self.modules = modules
        self.option_id = option_id
        self.poisson = float(poisson)
        self.burgers_m = float(burgers_m)
        self.audit = registry_audit | {
            "constitutive_source_hashes": hashes,
            "signed_kernel_family_path": str(paths["signed_kernel_family_json"]),
            "port_mode": "direct_immutable_audited_executable_closure",
            "emission_model": "audited_v10221_aggregate_persistent_emission",
            "emission_stochastic": False,
            "cleavage_first_passage_stochastic": True,
            "fatigue_endpoint": "stable_crack_birth",
        }

    def copy(self):
        return copy.deepcopy(self)

    _CAPSULE_ARRAYS = (
        "mobile_positive", "mobile_negative", "retained_positive", "retained_negative",
        "accumulated_slip_positive", "accumulated_slip_negative",
        "wake_mobile_positive", "wake_mobile_negative",
        "wake_retained_positive", "wake_retained_negative",
        "wake_slip_positive", "wake_slip_negative",
    )
    _CAPSULE_SCALARS = (
        "emitted_total", "escaped_total", "recovered_total", "advance_total_m",
        "wake_discarded_mobile_total", "wake_discarded_retained_total",
        "wake_discarded_slip_total", "time_s", "signed_source_activations_total",
        "signed_line_content_emitted_total",
    )

    def capsule(self):
        return {
            "schema": "V9_SIGNED_MPZ_PREBIRTH_CAPSULE_1",
            "option_id": self.option_id,
            "audit": copy.deepcopy(self.audit),
            "arrays": {name: np.asarray(getattr(self.state, name)).copy() for name in self._CAPSULE_ARRAYS},
            "scalars": {name: float(getattr(self.state, name)) for name in self._CAPSULE_SCALARS},
        }

    def restore_capsule(self, capsule):
        if capsule.get("schema") != "V9_SIGNED_MPZ_PREBIRTH_CAPSULE_1":
            raise RuntimeError("unsupported signed-MPZ capsule schema")
        if capsule.get("option_id") != self.option_id or capsule.get("audit") != self.audit:
            raise RuntimeError("signed-MPZ capsule provenance/request mismatch")
        for name in self._CAPSULE_ARRAYS:
            target = np.asarray(getattr(self.state, name))
            value = np.asarray(capsule["arrays"][name], float)
            if value.shape != target.shape or np.any(~np.isfinite(value)) or np.any(value < 0.0):
                raise RuntimeError(f"invalid signed-MPZ capsule array: {name}")
            setattr(self.state, name, value.copy())
        for name in self._CAPSULE_SCALARS:
            value = float(capsule["scalars"][name])
            if not np.isfinite(value) or value < 0.0:
                raise RuntimeError(f"invalid signed-MPZ capsule scalar: {name}")
            setattr(self.state, name, value)
        self.modules["signed_burgers_shared_v1025"]._sync_active(self.state)
        self.modules["signed_burgers_shared_v1025"]._sync_wake(self.state)

    def cleavage_log_rate_s(self, opening_stress_Pa: float, T_K: float) -> float:
        barrier = self.state.manifest.cleavage
        G = float(np.asarray(barrier.values_eV(opening_stress_Pa, T_K)))
        return float(np.log(barrier.attempt_frequency_s) - G * EV_TO_J / (KB * float(T_K)))

    def emission_log_rate_per_site_s(self, stress_Pa: float, T_K: float) -> float:
        barrier = self.state.manifest.emission
        G = float(np.asarray(barrier.values_eV(stress_Pa, T_K)))
        return float(np.log(barrier.attempt_frequency_s) - G * EV_TO_J / (KB * float(T_K)))

    def advance(self, dt_s: float, T_K: float, opening_stress_Pa: float, signed_shear_Pa):
        signed = np.asarray(signed_shear_Pa, float).reshape(-1)
        if signed.size != self.state.n_systems:
            raise ValueError("signed shear history must supply exactly two channels")
        self.state._anisotropic_tau_signed_Pa = signed.copy()
        scale = max(float(opening_stress_Pa), 1.0)
        self.state._anisotropic_drive_factors = np.abs(signed) / scale
        diag = self.state.evolve(dt_s, T_K, opening_stress_Pa, self.burgers_m)
        return diag | self.summary()

    def summary(self):
        s = self.state
        campaign = self.modules["campaign_calibrated_tip"]
        rho, tau, sigma = campaign._campaign_backstress(s)
        return {
            "mobile_positive": np.asarray(s.mobile_positive).copy(),
            "mobile_negative": np.asarray(s.mobile_negative).copy(),
            "retained_positive": np.asarray(s.retained_positive).copy(),
            "retained_negative": np.asarray(s.retained_negative).copy(),
            "accumulated_slip_positive": np.asarray(s.accumulated_slip_positive).copy(),
            "accumulated_slip_negative": np.asarray(s.accumulated_slip_negative).copy(),
            "rho_back_by_system_m2": np.asarray(rho).copy(),
            "tau_back_by_system_Pa": np.asarray(tau).copy(),
            "sigma_back_by_system_Pa": np.asarray(sigma).copy(),
            "signed_active_K_shield_Pa_sqrt_m": float(s.active_K_shielding()),
            "tip_radius_m": float(s.blunted_radius(s._persistent_r0_m, self.burgers_m)),
            "emitted_total": float(s.emitted_total),
            "escaped_total": float(s.escaped_total),
            "time_s": float(s.time_s),
        }
