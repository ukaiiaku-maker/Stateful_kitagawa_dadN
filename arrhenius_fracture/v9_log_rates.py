"""Log-domain physical rates for the v9 large-N solver.

No representability limit in this module is a physical rate floor.  Exact zero
is reserved for an explicit constitutive switch or a mathematically zero gate.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
import sys
from typing import Iterable


LOG_MAX_FLOAT = math.log(sys.float_info.max)
LOG_MIN_SUBNORMAL = math.log(float.fromhex("0x0.0000000000001p-1022"))
KB_EV_PER_K = 8.617333262145e-5


@dataclass(frozen=True)
class LogPhysicalRate:
    name: str
    units: str
    log_rate: float
    constitutive_zero: bool = False
    extinction_reason: str | None = None
    cap_floor_activations: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.constitutive_zero:
            object.__setattr__(self, "log_rate", -math.inf)
            if not self.extinction_reason:
                raise ValueError("constitutive zero requires a physical extinction reason")
        elif not math.isfinite(self.log_rate):
            raise ValueError("nonzero physical log-rate must be finite")

    @property
    def numerical_representability(self) -> str:
        if self.constitutive_zero:
            return "exact_constitutive_zero"
        if self.log_rate > LOG_MAX_FLOAT:
            return "overflow_in_float64_rate_space"
        if self.log_rate < LOG_MIN_SUBNORMAL:
            return "underflow_in_float64_rate_space"
        return "representable_float64"

    @property
    def physical_rate(self) -> float | None:
        if self.constitutive_zero:
            return 0.0
        if self.numerical_representability != "representable_float64":
            return None
        return math.exp(self.log_rate)

    def diagnostic(self) -> dict[str, object]:
        return asdict(self) | {
            "physical_rate": self.physical_rate,
            "numerical_representability": self.numerical_representability,
            "numerical_floor_used_as_physics": False,
        }

    def scaled(self, factor: float, *, name: str | None = None, reason: str | None = None) -> "LogPhysicalRate":
        if factor < 0.0 or not math.isfinite(factor):
            raise ValueError("rate scale must be finite and nonnegative")
        if factor == 0.0:
            return LogPhysicalRate(name or self.name, self.units, -math.inf, True, reason or "exact_zero_scale")
        if self.constitutive_zero:
            return LogPhysicalRate(name or self.name, self.units, -math.inf, True, self.extinction_reason)
        return LogPhysicalRate(name or self.name, self.units, self.log_rate + math.log(factor), cap_floor_activations=self.cap_floor_activations)


def arrhenius_log_rate(
    name: str,
    prefactor: float,
    barrier_eV: float,
    temperature_K: float,
    *,
    units: str = "s^-1",
    constitutive_enabled: bool = True,
    activations: Iterable[str] = (),
) -> LogPhysicalRate:
    if not constitutive_enabled:
        return LogPhysicalRate(name, units, -math.inf, True, "mechanism_disabled")
    if prefactor <= 0.0:
        return LogPhysicalRate(name, units, -math.inf, True, "nonpositive_physical_prefactor")
    if barrier_eV < 0.0 or temperature_K <= 0.0:
        raise ValueError("Arrhenius barrier must be nonnegative and temperature positive")
    return LogPhysicalRate(
        name,
        units,
        math.log(prefactor) - barrier_eV / (KB_EV_PER_K * temperature_K),
        cap_floor_activations=tuple(activations),
    )


def series_log_rate(name: str, rates: Iterable[LogPhysicalRate]) -> LogPhysicalRate:
    """Series residence-time rate: ``1/r = sum(1/r_i)`` in log space."""
    items = tuple(rates)
    if not items:
        raise ValueError("series rate requires at least one mechanism")
    zero = next((rate for rate in items if rate.constitutive_zero), None)
    if zero is not None:
        return LogPhysicalRate(name, items[0].units, -math.inf, True, f"series_blocked_by:{zero.name}")
    residence_logs = [-rate.log_rate for rate in items]
    pivot = max(residence_logs)
    log_residence = pivot + math.log(sum(math.exp(value - pivot) for value in residence_logs))
    activations = tuple(sorted({event for rate in items for event in rate.cap_floor_activations}))
    return LogPhysicalRate(name, items[0].units, -log_residence, cap_floor_activations=activations)


def logistic_log_rate(
    name: str,
    prefactor_per_second: float,
    frequency_Hz: float,
    z: float,
    *,
    gate: float = 1.0,
    gate_reason: str = "exact_zero_state_gate",
) -> LogPhysicalRate:
    if gate < 0.0 or gate > 1.0:
        raise ValueError("probability/activity gate must lie in [0,1]")
    if gate == 0.0 or prefactor_per_second <= 0.0:
        return LogPhysicalRate(name, "cycle^-1", -math.inf, True, gate_reason)
    if frequency_Hz <= 0.0:
        raise ValueError("frequency must be positive")
    log_sigmoid = -math.log1p(math.exp(-z)) if z >= 0.0 else z - math.log1p(math.exp(z))
    return LogPhysicalRate(name, "cycle^-1", math.log(prefactor_per_second / frequency_Hz) + log_sigmoid + math.log(gate))


def completion_gated_log_rate(
    cleavage: LogPhysicalRate,
    completion_probability: float,
    *,
    name: str = "birth",
) -> LogPhysicalRate:
    if completion_probability < 0.0 or completion_probability > 1.0:
        raise ValueError("completion probability must lie in [0,1]")
    if completion_probability == 0.0:
        return LogPhysicalRate(name, cleavage.units, -math.inf, True, "exact_zero_completion")
    return cleavage.scaled(completion_probability, name=name)


def classify_exact_power_law_from_logs(log_cycles: Iterable[float], log_rates: Iterable[float]) -> tuple[str, float]:
    """Classify a declared exact power law without exponentiating its rates."""
    x = tuple(float(value) for value in log_cycles)
    y = tuple(float(value) for value in log_rates)
    if len(x) != len(y) or len(x) < 2 or any(not math.isfinite(v) for v in x + y):
        raise ValueError("finite equal-length log histories with at least two points are required")
    mx = sum(x) / len(x)
    my = sum(y) / len(y)
    denominator = sum((value - mx) ** 2 for value in x)
    if denominator == 0.0:
        raise ValueError("log cycles must span a nonzero range")
    slope = sum((a - mx) * (b - my) for a, b in zip(x, y)) / denominator
    exponent = -slope
    return ("endurance_supported" if exponent > 1.0 else "no_endurance", exponent)


def trace_completion_gated_rate_paths(
    *,
    plastic_barriers: tuple[tuple[str, float, float], ...],
    cleavage_prefactor: float,
    cleavage_barrier_eV: float,
    temperature_K: float,
    frequency_Hz: float,
    completion_probability: float,
    stabilization_prefactor_s: float,
    healing_prefactor_s: float,
    stabilization_z: float,
    growth_prefactor_s: float,
    growth_z: float,
    stable_activity: float,
    linkage_prefactor_s: float,
    linkage_z: float,
    directional_activity: float,
    activations: Iterable[str] = (),
) -> dict[str, LogPhysicalRate]:
    """Trace every v8.7 rate family into an unclipped v9 diagnostic bundle.

    ``plastic_barriers`` is ordered emission, Peierls and Taylor with tuples of
    ``(name, prefactor_s^-1, barrier_eV)``. Stabilization/healing/growth/linkage
    are logistic kinetic paths rather than Arrhenius barriers in v8.7, but are
    included so the complete event chain is explicit.
    """
    elementary = {
        name: arrhenius_log_rate(name, prefactor, barrier, temperature_K, activations=activations)
        for name, prefactor, barrier in plastic_barriers
    }
    required = ("plastic_emission", "plastic_peierls", "plastic_taylor")
    if any(name not in elementary for name in required):
        raise ValueError(f"plastic barriers must contain {required}")
    escape = series_log_rate("plastic_escape", (elementary["plastic_peierls"], elementary["plastic_taylor"]))
    delivery = series_log_rate("plastic_completed_delivery", (elementary["plastic_emission"], escape))
    cleavage = arrhenius_log_rate(
        "cleavage", cleavage_prefactor, cleavage_barrier_eV, temperature_K, activations=activations
    )
    birth = completion_gated_log_rate(cleavage, completion_probability)
    return elementary | {
        "plastic_escape": escape,
        "plastic_completed_delivery": delivery,
        "cleavage": cleavage,
        "birth": birth,
        "stabilization": logistic_log_rate("stabilization", stabilization_prefactor_s, frequency_Hz, stabilization_z),
        "healing": logistic_log_rate("healing", healing_prefactor_s, frequency_Hz, -stabilization_z),
        "stable_growth": logistic_log_rate(
            "stable_growth", growth_prefactor_s, frequency_Hz, growth_z,
            gate=stable_activity, gate_reason="no_stable_defect_activity",
        ),
        "linkage_front": logistic_log_rate(
            "linkage_front", linkage_prefactor_s, frequency_Hz, linkage_z,
            gate=directional_activity, gate_reason="no_directional_front_activity",
        ),
    }
