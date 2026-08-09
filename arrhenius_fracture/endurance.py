"""Conservative asymptotic classification of positive hazard histories.

This module diagnoses an observed tail; it does not turn a finite runout into
an endurance limit.  Finite-window regression produces empirical evidence and
an extrapolation.  A physics classification requires an external asymptotic
bound/proof, except in tests where the supplied generating law is explicitly
declared closed and exact.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Iterable

import numpy as np


CLASSIFICATIONS = ("endurance_supported", "no_endurance", "undetermined")


@dataclass(frozen=True)
class WindowFit:
    start_cycle: float
    end_cycle: float
    point_count: int
    model: str
    parameters: dict[str, float]
    r_squared: float
    extrapolated_remaining_integrated_hazard: float | None
    empirical_evidence: str


@dataclass(frozen=True)
class EnduranceResult:
    empirical_tail_evidence: str
    physics_asymptotic_classification: str
    classification: str
    reason_code: str
    reason: str
    tail_model: str | None
    tail_parameters: dict[str, float]
    fit_window: tuple[float, float] | None
    window_sensitivity: list[dict[str, object]]
    empirical_extrapolated_remaining_integrated_hazard: float | None
    empirical_extrapolated_survival_asymptote: float | None
    physics_remaining_integrated_hazard_upper_bound: float | None
    physics_survival_asymptote_lower_bound: float | None
    observed_cumulative_hazard: float

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class EndpointEnduranceDiagnostics:
    """Keep nucleation/birth and the physical failure endpoint distinct."""

    birth_endurance_diagnostic: EnduranceResult
    physical_handoff_endurance_classification: str
    physical_handoff_reason_code: str

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def classify_physical_handoff_endurance(
    birth: EnduranceResult,
    *,
    handoff_divergence_proven: bool = False,
) -> EndpointEnduranceDiagnostics:
    """Propagate only logically sufficient birth information to handoff.

    Finite birth hazard supports a nonzero no-birth population and therefore a
    nonzero no-handoff population.  Divergent birth hazard alone says nothing
    decisive about handoff because healing/stabilization/growth/linkage remain.
    """
    if birth.physics_asymptotic_classification == "endurance_supported":
        classification = "endurance_supported"
        reason = "finite_birth_hazard_implies_nonzero_no_handoff_population"
    elif handoff_divergence_proven:
        classification = "no_endurance"
        reason = "full_handoff_process_divergence_proven"
    else:
        classification = "undetermined"
        reason = "births_do_not_determine_handoff_progression"
    return EndpointEnduranceDiagnostics(birth, classification, reason)


def _linear_fit(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    slope, intercept = np.polyfit(x, y, 1)
    predicted = slope * x + intercept
    ss_res = float(np.sum((y - predicted) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r_squared = 1.0 if ss_tot == 0.0 and ss_res == 0.0 else 1.0 - ss_res / ss_tot
    return float(slope), float(intercept), float(r_squared)


def _fit_window(n: np.ndarray, h: np.ndarray) -> WindowFit:
    log_h = np.log(h)
    power_slope, power_intercept, power_r2 = _linear_fit(np.log(n), log_h)
    exp_slope, exp_intercept, exp_r2 = _linear_fit(n, log_h)

    relative_span = float(np.ptp(h) / np.mean(h))
    if relative_span <= 1.0e-8:
        model = "constant"
        rate = float(np.mean(h))
        params = {"hazard": rate}
        residual = None
        evidence = "divergent_tail_observed"
        quality = 1.0
    elif exp_slope < 0.0 and exp_r2 > power_r2 + 1.0e-6:
        model = "exponential"
        scale = -1.0 / exp_slope
        amplitude = math.exp(exp_intercept)
        residual = float(h[-1] * scale)
        params = {"amplitude": amplitude, "scale_cycles": scale}
        evidence = "integrable_tail_observed"
        quality = exp_r2
    else:
        model = "power_law"
        exponent = -power_slope
        amplitude = math.exp(power_intercept)
        params = {"amplitude": amplitude, "exponent": exponent}
        evidence = "integrable_tail_observed" if exponent > 1.0 + 1.0e-10 else "divergent_tail_observed"
        residual = (
            float(amplitude * n[-1] ** (1.0 - exponent) / (exponent - 1.0))
            if exponent > 1.0
            else None
        )
        quality = power_r2

    return WindowFit(
        start_cycle=float(n[0]),
        end_cycle=float(n[-1]),
        point_count=len(n),
        model=model,
        parameters=params,
        r_squared=quality,
        extrapolated_remaining_integrated_hazard=residual,
        empirical_evidence=evidence,
    )


def analyze_endurance_tail(
    cycles: Iterable[float],
    hazard: Iterable[float],
    *,
    cumulative_hazard: Iterable[float] | None = None,
    window_fractions: tuple[float, ...] = (0.50, 0.65, 0.80),
    min_points: int = 30,
    min_r_squared: float = 0.98,
    exponent_guard: float = 0.08,
    min_tail_dynamic_range: float = 10.0,
    tail_model_is_exact: bool = False,
    physics_remaining_integrated_hazard_upper_bound: float | None = None,
    physics_divergence_proven: bool = False,
) -> EnduranceResult:
    """Classify a positive hazard tail using nested-window sensitivity.

    ``classification`` is a compatibility alias for
    ``physics_asymptotic_classification``.  Empirical data remain physically
    undetermined unless a model-derived bound/proof is supplied.  Set
    ``tail_model_is_exact`` only when the analytical generating law is known,
    as in closed synthetic tests.
    """
    n = np.asarray(tuple(cycles), dtype=float)
    h = np.asarray(tuple(hazard), dtype=float)
    if n.ndim != 1 or h.ndim != 1 or len(n) != len(h):
        raise ValueError("cycles and hazard must be one-dimensional and equal length")
    if len(n) < min_points or not np.all(np.isfinite(n)) or not np.all(np.isfinite(h)):
        return _undetermined("insufficient_or_nonfinite_tail", n, h)
    if np.any(n <= 0.0) or np.any(np.diff(n) <= 0.0) or np.any(h <= 0.0):
        return _undetermined("invalid_or_clipped_hazard", n, h)

    if cumulative_hazard is None:
        observed = float(np.trapezoid(h, n))
    else:
        cumulative = np.asarray(tuple(cumulative_hazard), dtype=float)
        if cumulative.shape != n.shape or np.any(np.diff(cumulative) < 0.0):
            raise ValueError("cumulative_hazard must match cycles and be nondecreasing")
        observed = float(cumulative[-1])

    fits = []
    for fraction in window_fractions:
        start = int(math.floor(len(n) * fraction))
        if len(n) - start >= max(8, min_points // 5):
            fits.append(_fit_window(n[start:], h[start:]))
    sensitivity = [asdict(fit) for fit in fits]
    if not fits or any(fit.r_squared < min_r_squared for fit in fits):
        return _undetermined("poor_tail_fit", n, h, observed, sensitivity)
    if len({fit.empirical_evidence for fit in fits}) != 1 or len({fit.model for fit in fits}) != 1:
        return _undetermined("fit_window_sensitive", n, h, observed, sensitivity)
    if fits[0].model != "constant" and h[int(math.floor(len(h) * window_fractions[0]))] / h[-1] < min_tail_dynamic_range:
        return _undetermined("insufficient_tail_dynamic_range", n, h, observed, sensitivity)
    if fits[0].model == "power_law" and any(
        abs(fit.parameters["exponent"] - 1.0) < exponent_guard for fit in fits
    ):
        # The exact harmonic tail is mathematically divergent.  Permit p=1
        # only when the fits identify it essentially exactly; noisy estimates
        # close to the boundary remain unresolved.
        if not all(abs(fit.parameters["exponent"] - 1.0) < 1.0e-8 for fit in fits):
            return _undetermined("integrability_boundary", n, h, observed, sensitivity)

    selected = fits[0]
    residuals = [fit.extrapolated_remaining_integrated_hazard for fit in fits]
    residual = max(value for value in residuals if value is not None) if all(
        value is not None for value in residuals
    ) else None
    empirical_survival = math.exp(-(observed + residual)) if residual is not None else None
    evidence = selected.empirical_evidence
    if physics_remaining_integrated_hazard_upper_bound is not None and physics_divergence_proven:
        raise ValueError("finite physics bound and divergence proof are mutually exclusive")
    if physics_remaining_integrated_hazard_upper_bound is not None:
        bound = float(physics_remaining_integrated_hazard_upper_bound)
        if not math.isfinite(bound) or bound < 0.0:
            raise ValueError("physics remaining-hazard upper bound must be finite and nonnegative")
        classification = "endurance_supported"
        reason_code = "physics_finite_hazard_bound"
        reason = "A model-derived finite upper bound supports nonzero asymptotic survival."
        physics_survival = math.exp(-(observed + bound))
    elif physics_divergence_proven:
        bound = None
        classification = "no_endurance"
        reason_code = "physics_hazard_divergence_proven"
        reason = "The physical model proves divergent cumulative hazard."
        physics_survival = 0.0
    elif tail_model_is_exact:
        bound = residual
        classification = "endurance_supported" if evidence == "integrable_tail_observed" else "no_endurance"
        reason_code = "closed_exact_integrable_model" if classification == "endurance_supported" else "closed_exact_divergent_model"
        reason = "The analytical generating law is declared exact and closed."
        physics_survival = empirical_survival if classification == "endurance_supported" else 0.0
    else:
        bound = None
        classification = "undetermined"
        reason_code = "empirical_fit_without_physics_bound"
        reason = "Finite-window evidence cannot exclude a hidden floor or late crossover."
        physics_survival = None
    return EnduranceResult(
        empirical_tail_evidence=evidence,
        physics_asymptotic_classification=classification,
        classification=classification,
        reason_code=reason_code,
        reason=reason,
        tail_model=selected.model,
        tail_parameters=selected.parameters,
        fit_window=(selected.start_cycle, selected.end_cycle),
        window_sensitivity=sensitivity,
        empirical_extrapolated_remaining_integrated_hazard=residual,
        empirical_extrapolated_survival_asymptote=empirical_survival,
        physics_remaining_integrated_hazard_upper_bound=bound,
        physics_survival_asymptote_lower_bound=physics_survival,
        observed_cumulative_hazard=observed,
    )


def _undetermined(
    reason_code: str,
    n: np.ndarray,
    h: np.ndarray,
    observed: float = 0.0,
    sensitivity: list[dict[str, object]] | None = None,
) -> EnduranceResult:
    return EnduranceResult(
        empirical_tail_evidence="ambiguous_tail_observed",
        physics_asymptotic_classification="undetermined",
        classification="undetermined",
        reason_code=reason_code,
        reason="The observed tail does not support a window-stable asymptotic bound.",
        tail_model=None,
        tail_parameters={},
        fit_window=(float(n[0]), float(n[-1])) if len(n) else None,
        window_sensitivity=sensitivity or [],
        empirical_extrapolated_remaining_integrated_hazard=None,
        empirical_extrapolated_survival_asymptote=None,
        physics_remaining_integrated_hazard_upper_bound=None,
        physics_survival_asymptote_lower_bound=None,
        observed_cumulative_hazard=observed,
    )
