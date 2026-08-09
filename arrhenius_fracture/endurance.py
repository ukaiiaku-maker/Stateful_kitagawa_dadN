"""Conservative asymptotic classification of positive hazard histories.

This module diagnoses an observed tail; it does not turn a finite runout into
an endurance limit.  A classification is supported only when several nested
fit windows agree on an integrable or divergent analytical tail.
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
    residual_hazard: float | None
    classification: str


@dataclass(frozen=True)
class EnduranceResult:
    classification: str
    reason_code: str
    reason: str
    tail_model: str | None
    tail_parameters: dict[str, float]
    fit_window: tuple[float, float] | None
    window_sensitivity: list[dict[str, object]]
    residual_integrated_hazard: float | None
    survival_asymptote: float | None
    observed_cumulative_hazard: float

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


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
        classification = "no_endurance"
        quality = 1.0
    elif exp_slope < 0.0 and exp_r2 > power_r2 + 1.0e-6:
        model = "exponential"
        scale = -1.0 / exp_slope
        amplitude = math.exp(exp_intercept)
        residual = float(h[-1] * scale)
        params = {"amplitude": amplitude, "scale_cycles": scale}
        classification = "endurance_supported"
        quality = exp_r2
    else:
        model = "power_law"
        exponent = -power_slope
        amplitude = math.exp(power_intercept)
        params = {"amplitude": amplitude, "exponent": exponent}
        classification = "endurance_supported" if exponent > 1.0 + 1.0e-10 else "no_endurance"
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
        residual_hazard=residual,
        classification=classification,
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
) -> EnduranceResult:
    """Classify a positive hazard tail using nested-window sensitivity.

    ``undetermined`` is returned for insufficient, non-finite, non-positive,
    poorly fitted, window-sensitive, or borderline-integrability histories.
    Zero hazards are rejected because numerical clipping must not manufacture
    an endurance population.
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
    if len({fit.classification for fit in fits}) != 1 or len({fit.model for fit in fits}) != 1:
        return _undetermined("fit_window_sensitive", n, h, observed, sensitivity)
    if fits[0].model == "power_law" and any(
        abs(fit.parameters["exponent"] - 1.0) < exponent_guard for fit in fits
    ):
        # The exact harmonic tail is mathematically divergent.  Permit p=1
        # only when the fits identify it essentially exactly; noisy estimates
        # close to the boundary remain unresolved.
        if not all(abs(fit.parameters["exponent"] - 1.0) < 1.0e-8 for fit in fits):
            return _undetermined("integrability_boundary", n, h, observed, sensitivity)

    selected = fits[0]
    residuals = [fit.residual_hazard for fit in fits]
    residual = max(value for value in residuals if value is not None) if all(
        value is not None for value in residuals
    ) else None
    survival = math.exp(-(observed + residual)) if residual is not None else None
    classification = selected.classification
    return EnduranceResult(
        classification=classification,
        reason_code="stable_integrable_tail" if classification == "endurance_supported" else "stable_divergent_tail",
        reason="Nested tail windows agree on an integrable model." if classification == "endurance_supported" else "Nested tail windows agree on a divergent model.",
        tail_model=selected.model,
        tail_parameters=selected.parameters,
        fit_window=(selected.start_cycle, selected.end_cycle),
        window_sensitivity=sensitivity,
        residual_integrated_hazard=residual,
        survival_asymptote=survival,
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
        classification="undetermined",
        reason_code=reason_code,
        reason="The observed tail does not support a window-stable asymptotic bound.",
        tail_model=None,
        tail_parameters={},
        fit_window=(float(n[0]), float(n[-1])) if len(n) else None,
        window_sensitivity=sensitivity or [],
        residual_integrated_hazard=None,
        survival_asymptote=None,
        observed_cumulative_hazard=observed,
    )
