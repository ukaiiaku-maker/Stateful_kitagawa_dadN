"""Numerical primitives for the real v9 completion-gated physical solver.

The routines here deliberately retain logarithms of rates, memories and
cumulative hazards.  A float underflow is never interpreted as a constitutive
zero.  Persistent competing-risk clocks are sampled once and advanced by
integrated hazard, making their event identity independent of block layout.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .config import EV_TO_J, KB


KB_EV_PER_K = KB / EV_TO_J


NEG_INF = -math.inf
_GL_X, _GL_W = np.polynomial.legendre.leggauss(12)


def log1mexp(x: np.ndarray | float) -> np.ndarray:
    """Return log(1-exp(x)) for x<=0 without cancellation."""
    value = np.asarray(x, dtype=float)
    if np.any(value > 0.0):
        raise ValueError("log1mexp requires x <= 0")
    split = -math.log(2.0)
    return np.where(
        value < split,
        np.log1p(-np.exp(value)),
        np.log(-np.expm1(value)),
    )


def logaddexp_accumulate(log_total: np.ndarray, log_increment: np.ndarray) -> np.ndarray:
    total = np.asarray(log_total, dtype=float)
    increment = np.asarray(log_increment, dtype=float)
    if total.shape != increment.shape:
        raise ValueError("log totals and increments must have equal shapes")
    return np.logaddexp(total, increment)


def log_completion_k2_from_log_memory(log_memory: np.ndarray | float) -> np.ndarray:
    """Log of ``Q(2,Lambda)=1-exp(-Lambda)(1+Lambda)``.

    Below ``Lambda=1e-4`` a factored alternating series avoids cancellation;
    the calculation therefore remains meaningful when Lambda itself cannot be
    represented in float64 rate space.
    """
    log_lam = np.asarray(log_memory, dtype=float)
    out = np.full(log_lam.shape, NEG_INF, dtype=float)
    finite = np.isfinite(log_lam)
    small = finite & (log_lam < math.log(1.0e-4))
    if np.any(small):
        # Q(2,L)=L^2/2 * (1 - 2L/3 + L^2/4 - L^3/15 + ...)
        ls = log_lam[small]
        # Terms beyond the leading term are irrelevant below exp(-350), but
        # keeping the correction where representable gives a smooth splice.
        lam = np.exp(np.maximum(ls, -745.0))
        correction = 1.0 - (2.0 / 3.0) * lam + 0.25 * lam * lam
        out[small] = 2.0 * ls - math.log(2.0) + np.log(correction)
    regular = finite & ~small
    if np.any(regular):
        lam = np.exp(log_lam[regular])
        # log(exp(-L)*(1+L)) = -L + log1p(L)
        out[regular] = log1mexp(-lam + np.log1p(lam))
    return out


def advance_constant_memory_log(
    log_memory0: np.ndarray,
    log_delivery_rate_s: np.ndarray,
    dt_s: float,
    tau_s: float,
) -> np.ndarray:
    """Exact positive memory update for a constant delivery interval."""
    if dt_s < 0.0 or tau_s <= 0.0:
        raise ValueError("dt must be nonnegative and tau positive")
    m0 = np.asarray(log_memory0, dtype=float)
    rate = np.asarray(log_delivery_rate_s, dtype=float)
    if m0.shape != rate.shape:
        raise ValueError("memory and delivery-rate shapes differ")
    if dt_s == 0.0:
        return m0.copy()
    decay_log = -dt_s / tau_s
    source_factor_log = math.log(tau_s) + float(log1mexp(decay_log))
    return np.logaddexp(m0 + decay_log, rate + source_factor_log)


def _logsumexp_axis0(values: np.ndarray) -> np.ndarray:
    pivot = np.max(values, axis=0)
    out = np.full(pivot.shape, NEG_INF)
    finite = np.isfinite(pivot)
    if np.any(finite):
        out[finite] = pivot[finite] + np.log(
            np.sum(np.exp(values[:, finite] - pivot[None, finite]), axis=0)
        )
    return out


def advance_constant_phase_log(
    log_memory0: np.ndarray,
    log_delivery_rate_s: np.ndarray,
    log_cleavage_rate_s: np.ndarray,
    dt_s: float,
    tau_s: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Advance one constant phase and integrate completion-gated K=2 hazard."""
    memory0 = np.asarray(log_memory0, dtype=float)
    delivery = np.asarray(log_delivery_rate_s, dtype=float)
    cleavage = np.asarray(log_cleavage_rate_s, dtype=float)
    if memory0.shape != delivery.shape or memory0.shape != cleavage.shape:
        raise ValueError("phase memory/delivery/cleavage shapes differ")
    if dt_s <= 0.0:
        return memory0.copy(), np.full(memory0.shape, NEG_INF)
    times = 0.5 * dt_s * (_GL_X + 1.0)
    log_integrands = []
    for time, weight in zip(times, _GL_W):
        memory = advance_constant_memory_log(memory0, delivery, float(time), tau_s)
        log_integrands.append(
            cleavage + log_completion_k2_from_log_memory(memory) + math.log(float(weight))
        )
    log_hazard = _logsumexp_axis0(np.stack(log_integrands)) + math.log(0.5 * dt_s)
    return advance_constant_memory_log(memory0, delivery, dt_s, tau_s), log_hazard


def advance_phase_history_log(
    log_memory0: np.ndarray,
    log_delivery_rate_phase_s: np.ndarray,
    log_cleavage_rate_phase_s: np.ndarray,
    phase_dt_s: float,
    tau_s: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Advance an ordered phase history, accumulating hazard in log space."""
    delivery = np.asarray(log_delivery_rate_phase_s, dtype=float)
    cleavage = np.asarray(log_cleavage_rate_phase_s, dtype=float)
    if delivery.ndim != 2 or delivery.shape != cleavage.shape:
        raise ValueError("phase rates must be aligned (n_phase,n_site) arrays")
    memory = np.asarray(log_memory0, dtype=float).copy()
    if memory.shape != delivery.shape[1:]:
        raise ValueError("memory does not match phase-site shape")
    total = np.full(memory.shape, NEG_INF)
    for phase_delivery, phase_cleavage in zip(delivery, cleavage):
        memory, increment = advance_constant_phase_log(
            memory, phase_delivery, phase_cleavage, phase_dt_s, tau_s
        )
        total = np.logaddexp(total, increment)
    return memory, total


def constant_log_hazard_increment(log_rate_per_cycle: np.ndarray, dN: float) -> np.ndarray:
    if dN < 0.0:
        raise ValueError("cycle increment must be nonnegative")
    rate = np.asarray(log_rate_per_cycle, dtype=float)
    if dN == 0.0:
        return np.full(rate.shape, NEG_INF)
    return rate + math.log(dN)


def log_series_rate_arrays(*log_rates: np.ndarray) -> np.ndarray:
    """Series residence-time rate for aligned arrays, entirely in log space."""
    if not log_rates:
        raise ValueError("at least one series rate is required")
    aligned = np.broadcast_arrays(*(np.asarray(rate, dtype=float) for rate in log_rates))
    residence = np.stack([-rate for rate in aligned])
    return -_logsumexp_axis0(residence)


def plastic_chain_log_rates(chain, sigma_eq_Pa, rho_m2, temperature_K: float) -> dict[str, np.ndarray]:
    """Evaluate the frozen EXP-floor barrier surfaces without exponent clipping."""
    if temperature_K <= 0.0:
        raise ValueError("temperature must be positive")
    sigma = np.maximum(np.asarray(sigma_eq_Pa, dtype=float), 0.0)
    phi = chain.taylor_phi(rho_m2)

    def elementary(barrier, stress):
        prefactor = float(barrier.rate_prefactor)
        if prefactor <= 0.0:
            return np.full(np.broadcast_shapes(np.shape(stress), np.shape(phi)), NEG_INF)
        energy = np.asarray(barrier.deltaG_eV(stress, temperature_K), dtype=float)
        return math.log(prefactor) - energy / (KB_EV_PER_K * temperature_K)

    emission = elementary(chain.emit, sigma)
    peierls = elementary(chain.peierls, sigma)
    taylor = elementary(chain.taylor, phi * sigma)
    escape = log_series_rate_arrays(peierls, taylor)
    flow = log_series_rate_arrays(emission, escape)
    return {
        "lambda_emit": emission,
        "lambda_peierls": peierls,
        "lambda_taylor": taylor,
        "lambda_escape": escape,
        "lambda_flow": flow,
        "dot_ep": flow + math.log(max(float(chain.plastic_event_strain), 1e-300)),
        "phi_taylor": np.asarray(phi, dtype=float),
    }


def cleavage_log_rate(barrier, sigma_Pa, temperature_K: float, state_shift_eV=0.0) -> np.ndarray:
    """Frozen cleavage free-energy surface with an unclipped log rate."""
    energy = np.asarray(barrier.deltaG_eV(sigma_Pa, temperature_K), dtype=float)
    energy = np.maximum(energy + np.asarray(state_shift_eV, dtype=float), 1e-12)
    prefactor = float(barrier.rate_prefactor)
    shape = energy.shape
    if prefactor <= 0.0:
        return np.full(shape, NEG_INF)
    return math.log(prefactor) - energy / (KB_EV_PER_K * temperature_K)


@dataclass(frozen=True)
class PersistentCompetingClockLedger:
    """Per-site persistent cumulative-hazard clocks for one active stage."""

    threshold: np.ndarray
    log_cumulative_hazard: np.ndarray
    active: np.ndarray

    def __post_init__(self) -> None:
        threshold = np.asarray(self.threshold, dtype=float)
        cumulative = np.asarray(self.log_cumulative_hazard, dtype=float)
        active = np.asarray(self.active, dtype=bool)
        if threshold.shape != cumulative.shape or threshold.shape != active.shape:
            raise ValueError("persistent clock arrays must share a shape")
        if np.any(~np.isfinite(threshold)) or np.any(threshold <= 0.0):
            raise ValueError("persistent thresholds must be positive and finite")
        if np.any(np.isnan(cumulative)):
            raise ValueError("cumulative log hazards cannot be NaN")

    @classmethod
    def sampled(cls, count: int, rng: np.random.Generator) -> "PersistentCompetingClockLedger":
        if count < 0:
            raise ValueError("clock count must be nonnegative")
        return cls(
            threshold=rng.exponential(1.0, count),
            log_cumulative_hazard=np.full(count, NEG_INF),
            active=np.ones(count, dtype=bool),
        )

    def advance_constant(
        self,
        log_total_rate_per_cycle: np.ndarray,
        dN: float,
    ) -> tuple["PersistentCompetingClockLedger", np.ndarray, np.ndarray]:
        """Advance without RNG draws and return crossing mask/fractions."""
        rates = np.asarray(log_total_rate_per_cycle, dtype=float)
        if rates.shape != self.threshold.shape:
            raise ValueError("one total rate is required per clock")
        log_inc = constant_log_hazard_increment(rates, dN)
        old = np.asarray(self.log_cumulative_hazard, dtype=float)
        new = np.where(self.active, np.logaddexp(old, log_inc), old)
        log_threshold = np.log(self.threshold)
        crossed = self.active & (new >= log_threshold)
        fractions = np.full(self.threshold.shape, np.nan)
        if np.any(crossed):
            # At a crossing the required linear-hazard increment is safely
            # representable relative to the O(1) exponential threshold.
            old_linear = np.exp(old[crossed])
            needed = self.threshold[crossed] - old_linear
            increment = np.exp(log_inc[crossed])
            fractions[crossed] = np.clip(needed / increment, 0.0, 1.0)
            new[crossed] = log_threshold[crossed]
        return (
            PersistentCompetingClockLedger(self.threshold.copy(), new, self.active & ~crossed),
            crossed,
            fractions,
        )

    def choose_competing_event(
        self,
        crossed: np.ndarray,
        log_rate_a: np.ndarray,
        log_rate_b: np.ndarray,
        uniforms: np.ndarray,
    ) -> np.ndarray:
        """Choose A/B at localized crossings using preassigned uniforms."""
        mask = np.asarray(crossed, dtype=bool)
        a = np.asarray(log_rate_a, dtype=float)
        b = np.asarray(log_rate_b, dtype=float)
        u = np.asarray(uniforms, dtype=float)
        if any(value.shape != self.threshold.shape for value in (mask, a, b, u)):
            raise ValueError("competing event arrays must share the ledger shape")
        log_total = np.logaddexp(a, b)
        probability_a = np.exp(a[mask] - log_total[mask])
        selected_a = np.zeros(mask.shape, dtype=bool)
        selected_a[mask] = u[mask] < probability_a
        return selected_a
