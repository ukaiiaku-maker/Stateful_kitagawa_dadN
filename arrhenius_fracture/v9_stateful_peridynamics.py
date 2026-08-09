"""V9 transactional extensions to the frozen v8.7 Stateful-PD patch."""

from __future__ import annotations

import math

import numpy as np

from .config import EV_TO_J, KB
from .stateful_peridynamics_v8_7_local_front_spacing import StatefulPDPatch
from .v9_physical_integrator import advance_constant_phase_log


_SITE_AVAILABLE = 0
_SITE_EMBRYO = 1
_SITE_STABLE = 2
_SITE_INACTIVE = 3


class V9StatefulPDPatch(StatefulPDPatch):
    """Preserve v8.7 spatial physics while making transition clocks persistent."""

    def initial_state(self):
        state = super().initial_state()
        count = len(state.site_node_index)
        # These draws occur once at site creation. They remain latent until a
        # site enters the embryo stage and are never redrawn by block layout.
        state.site_transition_threshold = self._event_rng.exponential(1.0, count)
        state.site_transition_cumulative_hazard = np.zeros(count)
        state.site_transition_outcome_uniform = self._event_rng.random(count)
        state.log_delivery_memory = np.full(len(self.xy), -math.inf)
        state.log_birth_cumulative_hazard = np.full(len(self.xy), -math.inf)
        return state

    @staticmethod
    def v9_extra_state_arrays(state):
        required = (
            "site_transition_threshold",
            "site_transition_cumulative_hazard",
            "site_transition_outcome_uniform",
            "log_delivery_memory",
            "log_birth_cumulative_hazard",
        )
        missing = [name for name in required if not hasattr(state, name)]
        if missing:
            raise RuntimeError("missing v9 persistent transition arrays: " + ", ".join(missing))
        return {name: np.asarray(getattr(state, name)) for name in required}

    def _phase_resolved_delivery_nucleation(
        self, delivery_rate_phase_s, nucleation_rate_phase_s,
        frequency_Hz, memory0, dN=None,
    ):
        if dN is None or not hasattr(self, "_v9_log_memory_context"):
            return super()._phase_resolved_delivery_nucleation(
                delivery_rate_phase_s, nucleation_rate_phase_s,
                frequency_Hz, memory0, dN=dN,
            )
        delivery = np.maximum(np.asarray(delivery_rate_phase_s, float), 0.0)
        cleavage = np.maximum(np.asarray(nucleation_rate_phase_s, float), 0.0)
        log_delivery = np.asarray(
            getattr(self, "_v9_log_delivery_phase_context", np.where(delivery > 0.0, np.log(delivery), -math.inf)),
            float,
        )
        log_cleavage = np.asarray(
            getattr(self, "_v9_log_cleavage_phase_context", np.where(cleavage > 0.0, np.log(cleavage), -math.inf)),
            float,
        )
        nphase = delivery.shape[0]
        frequency = float(frequency_Hz)
        phase_cycles = 1.0 / nphase
        phase_dt_s = phase_cycles / frequency
        remaining = float(dN)
        phase_position = float(self._v9_phase_cycle_context % 1.0)
        log_memory = np.asarray(self._v9_log_memory_context, float).copy()
        log_birth = np.full(log_memory.shape, -math.inf)
        max_memory = np.where(np.isfinite(log_memory), np.exp(np.maximum(log_memory, -745.0)), 0.0)
        segments = 0
        def advance_span(memory, position, span):
            local_hazard = np.full(memory.shape, -math.inf)
            local_max = np.where(np.isfinite(memory), np.exp(np.maximum(memory, -745.0)), 0.0)
            local_segments = 0
            left = float(span)
            while left > 2e-15:
                phase_scaled = position * nphase
                phase_floor = math.floor(phase_scaled + 1e-14)
                index = int(phase_floor) % nphase
                within = max(phase_scaled - phase_floor, 0.0)
                available_cycles = (1.0 - within) * phase_cycles
                step_cycles = min(left, available_cycles)
                memory, increment = advance_constant_phase_log(
                    memory, log_delivery[index], log_cleavage[index],
                    step_cycles / frequency, self.cfg.delivery_memory_s,
                )
                local_hazard = np.logaddexp(local_hazard, increment)
                represented = np.where(
                    np.isfinite(memory), np.exp(np.maximum(memory, -745.0)), 0.0
                )
                local_max = np.maximum(local_max, represented)
                left -= step_cycles
                position = (position + step_cycles) % 1.0
                local_segments += 1
            return memory, position, local_hazard, local_max, local_segments

        # The phase forcing is periodic and the memory map contracts exactly
        # by exp[-period/tau]. Iterate real cycles until both memory and gated
        # hazard reach a certified periodic orbit, then sum identical cycles in
        # log space. No 512-cycle cutoff or fixed physical rate floor is used.
        # With tau=1 ms and f=1 kHz the exact memory map contracts by e^-1
        # each cycle.  Waiting 1024 cycles before using the certified periodic
        # map forced tens of thousands of unnecessary phase updates around a
        # localized event.  Thirty-two cycles already bound the unrepresented
        # initial memory by exp(-32) < 1.3e-14.
        if remaining > 32.0:
            previous_hazard = None
            converged = False
            contraction = math.exp(-(1.0 / frequency) / self.cfg.delivery_memory_s)
            for _ in range(256):
                old_memory = log_memory.copy()
                log_memory, phase_position, cycle_hazard, cycle_max, used = advance_span(
                    log_memory, phase_position, 1.0
                )
                log_birth = np.logaddexp(log_birth, cycle_hazard)
                max_memory = np.maximum(max_memory, cycle_max)
                remaining -= 1.0
                segments += used
                old_linear = np.where(np.isfinite(old_memory), np.exp(np.maximum(old_memory, -745.0)), 0.0)
                new_linear = np.where(np.isfinite(log_memory), np.exp(np.maximum(log_memory, -745.0)), 0.0)
                memory_error = float(np.max(np.abs(new_linear - old_linear) / np.maximum(new_linear, 1e-300)))
                hazard_error = math.inf if previous_hazard is None else float(
                    np.max(np.abs(cycle_hazard - previous_hazard))
                )
                previous_hazard = cycle_hazard
                if memory_error <= 1e-12 and hazard_error <= 1e-12:
                    converged = True
                    break
            if not converged:
                raise RuntimeError("v9 periodic memory/hazard orbit did not converge")
            skipped = int(math.floor(remaining + 1e-12))
            if skipped > 0:
                log_birth = np.logaddexp(log_birth, previous_hazard + math.log(skipped))
                remaining -= skipped
                # Conservative geometric relative memory remainder diagnostic.
                self._v9_periodic_remainder_bound = memory_error / max(1.0 - contraction, 1e-15)

        if remaining > 2e-15:
            log_memory, phase_position, increment, local_max, used = advance_span(
                log_memory, phase_position, remaining
            )
            log_birth = np.logaddexp(log_birth, increment)
            max_memory = np.maximum(max_memory, local_max)
            segments += used
        self._v9_log_memory_result = log_memory
        self._v9_log_birth_increment = log_birth
        memory_end = np.where(log_memory >= -745.0, np.exp(log_memory), 0.0)
        birth = np.where(log_birth >= -745.0, np.exp(log_birth), 0.0)
        from scipy.special import gammainc
        completion = gammainc(self.cfg.delivery_hit_count, max_memory)
        period = 1.0 / frequency
        return {
            "delivery_events_cycle": np.sum(delivery, axis=0) * period / nphase,
            "delivery_rate_peak_s": np.max(delivery, axis=0),
            "nucleation_hazard_cycle": np.sum(cleavage, axis=0) * period / nphase,
            "nucleation_rate_peak_s": np.max(cleavage, axis=0),
            "memory_end": memory_end,
            "memory_max_block": max_memory,
            "completion_max_block": completion,
            "birth_hazard_block": birth,
            "birth_hazard_mean_per_cycle": birth / max(float(dN), 1e-300),
            "transient_cycles_used": segments,
        }

    def preview_rates(
        self, state, crack_barrier, sigma_hist_global,
        delivery_rate_phase_global, T_K, frequency_Hz,
        state_shift_eV_global, sigma_back_global, chi,
        plastic_state_global=None, point_amp=None, bond_amp=None,
    ):
        rates = super().preview_rates(
            state, crack_barrier, sigma_hist_global,
            delivery_rate_phase_global, T_K, frequency_Hz,
            state_shift_eV_global, sigma_back_global, chi,
            plastic_state_global, point_amp, bond_amp,
        )
        if hasattr(self, "_v9_log_memory_context"):
            if point_amp is None:
                point_amp = np.ones(len(self.xy))
            s1, _, _ = self._point_drivers(sigma_hist_global, point_amp)
            shift = np.asarray(state_shift_eV_global, float)[self.global_nodes]
            back = np.asarray(sigma_back_global, float)[self.global_nodes]
            sig_open = np.maximum(s1 - float(chi) * back[None, :], 0.0)
            energy = np.maximum(
                crack_barrier.deltaG_eV(sig_open, T_K) + shift[None, :], 1e-12
            )
            if crack_barrier.rate_prefactor <= 0.0:
                self._v9_log_cleavage_phase_context = np.full_like(energy, -math.inf)
            else:
                self._v9_log_cleavage_phase_context = (
                    math.log(float(crack_barrier.rate_prefactor))
                    - energy / ((KB / EV_TO_J) * float(T_K))
                )
        return rates

    def update(self, state, *args, **kwargs):
        # Parent positional convention: dN index 5, cycles_old index 6 after state.
        log_delivery_global = kwargs.pop("log_delivery_rate_phase_global", None)
        dN = kwargs.get("dN", args[5] if len(args) > 5 else None)
        cycles_old = kwargs.get("cycles_old", args[6] if len(args) > 6 else None)
        self._v9_log_memory_context = np.asarray(state.log_delivery_memory, float)
        self._v9_phase_cycle_context = float(cycles_old)
        if log_delivery_global is not None:
            self._v9_log_delivery_phase_context = np.asarray(log_delivery_global, float)[:, self.global_nodes]
        try:
            diagnostics = super().update(state, *args, **kwargs)
            state.log_delivery_memory = np.asarray(self._v9_log_memory_result, float).copy()
            state.log_birth_cumulative_hazard = np.logaddexp(
                np.asarray(state.log_birth_cumulative_hazard, float),
                np.asarray(self._v9_log_birth_increment, float),
            )
            return diagnostics
        finally:
            for name in (
                "_v9_log_memory_context", "_v9_phase_cycle_context",
                "_v9_log_memory_result", "_v9_log_birth_increment",
                "_v9_log_delivery_phase_context", "_v9_log_cleavage_phase_context",
                "_v9_periodic_remainder_bound",
            ):
                if hasattr(self, name):
                    delattr(self, name)

    def _advance_discrete_embryo_transitions(self, state, mu_stab, mu_heal, cycles_old, dN):
        """Integrated persistent competing-risk clocks; no block-local redraw."""
        self._sync_site_ledger(state)
        extra = self.v9_extra_state_arrays(state)
        threshold = extra["site_transition_threshold"]
        cumulative = extra["site_transition_cumulative_hazard"]
        outcome_uniform = extra["site_transition_outcome_uniform"]
        mu_s = np.maximum(np.asarray(mu_stab, float), 0.0)
        mu_h = np.maximum(np.asarray(mu_heal, float), 0.0)
        nodes = np.asarray(state.site_node_index, dtype=np.int64)
        embryo_ids = np.where(np.asarray(state.site_status, dtype=np.uint8) == _SITE_EMBRYO)[0]
        stabilized_counts = np.zeros(len(self.xy), dtype=np.int64)
        healed_counts = np.zeros(len(self.xy), dtype=np.int64)
        stable_cycles = []
        cycles_new = float(cycles_old) + float(dN)
        for site in embryo_ids:
            node = int(nodes[site])
            rate_s = float(mu_s[node])
            rate_h = float(mu_h[node])
            rate = rate_s + rate_h
            if rate <= 0.0:
                continue
            entered = float(state.site_birth_cycle[site])
            start = max(float(cycles_old), entered if np.isfinite(entered) else float(cycles_old))
            exposure = max(cycles_new - start, 0.0)
            if exposure <= 0.0:
                continue
            old_hazard = float(cumulative[site])
            increment = rate * exposure
            new_hazard = old_hazard + increment
            crossing_tol = 8.0 * np.finfo(float).eps * max(1.0, abs(float(threshold[site])))
            if new_hazard + crossing_tol < float(threshold[site]):
                cumulative[site] = new_hazard
                continue
            fraction = (float(threshold[site]) - old_hazard) / increment
            transition_cycle = start + float(np.clip(fraction, 0.0, 1.0)) * exposure
            cumulative[site] = threshold[site]
            stabilize_probability = rate_s / rate
            if float(outcome_uniform[site]) < stabilize_probability:
                state.site_status[site] = _SITE_STABLE
                state.site_stable_cycle[site] = transition_cycle
                stabilized_counts[node] += 1
                stable_cycles.append(transition_cycle)
            else:
                healed_counts[node] += 1
                state.healed_sites_cumulative[node] += 1
                if self._event_rng.random() < float(np.clip(self.cfg.heal_return_fraction, 0.0, 1.0)):
                    state.site_status[site] = _SITE_AVAILABLE
                    H = float(state.birth_cumulative_hazard[node])
                    state.site_birth_threshold[site] = H + float(self._event_rng.exponential(1.0))
                    state.site_birth_cycle[site] = np.nan
                    state.site_stable_cycle[site] = np.nan
                    # A new embryo episode receives a new persistent competing
                    # clock now, not once per future integration block.
                    threshold[site] = float(self._event_rng.exponential(1.0))
                    cumulative[site] = 0.0
                    outcome_uniform[site] = float(self._event_rng.random())
                else:
                    state.site_status[site] = _SITE_INACTIVE
        self._refresh_node_counts_from_site_ledger(state)
        return stabilized_counts, healed_counts, (min(stable_cycles) if stable_cycles else None)

    def next_embryo_transition_wait_cycles(self, state, mu_stab, mu_heal):
        """Exact wait to the next persistent embryo transition at frozen rates."""
        self._sync_site_ledger(state)
        extra = self.v9_extra_state_arrays(state)
        threshold = np.asarray(extra["site_transition_threshold"], float)
        cumulative = np.asarray(extra["site_transition_cumulative_hazard"], float)
        nodes = np.asarray(state.site_node_index, dtype=np.int64)
        embryo = np.asarray(state.site_status, dtype=np.uint8) == _SITE_EMBRYO
        if not np.any(embryo):
            return float("inf")
        rates = (
            np.maximum(np.asarray(mu_stab, float), 0.0)
            + np.maximum(np.asarray(mu_heal, float), 0.0)
        )[nodes[embryo]]
        residual = np.maximum(threshold[embryo] - cumulative[embryo], 0.0)
        waits = np.divide(residual, rates, out=np.full_like(residual, np.inf), where=rates > 0.0)
        return float(np.min(waits)) if waits.size else float("inf")
