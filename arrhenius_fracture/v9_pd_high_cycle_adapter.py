"""Spatial Stateful-PD adapter for :mod:`v9_pd_high_cycle`.

The driver supplies ``cycle_evaluator`` because it owns the accepted FEM phase
ordering.  The evaluator receives a private adapter clone and must return phase
data plus the endpoint active arrays.  This class owns all state separation and
commit/protection rules, so an evaluator cannot accidentally advance clocks.
"""
from __future__ import annotations

from copy import copy, deepcopy
from dataclasses import fields
import math
from typing import Any, Callable

import numpy as np

from .v9_pd_high_cycle import ActiveState, CycleEvaluation, ProtectedSignatures, _digest


ACTIVE_NAMES = ("ep_gp", "rho_gp", "epsp_acc_gp", "u", "log_delivery_memory",
                "available", "embryo", "stable", "inactive", "completion")
FIELD_RULES = {
    "ep_gp": {"transform": "scaled_linear", "scale": 1.0e-4, "rtol": 2.0e-7, "lower": -0.25, "upper": 0.25},
    "rho_gp": {"transform": "log_positive", "scale": 1.0e12, "rtol": 2.0e-7, "lower": 1.0e6, "upper": 1.0e20},
    "epsp_acc_gp": {"transform": "log1p_nonnegative", "scale": 1.0e-5, "rtol": 2.0e-7, "lower": 0.0, "upper": 10.0},
    "u": {"transform": "scaled_linear", "scale": 1.0e-6, "rtol": 2.0e-8, "lower": -0.1, "upper": 0.1},
    "log_delivery_memory": {"transform": "log1p_from_log", "scale": 1.0, "rtol": 2.0e-7, "lower": -math.inf, "upper": 709.0},
    "available": {"transform": "scaled_linear", "scale": 1.0, "rtol": 2.0e-7, "lower": 0.0, "upper": 1.0},
    "embryo": {"transform": "scaled_linear", "scale": 1.0, "rtol": 2.0e-7, "lower": 0.0, "upper": 1.0},
    "stable": {"transform": "scaled_linear", "scale": 1.0, "rtol": 2.0e-7, "lower": 0.0, "upper": 1.0},
    "inactive": {"transform": "scaled_linear", "scale": 1.0, "rtol": 2.0e-7, "lower": 0.0, "upper": 1.0},
    "completion": {"transform": "scaled_linear", "scale": 1.0, "rtol": 2.0e-7, "lower": 0.0, "upper": 1.0},
}
LEDGER_NAMES = (
    "born_cumulative", "healed_cumulative", "born_sites_cumulative",
    "healed_sites_cumulative", "log_birth_cumulative_hazard",
)
STOCHASTIC_NAMES = (
    "site_birth_threshold", "birth_cumulative_hazard",
    "site_transition_threshold", "site_transition_cumulative_hazard",
    "site_transition_outcome_uniform",
)
TOPOLOGY_NAMES = (
    "site_status", "bond_damage", "primary_seed_node", "active_front",
    "active_front_bonds", "front_backbone_bonds", "front_wake_bonds",
    "front_process_bonds", "active_front_path_xy",
)


class SpatialPDDormantAdapter:
    def __init__(self, *, patch, pd_state, mesh, ep_gp, rho_gp, epsp_acc_gp, u,
                 cycles: float, plastic_work: float,
                 cycle_evaluator: Callable[["SpatialPDDormantAdapter"], dict[str, Any]]):
        self.patch = patch
        self.pd_state = pd_state
        self.mesh = mesh
        self.ep_gp = np.asarray(ep_gp, float).copy()
        self.rho_gp = np.asarray(rho_gp, float).copy()
        self.epsp_acc_gp = np.asarray(epsp_acc_gp, float).copy()
        self.u = np.asarray(u, float).copy()
        self.cycles = float(cycles)
        self.plastic_work = float(plastic_work)
        self.external_ledgers = {"born_expectation": 0.0, "healed_expectation": 0.0}
        self.cycle_evaluator = cycle_evaluator

    def dormant_eligibility(self):
        s = self.pd_state
        checks = (
            (not np.any(np.asarray(s.site_status) == 1), "active_embryo"),
            (not np.any(np.asarray(s.site_status) == 2), "stabilized_site"),
            (not np.any(np.asarray(s.bond_damage) > 0.0), "bond_damage"),
            (not bool(s.active_front), "captured_front"),
            (int(s.primary_seed_node) < 0, "primary_seed"),
        )
        for passed, reason in checks:
            if not passed: return False, reason
        return True, "dormant_fixed_topology"

    def _active_arrays(self):
        return (self.ep_gp, self.rho_gp, self.epsp_acc_gp, self.u,
                np.asarray(self.pd_state.log_delivery_memory, float),
                np.asarray(self.pd_state.available,float),np.asarray(self.pd_state.embryo,float),
                np.asarray(self.pd_state.stable,float),np.asarray(self.pd_state.inactive,float),
                np.asarray(self.pd_state.completion,float))

    def active_state(self):
        arrays = self._active_arrays()
        transformed = []
        for name, a in zip(ACTIVE_NAMES, arrays):
            rule = FIELD_RULES[name]; scale = rule["scale"]
            if rule["transform"] == "scaled_linear": value = a / scale
            elif rule["transform"] == "log_positive": value = np.log(a / scale)
            elif rule["transform"] == "log1p_nonnegative": value = np.log1p(a / scale)
            elif rule["transform"] == "log_positive_offset": value = np.log((a + rule["offset"]) / scale)
            else:
                linear = np.where(np.isfinite(a), np.exp(np.minimum(a, 709.0)), 0.0)
                value = np.log1p(linear / scale)
            transformed.append(value)
        spec = tuple((name, a.shape, str(a.dtype)) for name, a in zip(ACTIVE_NAMES, arrays))
        vector = np.concatenate([a.ravel() for a in transformed]) if arrays else np.empty(0)
        return ActiveState(vector, spec)

    def restore_active_state(self, snapshot, vector):
        x = np.asarray(vector, float)
        offset, restored = 0, []
        names = []
        for name, shape, _ in snapshot.specification:
            size = int(np.prod(shape, dtype=int))
            value = x[offset:offset + size].reshape(shape).copy(); offset += size
            restored.append(value); names.append(name)
        if offset != x.size or len(restored) != len(ACTIVE_NAMES):
            raise ValueError("active state layout mismatch")
        physical = []
        for name, value in zip(names, restored):
            rule = FIELD_RULES[name]; scale = rule["scale"]
            if rule["transform"] == "scaled_linear": field = value * scale
            elif rule["transform"] == "log_positive": field = scale * np.exp(value)
            elif rule["transform"] == "log1p_nonnegative": field = scale * np.expm1(value)
            elif rule["transform"] == "log_positive_offset": field = scale * np.exp(value) - rule["offset"]
            else:
                linear = scale * np.expm1(value)
                field = np.full_like(linear, -math.inf)
                positive = linear > 0.0
                field[positive] = np.log(linear[positive])
            below = field < rule["lower"]
            above = field > rule["upper"]
            nan = np.isnan(field)
            if np.any(below) or np.any(above) or np.any(nan):
                raise ValueError(
                    f"projected {name} violates its constitutive domain "
                    f"(coordinate=[{float(np.min(value)):.17g}, {float(np.max(value)):.17g}], "
                    f"physical=[{float(np.nanmin(field)):.17g}, {float(np.nanmax(field)):.17g}], "
                    f"below={int(np.count_nonzero(below))}, above={int(np.count_nonzero(above))}, "
                    f"nan={int(np.count_nonzero(nan))})"
                )
            physical.append(field)
        (self.ep_gp,self.rho_gp,self.epsp_acc_gp,self.u,logmem,available,embryo,
         stable,inactive,completion)=physical
        self.pd_state.log_delivery_memory = logmem
        self.pd_state.delivery_memory = np.where(np.isfinite(logmem),np.exp(np.minimum(logmem,709.0)),0.0)
        self.pd_state.available=available; self.pd_state.embryo=embryo; self.pd_state.stable=stable
        self.pd_state.inactive=inactive; self.pd_state.completion=completion

    def active_residual(self, a, b):
        offset = 0; by_field = {}
        for name, shape, _ in a.specification:
            size = int(np.prod(shape, dtype=int)); av=a.vector[offset:offset+size]; bv=b.vector[offset:offset+size]; offset += size
            by_field[name] = float(np.max(np.abs(av-bv) / np.maximum.reduce([np.abs(av),np.abs(bv),np.ones_like(av)]))) if size else 0.0
        return max(by_field.values(), default=0.0), by_field

    def project_active_state(self, start, first, second, horizon):
        """Project smooth population recurrences without crossing their bounds.

        The expected embryo/stable/inactive and completion fields have additive
        sources, so logarithmic tangent projection is singular when a component
        first leaves exact zero.  Their exact fixed-rate recurrence is affine.
        A three-point geometric sum is therefore the appropriate local model;
        exact midpoint/end cycle maps still qualify every proposed segment.
        Other constitutive coordinates retain ordinary tangent projection.
        """
        h = int(horizon)
        out = start.vector + h * (first.vector - start.vector)
        offset = 0
        geometric_names = {"embryo", "stable", "inactive", "completion"}
        for name, shape, _ in start.specification:
            size = int(np.prod(shape, dtype=int))
            sl = slice(offset, offset + size); offset += size
            if name == "u":
                # ``u`` is the converged quasistatic solve used only as the
                # next nonlinear-solver warm start. It is not constitutive
                # memory and must not acquire a projective tangent drift.
                out[sl] = second.vector[sl]
                continue
            if name not in geometric_names:
                continue
            x0, x1, x2 = start.vector[sl], first.vector[sl], second.vector[sl]
            if name == "completion":
                # Completion is the cycle maximum reconstructed from the
                # finite delivery-memory state, not an accumulated inventory.
                # Hold the last trained value; exact endpoint maps validate it.
                out[sl] = x2
                continue
            d1, d2 = x1 - x0, x2 - x1
            projected = x0 + h * d1
            informative = np.abs(d1) > 32.0 * np.finfo(float).eps * np.maximum(1.0, np.abs(x0))
            ratio = np.zeros_like(d1)
            np.divide(d2, d1, out=ratio, where=informative)
            geometric = informative & (ratio >= 0.0) & (ratio < 1.0 - 1e-12)
            if np.any(geometric):
                r = ratio[geometric]
                factor = np.ones_like(r)
                positive_r = r > 0.0
                factor[positive_r] = (-np.expm1(h * np.log(r[positive_r]))) / (1.0 - r[positive_r])
                projected[geometric] = x0[geometric] + d1[geometric] * factor
            # The closed affine recurrence is nonnegative analytically.  Remove
            # only cancellation at the exact constitutive boundary (tens of
            # machine eps), leaving any material negative proposal to fail.
            roundoff = 64.0 * np.finfo(float).eps * np.maximum.reduce(
                [np.ones_like(projected), np.abs(x0), np.abs(x1), np.abs(x2)]
            )
            projected[(projected < 0.0) & (projected >= -roundoff)] = 0.0
            out[sl] = projected
        return out

    def projective_curvature(self, start, first, second):
        """Second-difference metric on physical memory, excluding diagnostics.

        Quasistatic ``u`` is a solver warm start and completion is reconstructed
        from delivery memory; neither is an evolving constitutive coordinate.
        """
        offset = 0
        values = []
        for name, shape, _ in start.specification:
            size = int(np.prod(shape, dtype=int)); sl = slice(offset, offset + size); offset += size
            if name in {"u", "completion"} or size == 0:
                continue
            d1 = first.vector[sl] - start.vector[sl]
            d2 = second.vector[sl] - first.vector[sl]
            scale = np.maximum.reduce([np.abs(d1), np.abs(d2), np.ones_like(d1)])
            values.append(float(np.max(np.abs(d2 - d1) / scale)))
        return max(values, default=0.0)

    @staticmethod
    def _active_field(snapshot, vector, target):
        offset = 0
        for name, shape, _ in snapshot.specification:
            size = int(np.prod(shape, dtype=int))
            if name == target:
                return np.asarray(vector[offset:offset + size], float).reshape(shape)
            offset += size
        raise KeyError(target)

    def conservative_population_ledgers(self, start, end_vector, ledgers):
        """Close expected birth/healing ledgers from population conservation."""
        result = dict(ledgers)
        available0 = self._active_field(start, start.vector, "available")
        available1 = self._active_field(start, end_vector, "available")
        inactive0 = self._active_field(start, start.vector, "inactive")
        inactive1 = self._active_field(start, end_vector, "inactive")
        returned = float(self.patch.cfg.heal_return_fraction)
        if not (0.0 <= returned < 1.0):
            raise ValueError("population ledger closure requires heal_return_fraction in [0,1)")
        healed = (inactive1 - inactive0) / (1.0 - returned)
        born = -(available1 - available0) + returned * healed
        scale = np.maximum.reduce([np.ones_like(born), np.abs(available0), np.abs(available1)])
        tolerance = 128.0 * np.finfo(float).eps * scale
        for name, value in (("healed_cumulative", healed), ("born_cumulative", born)):
            clean = np.asarray(value, float).copy()
            clean[(clean < 0.0) & (clean >= -tolerance)] = 0.0
            if np.any(clean < 0.0):
                raise ValueError(f"projected {name} violates population conservation")
            result[name] = clean
        return result, {"healed_cumulative", "born_cumulative"}

    def protected_signatures(self):
        s = self.pd_state
        ledgers = {name: getattr(s, name) for name in LEDGER_NAMES}
        ledgers["plastic_work"] = self.plastic_work
        ledgers["external_ledgers"] = self.external_ledgers
        stochastic = {name: getattr(s, name) for name in STOCHASTIC_NAMES}
        stochastic["candidate_rng"] = self.patch._candidate_rng.bit_generator.state
        stochastic["event_rng"] = self.patch._event_rng.bit_generator.state
        topology = {name: getattr(s, name) for name in TOPOLOGY_NAMES}
        topology["mesh_nodes"] = np.asarray(self.mesh.nodes)
        topology["bonds"] = np.asarray(self.patch.bonds)
        return ProtectedSignatures(_digest(ledgers), _digest(stochastic), _digest(topology))

    def exact_private_cycle(self):
        # Geometry, patch operators, and the FEM cache are immutable in the
        # eligible regime and may be shared. Deep-copying them dominated real
        # 48x96 private-map wall time without adding isolation.
        clone = copy(self)
        clone.pd_state = deepcopy(self.pd_state)
        clone.ep_gp = self.ep_gp.copy(); clone.rho_gp = self.rho_gp.copy()
        clone.epsp_acc_gp = self.epsp_acc_gp.copy(); clone.u = self.u.copy()
        clone.external_ledgers = deepcopy(self.external_ledgers)
        before = clone.protected_signatures()
        start = clone.active_state()
        payload = clone.cycle_evaluator(clone)
        end = clone.active_state()
        after = clone.protected_signatures()
        if before != after or clone.cycles != self.cycles:
            raise RuntimeError("cycle evaluator mutated protected state or physical cycles")
        log_action = np.asarray(payload["log_birth_action"], float)
        return CycleEvaluation(start, end, log_action,
            dict(payload.get("ledger_increments", {})),
            np.asarray(payload.get("phase", []), float),
            np.asarray(payload.get("phase_log_birth_rate", []), float),
            dict(payload.get("diagnostics", {})),
            str(payload.get("transition_signature", "dormant")), before.topology)

    def commit_private_cycle(self, evaluation):
        self.restore_active_state(evaluation.state_end, evaluation.state_end.vector)
        self.commit_ledger_increments(evaluation.ledger_increments)

    def commit_ledger_increments(self, increments):
        for name, increment in increments.items():
            if name == "plastic_work": self.plastic_work += float(increment)
            elif hasattr(self.pd_state, name):
                setattr(self.pd_state, name, getattr(self.pd_state, name) + increment)
            else:
                self.external_ledgers[name] = self.external_ledgers.get(name, 0.0) + float(increment)

    def remaining_birth_actions(self):
        s = self.pd_state
        nodes = np.asarray(s.site_node_index, dtype=np.int64)
        available = np.asarray(s.site_status) == 0
        site_residual = np.full(len(nodes), math.inf)
        site_residual[available] = np.maximum(
            np.asarray(s.site_birth_threshold)[available]
            - np.asarray(s.birth_cumulative_hazard)[nodes[available]], 0.0)
        residual = np.full_like(s.birth_cumulative_hazard, math.inf, dtype=float)
        np.minimum.at(residual, nodes, site_residual)
        return residual

    def commit_birth_action(self, increment, cycles):
        inc = np.asarray(increment, float)
        s = self.pd_state
        if inc.shape == np.asarray(s.birth_cumulative_hazard).shape:
            node_inc = inc
        elif inc.shape == np.asarray(s.site_birth_threshold).shape:
            node_inc = np.zeros_like(s.birth_cumulative_hazard, dtype=float)
            np.maximum.at(node_inc, np.asarray(s.site_node_index, int), inc)
        else:
            raise ValueError("birth action shape does not match PD nodes or sites")
        s.birth_cumulative_hazard = np.asarray(s.birth_cumulative_hazard) + node_inc
        log_inc = np.where(node_inc > 0.0, np.log(node_inc), -math.inf)
        s.log_birth_cumulative_hazard = np.logaddexp(s.log_birth_cumulative_hazard, log_inc)

    def commit_log_birth_action(self, log_increment, cycles):
        log_inc = np.asarray(log_increment, float)
        s = self.pd_state
        if log_inc.shape == np.asarray(s.birth_cumulative_hazard).shape:
            node_log_inc = log_inc
        elif log_inc.shape == np.asarray(s.site_birth_threshold).shape:
            node_log_inc = np.full_like(s.birth_cumulative_hazard, -math.inf, dtype=float)
            nodes = np.asarray(s.site_node_index, int)
            for node in np.unique(nodes):
                node_log_inc[node] = np.max(log_inc[nodes == node])
        else:
            raise ValueError("log birth action shape does not match PD nodes or sites")
        total_log = np.logaddexp(np.asarray(s.log_birth_cumulative_hazard, float), node_log_inc)
        s.log_birth_cumulative_hazard = total_log
        representable = total_log >= math.log(np.nextafter(0.0, 1.0))
        linear = np.zeros_like(total_log)
        linear[representable] = np.exp(total_log[representable])
        s.birth_cumulative_hazard = linear

    def physical_cycles(self): return self.cycles
    def set_physical_cycles(self, cycles): self.cycles = float(cycles)


__all__ = ["SpatialPDDormantAdapter", "ACTIVE_NAMES", "FIELD_RULES", "LEDGER_NAMES",
           "STOCHASTIC_NAMES", "TOPOLOGY_NAMES"]
