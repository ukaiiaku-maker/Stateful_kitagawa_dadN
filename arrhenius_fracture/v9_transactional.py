"""Versioned transactional event engine and complete v9 state capsule."""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Protocol
import uuid


SOLVER_VERSION = "STATEFUL_PD_V9_TRANSACTIONAL_1"


def canonical_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


@dataclass(frozen=True)
class PersistentClock:
    identity: str
    event_kind: str
    threshold: float
    cumulative_hazard: float = 0.0
    active: bool = True


@dataclass(frozen=True)
class V9StateCapsule:
    cycle: float
    fem_plastic_state_json: str
    rho_backstress_shielding_json: str
    delivery_memory_json: str
    completion_json: str
    candidate_sites_json: str
    clocks: tuple[PersistentClock, ...]
    embryo_stable_transition_json: str
    pd_damage_front_json: str
    rng_states_json: str
    rng_digests_json: str
    geometry_identity_state_json: str
    source_hashes_json: str
    request_hash: str
    solver_version: str = SOLVER_VERSION
    aggregate_cumulative_hazard: float = 0.0
    event_sequence: tuple[tuple[str, str, float], ...] = ()
    cap_floor_diagnostics: tuple[str, ...] = ()
    commit_count: int = 0
    committed_cycle_measure: float = 0.0
    controller_next_block_cycles: float = 0.0

    def __post_init__(self) -> None:
        if self.cycle < 0.0 or self.aggregate_cumulative_hazard < 0.0:
            raise ValueError("cycle and cumulative hazard must be nonnegative")
        for field in (
            self.fem_plastic_state_json, self.rho_backstress_shielding_json,
            self.delivery_memory_json, self.completion_json, self.candidate_sites_json,
            self.embryo_stable_transition_json, self.pd_damage_front_json,
            self.rng_states_json, self.rng_digests_json,
            self.geometry_identity_state_json, self.source_hashes_json,
        ):
            json.loads(field)

    @property
    def identity_digest(self) -> str:
        return hashlib.sha256(canonical_json(self.to_dict()).encode()).hexdigest()

    def to_dict(self) -> dict[str, object]:
        return asdict(self)

    @classmethod
    def from_dict(cls, value: dict[str, object]) -> "V9StateCapsule":
        data = dict(value)
        data["clocks"] = tuple(PersistentClock(**clock) for clock in data["clocks"])
        data["event_sequence"] = tuple(tuple(event) for event in data["event_sequence"])
        data["cap_floor_diagnostics"] = tuple(data["cap_floor_diagnostics"])
        return cls(**data)


def empty_capsule(*, clocks: tuple[PersistentClock, ...], request_hash: str = "synthetic") -> V9StateCapsule:
    empty = canonical_json({})
    return V9StateCapsule(
        cycle=0.0,
        fem_plastic_state_json=canonical_json({"x": 0.0}),
        rho_backstress_shielding_json=empty,
        delivery_memory_json=empty,
        completion_json=empty,
        candidate_sites_json=empty,
        clocks=clocks,
        embryo_stable_transition_json=empty,
        pd_damage_front_json=empty,
        rng_states_json=empty,
        rng_digests_json=empty,
        geometry_identity_state_json=canonical_json({"kind": "synthetic"}),
        source_hashes_json=canonical_json({"v9_transactional": "synthetic"}),
        request_hash=request_hash,
    )


@dataclass(frozen=True)
class Proposal:
    dN: float
    fem_plastic_state_json: str
    clock_hazard_increments: tuple[float, ...]
    aggregate_hazard_increment: float
    error_estimate: float
    cap_floor_activations: tuple[str, ...] = ()


class TransactionModel(Protocol):
    def propose(self, state: V9StateCapsule, dN: float) -> Proposal: ...
    def crossing_fraction(self, state: V9StateCapsule, clock_index: int, target_increment: float, dN: float) -> float: ...


@dataclass(frozen=True)
class EngineConfig:
    initial_block_cycles: float = 1.0
    quiet_growth_factor: float = 10.0
    error_tolerance: float = 1.0e-10
    minimum_block_cycles: float = 1.0e-12
    maximum_rejections: int = 128


class TransactionalEngine:
    def __init__(self, model: TransactionModel, config: EngineConfig = EngineConfig()):
        self.model = model
        self.config = config

    def advance(self, initial: V9StateCapsule, cycle_end: float) -> V9StateCapsule:
        state = initial
        next_block = state.controller_next_block_cycles or self.config.initial_block_cycles
        candidate = min(next_block, cycle_end - state.cycle)
        rejections = 0
        while state.cycle < cycle_end:
            dN = min(candidate, cycle_end - state.cycle)
            proposal = self.model.propose(state, dN)
            if proposal.error_estimate > self.config.error_tolerance:
                dN *= 0.5
                candidate = dN
                rejections += 1
                if dN < self.config.minimum_block_cycles or rejections > self.config.maximum_rejections:
                    raise RuntimeError("transactional block subdivision failed")
                continue

            earliest = 1.0
            for index, (clock, increment) in enumerate(zip(state.clocks, proposal.clock_hazard_increments)):
                remaining = clock.threshold - clock.cumulative_hazard
                if clock.active and increment > 0.0 and increment >= remaining:
                    fraction = self.model.crossing_fraction(state, index, remaining, dN)
                    earliest = min(earliest, max(0.0, min(1.0, fraction)))
            if earliest < 1.0 - 1.0e-14:
                dN *= earliest
                proposal = self.model.propose(state, dN)

            state = self._commit_once(state, proposal)
            rejections = 0
            candidate = max(self.config.minimum_block_cycles, dN * self.config.quiet_growth_factor)
            state = replace(state, controller_next_block_cycles=candidate)
        return state

    @staticmethod
    def _commit_once(state: V9StateCapsule, proposal: Proposal) -> V9StateCapsule:
        new_cycle = state.cycle + proposal.dN
        clocks = []
        events = list(state.event_sequence)
        for clock, increment in zip(state.clocks, proposal.clock_hazard_increments):
            effective_increment = increment if clock.active else 0.0
            cumulative = clock.cumulative_hazard + effective_increment
            crossed = clock.active and cumulative >= clock.threshold - 64.0 * math.ulp(max(1.0, clock.threshold))
            if crossed:
                cumulative = clock.threshold
                events.append((clock.event_kind, clock.identity, new_cycle))
            clocks.append(replace(clock, cumulative_hazard=cumulative, active=clock.active and not crossed))
        diagnostics = tuple(dict.fromkeys(state.cap_floor_diagnostics + proposal.cap_floor_activations))
        return replace(
            state,
            cycle=new_cycle,
            fem_plastic_state_json=proposal.fem_plastic_state_json,
            clocks=tuple(clocks),
            aggregate_cumulative_hazard=state.aggregate_cumulative_hazard + proposal.aggregate_hazard_increment,
            event_sequence=tuple(events),
            cap_floor_diagnostics=diagnostics,
            commit_count=state.commit_count + 1,
            committed_cycle_measure=state.committed_cycle_measure + proposal.dN,
        )


class AnalyticHazardModel:
    """Synthetic exact model used for numerical proofs, not production physics."""

    def __init__(self, hazards: tuple[tuple[float, float], ...], *, state_rate: float = 1.0, state_cap: float | None = None):
        self.hazards = hazards  # h(N) = intercept + slope*N for each clock
        self.state_rate = state_rate
        self.state_cap = state_cap

    def _integral(self, index: int, n0: float, n1: float) -> float:
        intercept, slope = self.hazards[index]
        return intercept * (n1 - n0) + 0.5 * slope * (n1 * n1 - n0 * n0)

    def propose(self, state: V9StateCapsule, dN: float) -> Proposal:
        n1 = state.cycle + dN
        increments = tuple(self._integral(i, state.cycle, n1) for i in range(len(state.clocks)))
        physical = json.loads(state.fem_plastic_state_json)
        x = float(physical.get("x", 0.0)) + self.state_rate * dN
        activations = ()
        if self.state_cap is not None and x > self.state_cap:
            x = self.state_cap
            activations = ("synthetic_state_cap",)
        return Proposal(dN, canonical_json({"x": x}), increments, sum(increments), 0.0, activations)

    def crossing_fraction(self, state: V9StateCapsule, clock_index: int, target_increment: float, dN: float) -> float:
        intercept, slope = self.hazards[clock_index]
        n0 = state.cycle
        rate0 = intercept + slope * n0
        if abs(slope) < 1.0e-30:
            delta = target_increment / rate0
        else:
            disc = rate0 * rate0 + 2.0 * slope * target_increment
            delta = (-rate0 + math.sqrt(max(disc, 0.0))) / slope
        return delta / dN


class AtomicGenerationStore:
    """Atomic checkpoint/summary/manifest files; manifest is published last."""

    def __init__(self, root: Path):
        self.root = Path(root)

    def write_generation(self, state: V9StateCapsule, summary: dict[str, object]) -> str:
        self.root.mkdir(parents=True, exist_ok=True)
        generation = f"g{state.commit_count:08d}-{uuid.uuid4().hex[:12]}"
        checkpoint_name = f"checkpoint-{generation}.json"
        summary_name = f"summary-{generation}.json"
        checkpoint_bytes = (canonical_json(state.to_dict()) + "\n").encode()
        summary_bytes = (canonical_json(summary | {"generation": generation}) + "\n").encode()
        self._atomic_write(checkpoint_name, checkpoint_bytes)
        self._atomic_write(summary_name, summary_bytes)
        manifest = {
            "generation": generation,
            "solver_version": state.solver_version,
            "request_hash": state.request_hash,
            "state_digest": state.identity_digest,
            "checkpoint": checkpoint_name,
            "checkpoint_sha256": hashlib.sha256(checkpoint_bytes).hexdigest(),
            "summary": summary_name,
            "summary_sha256": hashlib.sha256(summary_bytes).hexdigest(),
        }
        self._atomic_write("manifest-latest.json", (canonical_json(manifest) + "\n").encode())
        return generation

    def load_latest(self) -> V9StateCapsule:
        manifest = json.loads((self.root / "manifest-latest.json").read_text())
        raw = (self.root / manifest["checkpoint"]).read_bytes()
        if hashlib.sha256(raw).hexdigest() != manifest["checkpoint_sha256"]:
            raise RuntimeError("checkpoint generation hash mismatch")
        state = V9StateCapsule.from_dict(json.loads(raw))
        if state.identity_digest != manifest["state_digest"]:
            raise RuntimeError("state digest mismatch")
        return state

    def _atomic_write(self, name: str, payload: bytes) -> None:
        temporary = self.root / f".{name}.{uuid.uuid4().hex}.tmp"
        with temporary.open("wb") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, self.root / name)
