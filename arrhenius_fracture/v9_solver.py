"""Public v9 solver facade over log rates and transactional state commits.

This first architecture gate intentionally supplies no v8.7 FEM/PD mutation
adapter.  A physical adapter must populate every state section explicitly; the
facade refuses implicit reconstruction of stochastic or geometry state.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib

from .v9_transactional import (
    AtomicGenerationStore,
    EngineConfig,
    PersistentClock,
    SOLVER_VERSION,
    TransactionModel,
    TransactionalEngine,
    V9StateCapsule,
    canonical_json,
)


BASELINE_MECHANISM = "completion_gated_independent_cleavage"


@dataclass(frozen=True)
class V9SolverRequest:
    mechanism: str = BASELINE_MECHANISM
    temperature_K: float = 300.0
    frequency_Hz: float = 1000.0
    source_hashes: tuple[tuple[str, str], ...] = ()
    configuration_version: str = SOLVER_VERSION

    @property
    def request_hash(self) -> str:
        return hashlib.sha256(canonical_json(asdict(self)).encode()).hexdigest()


def build_physical_capsule(
    request: V9SolverRequest,
    *,
    cycle: float,
    fem_plastic_state: object,
    rho_backstress_shielding: object,
    delivery_memory: object,
    completion: object,
    candidate_sites: object,
    clocks: tuple[PersistentClock, ...],
    embryo_stable_transition: object,
    pd_damage_front: object,
    rng_states: object,
    rng_digests: object,
    geometry_identity_state: object,
) -> V9StateCapsule:
    required = {
        "fem_plastic_state": fem_plastic_state,
        "rho_backstress_shielding": rho_backstress_shielding,
        "delivery_memory": delivery_memory,
        "completion": completion,
        "candidate_sites": candidate_sites,
        "embryo_stable_transition": embryo_stable_transition,
        "pd_damage_front": pd_damage_front,
        "rng_states": rng_states,
        "rng_digests": rng_digests,
        "geometry_identity_state": geometry_identity_state,
    }
    missing = [name for name, value in required.items() if value is None]
    if missing:
        raise ValueError("authoritative v9 state sections missing: " + ", ".join(missing))
    return V9StateCapsule(
        cycle=cycle,
        fem_plastic_state_json=canonical_json(fem_plastic_state),
        rho_backstress_shielding_json=canonical_json(rho_backstress_shielding),
        delivery_memory_json=canonical_json(delivery_memory),
        completion_json=canonical_json(completion),
        candidate_sites_json=canonical_json(candidate_sites),
        clocks=clocks,
        embryo_stable_transition_json=canonical_json(embryo_stable_transition),
        pd_damage_front_json=canonical_json(pd_damage_front),
        rng_states_json=canonical_json(rng_states),
        rng_digests_json=canonical_json(rng_digests),
        geometry_identity_state_json=canonical_json(geometry_identity_state),
        source_hashes_json=canonical_json(dict(request.source_hashes)),
        request_hash=request.request_hash,
        solver_version=request.configuration_version,
    )


def advance_v9(
    state: V9StateCapsule,
    model: TransactionModel,
    cycle_end: float,
    *,
    engine_config: EngineConfig = EngineConfig(),
    generation_store: AtomicGenerationStore | None = None,
) -> V9StateCapsule:
    final = TransactionalEngine(model, engine_config).advance(state, cycle_end)
    if generation_store is not None:
        generation_store.write_generation(final, {
            "status": "accepted_interval_endpoint",
            "cycle": final.cycle,
            "mechanism": BASELINE_MECHANISM,
            "birth_endurance_diagnostic": "not_evaluated_by_stepper",
            "physical_handoff_endurance_classification": "undetermined",
        })
    return final
