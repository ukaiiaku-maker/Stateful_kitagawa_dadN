# V9 numerical architecture audit — first gate

## Scope and immutable baseline

This gate introduces new v9 modules only:

```text
arrhenius_fracture/v9_log_rates.py
arrhenius_fracture/v9_transactional.py
arrhenius_fracture/v9_solver.py
```

The three v8.7 source files remain byte-identical. No production S-N run, 3x3x2 pilot, parameter optimization, or K360 rerun/continuation was started.

The baseline v9 mechanism name is:

```text
completion_gated_independent_cleavage
```

`persistent_direct_cleavage` remains only a proposed comparison name. It has not been implemented and carries no assumed endurance class.

## Phase 1: log-domain physical rates

`LogPhysicalRate` represents a physical rate as a natural logarithm without evaluating `exp(log_rate)` when float64 cannot represent the result. Diagnostics distinguish:

- `log_rate` and units;
- `physical_rate` when representable, otherwise `null` rather than a numerical zero/floor;
- `numerical_representability` (`representable`, float overflow, float underflow, or exact constitutive zero);
- exact `constitutive_zero` and its physical reason;
- cap/floor activation names;
- an explicit assertion that no numerical floor was used as physics.

Arrhenius rates use

```text
log r = log(nu0) - DeltaG/(kB T)
```

with no `[-700,0]` exponent clip. Series emission/mobility rates use log-domain residence-time addition, not `max(rate,1e-300)`. Completion gating adds `log Q`; `Q=0` is exact only when the constitutive completion law returns exact zero. Stabilization, healing, growth and linkage in v8.7 are logistic kinetic rates rather than Arrhenius barriers; v9 evaluates their log-sigmoid forms without underflow. Exact zero growth/link gates require an explicit zero state/activity gate.

### Complete rate-path trace

`trace_completion_gated_rate_paths` reports:

```text
plastic_emission
plastic_peierls
plastic_taylor
plastic_escape
plastic_completed_delivery
cleavage
birth
stabilization
healing
stable_growth
linkage_front
```

A shortened-real-case test constructs the actual v8.7 default parser, case-64-M1 barriers and plastic chain, evaluates their barriers at a finite stress/state, and traces every path into v9 diagnostics. This is read-only adaptation; v8.7 rate methods are not modified.

### Activity-extinction contract

The contract is intentionally strict:

1. `rate below epsilon` is never extinction.
2. Float underflow is never extinction and yields `physical_rate=null` plus a retained finite log-rate.
3. Exact extinction requires a constitutive statement: disabled mechanism, nonpositive physical prefactor, exact zero completion, or exact zero state/activity gate, with a reason string.
4. If the constitutive model stays positive, v9 retains the log-rate and endurance must be classified from an analytical/log-domain integrated-hazard bound.

`classify_exact_power_law_from_logs` verifies that endurance classification is invariant to arbitrary log-rate presentation offsets, including histories entirely below float64 rate range.

## Phase 2: transactional large-N stepping

`TransactionalEngine` implements:

```text
immutable capsule -> propose -> error/event checks -> reject/subdivide or truncate -> accept once
```

Candidate blocks grow geometrically in quiet regions. A proposed block contains predicted physical state, integrated hazard for every persistent clock, aggregate hazard, an error estimate and cap/floor activations. Accuracy comes from model error estimates and event truncation, not a fixed maximum `dN`.

For every proposal the engine:

1. rejects and halves the interval when the error estimate exceeds tolerance;
2. tests every active birth, embryo-transition and progression clock;
3. asks the model for the crossing fraction and selects the earliest event;
4. recomputes the proposal on the truncated interval;
5. commits state, hazards and events exactly once;
6. records `commit_count` and `committed_cycle_measure` for double-commit auditing;
7. retains the next adaptive candidate block in the capsule for exact restart scheduling.

The architecture supports arbitrary physical models through the `TransactionModel` protocol. `AnalyticHazardModel` is a proof fixture, not production physics. It provides exact constant/affine-in-cycle hazard integrals and inverses.

## Phase 3: complete state capsule and atomic generations

`V9StateCapsule` is frozen and versioned as `STATEFUL_PD_V9_TRANSACTIONAL_1`. It includes explicit sections for:

- cycle coordinate and adaptive-controller state;
- FEM/plastic state;
- density, back stress and shielding;
- delivery memory Lambda and completion Q;
- candidate identities/status;
- persistent clocks with thresholds and cumulative hazards;
- aggregate cumulative hazard;
- embryo/stable/transition state;
- PD damage/front state;
- complete RNG states and digests;
- geometry identity/state;
- source hashes, request hash and solver/configuration version;
- ordered events, commit accounting and cap/floor diagnostics.

State sections are canonical JSON strings inside the frozen capsule. This prevents the facade from partially guessing fields while remaining portable for later NumPy array codecs.

`build_physical_capsule` requires every authoritative section and fails closed when any section is missing. Request identity is a SHA-256 of canonical request/configuration content.

`AtomicGenerationStore` writes checkpoint and summary to temporary files, flushes/fsyncs and atomically replaces their generation names. A hash-addressed manifest is published last. Loading verifies checkpoint bytes and the canonical state digest. Summary/checkpoint/manifest generation IDs are identical.

## Phase 4: numerical proofs

Tests added under `tests/test_v9_log_rates.py` and `tests/test_v9_transactional.py` cover:

### A. Block-partition invariance

The same affine-hazard/state trajectory uses schedules beginning at `0.1` and `10` cycles with different geometric growth (more than 10x schedule separation). Endpoint cycle, physical state, aggregate hazard and persistent-clock hazard agree to test tolerance.

### B. Exact threshold crossing

- Constant `h=2`, threshold 5: event at `N=2.5`.
- Time-varying `h=2N`, threshold 9: event at `N=3`.
- Multiple clocks verify the earliest embryo event truncates before later birth/progression events.

### C. Restart equivalence

An uninterrupted trajectory is compared with interruption at an accepted boundary, atomic generation write/load and resume. Physical state, thresholds/cumulative hazards, event sequence, endpoint and RNG state section agree.

### D. No double state commit

For every trajectory, `committed_cycle_measure == final_cycle - initial_cycle`. Rejected proposals do not increment commit count or physical state. A forced-error model proves subdivision occurs before commit.

### E. Log-rate invariance

Arrhenius/logistic rates far below float64 range remain finite log-rates and nonzero physical mechanisms. Adding an arbitrary `+1200` log presentation offset does not change the exact power-law exponent or endurance class.

### F. Cap/floor activation

A synthetic state cap is explicitly emitted in capsule diagnostics. The capped and uncapped cases retain identical hazard integrals when the cap is not part of the hazard law, demonstrating that cap reporting is separate from asymptotic classification. The shortened v8.7 trace propagates constitutive floor activation names without converting them into rate floors.

## Preserved K360 continuation audit

Exact continuation is **not possible** from the committed K360 reference bundle. The directory contains only:

```text
README.md
K360_CENSOR_TAIL_AUDIT.json
K360_CENSOR_DERIVED_TAIL.csv
```

These files preserve derived aggregate history and fit evidence, not an authoritative solver checkpoint. Exact v9 continuation requires at least the following missing state:

- FEM integration-point plastic strain, accumulated plastic strain, density and residual fields;
- current mesh nodes, feature/root geometry and local PD geometry identity;
- delivery memory Lambda and completion state for every PD point;
- candidate-site identities/counts/statuses;
- persistent exponential thresholds and cumulative per-site hazards;
- embryo, stable, inactive, growth and transition-clock states;
- bond damage, active-front path/tip/contact, candidate modes and morphology counters;
- candidate and event RNG bit-generator states plus digests;
- adaptive block/controller state and exact accepted cycle coordinate;
- full run arguments, request/signature hash, model/source hashes and checkpoint version.

The derived cumulative expected-birth history cannot reconstruct these quantities uniquely. In particular, inventing thresholds or RNG streams would destroy restart equivalence and stochastic identity. Therefore Phase 5 stops at insufficiency reporting; no K360 state was silently reconstructed and no continuation was run.

## Departures and limitations

This gate establishes the numerical kernel and proofs, not the production physical adapter:

- The transactional protocol is exercised with exact synthetic hazards and a read-only shortened v8.7 rate trace. It does not yet advance full FEM/PD arrays.
- Error estimates are supplied by the model protocol. A future physical adapter must implement embedded state/hazard/memory error bounds and cap-boundary localization.
- Generic clocks support birth, embryo and progression events, but the v8.7 stochastic transition laws have not yet been migrated to persistent v9 transition clocks.
- Capsule array storage currently uses canonical JSON sections; production-scale arrays should use an atomic NPZ/Zarr-like generation codec while retaining the same canonical manifest identity.
- Physical-handoff endurance remains `undetermined`; the stepper does not infer it from birth events.
- A fixed `minimum_block_cycles` exists solely as failure protection against endless subdivision. There is no fixed maximum block controlling accuracy.

These limitations are deliberate stop conditions. Coupling full v8.7 physics before the numerical contracts pass would make failures hard to localize.

## Next reviewed gate

After review, the next implementation should add a read-only-v8.7/new-v9 physical adapter with:

1. log-domain phase integration for plastic delivery and cleavage;
2. log-domain completion/memory bounds for large blocks;
3. embedded FEM/plastic/PD state error estimators;
4. persistent embryo/progression clocks and conservative event localization;
5. production array generation codec;
6. shortened real FEM/PD block-partition and checkpoint/restart equivalence.

Only an authoritative checkpoint containing the complete capsule fields may enable K360 continuation. The 3x3x2 pilot remains prohibited.
