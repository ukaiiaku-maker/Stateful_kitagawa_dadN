# V9 autonomous S-N build progress

## Current branch state

- Branch: `codex/sn-endurance-limit-v9`
- Autonomous handoff merged from `origin/codex/sn-endurance-limit`.
- Frozen v8.7 driver, PD and geometry sources remain unchanged.
- No K360 or campaign condition has been launched before real-trajectory qualification.

## Completed in the current work package

The immutable v2.8 K360 requests are recoverable from the committed orchestrator. The failure and long-life anchor definitions include exact stresses, temperature, load ratio, frequency, geometry, mesh/PD resolution, random seeds, shielding parameters, and progression parameters; no physical parameter needs to be invented.

The first real numerical layer is now implemented in `v9_physical_integrator.py`:

- exact frozen EXP-floor barrier evaluation without exponent clipping;
- series residence-time plastic rates in log space;
- exact constant-delivery finite-memory evolution in log space;
- cancellation-safe `log Q(2,Lambda)` below float rate range;
- log-domain Gauss integration of completion-gated cleavage hazard;
- persistent integrated-hazard clocks with localized crossings;
- persistent competing-risk selection from preassigned uniforms;
- semigroup, restart, >10x partition and real-v8.7 rate-equivalence tests.

The v9 Boltzmann constant now uses the same exact repository `KB/EV_TO_J` value as the frozen barrier implementation. The prior rounded decimal caused approximately `7e-13` relative drift in otherwise identical representable rates.

## Physical migration decisions resolved from the frozen equations

- Birth remains a persistent exponential first-passage clock.
- Embryo stabilization/healing is an exact competing-risk process: one persistent total-hazard threshold and a preassigned outcome uniform replace block-local redraws.
- Growth and bond linkage are not converted into new random events. Their frozen updates are deterministic survival-state equations. V9 stores/adds their cumulative hazards (`H=-log(1-x)`) and localizes the existing growth/damage state boundaries. This preserves the intended model while removing partition dependence.
- Float underflow in a state increment is recorded as a representability condition, not physical extinction. Log hazard/memory state continues to advance.

## Qualification status

| Gate | Status |
|---|---|
| Log temporal primitives | passed focused tests |
| Persistent competing clocks | passed focused tests |
| Frozen representable-rate equivalence | passed |
| Full FEM/rho embedded transaction | in progress |
| Real shortened partition proof | not yet qualified |
| Real atomic restart proof | not yet qualified |
| K360 735.921 MPa | not started; correctly gated |
| K360 690.443 MPa | not started; correctly gated |

## Next automatic action

Implement the full FEM/rho/PD proposal-reject-commit bridge. It must use embedded full/half-step errors, localize rho/state/cap boundaries, advance deterministic progression hazards once, and serialize the new persistent transition state through atomic array generations. Then rerun the shortened real >10x partition and restart proofs before starting K360.
