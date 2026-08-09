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

The real v9 driver and PD bridge now additionally provide:

- side-effect-free Heun/Euler FEM, plastic-strain and rho proposals;
- componentwise embedded error rejection before any PD/RNG commit;
- rho cap/floor boundary rejection instead of accepted-state clipping;
- a persistent physical cycle-phase coordinate for fractional blocks;
- exact ordered log-domain phase continuation for Lambda and gated cleavage;
- converged periodic-orbit acceleration with a geometric memory remainder;
- schema-9 persistent embryo competing-risk thresholds/outcome uniforms;
- hash-verified atomic array generations used as the resume authority.

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
| Full FEM/rho embedded transaction | implemented; PD error coupling still being extended |
| Real shortened partition proof | passed at production K360 geometry for `dN_max=1` versus `0.05` |
| Real atomic restart proof | passed at production K360 geometry through schema-9 ACTIVE generation |
| K360 735.921 MPa | active; restartable generation at `N=9478.027455854328` |
| K360 690.443 MPa | not started; correctly gated |

## Next automatic action

Run the K360 735.921 MPa anchor from its immutable request with the now-qualified schema-9 driver. Monitor accepted/rejected blocks and periodic remainder behavior; preserve and resume the newest valid ACTIVE generation for any implementation interruption.

The K360 failure anchor has now started and published a valid ACTIVE generation at `N=9478.027455854328` (10 accepted blocks, no birth). Profiling showed repeated assembly/factorization of the unchanged fixed-geometry linear stiffness dominated wall time. The run was interrupted only after that atomic generation was present. A v9-only cached FEM implementation reuses the immutable stiffness and sparse factorization; a production-geometry A/B checkpoint at `N=0.1` agrees with the uncached implementation to about `1e-13` relative or better across FEM, Lambda and birth hazard. The exact predecessor source hashes are explicitly allowlisted for this verified cache-only migration so the K360 generation is resumed rather than restarted.

The resumed anchor reached a persistent birth/stable event near the historical region, captured a localized active front, and published another valid generation at `N=2.588725650997854e6` (block 510). At that boundary it had 38 broken bonds, all active-front stall counters remained clear, and no handoff was yet claimed. A second fixed-geometry optimization replaces repeated phase-by-phase linear solves with the exact affine displacement/reaction response; tests compare the complete phase stress histories against repeated solves. The small differences are sparse-solve roundoff (stress relative errors below about `2e-10`), while a production endpoint comparison for the preceding cache migration was approximately `1e-13` relative across physical and endurance state.

### Real shortened proof evidence

At `N=1` on the production K360 mesh, initial candidate maxima `1.0` and `0.05` cycles produced identical RNG states and no events. Maximum relative endpoint differences were approximately:

- FEM plastic strain: `2.48e-11`;
- rho: `1.22e-15`;
- log delivery memory: `1.55e-8` absolute (`8.78e-10` scaled);
- log cumulative birth hazard: `3.06e-7` absolute (`5.26e-9` scaled);
- linear cumulative birth hazard: `3.06e-7` relative.

The prior 77–99% Lambda/Q/hazard discrepancy was traced to v8.7 restarting phase zero at every fractional block. Persistent phase continuation removed it. An interruption at `N=0.1`, atomic reload, and continuation to `N=0.2` succeeded from the ACTIVE generation with all component hashes, schema, stochastic arrays and cycle identity verified.
