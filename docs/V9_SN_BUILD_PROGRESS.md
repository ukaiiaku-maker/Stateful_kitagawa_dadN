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
| K360 735.921 MPa | valid physical handoff at `N=3,474,029.481029375` |
| K360 690.443 MPa | valid right-censored long-life diagnostic at `N=1e8`; no realized birth |
| S001 804.139 MPa | valid physical handoff at `N=1,632,430.857` |
| S002 906.465 MPa | valid physical handoff at `N=1,530,050.877` |

## First physical anchor result

The K360 shielded `735.9214951373579 MPa` condition completed after 852 accepted
blocks at `N=3,474,029.481029375`.  The terminal classification is
`physical_handoff`, not a numerical or geometry termination.  The crack was
root-connected and active with 75 connected bonds, a `62.625 um` centerline,
slenderness `8.35`, width ratio `0.2395`, orientation coherence `0.9784`, no
off-front broken bonds, and `527.375 um` boundary clearance.  Every physical
handoff sub-gate and the production geometry audit passed.

The event differs from the historical v8.7 reference near `3.696e6` cycles by
about `-6.0%`; no parameter was fitted to reduce that discrepancy.  Persistent
birth occurred at `N=1,793,122.651`, front capture at `N=2,351,351.153`, and the
physical-handoff length was reached at the terminal accepted boundary.

The terminal schema-9 ACTIVE generation is an accepted physical boundary at
exactly the reported cycle coordinate.  It contains all 46 array components,
including log delivery memory, log cumulative birth hazard, persistent
transition thresholds/outcome uniforms, and both independent RNG states.

## Long-life anchor result and performance hardening

The exact K360 shielded `690.4432004940379 MPa` condition reached `N=1e8` with
no realized birth, no broken bond, valid fixed geometry, and status
`right_censored`.  Its cumulative expected births were `0.627593637`; the
terminal maximum birth intensity was `6.979e-13/cycle`, delivery memory
`4.986e-6`, and completion `Q=1.243e-11`.  Nested effective-hazard tail fits
over the last 50%, 35%, and 20% of points gave exponents `2.049`, `2.285`, and
`2.335` with extrapolated remaining hazards `0.00500`, `0.00355`, and `0.00337`.
This is empirical integrable-tail evidence, while the physics-asymptotic
classification remains explicitly `undetermined`; finite runout is not called
endurance.

The completed trajectory used 971 accepted blocks and a largest accepted block
of `560,889.6` cycles.  Profiling exposed two v9 controller defects: discarded
next-step estimates caused repeated FEM solves, and a near-zero `ep_gp`
component with `1e-18` absolute tolerance controlled quiet-tail steps.  The
controller candidate is now checkpointed and the strain absolute tolerance is
`1e-8` (about 4.1 kPa at 410 GPa), with active-component relative error, rho,
hazard, and event controls unchanged.  The production tolerance was compared
against `1e-9`; the complete 116-test regression suite passes.  Real accepted
blocks consequently grew to roughly `2e5` cycles through the late tail.

## Next automatic action

The first bootstrap condition completed with 75 connected bonds, a `62.625 um`
centerline, width ratio `0.2395`, orientation coherence `0.9784`, and every
handoff/geometry gate passing.  It used 182 accepted blocks and 28 rejected
proposals.  The next one-at-a-time choice uses the two finite lives and bounds
upper-stress extrapolation to 1.5 times their stress spacing, proposing
`906.4651000498077 MPa` with predicted `log10(N)=5.721`.  After it completes,
recompute the selector from all valid conditions.

S002 completed with first stable birth at `N=5,602.142`, front capture at
`N=705,449.704`, and physical handoff at `N=1,530,050.877`.  The terminal crack
had 43 broken bonds, `63.173 um` centerline, width ratio `0.2082`, orientation
coherence `0.9718`, and all handoff/geometry gates passed.  Five births and two
primary-seed reselections were preserved.  The weak life reduction from S001
despite much earlier nucleation indicates a front-progression-controlled
finite-life plateau under the frozen physics.  The bounded selector proposes
`1059.9543444710125 MPa` next to test whether shorter-life decades remain
physically accessible.

During S002, a hot-stack capture found exact phase integration below an
arbitrary 1024-cycle acceleration threshold.  At `tau=1 ms` and `f=1 kHz`, the
memory map contracts by `e^-1` per cycle, so the threshold is now 32 cycles
(`exp(-32)<1.3e-14`).  A 100-cycle accelerated calculation matches four direct
25-cycle spans within `2e-10` in log hazard.  S002 then completed in 287
accepted blocks with 30 rejected proposals.

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
