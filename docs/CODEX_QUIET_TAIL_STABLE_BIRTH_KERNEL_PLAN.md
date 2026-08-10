# Autonomous plan: quiet-tail acceleration and stable-crack-birth renewal kernel at 300 K

## Governing scope

Continue from the qualified model:

```text
v9_canonical_four_class_energy_gated_stable_birth_v2
```

300 K only.

The endpoint is the first **energy-admissible stable crack birth**. Cleavage first-passage attempts that are rejected by the existing post-passage fixed-opening energy balance are consumed nonpropagating attempts and do not terminate the S-N condition.

Do not reintroduce post-birth crack growth, handoff, or da/dN.

## Current status to preserve

Preserve the completed cleavage-attempt tails and the current endpoint implementation.

The existing low-stress 1e8 trajectories are now diagnostics for the deterministic pre-birth physical state and cleavage-attempt process. They must not be discarded or restarted merely because the stable-birth estimator is still pending.

The previous positive intrinsic cleavage floor remains a valid statement about cleavage attempts, not stable-crack birth.

## Primary technical objective

Do not brute-force production FEM cycle blocks from 1e8 to 1e14.

Build a **phase-resolved quiet-tail kernel** for the deterministic pre-birth state and use it to evaluate the stochastic marked-renewal stable-birth process efficiently and exactly enough to report survival at 1e10, 1e12, and 1e14.

The central physical question is what the post-first-passage energy gate does after cyclic shakedown.

For each class/stress, determine whether the late-time periodic state has:

1. an energy gate closed for every possible event phase and every threshold-scaled proposal -> endurance-capable / potentially strict endurance;
2. a nonzero-measure set of event phase + Xi combinations with admitted crack length > 0 -> eventual stable birth if cleavage attempts continue indefinitely;
3. an unresolved/crossing regime -> asymptotically undetermined.

## Phase A — deterministic late-cycle kernel

For each low-stress condition, represent the accepted pre-birth trajectory independently of stochastic threshold sequence.

At decade checkpoints and especially in the late tail, record a complete phase-resolved cycle map including:

```text
phase
root opening / stress tensor
signed resolved channel stresses
mobile positive/negative arrays
retained positive/negative arrays
wake state
accumulated slip
backstress
K_shield
blunted radius
raw cleavage log-rate
renewed cleavage log-rate
cumulative cleavage action within the cycle
```

Use the exact accepted FEM + signed-MPZ ordering already qualified.

Characterize convergence toward a periodic orbit / stationary cycle map. Do not assume stationarity from one scalar hazard value.

Provide convergence errors for every active state family and for the phase-resolved cleavage action.

## Phase B — qualify the quiet-tail acceleration

Implement the least invasive acceleration that preserves the accepted deterministic pre-birth trajectory. Possible methods include a converged periodic-orbit map, fixed-point cycle solve, analytically accelerated signed-state update, or another method justified by the executable equations.

Do not prescribe a particular method if another is more exact.

Qualification must compare accelerated continuation against direct production-FEM continuation over overlapping windows.

Require bounds for:

```text
full signed MPZ state
FEM root response
backstress
shielding
blunting
phase-resolved cleavage action
cycle-integrated H increment
energy-gate admissibility envelope
```

Use more than one starting checkpoint where practical so an accidentally good local match cannot qualify the acceleration.

No tolerance may be loosened merely to reach 1e14.

## Phase C — energy-admissibility envelope

At each late-tail checkpoint / stationary cycle, evaluate the existing v10.2.30 post-first-passage energy gate over the full relevant stochastic event space.

The threshold-scaled event proposal uses the same Xi as the waiting threshold and clips its event-length scale to the qualified bounded range. Therefore map the gate across:

- the complete physical cycle phase;
- Xi below the lower clip bound;
- Xi through the interior threshold-scaled range;
- Xi above the upper clip bound.

For each `(phase, Xi)` evaluate and retain:

```text
proposed_event_length
fixed_opening_released_energy
active resistance
admitted_event_length
admission_reason / rejection_reason
mesh-resolved admissibility
```

Do not replace this with a scalar K threshold.

Produce an explicit late-tail gate map:

```text
L_admit(phase, Xi; N)
```

and, where the state is stationary,

```text
L_admit_star(phase, Xi)
```

## Phase D — marked-renewal stable-birth solver

Because rejected attempts do not change crack geometry and must not change the deterministic pre-birth physical state, separate the expensive physical trajectory from the threshold-renewal process.

Implement a cheap marked-renewal solver driven by the deterministic trajectory/kernel:

- Xi_i are the exact iid Exp(1) canonical threshold increments;
- event time is obtained by crossing cumulative cleavage action;
- the same Xi_i determines the threshold-scaled event proposal;
- the exact event phase/state is recovered from the deterministic trajectory/kernel;
- the v10.2.30 energy gate accepts or rejects the attempt;
- rejected attempts consume Xi_i and continue;
- the first admitted attempt is absorbing stable crack birth.

Use a deterministic quadrature / renewal-operator method where possible for accurate very-small failure probabilities. A large cheap Monte Carlo ensemble may be used as an independent validation, but Monte Carlo sampling error must not define the final high-survival tail if an operator/quadrature solution is available.

Cross-check the renewal solver against explicit v2 threshold-by-threshold simulation on shortened real trajectories containing both rejected and admitted attempts.

## Stationary-tail asymptotic classification

Once a stationary periodic tail is rigorously established, exploit it for a model-derived classification.

### Endurance-supported route

If after some finite cycle `N_star` the gate is proven closed for **all** physically possible event phases and all Xi/event proposals, then no future cleavage attempt can create a stable crack while the system remains in that stationary state.

Then:

```text
S_stable(infinity) = S_stable(N_star) > 0
```

and strict stable-crack-birth endurance is supported.

The closure proof must include numerical margins and show the gate is not merely missed because of phase/Xi discretization.

### No-endurance route

If the stationary tail has:

1. infinitely continuing cleavage attempts (already true if the positive attempt floor remains active), and
2. a positive-probability / positive-measure set of `(phase, Xi)` values for which the energy gate admits a stable event,

then eventual stable birth has probability 1 under the stationary marked-renewal process.

Report:

```text
stable_crack_birth_asymptotic_classification = no_endurance
```

with the accepted-event probability/rate or renewal-operator spectral argument as the basis.

Do not infer this merely from one admitted event early in the transient.

### Undetermined route

If state stationarity is not established, the gate margin approaches zero, or admissibility is unresolved over a non-negligible region, retain `undetermined`.

## Finite VHCF horizons remain required

Regardless of strict asymptotic classification, report practical stable-birth survival at:

```text
1e8
1e10
1e12
1e14 cycles
```

where the validated acceleration/kernel covers those horizons.

Keep strict asymptotic endurance separate from finite-horizon survival.

Report at least:

```text
P_stable_birth_by_N
S_stable_birth(N)
expected nonpropagating attempt count
probability of >=1 cleavage attempt
conditional probability that an attempt is admitted
```

At very low attempt rates, retain log probabilities to avoid underflow.

## Important floor interpretation

Do not let the mathematically positive zero-drive attempt floor dominate the engineering interpretation of 1e14 cycles without quantifying it.

For each class report the contribution that the intrinsic zero-drive attempt floor alone would make to cumulative attempt hazard by 1e12 and 1e14. Keep this separate from the actual cyclic tail.

This distinguishes:

- strict N->infinity behavior;
- physically relevant VHCF behavior through 1e14.

## Class sequence

Use the current low-stress 300 K conditions first:

```text
Peak:     735.9214951373579 MPa
DBTT:     700 MPa
weak-T:   350 MPa
ceramic:  400 MPa
```

These are tail probes, not fitted endurance limits.

Proceed one class at a time:

```text
Peak -> DBTT -> weak-T -> ceramic
```

If one condition is so far from the stable-birth transition that it provides no useful admission information, choose at most one additional lower/higher stress for that class based on the completed kernel. Do not densify the S-N knee.

## Required outputs

Maintain atomically:

```text
runs/sn_v9_canonical_four_class/300K_quiet_tail_state.csv
runs/sn_v9_canonical_four_class/300K_energy_gate_envelope.csv
runs/sn_v9_canonical_four_class/300K_stable_birth_survival.csv
runs/sn_v9_canonical_four_class/300K_stable_birth_endurance_summary.csv
runs/sn_v9_canonical_four_class/endurance_diagnostics.json
runs/sn_v9_canonical_four_class/condition_registry.json
runs/sn_v9_canonical_four_class/campaign_manifest.json
```

Every class summary must contain both:

```text
cleavage_attempt_asymptotic_classification
stable_crack_birth_asymptotic_classification
```

and clearly distinguish finite-horizon survival from strict asymptotic endurance.

## Autonomous mode

Continue beyond short unit tests and isolated audits.

Commit/push useful validated milestones but keep proceeding automatically through:

1. deterministic quiet-tail kernel;
2. acceleration qualification;
3. energy-admission envelope;
4. marked-renewal survival solver;
5. Peak 1e14-capable result;
6. DBTT, weak-T, ceramic sequentially.

Stop only for a genuine constitutive/provenance ambiguity, failure of bounded acceleration/restart equivalence, or an energy-gate mapping that cannot be resolved without changing physics.
