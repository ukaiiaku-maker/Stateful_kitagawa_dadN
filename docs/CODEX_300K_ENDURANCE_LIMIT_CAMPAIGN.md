# Codex autonomous plan: 300 K low-stress VHCF/endurance campaign

## Scope correction

This document supersedes prior guidance to densify the steep S-N transition or to proceed into 900/1300 K fatigue calculations.

The present scientific objective is **only the 300 K asymptotic fatigue problem**:

> As nominal cyclic stress is reduced, does the cumulative canonical stable-crack-birth hazard converge to a finite value, producing nonzero survival at infinite life, or does it continue to diverge so that stable-crack birth remains inevitable at sufficiently large cycle count?

The response-class labels Peak, DBTT, weak-T/FCC-like, and ceramic-like refer to the parameterizations' **fracture-toughness temperature response**. They do not imply that fatigue S-N temperature dependence should follow the same qualitative labels. Do not use those labels as constraints on fatigue behavior.

Do not run additional temperature cases in this work package.

## Canonical model and rows

Retain the already-qualified full-FEM canonical stable-birth path and exact final 2-D paper/parity rows:

```text
Peak     = v913_paper_peak01_0242980_persistent_sites
DBTT     = v913_paper_dbtt01_0202500_persistent_sites
weak-T   = v913_paper_weakT01_0129902_persistent_sites
ceramic  = v913_paper_ceramic01_0077080_persistent_sites
```

The rejected historical controls `0257068` and `0189364` remain preserved but are excluded from canonical production.

The endpoint remains:

```text
N_stable_crack_birth = N_first_canonical_cleavage_first_passage
```

No PD/K=2 embryo machinery, stabilization/healing, crack growth, handoff, event length, or da/dN is active.

For each fixed stress/material condition, the deterministic pre-birth trajectory defines

```text
S_stable(N) = exp[-H_cleave(N)]
```

with the canonical cumulative cleavage hazard `H_cleave` evaluated in the qualified log-safe formulation.

## Do not spend compute resolving the steep S-N drop

A rapid, approximately logarithmic stress-life falloff is expected. The existing points are sufficient to establish that the curves have steep stress sensitivity.

Do **not** add many intermediate stresses simply to make the plotted knee smooth.

Instead, move to **lower stresses and much longer cycle horizons**.

The target horizons for this work package are:

```text
1e8
1e10
1e12
1e14 cycles
```

where computationally and numerically tractable.

A condition may continue beyond a median-life crossing only when useful for validating asymptotic hazard behavior; otherwise use a lower-stress condition for the endurance study.

## Primary scientific observables

For each material class and low-stress condition, record at logarithmic cycle checkpoints at least:

```text
N
H_cleave(N)
log_h_cleave(N)
S_stable(N) = exp[-H_cleave(N)]
DeltaH_per_decade
h_ratio_per_decade = h(10N)/h(N)
effective_tail_exponent p_eff
signed mobile state summaries
signed retained state summaries
wake state summaries
Taylor backstress
signed shielding
blunting
emission rate/activity
Peierls activity
Taylor activity
accepted/rejected block counts
largest accepted dN
wall time
all floor/cap/underflow-representation diagnostics
```

Use exact log-domain quantities when the linear rate is below floating-point range. Underflow is not physical extinction.

## Endurance criterion

Finite runout at `1e12` or `1e14` cycles is **not** by itself proof of endurance.

For the canonical stable-birth process:

```text
H_infinity = integral_0^infinity h_cleave(N) dN
```

True statistical endurance requires

```text
H_infinity < infinity
```

which gives

```text
S_infinity = exp(-H_infinity) > 0.
```

Classify separately:

```text
finite_horizon_observation
empirical_tail_evidence
physics_asymptotic_classification
```

Allowed physics classifications remain:

```text
endurance_supported
no_endurance
undetermined
```

Use `endurance_supported` only when the actual constitutive equations/state trajectory permit a defensible finite upper bound on the remaining hazard.

Use `no_endurance` only when divergence is established from the actual asymptotic law, for example a positive hazard floor or a tail no steeper than `C/N^p` with `p <= 1`.

Otherwise report `undetermined` even at `1e14` cycles.

## Tail diagnostics

At each decade, evaluate the incremental hazard

```text
DeltaH_k = H(10^(k+1)) - H(10^k)
```

and the local effective power-law slope of the instantaneous hazard:

```text
p_eff = - d log h / d log N.
```

Interpretation:

- persistent positive hazard floor -> no endurance;
- `p_eff < 1` asymptotically -> divergent hazard;
- `p_eff ~ 1` -> logarithmic divergence unless a later crossover is proven;
- stable `p_eff > 1` over a broad tail -> empirical integrable-tail evidence, not by itself proof;
- exponentially decaying hazard/state activity -> endurance-capable if a constitutive upper bound on the residual integral can be established.

Also test adversarial tail alternatives before claiming endurance:

- small positive hazard floor;
- late crossover from `p>1` to `p<=1`;
- insufficient dynamic range;
- numerical clipping/underflow masquerading as zero;
- a late state transition that reactivates emission or cleavage.

## Stress-selection strategy at 300 K

Proceed one material class at a time and one stress at a time.

Do not use a uniform stress grid.

Use the current 300 K results only to identify the low-stress side, then step downward until the trajectory is informative over VHCF horizons.

### Peak

Existing evidence includes a long-lived condition near `735.921 MPa` with `H_cleave ~ 0.0515` by `1e8` cycles and no N10 crossing. Use this as the first Peak endurance anchor and extend/reproduce it with the qualified canonical production registry if needed for exact current-state outputs. Then test a modestly lower stress selected to produce a clearly smaller hazard through `1e8` while still measurable in log space.

Do not add extra points between 735.9 and 770 MPa just to resolve the S-N drop.

### DBTT

The existing 770 MPa point is finite-life. Select a lower 300 K stress expected to remain unfailed or low-hazard through at least `1e8`, then continue toward `1e12-1e14` if the adaptive solver remains efficient.

### weak-T / 0129902

The existing 400 MPa condition has `N50 ~ 3.36e6`, so move substantially below 400 MPa rather than filling 400-450 MPa. Seek a stress giving low but measurable cumulative hazard at `1e8`, then extend it to the VHCF horizons.

### ceramic / 0077080

The existing 400 MPa condition has not reached N10 by `1e8` and is already a useful long-life anchor. Extend this state to later decades where possible and add one lower stress only if necessary to distinguish finite-tail convergence from a very slow divergent hazard.

## Computational policy

The purpose of the v9 large-N engine is to make `1e12-1e14` trajectories feasible without cycle-by-cycle stepping.

For quiet tails:

- grow candidate blocks geometrically;
- preserve the validated FEM/MPZ error controller;
- use periodic-state or closed-form acceleration only when its error bound is demonstrated;
- retain exact log hazard accumulation;
- checkpoint at decade milestones and useful controller-state changes;
- profile before accepting multi-day execution.

Performance optimization is allowed only when the accepted physical trajectory and cumulative hazard remain partition/restart invariant.

Do not loosen tolerances merely to reach `1e14` faster.

## Sequential campaign order

Suggested order, because useful low-stress anchors already exist:

```text
1. Peak 300 K low-stress tail
2. ceramic 300 K low-stress tail
3. weak-T 300 K low-stress tail
4. DBTT 300 K low-stress tail
```

This order is practical, not a scientific ranking. Continue one condition at a time.

For each class, prefer to finish one informative low-stress trajectory through as many decades as needed before adding another nearby stress.

## Required outputs

Maintain/update atomically:

```text
runs/sn_v9_canonical_four_class/campaign_manifest.json
runs/sn_v9_canonical_four_class/sn_stable_crack_birth_results.csv
runs/sn_v9_canonical_four_class/sn_quantile_results.csv
runs/sn_v9_canonical_four_class/survival_results.csv
runs/sn_v9_canonical_four_class/endurance_diagnostics.json
runs/sn_v9_canonical_four_class/condition_registry.json
runs/sn_v9_canonical_four_class/300K_endurance_tail_summary.csv
```

`300K_endurance_tail_summary.csv` should contain one row per class/stress/checkpoint decade with the tail diagnostics listed above.

Create plots only after the numerical data are trustworthy. Useful final plots include:

- 300 K stress versus `N10/N50/N90` where quantiles exist;
- 300 K `H_cleave(N)` on log-N axes for the low-stress anchors;
- instantaneous hazard versus N;
- `DeltaH_per_decade` versus cycle decade;
- `p_eff(N)` versus cycle count;
- inferred `S_stable(N)` through the longest trusted horizon.

Do not imply a horizontal fatigue limit from plotting convention alone.

## Autonomous work mode

Continue autonomously through the 300 K VHCF/endurance program.

Do not stop merely because one tail fit, one decade, one test, or one material class completes.

Commit/push validated milestones and continue to the next class unless a genuine blocker occurs.

Genuine stop conditions are:

1. a constitutive ambiguity not defined by the audited source;
2. a real partition/restart failure that cannot be fixed without changing physics;
3. numerical scaling that prevents trustworthy extension even after algorithmic optimization;
4. a late-N state transition or cap/floor activation whose physical interpretation requires a user decision;
5. evidence that the current asymptotic classification logic itself is invalid for the realized state evolution.

## Preferred stopping point

Continue until the 300 K data can answer, for each of Peak, DBTT, weak-T, and ceramic:

- whether the canonical stable-birth hazard is still materially accumulating at `1e12-1e14` cycles;
- whether the observed tail is empirically integrable, divergent, or ambiguous;
- whether a model-derived endurance/no-endurance classification can be justified;
- the corresponding nonzero-survival estimate or divergence evidence where justified.

Do not run 900/1300 K or any other temperature in this work package.
