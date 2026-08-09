# Codex handoff: complete the large-N S-N / endurance-limit model

## Repository and branch

Work in:

```text
https://github.com/ukaiiaku-maker/Stateful_kitagawa_dadN.git
```

Start from branch:

```text
codex/sn-endurance-limit
```

Create your own child branch before modifying source, e.g.:

```text
codex/sn-endurance-limit-v9
```

Do not work directly on `main`.

## Read first

Before changing code, read in this order:

1. `README.md`
2. `docs/CURRENT_STATE.md`
3. `docs/PHYSICS_TARGET.md`
4. `docs/KNOWN_ISSUES.md`
5. `docs/DEFINITION_OF_DONE.md`

Then inventory the imported current source before modifying it. Preserve the validated v8.7 behavior as a frozen baseline and use new versioned modules for physics changes.

## End goal

Produce an S-N framework that can span very large cycle counts and scientifically distinguish:

- a mechanism with a true endurance regime; and
- a mechanism with no true endurance limit but potentially extremely long lives.

The distinction must be made from the asymptotic hazard/state evolution, not from a finite runout convention.

## Core mathematical requirement

For the no-failure trajectory, track total cumulative failure/initiation hazard

```text
H_tot(N) = integral h_tot(N) dN
```

(or the exact persistent-site equivalent).

The classifier must distinguish:

```text
H_tot(infinity) finite    -> nonzero infinite-life survival is possible
H_tot(infinity) divergent -> eventual failure probability tends to one
```

Do not classify a run as endurance solely because it reaches `1e8`, `1e10`, or `1e12` cycles.

## Work plan

### 1. Audit before editing

Run the current tests and produce `docs/CODEX_INITIAL_AUDIT.md` containing:

- branch/HEAD;
- Python/environment information;
- exact test results;
- active source hashes;
- import/dependency map for the current S-N solver;
- a map of where persistent thresholds, cumulative hazards, state commits, checkpoint writes, and physical-handoff classification occur;
- any current numerical caps/floors that can change the large-N asymptote.

Do not launch a long physical run during this audit.

### 2. Implement a standalone hazard-tail classifier first

Create a small, testable module (suggested name `arrhenius_fracture/endurance.py`) that accepts cycle/hazard history and returns at least:

```text
classification: endurance_supported | no_endurance | undetermined
reason_code
tail_model
tail_parameters
estimated_or_bounded_remaining_hazard
S_infinity_estimate_or_bound
fit_window
fit_quality
```

Use conservative bounds. If the tail cannot be distinguished, return `undetermined`.

Synthetic tests are mandatory:

- constant hazard;
- exponential decay;
- power-law p=0.5, 1.0, and 1.5;
- noisy/short tail that must remain undetermined.

### 3. Make large-N integration event-driven and restart-safe

The solver must span at least `1e3` to `1e12` cycles without stepping every cycle.

Requirements:

- adaptive/logarithmic cycle blocks when far from any event;
- exact or conservatively interpolated first-passage crossing inside a block;
- block subdivision when hazard/state changes invalidate a long extrapolation;
- persistent candidate thresholds never resampled on resume;
- state committed exactly once per accepted physical interval;
- atomic checkpoint generation;
- restart equivalence tests;
- request-hash identity and cache reuse;
- no ambiguous `--skip-existing` behavior.

### 4. Preserve the current mechanism as a reference

Do not alter the meaning of the current `plastic_shielded_case64_M1` baseline in place.

Create a new explicit model/version for experimental physics.

### 5. Add an activity-conditioned initiation hypothesis

The current architecture contains independent cleavage nucleation. A small positive asymptotic hazard would mathematically imply no strict endurance even if a `1e8` run censors.

Add a separately named closure in which initiation requires physically measurable cyclic irreversible activity (for example plastic-event delivery, dissipated plastic work rate, or an equivalent already-resolved state). The gate must be dimensional, documented, output as a diagnostic, and derived from state—not an arbitrary stress threshold.

This model is a hypothesis to test. Do not force it to exhibit endurance. Let its computed hazard tail decide.

### 6. Check whether the current hardening is artificially one-sided

The K360 lower-stress reference hardens/shields strongly and the birth rate falls by orders of magnitude. Audit whether the long-N state has:

- a physical saturation;
- recovery;
- cyclic softening/damage;
- an irreversible degradation channel;
- or only monotonic hardening.

If a missing slow process is required for the no-endurance comparison, add it as a separately selectable, documented mechanism with a physically interpretable rate. Do not add an arbitrary per-cycle damage constant solely to make failure happen.

### 7. Separate initiation life from propagation life

The primary S-N result is life to PD physical handoff.

Do not silently append sharp-front life. If total life is later requested, retain:

```text
N_initiation
N_propagation
N_total
```

and transfer an explicit tip-state capsule into the propagation model.

### 8. Statistical S-N campaign

Implement a campaign driver that supports:

- stress grid;
- multiple seeds;
- common-random-number comparisons between mechanism variants;
- valid right censoring;
- failure probability/survival probability by N;
- median and quantile S-N summaries where identifiable;
- monotonicity diagnostics versus stress.

Do not fit invalid, geometry-limited, max-block, corrupt, or infrastructure-failed runs as censors or failures.

### 9. Outputs

At minimum generate:

```text
sn_results.csv
survival_results.csv
endurance_diagnostics.json
campaign_manifest.json
```

and plots for:

- S vs log10(N), marking failures and censors;
- survival/failure probability vs N;
- instantaneous and cumulative hazard tails;
- representative internal-state evolution.

### 10. Validation ladder: do not skip

Follow `docs/DEFINITION_OF_DONE.md` literally.

Do not launch a many-hour or multi-day sweep until:

1. unit/synthetic tail tests pass;
2. current K360 fixtures validate;
3. block-partition invariance passes;
4. real checkpoint/restart equivalence passes;
5. a `3 stress x 3 seed x 2 model` pilot passes.

If any infrastructure job fails, fail the campaign closed. Never rank a missing summary or solver error as physical candidate performance.

## Constraints

- Preserve expensive completed data.
- Do not modify historical reference data after it is imported.
- Do not change fundamental case-64 activation barriers merely to fit the S-N shape unless a separate calibration task is explicitly justified.
- Do not add an athermal fracture cutoff just to create a fatigue limit.
- Do not classify endurance from finite runout alone.
- Do not hide numerical floors/caps; report any that influence long-N behavior.
- Commit in small, reviewable steps with tests after each step.

## First deliverable from Codex

Do **not** start by running a production campaign.

First commit should contain only:

1. `docs/CODEX_INITIAL_AUDIT.md`;
2. the standalone endurance-tail module plus synthetic tests;
3. no changes to the frozen baseline physics.

Then report the test output and explain whether the existing K360 censor trajectory, if its observed tail were continued unchanged, is presently more consistent with finite or divergent cumulative hazard—or remains mathematically undetermined.
