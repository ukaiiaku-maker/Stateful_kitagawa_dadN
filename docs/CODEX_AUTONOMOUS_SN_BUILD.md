# Codex autonomous S-N build plan

## Purpose

The next development period should no longer stop after every short audit. Work independently toward a functioning S-N calculation workflow, **one physical condition at a time**, while preserving the fail-closed numerical standards already established.

The target is a working program that can cover approximately `1e3–1e12` cycles, produce stress-life and survival outputs, and distinguish an endurance-capable response from a very long but ultimately finite-life response.

This is an implementation-and-calculation campaign, not another sequence of stand-alone audits.

## Working branch and current state

Work in the existing local Git tree:

```text
/Volumes/Data/Data/Nanopillar_calculation/Stateful_kitagawa_dadN
```

Continue on:

```text
codex/sn-endurance-limit-v9
```

Before starting new work, preserve and push the completed physical-adapter milestone. Do not rewrite its history.

Then fetch and merge the current handoff branch:

```bash
git fetch origin codex/sn-endurance-limit
git switch codex/sn-endurance-limit-v9
git merge --no-edit origin/codex/sn-endurance-limit
```

Read this document together with:

```text
docs/V9_NUMERICAL_ARCHITECTURE_AUDIT.md
docs/V9_PHYSICAL_ADAPTER_AUDIT.md
docs/ENDURANCE_ASYMPTOTIC_AUDIT.md
docs/PHYSICS_TARGET.md
docs/DEFINITION_OF_DONE.md
```

The frozen v8.7 source remains immutable.

## Change in working mode

Do **not** stop merely because one subtest or implementation stage has completed.

Proceed autonomously through the sequence below until one of these genuine stop conditions occurs:

1. a physical-model decision is required that cannot be resolved from the existing model/handoff;
2. a fail-closed identity/restart/stochastic-state requirement cannot be satisfied;
3. the real v9 accepted trajectory cannot be made block-partition invariant without changing physics;
4. a run exposes a numerical or physical inconsistency whose repair would alter the intended model rather than its implementation;
5. data corruption or repository ambiguity prevents trustworthy continuation.

Otherwise, fix implementation defects, add tests, commit working milestones, and continue to the next step.

Do not ask for approval after every audit or short smoke.

## Immediate technical objective: replace the partition-sensitive v8.7 block trajectory

The previous real-physics comparison found that native v8.7 1-cycle versus 0.05-cycle blocks produced large differences in delivery memory, completion and cumulative birth hazard. Therefore native v8.7 block integration cannot define the accepted large-N v9 trajectory.

Build the actual v9 physical `TransactionModel` and iterate until the **accepted** trajectory is converged under block partitioning.

Required implementation features:

- log-domain plastic-delivery and cleavage rates;
- exact/controlled finite-memory `Lambda` evolution and completion `Q(K,Lambda)`;
- stable integrated birth hazards for ultra-small rates;
- FEM/plastic/rho/back-stress updates with embedded error estimates;
- reject/subdivide rather than accepted-state clipping at law/cap boundaries;
- persistent birth, embryo, stabilization/healing, growth and progression clocks;
- earliest-event localization inside a proposed block;
- one physical/stochastic commit per accepted interval;
- atomic restartable array generations;
- all floor/cap activations recorded.

Use synthetic tests as needed, but do not stop when they pass. Continue directly into the first real condition.

## First real working condition

Use the previously validated **K360 shielded failure-like condition** as the first end-to-end target.

Target stress amplitude:

```text
sigma_a ≈ 735.921 MPa
```

Use the same temperature, R ratio, geometry, material/barrier case, machine/load convention, and stochastic/request settings as the preserved v2.8 K360 failure gate. Extract these from the committed v2.8 request/orchestrator/reference material. **Do not invent missing parameters.** If one essential parameter cannot be recovered unambiguously, that is a genuine stop condition.

This case is preferred because the historical v8.7 calculation reached a clean physical handoff near `3.696e6` cycles, so it provides a finite-life target without requiring a 1e8-cycle run before the real v9 solver is proven.

### Acceptance for condition 1

The new v9 physical solver must produce a valid result with:

- no invalid geometry or morphology classification;
- persistent thresholds sampled once and preserved;
- valid connected/slender physical handoff;
- complete checkpoint/summary/manifest output;
- restart equivalence;
- block-partition invariance using at least two substantially different initial candidate-block schedules;
- no hidden numerical rate floor defining the event;
- complete hazard/state diagnostics.

Historical life `~3.696e6` cycles is a reference, not an exact-fitting requirement. Diagnose any material discrepancy rather than tuning parameters simply to reproduce that number.

Once condition 1 works, commit and push it, then continue automatically.

## Second real working condition: long-life/endurance candidate

Next run the preserved K360 shielded lower-stress condition:

```text
sigma_a ≈ 690.443 MPa
```

Again use the exact v2.8 request definition where recoverable.

Run this with the working v9 adaptive engine. The objective is not merely to reproduce `right_censored at 1e8`; it is to determine the evolving large-N hazard/state behavior efficiently.

The run should be checkpointed and resumable. Extend sequentially through useful logarithmic milestones such as

```text
1e5, 1e6, 1e7, 1e8, 1e9, ...
```

but do not force execution to `1e12` if the asymptotic evidence is already decisive or if another regime/cap activation appears.

At each decade or other adaptive diagnostic checkpoint, record:

- delivery log-rate and activity;
- `Lambda` and `Q(K,Lambda)`;
- cleavage log-rate;
- instantaneous and cumulative birth hazard;
- stabilization/healing/growth/linkage activity;
- rho/back-stress/shielding state;
- tail exponent/bound diagnostics;
- floor/cap activation counters;
- current endurance evidence/classification.

A finite runout alone is not endurance.

Once condition 2 is trustworthy, commit and push it, then continue automatically.

## Build an S-N skeleton sequentially, not as a broad sweep

After the two anchor conditions work, construct a **single-seed S-N skeleton one condition at a time**.

Do not launch a grid of stresses in parallel.

Maintain a machine-readable registry of completed conditions, for example:

```text
runs/sn_v9_single_seed_skeleton/campaign_manifest.json
runs/sn_v9_single_seed_skeleton/sn_results.csv
runs/sn_v9_single_seed_skeleton/endurance_diagnostics.json
```

Every condition must have an immutable request hash and be reused/resumed rather than rerun when already valid.

### Adaptive stress selection

Choose the next stress based on the completed results, with the goal of covering log-life decades rather than uniform stress spacing.

Aim initially for valid finite-life points roughly spanning:

```text
1e4, 1e5, 1e6, 1e7, 1e8+ cycles
```

These are coverage targets, not forced exact lives.

Use interpolation/bracketing in stress versus `log10(N_handoff)` to select the next condition. If a candidate becomes long-lived/endurance-like, retain it as a censored/asymptotic point and choose the next stress accordingly.

Prefer the next condition that fills the largest useful gap in `log10(N)`.

Do not add many nearby stresses that provide redundant information.

### Classification rules

For every completed condition classify separately:

```text
physical_handoff
right_censored / finite-horizon observation
endurance_supported
no_endurance
undetermined
invalid
```

`right_censored` and `endurance_supported` are not synonyms.

Invalid geometry, numerical failure, missing summary, maximum-subdivision failure, request mismatch or incomplete checkpoint are never censors.

## Efficiency requirement

The v9 solver must make the large-N program practical.

Track per-condition:

- accepted block count;
- rejected/subdivided block count;
- wall time;
- largest accepted `dN`;
- cycles advanced per wall-second by decade;
- fraction of time in FEM, PD, hazard integration and I/O.

If a quiet long-life condition is progressing too slowly, profile and optimize the **validated algorithm** before simply allowing a multi-day run.

Optimization is allowed when it preserves the physical trajectory and the block/restart invariance tests. Examples include cached deterministic structures, vectorization, log-domain closed-form memory updates, and less frequent nonessential plotting.

Do not obtain speed by loosening error tolerances until the endurance-relevant state becomes block dependent.

## Checkpoint/resume discipline

Longer conditions may run for hours. They must be restartable.

- Publish atomic checkpoint generations periodically.
- Resume from the newest fully validated generation.
- Never restart a completed valid condition from zero merely because orchestration changed.
- The runner must distinguish `complete`, `restartable`, `invalid`, and `not_started`.
- A process interruption is not a scientific result.

Create a simple monitor/status command that reports the active condition, cycle coordinate, wall time, accepted/rejected blocks, last checkpoint and current hazard/endurance diagnostics.

## Autonomous debugging policy

During this work period Codex is expected to diagnose and repair implementation defects without stopping after each one.

For an implementation failure:

1. preserve the failing minimal reproducer;
2. identify the violated invariant;
3. fix it in a versioned source file;
4. add a regression test;
5. rerun the relevant short proof;
6. resume/restart the same physical condition if its checkpoint remains valid;
7. continue toward the working S-N condition.

Do not respond to every bug by starting a new campaign directory or rerunning all previous conditions.

## When to add stochastic replication

Do not begin with 3 seeds at every stress.

First establish the single-seed skeleton across several life decades.

Once at least four scientifically valid stress conditions exist and the solver is stable, add replication strategically:

1. the apparent transition/endurance region;
2. one intermediate finite-life condition;
3. one higher-stress finite-life condition.

Start with 3 seeds at those selected conditions. Use common random numbers for mechanism comparisons where appropriate.

Only after scatter is quantified should the workflow expand to a larger statistical S-N campaign.

## Mechanism comparison comes after the baseline skeleton works

Keep the first S-N skeleton on:

```text
completion_gated_independent_cleavage
```

Do not interleave development of `persistent_direct_cleavage` while the baseline solver itself is still being made operational.

Once the baseline skeleton works, add the comparison mechanism as a separately selectable model and run it through the same one-condition-at-a-time workflow.

No mechanism is to be labeled `no_endurance` by construction; its asymptotic hazard must establish that result.

## Required continuously updated outputs

As soon as two or more valid conditions exist, maintain:

```text
sn_results.csv
survival_results.csv
endurance_diagnostics.json
campaign_manifest.json
```

Also maintain figures that can be regenerated cheaply from completed data:

- stress amplitude/range versus `log10(N_handoff)` with censor symbols;
- instantaneous birth hazard versus cycle number;
- cumulative hazard versus cycle number;
- survival probability versus cycle number for long-life conditions;
- representative state evolution (`rho`, back stress, delivery, `Lambda`, `Q`);
- wall-time/performance diagnostics.

Do not wait for a full campaign before checking whether the emerging S-N shape is physically sensible.

## Git discipline during autonomous work

Make and push small meaningful commits as working capability is achieved, for example:

```text
1. implement converged v9 physical TransactionModel
2. validate first K360 failure condition
3. validate first K360 long-life condition
4. add sequential S-N condition selector/registry
5. add next finite-life S-N condition(s)
6. add selected-seed replication
```

Do not create a new branch for every bug.

Update a running document, e.g.:

```text
docs/V9_SN_BUILD_PROGRESS.md
```

with completed conditions, current blockers, performance, and next automatically selected condition.

## Final stopping point for this autonomous work package

Continue independently until **one** of the following is achieved:

### Preferred successful stopping point

A functioning baseline S-N building program exists with:

- converged/restartable v9 physical integration;
- at least four scientifically valid single-seed stress conditions;
- at least three decades of finite-life coverage if physically accessible;
- at least one long-life/endurance diagnostic condition;
- automatically updated S-N/endurance outputs;
- a demonstrated rule for choosing the next stress one condition at a time;
- no unresolved block-partition or restart defect.

At that point report the emerging S-N curve, endurance evidence, performance, remaining gaps, and recommended next conditions/seeds.

### Earlier genuine stop

Stop earlier only for one of the genuine stop conditions listed at the top of this document. Report the exact blocker, the last valid condition/checkpoint, and the narrowest scientific decision needed from the user.

Do not stop merely because a short audit, unit-test set, or intermediate implementation milestone passed.