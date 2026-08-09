# Codex local workspace start: recover current Stateful-PD baseline and build large-N S-N capability

## Purpose

Use the local canonical reconstruction as a **read-only source artifact** and the GitHub clone as the **only writable development tree**. Recover a runnable, provenance-controlled baseline first. Then continue from the last validated scientific checkpoint toward S-N calculations spanning approximately `1e3` to `1e12` cycles and capable of distinguishing an endurance regime from merely very long finite life.

Do not start another parameter sweep or multi-hour/multi-day campaign until the gates below pass.

## Local paths

Read-only reconstruction supplied by the user:

```text
/Volumes/Data/Data/Nanopillar_calculation/Stateful_kitagawa_dadN_canonical_reconstruction
```

GitHub repository:

```text
https://github.com/ukaiiaku-maker/Stateful_kitagawa_dadN.git
```

Preferred local Git working tree:

```text
/Volumes/Data/Data/Nanopillar_calculation/Stateful_kitagawa_dadN
```

The user may open the reconstruction folder in VS Code. That does **not** make it the development tree. Treat it as immutable evidence. All source edits, tests, commits, and generated development artifacts must go into the Git working tree.

## Step 0 — establish the two-root workspace

Before touching source, print and record:

```bash
REFERENCE_ROOT=/Volumes/Data/Data/Nanopillar_calculation/Stateful_kitagawa_dadN_canonical_reconstruction
WORK_ROOT=/Volumes/Data/Data/Nanopillar_calculation/Stateful_kitagawa_dadN

printf 'REFERENCE_ROOT=%s\n' "$REFERENCE_ROOT"
printf 'WORK_ROOT=%s\n' "$WORK_ROOT"
find "$REFERENCE_ROOT" -maxdepth 2 -type f | sort
```

If `WORK_ROOT` does not exist, clone it:

```bash
cd /Volumes/Data/Data/Nanopillar_calculation
git clone --branch codex/sn-endurance-limit \
  https://github.com/ukaiiaku-maker/Stateful_kitagawa_dadN.git \
  Stateful_kitagawa_dadN
```

If it already exists, do not delete or overwrite it. Audit it first.

Then:

```bash
cd "$WORK_ROOT"
git fetch --all --prune
git switch codex/sn-endurance-limit
git pull --ff-only

git branch --show-current
git rev-parse HEAD
git status --short --branch
```

Create a child branch for implementation. Do not develop directly on the handoff branch:

```bash
git switch -c codex/sn-endurance-limit-v9
```

If that branch already exists, inspect it before deciding whether to reuse it.

## Step 1 — read the handoff before modifying physics

Read, in this order:

```text
docs/CANONICAL_SOURCE_PROVENANCE.md
docs/CURRENT_STATE.md
docs/KNOWN_ISSUES.md
docs/PHYSICS_TARGET.md
docs/DEFINITION_OF_DONE.md
docs/CODEX_HANDOFF.md
GitHub Issue #1
```

Also inspect every file in `REFERENCE_ROOT` and make an explicit inventory of what is present and what is missing.

## Canonical baseline identity

The latest reconstructed/validated solver lineage is v8.7 plus the v2.8.1 gate/orchestration correction. These three v8.7 hashes are immutable baseline identity during dependency recovery:

```text
09e993183efff4857819b4ede0d668c909e20810936c775638edbc9b2c5809d6
  arrhenius_fracture/sn_pd2d_stateful_v8_7_generalized_features.py

bb8341a3d8ed883f605bd99f6d74b4b011837d436d913cc1914a882a133f05b3
  arrhenius_fracture/stateful_peridynamics_v8_7_local_front_spacing.py

8ae3e8c454297fcbc61f179b747ffbd05129aa3c8598528113aa2edf8913396e
  arrhenius_fracture/sn_feature_geometry_v8_7.py
```

Do not alter those files while reconstructing imports/dependencies. If a dependency conflict appears, resolve the dependency closure around v8.7 rather than silently editing v8.7 to make imports pass.

## Step 2 — recover a runnable baseline into Git

The reconstruction is incomplete as a standalone Python package. Recover only the missing dependency closure needed by the v8.7 solver and its tests.

Preferred evidence order:

1. files already in `REFERENCE_ROOT`;
2. exact historical artifacts/manifests referenced by `docs/CANONICAL_SOURCE_PROVENANCE.md`;
3. historical branches/files from `ukaiiaku-maker/PF-fracture-fatigue` or other clearly identified source repositories;
4. only if necessary, a compatible dependency selected by explicit audit and documented as a reconstruction choice.

Do not mix dependency versions opportunistically.

Trace imports recursively from at least:

```text
arrhenius_fracture/sn_pd2d_stateful_v8_7_generalized_features.py
arrhenius_fracture/stateful_peridynamics_v8_7_local_front_spacing.py
arrhenius_fracture/sn_feature_geometry_v8_7.py
```

Build a dependency map. For every recovered file record:

- source repository/branch/artifact;
- original path;
- SHA-256;
- why that version was chosen;
- whether it is exact or a compatibility reconstruction.

Write this to:

```text
docs/CODEX_INITIAL_AUDIT.md
reference/BASELINE_DEPENDENCY_MANIFEST.json
```

Copy source from `REFERENCE_ROOT` into `WORK_ROOT` only after comparing hashes. Never modify `REFERENCE_ROOT`.

## Step 3 — baseline import/test gate

Before changing physics:

1. compile the package;
2. import the v8.7 driver and PD classes;
3. run all available unit/regression tests that do not launch a long campaign;
4. add focused tests for the recovered dependency interfaces if gaps exist;
5. run one very short deterministic smoke sufficient to instantiate the production geometry, construct the local PD patch, advance a tiny number of blocks, and write a real summary/checkpoint.

The smoke is an infrastructure test, not a physics calibration.

Record exact Python executable, environment, package versions, branch, HEAD, source hashes, commands, and results in `docs/CODEX_INITIAL_AUDIT.md`.

Do not proceed if the canonical v8.7 hashes changed or if the smoke creates invalid geometry/morphology.

## Where the previous work stopped

The last meaningful v8.7 gated result established:

- a K360 shielded failure endpoint with valid physical handoff at about `3.696e6` cycles;
- a lower K360 endpoint that was right-censored at `1e8` cycles;
- local front spacing about `7.5 um` with `delta/h = 4`;
- clean connected/slender handoff morphology;
- correct reuse after the v2.8.1 summary-discovery fix;
- strong low-stress hardening/shielding during the censor trajectory, with the instantaneous initiation hazard decreasing dramatically.

The unresolved scientific issue is that this strong hazard decay may represent either:

1. a true endurance regime with finite remaining integrated hazard, or
2. extremely long finite life with a hazard tail whose integral still diverges.

A `1e8` or `1e12` runout alone cannot distinguish these.

Previous optimization attempts also exposed orchestration failures (`--skip-existing`, summary discovery, false ranking of failed jobs). Do not resume those optimizers until the core solver and large-N diagnostics are proven.

## Primary scientific target

Construct statistically meaningful S-N predictions spanning roughly:

```text
1e3 to 1e12 cycles
```

and classify each material/mechanism/stress regime as:

```text
endurance_supported
no_endurance
undetermined
```

based on the asymptotic initiation/failure hazard, not a finite runout convention.

For total hazard `h(N)` and cumulative hazard

```text
H(N) = integral_0^N h(n) dn
```

an endurance population is possible only if the remaining integrated hazard is finite:

```text
H(infinity) < infinity
```

so that survival has a nonzero asymptote:

```text
S(infinity) = exp(-H(infinity)) > 0.
```

If `H(infinity)` diverges, eventual failure probability tends to one even if `h(N)` becomes extremely small.

Never create an endurance limit by:

- an arbitrary stress cutoff;
- declaring runout to be infinite life;
- clipping a physical hazard to zero;
- inserting a fitted horizontal S-N plateau;
- choosing a numerical underflow threshold that makes the rate vanish.

## Phase A — implement hazard-tail diagnostics before changing physics

Create a standalone versioned module, e.g.:

```text
arrhenius_fracture/endurance.py
```

It must accept time/cycle and hazard/cumulative-hazard histories and return:

- `endurance_supported`;
- `no_endurance`;
- `undetermined`.

It must report:

- fitted/estimated tail form;
- fit window(s);
- uncertainty/sensitivity to window selection;
- remaining integrated-hazard estimate or bound;
- residual survival estimate when finite;
- reason for classification.

Synthetic tests are mandatory:

```text
constant positive hazard         -> no_endurance
exponential decay                -> endurance_supported
power law p = 0.5                -> no_endurance
power law p = 1.0                -> no_endurance
power law p = 1.5                -> endurance_supported
short/noisy ambiguous tail       -> undetermined
```

Do this before modifying v8.7 physics.

## Phase B — large-N numerical architecture

Create a new solver version. Do not rewrite v8.7 in place.

The new implementation must be capable of reaching `1e12` cycles without stepping cycle-by-cycle. Use adaptive/event-driven/log-cycle blocks while preserving exact first-passage semantics.

Required invariants:

1. Persistent exponential thresholds are sampled once and never silently resampled.
2. Integrated hazards advance over the exact physical cycle interval.
3. When a threshold lies inside a proposed block, locate the crossing conservatively/exactly and stop at it.
4. Physical state is committed once per accepted interval. No doubled constitutive time or repeated state commit.
5. Results are invariant, within explicit tolerance, to reasonable block partition changes.
6. Checkpoints are atomic and contain every stochastic threshold, hazard accumulator, physical state variable, geometry state, RNG state/digest, cycle coordinate, and request identity required for exact continuation.
7. Restart must reproduce uninterrupted execution for the same case.
8. Floors/caps used only for numerical safety must be audited so they cannot manufacture an endurance limit.
9. Any asymptotic extrapolation must fail to `undetermined` when the tail cannot be bounded convincingly.

## Phase C — compare two physical initiation hypotheses

Preserve the current independent-cleavage baseline as a named mechanism class.

Add a separate, explicitly named activity-conditioned initiation hypothesis in which cleavage/initiation depends on resolved cyclic irreversible activity/delivery. Appropriate activity variables may include physically existing quantities such as:

- irreversible plastic event delivery;
- mobile/escaping defect activity;
- plastic dissipation/work increment;
- another already resolved state with a documented physical interpretation.

Do not simply multiply by an arbitrary stress gate.

The purpose is to test whether true shakedown can cause the future initiation hazard to become integrable, while the independent thermally activated cleavage model may retain a non-integrable tail.

Also audit the existing low-stress trajectory for one-sided hardening/shielding. Determine whether recovery, softening, or slow degradation physics is absent. Do not add a degradation variable merely to force eventual failure; any new slow process needs a physical state, rate law, units, diagnostics, and an ablation.

## Phase D — S-N statistical pilot

Only after unit, tail-analysis, block-invariance, and restart gates pass, run a small pilot:

```text
3 stresses x 3 seeds x 2 mechanism classes
```

Use common-random-number comparisons between mechanism classes where possible.

For each realization record at least:

- stress amplitude/range and R ratio;
- seed/random-field identity;
- first nucleation/embryo event cycle;
- stabilization cycle;
- root connection/front capture cycle;
- physical handoff cycle;
- right-censor cycle if no handoff;
- instantaneous and cumulative initiation hazard;
- relevant shielding/hardening/activity states;
- geometry/morphology audit fields;
- source and request hashes.

Do not interpret invalid geometry, max-block termination, orchestration failure, or missing summary as a censor.

Before expanding the pilot, require:

- valid failure/censor classifications;
- no monotonicity inversion beyond statistical uncertainty;
- qualitatively sensible stress dependence;
- tail classifications that are stable enough to be meaningful;
- clear difference, if present, between endurance and no-endurance mechanism classes;
- acceptable compute-time projection for the full campaign.

## Phase E — production S-N outputs

The final workflow should generate:

```text
sn_results.csv
survival_results.csv
endurance_diagnostics.json
campaign_manifest.json
```

and figures for:

1. stress versus `log10(N)` with failures and right-censored points;
2. survival/failure probability versus N at selected stresses;
3. instantaneous hazard and cumulative hazard versus N;
4. representative state evolution showing why the hazard decays, saturates, or re-accelerates;
5. comparison of endurance-supported versus no-endurance mechanism classes.

Do not hide stochastic scatter behind a single deterministic S-N line.

## Initiation versus propagation

Keep initiation life to **PD physical handoff** separate from subsequent sharp-crack `da/dN` propagation life.

The immediate S-N target is the initiation/handoff distribution and its endurance behavior. Do not let long-crack propagation consume the first development cycle.

Once initiation is robust, propagation can be added as a second conditional life component:

```text
N_total = N_initiation_to_handoff + N_propagation
```

with each component separately auditable.

## Git discipline

- Work only in `WORK_ROOT`.
- Keep `REFERENCE_ROOT` read-only.
- Never force-push.
- Do not rewrite the canonical baseline commits.
- Make small commits aligned with validation gates.
- Recommended commit sequence:

```text
1. recover exact runnable v8.7 dependency closure
2. add endurance-tail analysis and synthetic tests
3. add large-N adaptive/event-driven stepping tests
4. add atomic restart/checkpoint equivalence
5. add activity-conditioned initiation mechanism
6. add 3x3x2 S-N pilot orchestration and analysis
```

After each milestone, update `docs/CODEX_INITIAL_AUDIT.md` or a versioned progress document with:

- branch/HEAD;
- tests run;
- physical/numerical changes;
- known remaining issues;
- next proposed gate;
- whether any long run is being requested.

## First Codex stopping point

For the first pass, stop after all of the following are complete:

1. two-root workspace audited;
2. Git child branch created;
3. reconstruction inventory complete;
4. exact runnable dependency closure recovered;
5. canonical v8.7 hashes verified unchanged;
6. package compiles/imports;
7. short deterministic smoke passes;
8. `docs/CODEX_INITIAL_AUDIT.md` and `reference/BASELINE_DEPENDENCY_MANIFEST.json` written;
9. endurance-tail module and its synthetic mathematical tests are implemented and passing.

Do **not** launch a production S-N sweep at this point. Report the audit, commits, test results, any dependency ambiguities, and the proposed design for the large-N solver before proceeding.
