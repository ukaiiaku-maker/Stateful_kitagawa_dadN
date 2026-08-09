# Codex plan: stable-crack-birth S-N endpoint and canonical four-class parameterizations

## Precedence and correction

This document supersedes any earlier guidance that treated **physical handoff, front capture, or subsequent stable-crack growth** as the primary S-N endpoint.

The scientific target for the present S-N work is the **birth of a stable crack**.

Do **not** develop, retune, or replace post-stable growth/linkage/front-progression physics merely to extend the S-N life range. Those processes are outside the current S-N endpoint definition.

The high-stress ~million-cycle physical-handoff plateau found in the first v9 skeleton is therefore not a blocker for the current S-N objective. It arose after stable-crack birth and should be retained only as a diagnostic of the legacy handoff model, not used to drive new progression physics.

## Primary S-N endpoint

The failure/event coordinate for the present S-N program is:

```text
N_stable_crack_birth
```

Operationally, in the existing persistent-site model this is the first persistent embryo/site that undergoes the **stabilization branch** of the stabilization/healing competing process and enters the stable state.

The endpoint is **not**:

```text
first embryo birth
front capture
minimum crack extension
physical handoff
subsequent da/dN propagation
```

A transient embryo that heals is not a stable crack and must not terminate the run.

Once the first stable crack is born, the S-N condition may terminate immediately after writing a complete accepted-boundary checkpoint/summary/manifest. There is no need to spend additional compute on growth, linkage, front capture, or physical handoff for the present S-N curve.

## Required implementation change

Add an explicit versioned endpoint mode, for example:

```text
--fatigue-endpoint stable_crack_birth
```

Requirements:

1. Preserve the existing embryo birth and competing stabilization/healing physics unchanged.
2. Preserve persistent site thresholds, transition thresholds/outcome uniforms, cycle-phase state, RNG states, and all accepted physical state up to the stabilization event.
3. Localize the first stable transition inside the proposed block using the same transactional/event-localization contract as other v9 events.
4. Stop only after the stable transition is committed exactly once at an accepted physical boundary.
5. Write the final event identity, site identity, embryo birth cycle, stable cycle, state variables, hazards, source hashes and request hash.
6. A run reaching its observation horizon without a stable crack is right-censored for the stable-crack-birth endpoint.
7. Invalid geometry, numerical failure, missing state, max-subdivision failure, or request mismatch are never censors.

Add an equivalence regression proving that, for a short case where the legacy/full solver later proceeds to growth/handoff, `stable_crack_birth` mode stops at exactly the same stable-site event cycle and state that the full trajectory records before continuing.

## Reinterpret the already-completed v9 skeleton before running anything new

Do not rerun the completed 300 K v9 conditions merely because the endpoint has changed.

Extract from their existing histories/checkpoints the exact:

```text
N_first_embryo_birth
N_stable_crack_birth
N_front_capture          # diagnostic only
N_physical_handoff       # diagnostic only
```

Rebuild the current single-seed S-N table using `N_stable_crack_birth` as the primary life.

The previously completed high-stress cases already indicate that stable-crack-birth life spans far more decades than physical-handoff life: first-stable timing moved from approximately 1e5 cycles to 1e2 cycles across the upper v9 conditions while handoff remained near 1e6 cycles. Recover the exact values from the authoritative histories rather than copying approximate handoff text.

The 690.443 MPa / 300 K condition had no realized birth through 1e8 cycles and therefore remains a valid stable-crack-birth right-censored observation.

Create/update outputs such as:

```text
sn_stable_crack_birth_results.csv
stable_crack_birth_endurance_diagnostics.json
campaign_manifest.json
```

Keep any existing physical-handoff results in a separate diagnostic table; do not overwrite or discard them.

## Survival/endurance semantics for this endpoint

The endurance question must now be stated for **stable-crack birth**.

Keep distinct:

```text
embryo_birth_hazard / cumulative embryo-birth hazard
stabilization hazard
healing hazard
stable_crack_birth event
```

`exp(-H_embryo_birth)` is a no-embryo-birth survival quantity, not in general the complete stable-crack-birth survival probability after embryos exist.

For ensemble S-N/survival work, estimate stable-crack-birth survival from repeated realizations or from an explicitly derived multistage first-passage law. Do not relabel no-embryo survival as stable-crack survival.

For asymptotic logic, a finite total embryo-birth hazard is sufficient to establish a nonzero probability of no stable crack, but divergent embryo-birth hazard alone does not prove eventual stable-crack birth because embryos can heal.

Report separately, where applicable:

```text
embryo_birth_endurance_diagnostic
stable_crack_birth_endurance_classification
```

A finite observation horizon is never itself an endurance limit.

# Canonical four fracture parameterizations

The S-N work must use the same four audited Arrhenius fracture/Peierls–Taylor parameterizations used in the current canonical FEM/PF fracture-parity work.

The intent is **parameter transfer, not refitting**. All four classes use the same numerical/model architecture and differ through their audited material/barrier parameter rows.

The canonical set for the present work is exactly:

```text
Peak     = v913_paper_peak01_0242980_persistent_sites
DBTT     = v913_paper_dbtt01_0202500_persistent_sites
weak-T   = v913_paper_weakT01_0257068_persistent_sites
ceramic  = v913_paper_ceramic01_0189364_persistent_sites
```

Do not manually reconstruct these rows from remembered constants. Locate and load them through the exact audited registry/selection machinery used by the corresponding fracture/fatigue/FEM-CZM campaigns, and record source repository/branch, registry file, option identifier and hashes in the v9 campaign manifest.

## 1. Peak

```text
v913_paper_peak01_0242980_persistent_sites
```

Representative strongly nonmonotonic temperature-dependent response. The canonical qualified response has a low-temperature brittle shelf, strong transition-region toughening, an intermediate-temperature maximum near the ~900 K region in the theta=0/rate1x qualification, and a high-temperature decline. Cleavage and Peierls/Taylor activity are both important.

## 2. Classical DBTT

```text
v913_paper_dbtt01_0202500_persistent_sites
```

Classical ductile-to-brittle-transition response arising from cleavage/emission/PT/blunting/backstress competition. High-temperature states may have intense emission/PT activity and large exact-event workloads.

## 3. Weak-temperature / FCC-like control

```text
v913_paper_weakT01_0257068_persistent_sites
```

Control with deliberately weak temperature sensitivity relative to Peak and DBTT. The important science test is that the framework preserves this weaker temperature dependence rather than automatically producing a strong DBTT-like S-N shift.

The user specifically expects the weak-T S-N response may differ substantially from the DBTT response. Do not tune them toward the same shape.

## 4. Ceramic-like control

```text
v913_paper_ceramic01_0189364_persistent_sites
```

Comparatively brittle / weak-plasticity control used to preserve the low-plasticity regime and class ordering without class-specific numerical fracture rules.

## Common physics to preserve across all four classes

For the current canonical parity framework, preserve the shared architecture:

- theta = 0 degrees for the canonical comparison where that loading convention applies;
- rate1x loading matched by the same nominal K-dot in the fracture-parity calculations rather than by arbitrary displacement-rate matching;
- persistent source-site formulation;
- two reduced signed slip channels;
- separate mobile and retained populations;
- Peierls and Taylor Arrhenius kinetics;
- unsigned active mobile+retained content controlling Taylor backstress;
- signed retained content controlling direct crack-tip shielding;
- accumulated slip controlling tip blunting;
- stochastic cleavage first passage;
- stochastic emission;
- moving crack-tip MPZ/advection/escape where the fracture model uses it;
- no athermal fracture criterion;
- no class-specific numerical rules or independent retuning.

For the S-N stable-crack-birth solver, map only those pieces that are physically relevant before stable crack birth. Do not add post-birth crack-growth machinery merely because it exists in the fracture code.

## Parameter provenance requirement

For the current canonical four-class FEM/CZM parity work:

```text
Peak and DBTT: audited v10.2.25 paper parameter registry/selection
weak-T and ceramic: audited v10.2.26 weak-T/ceramic registry/selection
```

Audit the exact source manifests/registries rather than relying on prose alone.

There is a separate later v10.4.x bulk-plasticity registry with different weak-T/ceramic candidate rows:

```text
Peak     = v913_paper_peak01_0242980_persistent_sites
DBTT     = v913_paper_dbtt01_0202500_persistent_sites
weak-T   = v913_paper_weakT01_0129902_persistent_sites
ceramic  = v913_paper_ceramic01_0077080_persistent_sites
```

Those **are not** the canonical weak-T/ceramic rows for the present S-N/parity work. Use them only if explicitly reproducing the named v10.4.x bulk-plasticity campaign after auditing its exact manifest.

# Four-class S-N implementation sequence

Proceed autonomously, one working condition at a time.

## Phase A — endpoint conversion and reuse of existing v9 results

1. Implement and test `stable_crack_birth` endpoint mode.
2. Reanalyze every completed v9 300 K condition without rerunning it.
3. Produce the stable-crack-birth S-N skeleton and stage timing table.
4. Confirm that the high-stress physical-handoff plateau disappears from the primary S-N endpoint because post-stable progression is no longer counted.
5. Preserve the 690.443 MPa observation as a right-censored stable-crack-birth condition if no stable transition occurred.

Do not proceed until the extracted stable event cycles agree exactly with authoritative historical state/event records.

## Phase B — exact four-class parameter loader

1. Locate the audited registries/selections for the four option IDs above.
2. Implement a versioned loader in this repository that selects by exact option identifier, preferably by importing/copying the audited registry mechanism rather than transcribing constants.
3. Record the full resolved row and SHA/hash/provenance in every request/manifest.
4. Add tests that each identifier resolves to exactly one row and that no v10.4.x alternate weak-T/ceramic row can be selected accidentally under the canonical four-class campaign name.
5. Do not alter any row to make one class match another.

## Phase C — one-condition four-class qualification at the existing validated S-N loading state

Use the already-qualified v9 S-N geometry/loading conventions as the numerical/mechanical testbed initially. Start at the current 300 K loading state because the v9 engine is already qualified there.

For each class, one at a time:

```text
Peak
DBTT
weak-T
ceramic
```

run/resolve enough stress conditions to demonstrate a valid stable-crack-birth finite-life or right-censored result with complete restartable state and the correct parameter identity.

Do not launch the four classes in parallel. Finish, register and commit one before moving to the next.

Stress values may be selected adaptively per class. Do not assume that a stress appropriate for Peak is appropriate for weak-T or ceramic.

## Phase D — temperature dependence

The scientific objective is to determine how the S-N/stable-crack-birth response changes with temperature for the four audited classes, especially the expected contrast between DBTT and weak-T.

Do **not** invent a new temperature grid from this document. Recover the canonical temperature set from the audited four-class fracture/fatigue campaign manifests/runner definitions and use that same set, or document exactly why a reduced subset is required for the first qualification.

At each temperature/class combination:

1. preserve the exact audited option row;
2. preserve common numerical rules;
3. build a small sequential stable-crack-birth S-N skeleton using adaptive stress selection;
4. maintain finite-life, right-censored and endurance classifications separately;
5. do not force common slopes, endurance limits or transition temperatures across classes.

Key hypothesis to test rather than impose:

- Peak: strong/nonmonotonic temperature effect may appear in stable-crack-birth S-N behavior;
- DBTT: strong systematic temperature sensitivity expected from cleavage/plasticity competition;
- weak-T: substantially weaker temperature sensitivity should be preserved if the parameterization behaves as intended;
- ceramic: comparatively brittle/weak-plasticity response should remain distinct.

## Phase E — endurance/no-endurance science

Once the four canonical classes produce working stable-crack-birth S-N skeletons, evaluate whether any class/temperature/stress regime supports a true endurance population.

Use the hardened asymptotic contract already implemented:

```text
empirical tail evidence != physics asymptotic proof
right censor != endurance
```

Focus the large-N diagnostics on the stable-crack-birth process, including embryo births and stabilization/healing competition.

Do not add a persistent-direct-cleavage comparison mechanism until the four audited physical classes are working and their existing behavior is understood. The four canonical parameterizations are the first scientific comparison set.

# Required outputs

At minimum maintain:

```text
sn_stable_crack_birth_results.csv
stable_crack_birth_stage_results.csv
stable_crack_birth_survival_results.csv
stable_crack_birth_endurance_diagnostics.json
campaign_manifest.json
```

Every row must include:

- class label;
- exact option identifier;
- parameter-row provenance/hash;
- temperature;
- stress amplitude/range and R ratio;
- geometry/loading identity;
- seed and random-stream identity;
- first embryo cycle if realized;
- stable crack birth cycle if realized;
- observation horizon if censored;
- classification;
- relevant delivery/cleavage/stabilization/healing diagnostics;
- source/request hashes;
- restart/validity status.

Keep physical-handoff/front-growth outputs in separately named legacy diagnostic files only.

# Stop conditions

Continue autonomously through endpoint conversion, exact four-class registry loading, existing-result reanalysis and one-condition-at-a-time four-class qualification.

Stop for user input only if:

1. the exact audited parameter row/provenance cannot be recovered unambiguously;
2. the semantics of `stable` in the current state machine do not correspond to an irreversible/stable crack nucleus and require a new physical definition;
3. transferring the fracture parameter row into the S-N model requires a constitutive mapping not defined by the existing common physics;
4. a temperature/loading convention is ambiguous in the authoritative campaign manifests;
5. a numerical defect cannot be repaired without changing physical assumptions.

Do not stop merely because another audit or unit-test milestone is complete.
