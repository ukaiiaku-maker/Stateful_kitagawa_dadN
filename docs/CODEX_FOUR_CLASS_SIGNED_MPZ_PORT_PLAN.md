# Codex plan: exact four-class signed persistent-site MPZ port through stable-crack birth

## Physics decision

Choose architecture 1 from `V9_STABLE_CRACK_BIRTH_FOUR_CLASS_AUDIT.md`:

**Port the audited signed persistent-site MPZ/emission/transport state through stable-crack birth.**

Do not use the reduced scalar projection as the production four-class route.

Reason: the scientific requirement is transfer of the same audited parameter rows and local constitutive physics used by the canonical fracture/fatigue models, without refitting. The four rows actively use independent cleavage/emission surfaces, Peierls/Taylor kinetics, Taylor correlation, persistent source sites, and separate signed mobile/retained populations. Projecting these onto the existing scalar-density/scaled-common-barrier v9 state would discard active constitutive coordinates and could alter class ordering and temperature dependence. That would no longer be exact parameter transfer.

The weak-T versus DBTT contrast is a primary science target. Preserve it; do not tune the classes toward a common response.

## Scope boundary

The S-N endpoint remains:

```text
N_stable_crack_birth
```

Port only the physics required from the initial flaw/notch state through the first stable crack birth.

Do not port or develop post-stable growth, linkage, front propagation, physical handoff, crack-event length, or da/dN as part of this work package unless a pre-birth dependency truly requires a read-only quantity from that machinery.

A transient embryo that heals is not failure. The run terminates when the first persistent embryo takes the stabilization branch and enters the stable state at an accepted transactional boundary.

## Canonical parameter classes

Use exactly these audited identifiers:

```text
Peak     = v913_paper_peak01_0242980_persistent_sites
DBTT     = v913_paper_dbtt01_0202500_persistent_sites
weak-T   = v913_paper_weakT01_0257068_persistent_sites
ceramic  = v913_paper_ceramic01_0189364_persistent_sites
```

The later v10.4.x weak-T `0129902` and ceramic `0077080` rows are not permitted for this campaign.

Load all four through the audited registry/selection machinery already qualified in v9. Do not reconstruct parameter rows manually.

## Preserve the current v9 milestone

Before new work:

1. verify `codex/sn-endurance-limit-v9` is clean;
2. push the completed stable-crack endpoint and registry-provenance commits;
3. fetch and merge `origin/codex/sn-endurance-limit`;
4. confirm the frozen v8.7 hashes and current v9 endpoint results are unchanged.

Do not rewrite existing v9 result histories.

## Autonomous working mode

This is a longer implementation-and-calculation work package. Do not stop after each small audit or unit-test milestone.

Continue autonomously through the sequence below unless one of these genuine stop conditions occurs:

- an audited pre-birth state variable or rate law cannot be recovered unambiguously;
- exact parameter-row provenance cannot be verified;
- a necessary mapping would require a new constitutive assumption rather than a direct port;
- transactional partition/restart invariance cannot be achieved without changing the audited physics;
- a class/temperature run reveals a physical inconsistency that requires a user decision.

Otherwise fix implementation defects, add tests, commit validated milestones, and continue.

# Phase 1 — recover the exact pre-birth constitutive closure

Trace the exact canonical fracture/fatigue implementations that consume the four parameter rows. Build the minimum dependency closure required before stable-crack birth.

Recover and document at least:

- independent cleavage barrier surface;
- independent emission barrier surface;
- Peierls transport surface;
- Taylor trapping/release surface;
- Taylor-correlation/prefactor formulation;
- persistent source-site population and thresholds;
- two signed slip channels;
- mobile population for each signed channel;
- retained population for each signed channel;
- active unsigned mobile+retained quantity entering Taylor backstress;
- signed retained quantity entering direct shielding;
- accumulated slip/blunting state if active before stable birth;
- emission/transport/retention/escape/advection updates active before stable birth;
- stochastic cleavage first passage;
- stochastic emission and any persistent event identity needed pre-birth;
- temperature/rate scaling actually used by the audited rows.

For every recovered equation/state record:

```text
source repository
branch/commit or immutable artifact
source file
function/class
option/registry identifier
SHA-256
physical units
state variables read
state variables written
```

Write the dependency/provenance map into a versioned manifest. Do not rely on prose summaries when executable audited source exists.

# Phase 2 — define the v9 pre-birth state explicitly

Create new v9 modules for the four-class pre-birth physics. Do not overwrite the earlier scalar v9 baseline.

Suggested architecture:

```text
v9_four_class_state.py
v9_four_class_mpz.py
v9_four_class_rates.py
v9_four_class_adapter.py
```

Names may differ, but the separation should remain clear.

The state capsule must preserve, at minimum, all active pre-birth quantities from the audited constitutive closure, including signed mobile/retained fields and persistent source identities.

Do not collapse signed populations to one scalar density merely for convenience.

Derived scalar observables may be cached, but the authoritative state must retain the constitutive coordinates needed to reproduce the audited equations.

# Phase 3 — mechanics/constitutive coupling rule

The global/local mechanics may remain the v9 FEM/notch representation; exact PF/FEM displacement-field equality is not required.

However, when a local stress/strain/loading history is supplied to the material closure, the **constitutive response must match the audited fracture/fatigue implementation** for the same parameter row and state.

This is the desired division:

```text
mechanics representation may differ
constitutive parameter row and local kinetic/state laws do not
```

Therefore first validate constitutive parity under identical prescribed local histories before drawing conclusions from different mechanics.

# Phase 4 — constitutive parity tests before full conditions

For each of the four classes, construct prescribed local loading histories at representative temperatures and compare the canonical source engine with the new v9 port.

At identical initial state and identical stress/loading history, require parity for all relevant pre-birth observables, including as applicable:

- cleavage log hazard;
- emission log hazard;
- Peierls rate/hazard;
- Taylor rate/hazard;
- signed mobile update;
- signed retained update;
- escape/advection update;
- Taylor backstress;
- direct shielding;
- accumulated slip/blunting;
- persistent source thresholds/status;
- embryo/stabilization/healing clocks.

Use log-domain comparisons for ultra-small rates.

Do not accept parity based only on final scalar stress or one aggregate rate.

Where the source code is stochastic, use identical persistent thresholds/RNG state or compare the deterministic hazard/state kernel before sampling.

# Phase 5 — transactional invariance with the new state

Before a long class run, prove the ported signed-state trajectory satisfies v9 numerical contracts:

- >10x block-partition invariance;
- atomic restart equivalence;
- persistent site/transition thresholds unchanged;
- RNG/event sequence preserved;
- no state double commit;
- signed population conservation/consistency rules;
- cap/floor activations explicitly reported;
- log-domain ultra-small rates do not become numerical zeros;
- first stable-crack transition localized consistently.

Use one representative Peak condition first because it exercises both cleavage and Peierls/Taylor activity strongly.

# Phase 6 — first four-class 300 K qualification, one condition at a time

Once the port passes the above gates, qualify one working stable-crack-birth condition for each canonical class at 300 K:

```text
Peak
DBTT
weak-T
ceramic
```

Do not run them in parallel.

Do not assume the existing case64 stresses will produce convenient lives for every class. Use short no-state-change/preflight hazard evaluations and sequential bracketing to choose a physically useful first stress without launching a broad sweep.

For each class, obtain at least:

- one finite stable-crack-birth condition;
- one lower-stress long-life/right-censored or asymptotic-diagnostic condition when computationally practical.

Reuse/resume every valid condition through request hashes. Never rerun a completed condition merely because orchestration changes.

# Phase 7 — build 300 K stable-crack-birth S-N skeletons

For each class independently, add stresses sequentially to span useful decades in stable-crack-birth life.

Target useful coverage, not uniform stress spacing:

```text
~1e2 to 1e8+ cycles where physically accessible
```

Select each next stress from the completed results for that class. Do not force the same stress set on all four parameterizations if their responses differ substantially.

Maintain separate machine-readable tables for class, stress, temperature, seed, endpoint, censoring, and endurance diagnostics.

# Phase 8 — temperature dependence

Use the audited canonical temperature grid already recovered:

```text
300, 600, 800, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300 K
```

Do not launch the full grid immediately.

First use a compact diagnostic temperature set that can reveal the class character, for example low/intermediate/high anchors selected from the canonical grid. A defensible initial set is:

```text
300 K
900 K
1300 K
```

with additional points around the transition/peak region as needed from the canonical grid.

Proceed one class/temperature/condition at a time.

Expected qualitative hypotheses, not fitting targets:

- Peak may show a strongly nonmonotonic S-N/stable-birth response with an intermediate-temperature maximum/shift;
- DBTT may show a strong systematic transition with temperature;
- weak-T may remain much less temperature-sensitive than DBTT;
- ceramic should remain a comparatively brittle/weak-plasticity control.

A weak-T response that differs substantially from DBTT is desirable evidence that the transferred framework preserves material-class physics. Do not tune it away.

If the first compact temperature set preserves the expected class distinctions without numerical pathologies, fill the remaining canonical temperatures adaptively.

# Phase 9 — outputs

Maintain at least:

```text
sn_stable_crack_birth_results.csv
stable_crack_birth_stage_results.csv
stable_crack_birth_endurance_diagnostics.json
four_class_campaign_manifest.json
```

Each result row must include:

- exact option identifier;
- registry/source hashes;
- temperature;
- stress amplitude/range and R ratio;
- seed/source-site stochastic identity;
- first embryo cycle;
- first stable-crack cycle;
- censor horizon/status;
- cleavage/emission/PT diagnostic histories or summary measures;
- signed mobile/retained state summaries;
- backstress/shielding/blunting summaries;
- endurance evidence/classification;
- request/source hashes;
- wall time and accepted/rejected block statistics.

Legacy front-capture/physical-handoff outputs remain separate diagnostics only.

# Phase 10 — endurance/no-endurance science

The endurance endpoint is stable-crack birth.

Preserve the distinction between:

```text
no-embryo survival
stable-crack-birth survival
```

A finite total embryo-birth hazard is sufficient for a nonzero probability of no stable crack. Divergent embryo-birth hazard alone does not prove eventual stable crack because embryos can heal.

As four-class histories accumulate, evaluate how the signed mobile/retained state and PT kinetics cause the stable-crack-birth hazard to decay, plateau, or remain active.

Do not create endurance with a finite runout, stress cutoff, rate clipping, or numerical underflow.

# Preferred autonomous stopping point

Continue until a useful first four-class comparison exists, preferably including:

- exact constitutive parity tests for all four parameter rows;
- one validated 300 K stable-crack-birth S-N skeleton per class;
- at least an initial low/intermediate/high temperature comparison for Peak, DBTT, weak-T, and ceramic;
- explicit weak-T versus DBTT comparison;
- restartable campaign outputs and source provenance;
- no unresolved partition/restart defect.

Stop earlier only for a genuine constitutive/provenance ambiguity requiring a user decision.

Do not implement post-stable crack growth or da/dN as part of this campaign.