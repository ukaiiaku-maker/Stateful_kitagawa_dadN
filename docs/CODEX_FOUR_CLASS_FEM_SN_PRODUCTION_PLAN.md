# Codex autonomous plan: full-FEM canonical four-class stable-crack-birth S-N production

## Governing scope

This plan follows the qualified `v9_canonical_four_class_cleavage_stable_birth_v1` implementation.

The canonical four-class S-N endpoint is the **first canonical sharp-front cleavage first-passage event at the fixed root/front-local process zone**:

```text
N_stable_crack_birth = N_first_canonical_cleavage_first_passage
```

Do not reintroduce the legacy PD candidate-site/K=2/stabilization model into this branch. Do not run crack growth, front capture, physical handoff, crack-event length, or da/dN after the event.

The authoritative four classes remain:

```text
Peak     = v913_paper_peak01_0242980_persistent_sites
DBTT     = v913_paper_dbtt01_0202500_persistent_sites
weak-T   = v913_paper_weakT01_0257068_persistent_sites
ceramic  = v913_paper_ceramic01_0189364_persistent_sites
```

The authoritative pre-birth constitutive model is:

- audited v10.2.21 aggregate persistent emission (deterministic conditional evolution);
- exact signed mobile/retained/wake MPZ state;
- independent audited emission, Peierls and Taylor kinetics;
- Taylor backstress, signed shielding and blunting;
- canonical stochastic cleavage first-passage clock only;
- no fictional stochastic emission ledger.

The frozen audited rows and source hashes must remain unchanged.

## Important statistical consequence

For this canonical branch, the pre-birth signed-MPZ trajectory is deterministic conditional on survival and the only stochastic fracture variable is the single persistent cleavage threshold

```text
Xi ~ Exp(1).
```

Therefore, for a fixed stress/temperature/material condition, the deterministic cumulative canonical cleavage hazard `H(N)` defines the **exact stable-crack-birth survival function**

```text
S_stable(N) = exp[-H(N)].
```

This should be used directly rather than relying on large Monte Carlo seed populations to estimate the S-N survival curve.

Quantile lives may be obtained from the same deterministic no-event trajectory by solving

```text
H(N_q) = -log(S_q)
```

for desired survival probabilities. For example:

```text
P(failure)=0.10 -> H = -log(0.90)
P(failure)=0.50 -> H = log(2)
P(failure)=0.90 -> H = -log(0.10)
```

A single sampled threshold/event trajectory remains useful as an event-localization/restart validation, but repeated seeds are not required merely to reconstruct the survival distribution when the deterministic-path assumptions above hold.

If any additional stochastic pre-birth process is later introduced, this exact-survival shortcut must be re-audited.

# Phase 1 — connect the canonical endpoint to the full accepted FEM cyclic history

Do not stop after a prescribed-history unit test. Wire the qualified canonical state into the actual v9 FEM cyclic driver.

For every accepted cyclic interval:

1. obtain the same resolved root-local opening and signed shear/stress history from the accepted FEM state;
2. advance the signed MPZ using the qualified canonical plastic-half-step / midpoint-cleavage / plastic-half-step ordering;
3. transform root loading to the canonical K-space driving coordinate in the already-qualified order:
   - nominal root stress -> nominal K using the initial root radius;
   - subtract signed active shielding in K space;
   - apply the current blunted-radius factor;
4. evaluate the exact registry-selected cleavage barrier and canonical `m=3`, `tau_c=1e-6 s` renewal law;
5. integrate cumulative cleavage hazard in log-safe form;
6. localize a threshold crossing inside the transaction without committing unused time;
7. terminate immediately after the first canonical cleavage crossing when running event mode.

No PD candidate-site state is active in this path.

## Full-coupling qualification

Before production conditions, require:

- prescribed-history parity remains exact for all four rows;
- full-FEM root-history parity against the source engine for a short accepted trajectory;
- >10x real block-partition convergence of signed state and cumulative cleavage hazard;
- atomic restart equivalence including FEM, signed MPZ, cleavage threshold, cumulative hazard and RNG state;
- no post-event RNG draw or crack-growth state modification;
- stable source/request/registry hashes in every generation.

Do not stop merely because these pass. Continue directly into the first physical condition.

# Phase 2 — qualify Peak at 300 K as the first full-FEM working condition

Use the already qualified K360 fixed-notch/FEM cyclic mechanical setup as the initial mechanics anchor unless the exact four-class campaign manifest requires a different pre-birth mechanical convention. Do not invent a new geometry, R ratio, frequency, or waveform.

Start with the already-working nominal stress neighborhood around the prior K360 conditions only as a numerical bracket seed, **not as a calibration target**.

The first objective is to obtain a trustworthy deterministic hazard trajectory and a stable-birth event/quantile bracket for Peak at 300 K.

If the starting stress is obviously outside a useful range:

- event much earlier than the desired S-N range -> reduce stress;
- negligible cumulative hazard through the long-life horizon -> increase stress;
- otherwise retain the condition and use its `H(N)` trajectory to choose the next stress.

Change only stress when bracketing. Do not retune material/barrier parameters.

For every condition record at least:

```text
material_class
option_id
T_K
sigma_a_MPa
R
frequency_Hz
request_hash
registry_hashes
H_cleave(N)
log_h_cleave(N)
S_stable(N)
N_sampled_event if a sampled threshold is used
N10, N50, N90 where crossed
terminal H at observation horizon
endurance_tail_diagnostics
signed mobile/retained/wake state summaries
backstress
shielding
blunting
accepted/rejected block counts
restart identity
```

# Phase 3 — build the 300 K four-class S-N curves sequentially

After the Peak 300 K full-FEM path is qualified, run classes sequentially:

```text
Peak -> DBTT -> weak-T -> ceramic
```

Do not launch them in parallel.

For each class:

1. construct an adaptive stress skeleton using one condition at a time;
2. use cumulative-hazard/quantile information from completed conditions to select the next stress;
3. fill useful gaps in `log10(N50)` rather than uniform stress spacing;
4. retain long-life conditions as asymptotic/endurance diagnostics rather than forcing an event;
5. stop adding nearby stress points when they are redundant.

Initial coverage goal, where physically accessible, is several decades of median stable-birth life, approximately spanning the useful range between `1e3` and `1e8+` cycles. These are coverage goals, not fitting constraints.

Because the exact survival curve is available from `H(N)`, report probabilistic S-N bands such as `N10/N50/N90`, not only one sampled event life.

Maintain the four classes separately; do not force a common stress grid if their useful brackets differ.

# Phase 4 — temperature dependence

Once each class has a valid 300 K skeleton, begin temperature dependence with the compact diagnostic subset:

```text
300 K, 900 K, 1300 K
```

for each class, one condition at a time.

Use the same canonical audited option row at every temperature; only the model's audited temperature dependence may change the response.

The expected behaviors are hypotheses, not fitting targets:

- Peak may be strongly nonmonotonic with temperature;
- DBTT may show a strong transition;
- weak-T should be allowed to remain substantially less temperature-sensitive;
- ceramic should remain a distinct weak-plasticity/brittle control.

The weak-T versus DBTT contrast is a central science output. Never tune weak-T toward DBTT or impose a shared S-N shape.

After the 300/900/1300 K diagnostic matrix is trustworthy, fill additional audited temperatures adaptively from:

```text
300, 600, 800, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300 K
```

Prioritize temperatures where the hazard/S-N behavior changes most rapidly rather than simply running the full grid blindly.

# Phase 5 — endurance classification

For every material/stress/temperature condition, retain the exact distinction between finite observation and asymptotic endurance.

The relevant quantity is now the canonical cleavage cumulative hazard:

```text
H_cleave(N) = integral h_cleave(N) dN.
```

True stable-crack-birth endurance requires

```text
H_cleave(infinity) < infinity
```

so that

```text
S_stable(infinity) = exp[-H_cleave(infinity)] > 0.
```

A finite `1e8` or `1e12` runout alone is not endurance.

Use empirical tail fits only as evidence. A physics classification must remain `undetermined` unless the signed-MPZ/cleavage equations support a finite upper bound or divergence proof.

Track whether the long-N behavior is controlled by:

- emission/activity extinction or persistence;
- Peierls/Taylor transport;
- retained signed shielding;
- Taylor backstress;
- blunting;
- intrinsic cleavage barrier behavior.

# Required outputs

Maintain continuously and atomically:

```text
runs/sn_v9_canonical_four_class/campaign_manifest.json
runs/sn_v9_canonical_four_class/sn_stable_crack_birth_results.csv
runs/sn_v9_canonical_four_class/sn_quantile_results.csv
runs/sn_v9_canonical_four_class/survival_results.csv
runs/sn_v9_canonical_four_class/endurance_diagnostics.json
runs/sn_v9_canonical_four_class/condition_registry.json
```

`survival_results.csv` must contain canonical **stable-crack-birth survival** from `exp(-H_cleave)`, clearly identified as such.

Do not mix legacy PD/K=2 or physical-handoff results into these tables.

# Autonomous work mode

Continue through implementation and sequential calculations without stopping after each audit, unit test, or smoke.

Commit and push validated milestones as useful, but continue automatically unless a genuine stop condition is reached.

Genuine stop conditions are limited to:

1. an unresolved physical mapping not defined by the audited source;
2. failure of real block-partition/restart equivalence that cannot be repaired without changing constitutive physics;
3. missing/ambiguous audited parameter provenance;
4. a new result that demonstrates an actual constitutive inconsistency requiring a user decision;
5. repository/data corruption that prevents trustworthy continuation.

Performance problems are not automatic stop conditions: profile and optimize the validated algorithm while preserving accepted-trajectory equivalence.

# Preferred stopping point for this work package

Continue until there is a useful first canonical four-class result set, preferably including:

- full-FEM canonical endpoint qualification;
- a 300 K probabilistic stable-birth S-N skeleton for all four classes;
- at least low/intermediate/high-temperature (`300/900/1300 K`) comparisons where computationally practical;
- explicit Peak/DBTT/weak-T/ceramic contrasts;
- stable-birth quantile curves from canonical cumulative hazard;
- initial endurance diagnostics;
- no unresolved partition/restart defect.

Do not perform post-birth crack growth or da/dN in this work package.
