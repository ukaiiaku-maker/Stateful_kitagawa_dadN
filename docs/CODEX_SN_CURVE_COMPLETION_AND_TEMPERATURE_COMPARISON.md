# Codex plan: complete canonical S-N curves and compare temperature response

## Governing state

The final canonical four-class rows are now:

```text
Peak     = v913_paper_peak01_0242980_persistent_sites
DBTT     = v913_paper_dbtt01_0202500_persistent_sites
weak-T   = v913_paper_weakT01_0129902_persistent_sites
ceramic  = v913_paper_ceramic01_0077080_persistent_sites
```

The earlier weak-T `0257068` and ceramic `0189364` results are preserved as rejected-transfer controls only.

Do not alter the common canonical cleavage/signed-MPZ physics.

The endpoint remains the first canonical cleavage first passage, with

```text
S_stable(N) = exp[-H_cleave(N)].
```

The current full-FEM coupling, block-partition, restart, and provenance gates are accepted.

## Scientific interpretation of the current results

The 300 K calculations show very steep stress-life transitions rather than broad nearly linear S-N trends. Existing examples include:

```text
Peak 804.139 MPa: N50 ~ 4.71e2 cycles
Peak 770 MPa:     N50 ~ 5.20e4 cycles
Peak 735.921 MPa: N50 not crossed by 1e8

weak-T 450 MPa:   N50 ~ 4.56e2 cycles
weak-T 400 MPa:   N50 ~ 3.36e6 cycles

ceramic 500 MPa:  N50 ~ 1.08e1 cycles
ceramic 450 MPa:  N50 ~ 1.94e4 cycles
ceramic 400 MPa:  N10 not crossed by 1e8

DBTT 770 MPa:     N50 ~ 1.14e4 cycles
```

These sharp transitions mean the next task is to resolve each class's transition region adaptively, not add broad coarse stress points.

## Phase 1 — complete the 300 K S-N skeletons

Proceed one condition at a time.

For each class, choose the next stress to fill the largest useful gap in `log10(N50)` or to bracket the long-life transition.

Do not use a common stress grid across classes.

### Peak

Existing points already bracket a very large gap between `770 MPa` and `735.921 MPa`. Add only as many intermediate stresses as needed to resolve median lives in the approximate `1e5–1e8+` range and the transition to long-life behavior.

Also add one intermediate point between `770` and `804.139 MPa` only if needed to resolve the short-life slope.

### DBTT

Only one useful 300 K median point currently exists. Build both a lower-stress long-life bracket and a higher-stress short-life bracket, then fill the largest `log10(N50)` gaps adaptively.

### weak-T

Resolve the very sharp transition between `400` and `450 MPa`. Use adaptive bracketing to obtain useful median lives across several decades rather than adding uniformly spaced stresses.

### ceramic

Resolve the transition between `400` and `450 MPa`, and retain `500 MPa` only as the short-life anchor. Do not run higher stress unless it adds genuinely new information.

For every accepted condition record `N10`, `N50`, `N90`, terminal `H`, survival, signed-MPZ state summaries, and endurance-tail diagnostics.

## Phase 2 — compare temperature using full S-N / iso-life quantities

The current one-cycle temperature preflight is a bracketing diagnostic only. It is not the primary material comparison because different classes use different anchor stresses.

The 900/1300 K preflight currently shows:

```text
Peak 770 MPa:
  H_1cycle(900 K)  ~ 1.10
  H_1cycle(1300 K) ~ 92.2

DBTT 770 MPa:
  H_1cycle(900 K)  ~ 0.340
  H_1cycle(1300 K) ~ 119

weak-T 400 MPa:
  H_1cycle(900 K)  ~ 1.55e-8
  H_1cycle(1300 K) ~ 4.29e-8

ceramic 450 MPa:
  H_1cycle(900 K)  ~ 21.2
  H_1cycle(1300 K) ~ 338
```

Use those values only to choose initial stress brackets at each temperature.

For quantitative temperature comparison, construct either full S-N curves or iso-life stress metrics.

Preferred comparison outputs are stresses required to reach selected median lives, for example:

```text
sigma_a(N50 = 1e3)
sigma_a(N50 = 1e5)
sigma_a(N50 = 1e7)
```

where physically accessible.

These iso-life stresses provide a direct comparison of Peak, DBTT, weak-T and ceramic temperature sensitivity without forcing a common stress grid.

## Phase 3 — complete 300 / 900 / 1300 K comparison

After the 300 K transition regions are adequately resolved, proceed sequentially through `900 K` and `1300 K`.

At each new material/temperature pair:

1. use the one-cycle diagnostic only to seed a safe stress bracket;
2. run one full deterministic hazard trajectory;
3. update the bracket from the resulting `N50` or long-life observation;
4. continue until the useful S-N transition is resolved;
5. stop adding points once interpolation in `sigma_a` versus `log10(N50)` is well constrained over the target range.

Do not run a high-stress case merely because the one-cycle hazard is large. Reduce stress to move the full median life into an informative range.

## Weak-T versus DBTT interpretation

Do not describe weak-T as having weak fatigue-temperature dependence merely because its original monotonic fracture label is weak-T.

The current full weak-T trajectories show:

```text
400 MPa, 300 K: N50 ~ 3.36e6
400 MPa, 900 K: N10 ~ 1.00e7; N50 > 1e8
400 MPa, 1300 K: N10 ~ 2.91e6; N50 > 1e8
```

This already demonstrates a response qualitatively distinct from DBTT, but the magnitude and sign of the fatigue-temperature shift must be determined from completed S-N/iso-life curves, not from the class name.

Preserve any nonmonotonicity that emerges. Do not tune weak-T toward DBTT.

## Endurance logic

For every long-life point, retain the distinction:

```text
finite horizon
empirical integrable-tail evidence
physics asymptotic classification
```

Do not label a finite `1e8` or `1e12` trajectory endurance by itself.

Use the canonical cumulative cleavage hazard for classification:

```text
H_cleave(infinity) < infinity  -> endurance-capable survival asymptote
H_cleave(infinity) = infinity  -> eventual stable crack birth with probability 1
```

Track which signed-MPZ state controls the tail: emission activity, Peierls/Taylor transport, retained shielding, Taylor backstress, blunting, or intrinsic cleavage barrier.

## Required campaign outputs

Continue maintaining atomically:

```text
runs/sn_v9_canonical_four_class/campaign_manifest.json
runs/sn_v9_canonical_four_class/sn_stable_crack_birth_results.csv
runs/sn_v9_canonical_four_class/sn_quantile_results.csv
runs/sn_v9_canonical_four_class/survival_results.csv
runs/sn_v9_canonical_four_class/endurance_diagnostics.json
runs/sn_v9_canonical_four_class/condition_registry.json
```

Add an iso-life temperature table, for example:

```text
runs/sn_v9_canonical_four_class/iso_life_temperature_results.csv
```

with fields at least:

```text
material_class
T_K
target_N50
sigma_a_MPa
bracket_low_MPa
bracket_high_MPa
interpolation_method
source_condition_ids
```

## Autonomous work mode

Do not stop after another single audit or stress point. Continue through:

1. completed 300 K transition-resolution skeletons for all four final classes;
2. initial 900 K S-N/iso-life results;
3. initial 1300 K S-N/iso-life results;
4. explicit Peak/DBTT/weak-T/ceramic temperature comparison;
5. endurance-tail diagnostics for the long-life branches.

Stop only for a genuine constitutive/provenance/partition/restart blocker.