# Endurance asymptotic hardening gate after initial v9 audit

## Context

The first Codex baseline-recovery milestone on local branch `codex/sn-endurance-limit-v9` reported two commits:

- `9352853d026c018bf017f24aaca7ee0d9eb9f558` — audit, dependency manifest, endurance analyzer, synthetic tests
- `0eda22789d080fba3878d029fa214f3a1f8aedb7` — exact runnable v8.7 closure and historical regression fixtures

Keep those commits. Do **not** start a production S-N campaign or the adaptive `1e12`-cycle solver yet.

The next gate is to harden the endurance classification and analyze the preserved K360 no-birth censor trajectory before changing baseline physics.

## 1. Correct the interpretation of finite-window tail fitting

A finite observed tail that is fit well by an exponential or by a power law with `p > 1` is evidence for integrability, but it is not proof of a true endurance population.

For example,

```text
h(N) = A exp(-N/tau) + h_floor
```

with `h_floor > 0` has

```text
integral_0^infinity h(N) dN = infinity
```

and therefore has no strict endurance population even if the observed finite interval is almost perfectly exponential.

Likewise an observed `p > 1` power-law segment can cross at later N to a positive floor or to `p <= 1`.

Revise the analysis architecture to distinguish at least:

```text
empirical_tail_evidence:
    integrable_tail_observed
    divergent_tail_observed
    ambiguous_tail_observed

physics_asymptotic_classification:
    endurance_supported
    no_endurance
    undetermined
```

A finite-window regression alone may set `empirical_tail_evidence = integrable_tail_observed`, but it must not set `physics_asymptotic_classification = endurance_supported` unless a model-derived asymptotic upper bound on remaining integrated hazard is also available.

## 2. Mandatory adversarial tests

Extend `tests/test_endurance.py` with at least:

```text
exponential + positive hazard floor
    -> not endurance_supported

p = 1.5 power law + positive hazard floor
    -> not endurance_supported

p > 1 over the observed interval, then late crossover to p = 0.8
    -> not endurance_supported when the crossover is represented in the supplied history

short integrable-looking tail with insufficient dynamic range
    -> undetermined

zero/clipped numerical tail
    -> undetermined
```

The current mandatory idealized tests remain required:

```text
constant positive hazard -> no_endurance
pure exponential decay   -> endurance_supported in the synthetic model-closed test
p = 0.5                  -> no_endurance
p = 1.0                  -> no_endurance
p = 1.5                  -> endurance_supported in the synthetic model-closed test
short/noisy tail         -> undetermined
```

For the pure synthetic exponential/power-law tests it is acceptable to declare endurance because the generating law is known by construction. For an empirical finite trajectory, fit quality alone is not sufficient.

Rename any field called `residual_hazard` when it actually represents

```text
integral_N^infinity h(n) dn
```

to `residual_integrated_hazard`.

## 3. Preserved K360 v2.8 reference data

Read-only reference data are now committed under:

```text
reference/k360_v2_8_tail/
```

Important files:

```text
K360_CENSOR_DERIVED_TAIL.csv
K360_CENSOR_TAIL_AUDIT.json
README.md
```

Do not rerun this expensive K360 censor simply to recover its tail. Use the preserved derived history. The full July runtime archive remains provenance outside the Git repository.

The censor condition is approximately:

```text
sigma_a = 690.443200 MPa
N_end = 1.0e8 cycles
realized births = 0
cumulative expected births at runout = 0.622315942
```

`K360_CENSOR_DERIVED_TAIL.csv` defines an interval-average aggregate birth hazard from the preserved no-birth trajectory as

```text
h_birth,agg = Delta(pd_expected_births_cumulative) / Delta N.
```

Nested 50/65/80 percent late-tail power-law fits give approximately:

```text
p = 2.174
p = 2.249
p = 2.296
```

with R^2 > 0.999 over those windows. A naive extrapolation of that observed power law gives a survival asymptote near `0.534`.

Treat this as strong empirical evidence of an integrable tail over the observed interval, **not as proof of a true endurance limit**.

Reproduce these numbers independently from the committed derived tail data and add a regression test or small analysis script so the derivation is transparent and repeatable.

## 4. Audit the actual baseline birth physics analytically

The v8.7 lineage inherits the completion-gated independent-cleavage structure:

```text
dLambda/dt = r_delivery(t) - Lambda/tau_delivery
h_birth(t) = h_cleave(t) Q(K, Lambda(t))
```

with persistent candidate-site Exp(1) thresholds in cumulative birth hazard.

For the production hit count `K = 2`, use the small-memory asymptotic

```text
Q(2, Lambda) ~ Lambda^2 / 2,    Lambda -> 0.
```

Derive sufficient conditions on the asymptotic cyclic delivery activity for the total birth hazard

```text
H_birth(infinity) = integral_0^infinity h_birth(N) dN
```

to be finite or divergent.

Examples to analyze explicitly:

1. If `h_cleave -> c > 0` and `r_delivery -> r_inf > 0`, then `Lambda` approaches a positive periodic/steady value and the birth hazard generally retains a positive asymptote -> no birth endurance.
2. If `r_delivery -> 0` and `Lambda` follows its stable finite-memory decay, then `Q(2,Lambda)` can decay quadratically in the activity/memory. Determine when that makes the birth hazard integrable.
3. If activity decays algebraically, derive the corresponding exponent condition for K=2.
4. Audit whether numerical floors/caps can prevent the physical delivery state from actually reaching its asymptotic regime.

Write the derivation and code-level mapping in:

```text
docs/ENDURANCE_ASYMPTOTIC_AUDIT.md
```

## 5. Birth endurance is not identical to physical-handoff endurance

The S-N failure endpoint is PD `physical_handoff`, not first birth.

Report at least two distinct diagnostics:

```text
birth_endurance_diagnostic
physical_handoff_endurance_classification
```

A finite total birth hazard is a sufficient condition for a nonzero infinite-life handoff population because some realizations never nucleate.

However, a divergent birth hazard does not by itself prove eventual physical handoff because embryos can heal and stabilization/growth/linkage probabilities may evolve.

Therefore v9 diagnostics should preserve and expose separately:

```text
cleavage rate / hazard
cyclic delivery rate
Lambda delivery memory
completion Q(K,Lambda)
aggregate birth hazard and cumulative hazard
persistent site thresholds and residual threshold distances
stabilization hazard
healing hazard
stable-defect growth activity
linkage/front activity
physical handoff state
```

## 6. Refine the mechanism comparison

Do not describe the current baseline simply as "activity independent". It is more accurately:

```text
completion-gated independent cleavage
```

because cleavage is a separate Arrhenius channel but birth is gated by multi-hit cyclic delivery completion.

For a contrasting clearly no-endurance hypothesis, consider a separately named model such as:

```text
persistent_direct_cleavage
```

where a physically justified nonzero direct initiation channel remains after cyclic delivery has shaken down. Do not insert an arbitrary stress cutoff or damage-per-cycle constant.

A useful scientific comparison is then:

```text
A. completion-gated cleavage capable of true shakedown/endurance
B. persistent direct thermally activated initiation capable of very long but ultimately finite life
```

Any new mechanism must be independently selectable and diagnostic.

## 7. Numerical/cap audit before large-N development

Before claiming endurance, explicitly audit at least:

```text
crack_floor_frac
Arrhenius exponent clipping / underflow handling
density caps
back-stress caps
stress-amplification caps
state-shift clipping
plastic increment caps
any hazard/rate lower bounds or zeroing thresholds
```

Classify each as:

```text
physical constitutive assumption
numerical safety bound
finite-domain/model-resolution bound
```

and determine whether it can change `H(infinity) < infinity` versus `H(infinity) = infinity`.

No numerical floor, underflow, clipping, or finite runout may manufacture endurance.

## 8. Stop condition for this gate

Do not start a long S-N or optimization campaign during this milestone.

Stop after all of the following pass:

1. hardened empirical-vs-physics endurance classifier;
2. adversarial tests above;
3. reproducible K360 preserved-tail analysis;
4. `docs/ENDURANCE_ASYMPTOTIC_AUDIT.md`;
5. explicit sufficient-condition derivation for v8.7 K=2 completion gating;
6. birth-vs-handoff endurance distinction in the data model;
7. caps/floors asymptotic audit;
8. all existing baseline/regression tests still pass.

Then report:

- branch and HEAD;
- commits created;
- tests and exact results;
- reproduced K360 tail metrics;
- whether the v8.7 baseline is structurally capable of a true endurance population;
- what remains only empirical evidence versus mathematically/physically bounded;
- proposed v9 large-N transactional stepping architecture;
- proposed no-endurance comparison mechanism.

Do not launch the `3 stresses x 3 seeds x 2 mechanisms` pilot until this gate is reviewed.
