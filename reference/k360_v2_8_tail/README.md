# K360 v2.8 tail reference

This directory preserves the derived late-tail evidence from the July v2.8 K360 censor trajectory so the large-N work does not need to rerun that expensive validation case merely to recover the observed hazard decay.

## Provenance

The source runtime archive was:

```text
stateful_pd_kitagawa_v2_8_gated_pilot(2).zip
```

A separate extracted reference bundle was generated with SHA-256:

```text
f63fae1d03d37d2b2cd01068f5d3c04d8846a775250290003f3fcf3f4254b260
```

The full runtime archive/checkpoint is not committed here. These committed files are read-only derived evidence for the endurance analysis.

## Censor record

- sigma_a = 690.443200 MPa
- N_end = 100000000 cycles
- cumulative expected births at runout = 0.622315942
- realized births = 0 in the preserved censor trajectory

`K360_CENSOR_DERIVED_TAIL.csv` computes an interval-average aggregate birth hazard as

```text
d(pd_expected_births_cumulative) / dN
```

on the no-birth censor trajectory.

Nested 50/65/80% tail fits give power-law exponents approximately:

- 2.174
- 2.249
- 2.296

All are > 1 and have R^2 > 0.999 in this preserved interval. A naive extrapolation gives a survival asymptote near 0.534.

## Critical caveat

This is **not proof of a true endurance limit**. A finite observed interval can look exponentially or super-harmonically decaying while hiding a later positive hazard floor or a crossover to p <= 1. The v9 analysis must combine empirical tail fits with structural/physics-based asymptotic bounds.

Do not reinterpret `right_censored` at 1e8 cycles as endurance by itself.

See also:

```text
docs/ENDURANCE_ASYMPTOTIC_NEXT_STEPS.md
```
