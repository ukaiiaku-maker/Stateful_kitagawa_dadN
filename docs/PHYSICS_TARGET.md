# Physics target: large-N S-N response and fatigue endurance

## Primary question

Can the stateful model distinguish a material/system that has a true fatigue-endurance regime from one that merely has a very long fatigue life?

The distinction must arise from the evolving physical hazard and state, not from an arbitrary runout cycle or a hard-coded stress cutoff.

## Statistical definition

Let `h_i(N)` be the failure/initiation hazard for persistent candidate site `i` along a no-failure trajectory and

```text
H_i(N) = integral_0^N h_i(n) dn.
```

For exponential first-passage thresholds, the no-failure survival probability is controlled by the accumulated hazard. A useful system-level total is

```text
H_tot(N) = sum_i H_i(N)
S(N) = exp[-H_tot(N)]
```

when the site hazards are conditionally independent along the survival path.

The important infinite-cycle distinction is:

### Endurance-capable regime

```text
H_tot(infinity) < infinity
```

so

```text
S(infinity) > 0.
```

A nonzero fraction of realizations survive indefinitely.

### No-endurance regime

```text
H_tot(infinity) = infinity
```

so

```text
S(infinity) = 0.
```

Eventual failure probability tends to one even if the instantaneous hazard becomes extremely small.

Examples for tail diagnostics:

- `h ~ constant > 0`: divergent, no endurance.
- `h ~ C/N^p` with `p <= 1`: divergent, no endurance.
- `h ~ C/N^p` with `p > 1`: integrable, endurance-capable.
- `h ~ C exp(-N/N0)`: integrable, endurance-capable.

A finite runout at `1e8`, `1e10`, or even `1e12` cycles cannot by itself establish which class applies.

## Mechanistic comparison to implement

Preserve the current independent-cleavage closure as one reference model. Then add a separately named, physically motivated closure in which initiation hazard requires continuing cyclic irreversible activity / delivery. The intent is to test two hypotheses, not to force a desired answer:

1. **Independent thermally activated pathway**: a residual cleavage/initiation hazard can remain after plastic shakedown. If the long-N hazard has a positive or non-integrable tail, there is no strict endurance limit.
2. **Activity-conditioned pathway**: initiation requires cyclic irreversible activity. If the material shakes down and the relevant activity decays sufficiently rapidly, the failure hazard can become integrable and an endurance regime can emerge.

Do not create endurance by:

- setting hazard to zero below an arbitrary stress;
- declaring `N > N_runout` to mean infinite life;
- clipping a small positive rate to zero for numerical convenience;
- fitting a horizontal S-N line and calling it mechanistic.

Any new activity gate must be expressed in terms of an existing physical state (plastic increment, dissipated work, mobile emission/delivery activity, etc.), have units documented, and be independently inspectable in output.

## Required observables

For every stress and seed, retain at minimum:

- `N` and log10(N);
- status: physical handoff / right censor / invalid;
- instantaneous total birth/initiation hazard;
- accumulated total hazard;
- maximum site hazard and accumulated site-hazard margin to threshold;
- plastic activity and plastic work rate;
- dislocation-density / back-stress / shielding state;
- active candidate count;
- embryo, stable-defect, softening, front-capture, and handoff event times;
- geometry/morphology validity metrics.

For the S-N campaign, report failure probability or survival probability versus N in addition to individual deterministic seed lives.

## Initiation versus propagation

The first milestone is an S-N curve to physical handoff. Crack propagation / `da/dN` after handoff is a separate observable. Do not combine initiation life and propagation life implicitly. A later state-preserving handoff to the sharp-front model can supply total-life calculations, but it must remain separable in the data model.
