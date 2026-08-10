# V9 PD production progress

## Current capability

- Synthetic milestone `f7c8fb1`: protected-state dormant engine, periodic and
  projective modes, event guard, atomic mode diagnostics.
- Production repair in progress: monotone ledgers now advance in stationary
  and projective segments; event calculations remain in log space below
  Float64 linear representability; spatial active coordinates now use declared
  field-specific transforms and constitutive bounds.
- Focused kernel qualification: 8 tests passing, including a formal `1e300`
  skip with `log(q)=-1000`.
- A real cached-FEM/spatial-PD exact-cycle callback now executes through the
  production phase ordering. A 10-cycle direct/callback comparison agrees in
  FEM state to roundoff, PD memory/population state to about `1e-12` relative,
  and cumulative log action to about `2e-13` relative. No threshold, RNG,
  site-status, bond, topology, geometry, or physical-cycle state changes during
  private evaluation.
- The real high-cycle orchestration path is opt-in via `--pd-high-cycle`, writes
  mode/controller histories, returns to the direct loop at an event guard, and
  uses no PD images.

## Real production status

The historical finite-life/front control and long-horizon production
qualification remain in progress. No synthetic result is classified as a
material calculation. Existing v3 and historical PD trajectories remain
unchanged.

## Historical K360 real replay

The exact archived 735.921495 MPa request was replayed with `pd_image_policy =
none`. The integrated path used high-cycle qualification attempts with
transactional fallback through curved transients and reproduced the ordered
sequence at node 23:

| event | preserved control N | replay N | relative shift |
|---|---:|---:|---:|
| first embryo | 1,793,122.651 | 1,805,851.926 | 0.71% |
| stable site | 1,793,124.874 | 1,805,854.149 | 0.71% |
| first softening | 1,817,260.678 | 1,821,610.036 | 0.24% |
| root connection | 2,300,120.914 | 2,367,711.399 | 2.94% |
| front capture | 2,351,351.153 | <=2,400,000 | <=2.07% |

The replay initially exposed a sub-ULP event-time stall and loss of controller
scale across localization. Both are repaired. The corrected trajectory retains
primary seed 23 with zero reselections, reaches maximum bond damage 0.974, 18
broken bonds, root connection, and active-front capture.

## Real low-stress high-cycle qualification

The archived shielded K360 690.443200494 MPa trajectory has been preserved and
continued atomically beyond `N=1.60e7` without a realized embryo, stable site,
or topology change. Real private-map qualification exposed and repaired three
numerical-coordinate defects: spatial ledger variation was confused with time
variation; structural-zero expected populations were extrapolated as an
exponential source; and the quasistatic FEM displacement warm start was treated
as constitutive memory. Expected birth/healing ledgers are now closed by the
validated population conservation identities rather than an inaccurate
independent extrapolation.

The real engine has accepted and persisted projective segments (38 cycles in
the first fully conservation-closed qualification, with state error at most
`3.53e-6` and log-hazard validation error `1.73e-4`). This proves adoption of a
real projected physical state, but one-cycle training still limits useful
segment growth. The next numerical task is a multi-cycle exact training burst;
the direct embedded macro-stepper remains the accepted fallback and no
accelerated low-stress material result is yet claimed.

The same archived 690.443200494 MPa trajectory subsequently completed an
atomic `N=1e8` boundary (generation
`generation_3655f6ccf3044ebf8a826261d9987b23`). It remains dormant: 4,349
available realized sites, zero realized births/stable sites, zero bond damage,
and no seed/front/topology transition. The maximum nodal cumulative birth
action is `0.140255573333055`, the maximum final per-cycle birth intensity is
`6.978364374432825e-13`, and the minimum remaining persistent-site threshold
action is `9.647106699392403e-6`. This is a finite runout/checkpoint, not an
endurance classification. Expected-population diagnostics reach about 0.627
births/stable sites and remain distinct from the realized persistent clocks.

Production fallback was bounded with a private-map efficiency budget. Once
the projective fit ceased providing useful cycles per exact map, control
returned to the authoritative embedded macro-stepper and deferred another
qualification until the next logarithmic checkpoint. Mode-history writes are
append-only across subsequent atomic resumes.

An extension of the same condition reached atomic generation
`generation_c18519735d8640a8b6cad153c26158ff` at
`N=168634947.0289538`, again with zero realized births, stable sites, damage,
or topology change. The maximum cumulative action increased only from
`0.140255573333055` to `0.140256540348356` over the additional 6.86e7 cycles;
the nearest remaining threshold changed from `9.647106699392403e-6` to
`9.645864666437046e-6`. The actual final maximum birth rate is
`1.4023253737703178e-13` per cycle and the direct controller estimates a
`2.00129794714e10`-cycle next wait. This strong rate decay is scientifically
important but is not by itself a proof of endurance.

The current single-cycle projective fit is not production-efficient on this
tail: exact state curvature is below `1e-9`, but two-cycle log-hazard prediction
error remains about `9.4e-4`, above the `2e-4` admission tolerance. A longer
exact training window/rate-separated fit is therefore required before a
credible 1e10--1e12 continuation; brute-force continuation was stopped at the
hash-verified atomic boundary.
