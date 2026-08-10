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
