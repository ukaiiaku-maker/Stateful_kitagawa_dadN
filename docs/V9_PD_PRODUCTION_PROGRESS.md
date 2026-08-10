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
