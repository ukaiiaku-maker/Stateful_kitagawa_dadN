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
