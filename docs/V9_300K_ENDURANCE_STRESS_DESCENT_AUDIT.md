# v9 300 K Endurance Stress-Descent Audit

## Scope

The accepted canonical elastic-FEM/signed-MPZ v3 model was used unchanged.
No temperature cases, post-birth propagation, parameter changes, dense knee
sweep, or pointwise-stationarity campaign was performed.

Class order was Peak, DBTT, weak-T, ceramic.  Each reported point is an actual
v3 trajectory with atomic checkpoints at `1e8`, `1e10`, `1e12`, and `1e14`.

## Adaptive selection

Sparse initial-cycle predictors selected one materially lower stress per
class.  Each actual `1e8` result predicted an informative `1e14` action, so the
condition was continued.  A second, lower stress was then selected from the
actual two-point log-hazard response to bracket `S(1e14)=0.9`.  This produced
eight descent conditions rather than a uniform sweep.

## Actual 1e14 brackets

| class | lower MPa | H | S_stable | upper MPa | H | S_stable |
|---|---:|---:|---:|---:|---:|---:|
| Peak | 631 | 0.07468102608 | 0.92803946 | 650 | 3.080632180 | 0.04593021 |
| DBTT | 648 | 0.1035753339 | 0.90160810 | 660 | 1.337240038 | 0.26256935 |
| weak-T | 310 | 0.08488388614 | 0.91861894 | 329 | 7.645746665 | 0.0004780732 |
| ceramic | 332 | 0.07439511296 | 0.92830484 | 342 | 0.9250680904 | 0.39650442 |

At all 32 descent checkpoints, the deterministic marked-renewal quadrature
gave attempt-admission probability one to numerical precision.  The evaluated
32-phase by 16-Xi envelope rejected no points.  Admitted lengths span
`2.2973400956248734e-6` to `1.8378720764998988e-5 m`.

`300K_endurance_stress_descent.csv` records at every checkpoint:

- cleavage action and both survival quantities;
- conditional admission probability;
- rejected phase/Xi fraction;
- minimum and maximum admitted length;
- minimum exact trial energy residual in J/m.

The trial-residual minimum may be negative even when the final proposal is
admitted because it includes all exact subgrid/topology trial lengths.  It is
not substituted for the gate's mesh-resolved committed-length decision.

## Practical 1e14 survival stresses

Bracketed secant interpolation is linear in stress versus log stable-birth
action.  These are practical finite-horizon estimates, not strict endurance
limits.

| class | S=0.9 MPa | S=0.5 MPa | S=0.1 MPa |
|---|---:|---:|---:|
| Peak | 632.7580 | 642.3807 | 648.5131 |
| DBTT | 648.0802 | 656.9174 | 662.5895 |
| weak-T | 310.9123 | 318.8653 | 323.9335 |
| ceramic | 333.3807 | 340.8549 | 346.0884 |

The `S=0.1` DBTT and ceramic estimates use their completed higher-stress v3
anchors as the upper bracket; all other estimates use the two new descent
conditions.

## Strict classification

The cleavage-attempt classification remains `no_endurance` for all classes
from the previously accepted positive intrinsic renewal floor.

Energy-admissible stable birth remains `undetermined` for all classes.  The
finite envelopes are open, but an admitted transient is not a positive
asymptotic lower bound on accepted-event measure.  Conversely, no permanently
closed gate was found or proven.  No infinite-life classification is inferred
from the `1e14` results.
