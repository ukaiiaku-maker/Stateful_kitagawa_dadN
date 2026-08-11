# Peak shared-root m1 numerical/probabilistic qualification

Status: numerical production gate qualified at the existing Peak 1500/2000 MPa
controls. This document does not select or report a new physical stress.

The first S-N integration policy is now fixed at 0.6103515625 cycles for
pre-attempt localization and reversible-embryo transitions, and 625 cycles for
post-stable topology. The measured residual attempt-coordinate uncertainty is
approximately 0.027%; representative front-coordinate discretization is
approximately 0.1--0.15%. A 156.25-cycle morphology check remains required at
one lower-stress production condition before publication.

The high-cycle action guard treats canonical phase rates as rates per second and
multiplies the requested guard cycles by the explicit seconds per cycle. Its
upper bound also includes the measured state and log-hazard validation
residuals. Tests at 20, 1000, and 2500 Hz prevent a hidden unit cancellation.

## Corrected local-stress reference

The immutable first audit remains in `fem_pd_local_tip_reference`. The corrected
version is `fem_pd_local_tip_reference_v2` and distinguishes:

- `Kt_root_opening`: clock-root opening maximum / remote maximum;
- `Kt_root_principal`: principal maximum at the exact clock node / remote maximum;
- `Kt_hotspot_principal`: field-wide principal maximum / remote maximum.

For blunt PD, node 11 is the clock root: `Kt_root_opening=1.68315` and
`Kt_root_principal=1.70697`. Node 124 is the nearby principal hotspot,
14.3613 um away, with `Kt_hotspot_principal=2.03369`. K360 has
`Kt_root_opening=5.56548`; the opening-transfer ratio is 3.307. The prior 9%
PD-2000 extrapolation is only a geometry-transfer consistency check. Exact
same-tensor parity is the constitutive validation.

The 8x8 saved matrix cross-applies the K360 631/650/735.921 MPa capsules at
1e8/1e14 and PD 1500/2000 pre-attempt capsules. All 64 combinations pass:
raw m1 action is identical in log-domain evaluation, maximum m3 relative action
error is 1.42e-14, and signed-state error is zero.

## Event boundaries

`PD_EVENT_BOUNDARY_STATES.csv/json` extracts initial, last pre-attempt, exact
attempt, stable seed, and terminal/front-capture boundaries without rerunning.
Only pre-attempt and exact attempt capsules are hazard-comparison states.
Stable/terminal states are post-seed MPZ/topology diagnostics.

## Transition and macro integration

The internal ceiling is now explicit and can be rate-separated. Peak 2000 MPa
transition-only refinement used 39.0625, 9.765625, 2.44140625, and
0.6103515625 cycles. The finest-pair attempt difference is 0.153664 cycles
(2.716e-4 relative). Stable-minus-attempt delay is exactly partition invariant
over the matrix. Production therefore uses 0.6103515625 cycles while the
reversible embryo transition clock is active, and a separately reported
post-stable topology ceiling. The residual sub-cycle attempt uncertainty is
reported rather than hidden.

A persisted restart from the exact stable boundary reproduces attempt, stable,
softening, root connection, front capture, selected mark, global action, and
global threshold endpoints. Continuous reconstructed arrays differ only at
floating solve/encoding tolerance; the largest recorded difference is
1.79e-10 in a dormant legacy log-birth ledger.

## Shared-root high-cycle path

The high-cycle active inventory now includes FEM/PD smooth state and the complete
authoritative signed-MPZ constitutive arrays/scalars. Protected state includes
global action/threshold, hazard RNG, mark RNG, attempt identity/count, PD
transition thresholds/actions/outcomes, and topology. Candidate-local PD fields
do not enter the conserved root hazard.

At Peak 1500 MPa:

- 64-cycle overlap: 21.3 cycles/exact map, H error 9.43e-13 relative;
- 100,000-cycle overlap: 16,667 cycles/exact map, H error 1.30e-10;
- 1e8-cycle overlap: 8.33e6 cycles/exact map, H error 2.97e-10;
- restart from 1e8 to 1.1e8: H error 1.73e-13 and exact threshold identity;
- accelerated event guard localized the known attempt within 0.00528 cycle.

Acceleration is claimed only when `campaign_accepted_projected_cycles>0`.
At the strongly evolving early 2000 MPa transient the controller rejected
projection and direct integration was retained.

## Probability semantics

`analytical_attempt_quantiles.csv` reports global first-attempt N10/N50/N90
from `S_no_attempt=exp(-H_attempt)`. These are not stable-front quantiles.
Stabilization/healing and front-capture probabilities remain conditional mark /
transition quantities and are reported separately from the canonical attempt
survival.

The conditional 1500 MPa ensemble now has five mark/transition seeds. All five
stabilized; three captured a front within their individual observation windows.
Seed 46 is a rate-separated result and was right-censored after exactly 1e6
post-stable cycles without root connection/front capture. Because censor windows
are unequal, the observed 3/5 capture fraction is a compact diagnostic, not a
converged stable-front survival estimate.
