# Parallel 300 K four-class Stateful-PD audit

## Scope

This milestone restores Stateful-PD as a treatment parallel to, and explicitly
separate from, the accepted
`elastic_FEM_signed_MPZ` v3 campaign. Nothing under
`runs/sn_v9_canonical_four_class_elastic_v3` is modified or reinterpreted.

The new driver/model contract is
`v9_four_class_stateful_PD_blunt_notch_stable_spatial_birth_v1`. Its primary
endpoint is the first captured spatial front after a stabilized seed. The
following coordinates remain separately reported by the inherited state:

- `cycles_first_embryo`;
- `cycles_first_stable`;
- `cycles_root_connected`;
- `cycles_front_capture`; and
- `cycles_connected` (later handoff diagnostic only).

A transient embryo and one damaged bond are not stable spatial crack birth.
The driver stops at front capture and does not intentionally continue to
physical handoff.

## Historical geometry provenance

`sn_feature_geometry_v8_7.BluntNotchGeometry` retains the original S-N
defaults:

- plate: 2.0 mm by 4.0 mm;
- half-elliptical edge notch;
- depth `a = 150 um`;
- half-height `b = 300 um`; and
- analytic root radius `b^2/a = 600 um`.

The four-class wrapper fixes those values as its geometry contract. The first
real h15 construction resolved a radius of 598.745 um and had
`rho_mesh/h_PD = 49.9`; the production-resolution preflight passed with
`delta/h = 3.75` and 25 feature-surface nodes. No intermediate radius was
introduced.

## Existing Stateful-PD physics retained

The implementation starts from `V9StatefulPDPatch`, not a new damage toy
model. It retains:

- spatial candidate sites and persistent per-site first-passage thresholds;
- log-domain finite-memory K=2 delivery/completion;
- reversible embryos and persistent stabilization/healing clocks;
- directional cohesive-bond softening and damage;
- local PD spring re-equilibration and stress redistribution;
- stable-seed selection and competition;
- root-connected topology, active-front capture, and morphology audits; and
- atomic checkpoint state including candidate/event RNG streams and the v9
  transition arrays.

Existing restart and block-layout invariance tests remain applicable to these
inherited state variables. The new adapter adds no stochastic stream.

## Four-class transfer contract

The wrapper fails closed through `select_canonical_option` and records the
registry, selection, active-row fingerprint, exact row, and executable source
hashes. The only accepted rows are:

- Peak `v913_paper_peak01_0242980_persistent_sites`;
- DBTT `v913_paper_dbtt01_0202500_persistent_sites`;
- weak-T `v913_paper_weakT01_0129902_persistent_sites`; and
- ceramic `v913_paper_ceramic01_0077080_persistent_sites`.

The audited cleavage and emission EXP-floor surfaces and their attempt
frequencies are transferred exactly. Peierls and Taylor transport surfaces are
constructed with the authoritative `TransportBarrier.as_surface(emission)`
definition. Unit tests compare the local executable EXP-floor formula and the
PD cleavage adapter at 300 K.

This does **not** make the complete PD state identical to the signed front-local
MPZ. The following remain model-distinguishing PD closure choices and are
declared as such in every contract:

- spatial candidate density and K=2 delivery;
- scalar FEM plastic/rho state in the legacy PD mechanics scaffold;
- stabilization/healing kinetics;
- bond growth/linkage and front-capture criteria; and
- absence of signed mobile/retained MPZ populations at every PD point.

No class-specific fitting was introduced. In particular, source density,
stabilization, healing, growth, and linkage were not tuned to move the PD curve
toward FEM-v3.

## Storage policy

The new production wrapper defaults to `pd_image_policy = none`. The inherited
driver now supports:

- `none`: no PD or FEM images;
- `event_only`: terminal images only when a stable/front event exists; and
- `selected`: the historical periodic plus terminal diagnostic images.

Atomic restart arrays, numerical histories, summaries, and topology audits are
still written. A completed 100,000-cycle Peak probe produced no PNG files.

## Qualification probes and blocker

Two deliberately sequential Peak probes were used to examine the physical
mapping; neither is a production S-N point.

1. 700 MPa, 100,000 cycles, initial exact-surface adapter with completed-flow
   delivery. The mesh/root audit passed and image policy `none` produced no
   images. The local FEM maximum principal stress was 3.164 GPa, but the
   completed-flow delivery maximum was `9.39e-44 s^-1`, maximum K=2 completion
   was about `3.49e-95`, and maximum cumulative birth hazard was
   `9.86e-111`. It right-censored without an embryo.

2. 1800 MPa, emission-delivery correction, advanced atomically beyond 80,000
   cycles before diagnostic termination. Even with a much higher applied
   stress, delivery memory was about `5.89e-22`, K=2 completion about
   `1.73e-43`, and no embryo or damaged bond appeared.

The second probe reveals a genuine constitutive mapping ambiguity, not a
numerical failure. In the authoritative signed-MPZ model, aggregate emission is
continuous and obeys `N = M lambda(sigma_drive - sigma_back(N),T) dt`.
The legacy PD birth law instead treats a local emission/delivery rate as a
finite-memory K=2 prerequisite for cleavage. For the audited Peak row at 300 K
and a genuinely blunt notch, this extra K=2 completion factor suppresses the
birth hazard by tens of decades even at an applied stress far above a useful
fatigue regime.

It would be scientifically improper to obtain a finite-life condition by
raising stress further, multiplying delivery by a fitted factor, changing the
site density, removing K=2 silently, or retuning the four-class row.

## Decision required before production

The code/provenance/storage/endpoint separation is ready, but a four-class PD
S-N skeleton is blocked until the relationship between audited aggregate
emission and the PD embryo prerequisite is chosen. The narrow alternatives are:

1. retain K=2 and derive its delivery measure/multiplicity from an audited
   physical area/source-count mapping;
2. define canonical cleavage attempts as embryo births in the PD patch and use
   stabilization/healing plus spatial topology as the additional PD physics;
   or
3. provide another audited delivery observable from the signed-MPZ state.

Until that decision is made, the current probes are mapping diagnostics only,
not four-class production points, and no adaptive PD stress campaign should be
launched.
