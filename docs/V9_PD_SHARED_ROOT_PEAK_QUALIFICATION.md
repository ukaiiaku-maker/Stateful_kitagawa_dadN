# v9 shared-root marked-cleavage Peak qualification

Model: `v9_four_class_stateful_PD_shared_root_MPZ_marked_cleavage_v1`

This is an active production qualification, not a completed four-class S-N
campaign.  All conditions are Peak, 300 K, blunt half ellipse (`a=150 um`,
`b=300 um`, nominal `rho=600 um`) with `pd_image_policy=none`.

## Coupling evidence

- The global canonical action is computed only from the exact root FEM phase
  tensor and one authoritative signed MPZ.
- Local PD cleavage propensity and initiation fields enter only the normalized
  spatial mark.
- Doubling candidate density changed the realized population from 3,804 to
  7,468 sites but left the evolved five-cycle action, action increment, and
  threshold bit-identical (`H=4.661695343009272e-56`).
- The ledger identity `H_new=H_old+dH` closes to reported absolute error zero
  in the completed large-N transactions.
- The earlier `f0b5a71` attempt coordinates and second attempt are retained as
  provisional controls only.  They used scaled-average localization and the
  pre-policy repeated-attempt state machine.
- The corrected ordered run accepts 19,949 complete cycles, then localizes the
  first phase crossing at `N=19949.00251302945` (phase 0, within-phase fraction
  `0.0804169424`).  Site 801 / PD node 124 is marked from that exact phase.
- The birth row has one realized embryo, zero smooth stable exposure and zero
  realized stable sites.  The clock is paused through the reversible embryo,
  which stabilizes at `N=19951.77899388409`; global action remains fixed at
  `0.037363528078493934` and no second attempt is permitted.
- The same corrected trajectory first softens at `N=33811.82370133839`, first
  connects to the root at `N=536759.639157665`, and reaches first front capture
  at `N=586759.639157665`.  It stops there without post-capture propagation;
  physical handoff is false at the endpoint.

## Physical Peak bracket

| sigma_a (MPa) | actual checkpoint N | H_attempt | S_no_attempt | tail dH/dN | spatial status |
|---:|---:|---:|---:|---:|---|
| 1000 | 20 | 1.86468e-55 | ~1 | 9.32339e-57 | no attempt; diagnostic smoke |
| 2000 | 1e14 | 3.39552e-4 | 0.999660506 | 3.39552e-18 | no attempt/front; censored |
| 2150 | 1e10 | 3.21132e-4 | 0.999678919 | 3.21132e-14 | no attempt/front; censored |
| 2150 | 1e12 | 3.21132e-2 | 0.968396914 | 3.21132e-14 | no attempt/front; censored |
| 2500 | 5.8676e5 | 3.73635e-2 | 0.963325 | clock stopped after stabilization | one attempt; front captured |

The 2000/2150/2500 stresses are new shared-root production qualifications;
the historical 4 GPa and 12 GPa trajectories remain controls only.
Maximum root-cycle stresses are about 9.04, 9.72, and 11.30 GPa respectively,
with negligible continuum accumulated plastic strain in the reported tail
checkpoints.  Geometry resolution remained valid (`rho_mesh/h ~= 49.9`).

## Current physical interpretation

Replacing independent multiplicity-sensitive K2 site clocks with the conserved
root clock materially increases and regularizes the predicted attempt rate.
At 2500 MPa it produces canonical attempts and stable spatial seeds on a
roughly 1e4--1e5-cycle scale, where the legacy K2 bridge remains suppressed by
its extremely small delivery completion.  At 2000 MPa, however, the actual
integrated attempt action remains only `3.40e-4` at `1e14` cycles.  Thus the new
formulation predicts a very sharp low-stress life increase, not a resolved
finite-horizon cutoff.

No strict endurance limit is established.  The canonical cleavage model has a
positive intrinsic attempt floor, so absence of an attempt through `1e14` is
censoring.  More importantly, stable-front survival is not `exp(-H_attempt)`:
it also depends on marked embryo healing/stabilization and subsequent topology.
The corrected 2500 MPa seed reaches a captured front, but one seed is not a
stable-front survival estimate.  The endurance question remains open pending a
mark/transition-seed ensemble and lower-stress front statistics.

## Event-boundary and acceleration audit

- Shared-root mode now fails closed when `--pd-high-cycle` selects the legacy
  adapter.  That adapter does not protect the signed MPZ, global action and
  threshold, or the hazard/mark RNG streams.
- The completed 2000 MPa (`1e14`) and 2150 MPa (`1e10`, `1e12`) checkpoints all
  report `pd_high_cycle=false` and have no accepted high-cycle mode history;
  their finite-horizon action was produced by direct macro-transactions.
- The shared MPZ uses the mesh-resolved initial radius
  `598.7445183220907 um`; nominal `600 um` remains a separate geometry field.
- A deliberately repartitioned split at `N=19950` preserved attempt identity,
  clock/RNG state and endpoint class but produced small smooth-field differences
  downstream (largest observed bond-damage difference `2.42e-9`).  It is a
  partition-convergence result, not bitwise restart identity.  Exact persisted
  restart identity must use the same accepted partition and remains an explicit
  qualification item.

## Remaining production work

1. Complete exact restart identity after the expanded authoritative MPZ
   capsule qualification and retain a machine-readable comparison.
2. Qualify shared-root high-cycle projection with protected global-clock and
   mark state before using it for event-producing `1e12--1e14` trajectories.
3. Estimate stable-front survival across mark/transition seeds at adaptively
   selected stresses.  Do not infer it from no-attempt survival.
