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
- The first 2500 MPa global attempt localized at
  `N=19949.04475728977`; one marked embryo was injected and stabilized at
  `N=19951.82123814441`.  A second attempt was present by approximately
  `N=5.42e4`.  No legacy per-site K2 clock created either embryo.
- The 2500 MPa continuation is checkpointed but remains censored before front
  capture.  It must not be reported as a front-capture lifetime.

## Physical Peak bracket

| sigma_a (MPa) | actual checkpoint N | H_attempt | S_no_attempt | tail dH/dN | spatial status |
|---:|---:|---:|---:|---:|---|
| 1000 | 20 | 1.86468e-55 | ~1 | 9.32339e-57 | no attempt; diagnostic smoke |
| 2000 | 1e14 | 3.39552e-4 | 0.999660506 | 3.39552e-18 | no attempt/front; censored |
| 2150 | 1e10 | 3.21132e-4 | 0.999678919 | 3.21132e-14 | no attempt/front; censored |
| 2150 | 1e12 | 3.21132e-2 | 0.968396914 | 3.21132e-14 | no attempt/front; censored |
| 2500 | 5.7e4 checkpoint | ~1.1e-1 | attempt renewal active | ~1.87295e-6 | two attempts, stable sites, no front yet |

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
The available 2500 MPa trajectory has stable sites but no captured front, so
the stable-front endurance question remains open pending completion and a
mark/transition-seed ensemble.

## Remaining production work

1. Finish the checkpointed 2500 MPa trajectory to first front capture or a
   defined censoring boundary.
2. Complete exact restart identity after the expanded authoritative MPZ
   capsule qualification and retain a machine-readable comparison.
3. Qualify shared-root high-cycle projection with protected global-clock and
   mark state before using it for event-producing `1e12--1e14` trajectories.
4. Estimate stable-front survival across mark/transition seeds at adaptively
   selected stresses.  Do not infer it from no-attempt survival.

