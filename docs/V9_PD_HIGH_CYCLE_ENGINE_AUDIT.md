# V9 Stateful-PD high-cycle engine audit

> **Historical/superseded audit.** This document records the synthetic-engine
> boundary at its original milestone and is preserved without rewriting those
> claims. Real callback extraction, physical qualification, subsequent repairs,
> and current production results are governed by
> `docs/V9_PD_PRODUCTION_PROGRESS.md`.

## Scope

`v9_stateful_pd_dormant_high_cycle_event_engine_v1` is a separately versioned,
event-to-event engine for the dormant, fixed-topology interval.  It is
fail-closed when any embryo, stable site, bond damage, primary seed, captured
front, or topology change exists.  It does not change the frozen v8.7 spatial
physics.

The explicit inventory is implemented in `state_inventory()`.  Only FEM
`ep_gp`, `rho_gp`, `epsp_acc_gp`, `u`, and PD `log_delivery_memory` are active
coordinates.  Cumulative ledgers, all site birth/transition actions and
thresholds, both RNG states, site status, bond/front state, mesh, and bond
connectivity are protected rather than projected.

## Numerical contract

Every private cycle, periodic iteration, and midpoint/endpoint projective
validation compares SHA-256 signatures of ledgers, stochastic state, and
topology and verifies that physical cycle count is unchanged.  The projective
path uses independent exact maps at midpoint and endpoint, Simpson action
quadrature, a start/mid/end upper action guard, physical-bound checks, and
reject/halve semantics.  Stationary propagation requires both a converged map
and strict current-to-fixed-point admission distance.  It advances only whole
cycles before an exact event guard.

Atomic controller/mode sidecars record physical cycle, topology signature,
map counts, accepted projected cycles, and mode history.  Complete physical
restart state remains in the existing hash-verified v9 generation checkpoint.

## Qualification

`tests/test_v9_pd_high_cycle.py` covers private-state invariance, stationary
right-censors at 1e12 and 1e14, first-passage guarding near 1e12, partition
equivalence, projective linear drift, topology invalidation, and stationary /
projective restart equivalence.  The complete regression result is 168 passed
plus 7 subtests.

`scripts/qualify_v9_pd_high_cycle.py` writes numerical comparison and mode
history files plus a compact map-count plot.  These are synthetic numerical
engine qualifications, not material predictions.

## Current production boundary

The spatial adapter enforces the real PD state separation, but the monolithic
FEM/PD driver has not yet supplied its exact one-cycle evaluator callback.
Consequently the previously completed K360 1e8 trajectory remains a real
transactional-driver result, not a high-cycle-engine result, and no claim of a
real 1e12/1e14 PD trajectory is made here.  Production use must remain disabled
until the exact driver phase ordering is extracted into that callback and is
qualified against the historical finite-life/front trajectory.  This is the
remaining implementation boundary, not a physics or provenance ambiguity.
