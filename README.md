# Stateful Kitagawa / da/dN

This repository is the consolidated working tree for the stateful fatigue model developed for blunt-defect S-N / Kitagawa calculations and eventual crack-growth (`da/dN`) continuation.

## Immediate scientific objective

Produce statistically meaningful S-N calculations spanning roughly `10^3` through `10^12` cycles and distinguish:

- a system with a genuine fatigue-endurance regime, for which the infinite-cycle failure probability remains below one; and
- a system with no true endurance limit, for which the cumulative failure hazard diverges and eventual failure probability tends to one even when the finite-cycle rate becomes extremely small.

A finite `right_censored` run at `10^8` cycles is **not** by itself evidence for a fatigue limit.

## Current validated baseline

The current finite-feature architecture is the stateful FEM + local peridynamics solver developed through v8.7. The July gated K360 pilot established a real failure endpoint (`physical_handoff` at about `3.696e6` cycles at `sigma_a = 735.921495 MPa`) and a lower-stress `right_censored` endpoint at `1e8` cycles (`sigma_a = 690.443200 MPa`). These are regression anchors, not a completed S-N calibration.

Read `docs/CODEX_HANDOFF.md` before making physics changes.

## Development rule

Do not launch large sweeps until the validation ladder in `docs/DEFINITION_OF_DONE.md` passes. Previous campaign wrappers failed in ways that consumed time without producing physical results; fail-closed orchestration, restart equivalence, and asymptotic-hazard diagnostics are now core requirements.
