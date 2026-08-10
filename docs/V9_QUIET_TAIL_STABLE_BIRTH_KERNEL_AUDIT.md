# v9 Quiet-Tail Stable-Birth Kernel Audit

## Scope and endpoint

This milestone implements the 300 K, pre-birth quiet-tail path required by
`CODEX_QUIET_TAIL_STABLE_BIRTH_KERNEL_PLAN.md`.  It does not gate cleavage
hazard, propagate a crack, change any material row, or change the v10.2.30
fixed-opening energy gate.  Cleavage attempts and energy-admissible stable
birth remain distinct endpoints.

The work stopped during Peak, before DBTT, weak-T, or ceramic production, for
the fail-closed constitutive boundary described below.

## Implemented kernel and renewal machinery

`v9_quiet_tail_kernel.py` records a complete phase-resolved cycle: FEM root
stress and opening, signed channel stress, signed mobile/retained/wake state,
Taylor backstress, signed shielding, blunting, cleavage log rate, and
within-cycle cleavage action.  Checkpoints use hash-verified atomic array
generations and can restore an explicitly named generation.

The exact v10.2.30 energy gate is evaluated over phase and the prescribed Xi
grid.  Repeated topology solves are cached by phase and topology without
changing the returned result.  The deterministic marked-renewal operator keeps
the same Xi for attempt timing and the threshold-scaled proposal.  Its
very-small-cycle-hazard limit uses uniform cumulative-action phase with a
reported total-variation bound of `h_cycle/2`.

The production block controller now persists its proposed next block across
calls.  This removes repeated re-subdivision while retaining rejected-block
rollback and restart state.

## Peak 300 K direct tail

Stress: `735.9214951373579 MPa`; option:
`v913_paper_peak01_0242980_persistent_sites`.

| N | H_cleave | per-cycle hazard | stationary vs prior checkpoint |
|---:|---:|---:|:---|
| 1e6 | 4.129489053787163e-2 | 1.563617181695305e-8 | no |
| 1e7 | 5.152281158061033e-2 | 4.418690911433176e-13 | no |
| 1e8 | 5.152334806131571e-2 | 1.392886106600335e-26 | no |

At 1e8 the 32-phase by 16-Xi map contains 512 evaluated gate points; every
sampled point has nonzero admitted length.  This finite sampled result is not
used as an asymptotic no-endurance proof.  The piecewise checkpoint estimate is
`S_stable_birth(1e8) = 0.9497814741472245`; because the tail is not stationary
and the full continuous mark domain is not yet certified, it remains a
diagnostic rather than a production 1e14 result.

The original checkpoint label formatted multiple nearby cycle counts as
`N_1e+08`.  The exact 1e8 generation
`generation_9dff58c4beb5493c8bb3cd3d69d6f284` was recovered by hash and copied
to the unique atomic path `N_100000000`.  New writers use 17-digit labels.

## Acceleration qualification

A frozen periodic kernel failed: from 1e8 to 1.01e8, root-tensor relative
error was `8.8568e-2`, phase-log-rate error `2.9215e-1`, and cycle-hazard
relative error `1.5728e-1`.  An affine macro-cycle DMD prototype also failed
the phase/root/hazard acceptance bounds and is retained only as negative
evidence.

A one-way high-order DOP853 FEM integration in log-cycle coordinates passed
direct production-FEM overlaps:

| overlap | signed-state rel. error | root rel. error | phase log-rate abs. error | hazard rel. error |
|---:|---:|---:|---:|---:|
| 1e6 | 1.19e-8 | 1.748e-5 | 1.84e-8 | 1.10e-8 |
| 1e7 | 1.083e-5 | 7.02e-5 | 2.64e-6 | 1.59e-6 |

All are below the fixed `1e-3` qualification limits.  This validates the
accelerator over the overlapping windows; it does not authorize crossing an
existing constitutive boundary.

## Genuine blocker: FEM density cap before 1e10

At the active 1.10e8 diagnostic state, the maximum FEM density is
`3.814872993990296e15 m^-2`, its current positive increment is
`2.43506515e7 m^-2/cycle`, and the frozen production cap is `1e17 m^-2`.
The local boundary estimate is approximately `3.950002200e9` additional
cycles, placing the boundary near `4.06e9`, below the requested 1e10 horizon.

The production transaction explicitly defines this cap as an integration
boundary, not an accepted-state clip.  The constitutive RHS used to calculate
rates clips density internally, so an unconstrained adaptive solver stalls in
vanishing steps as it approaches the cap.  The accelerator now detects this
boundary fail-closed and must not extrapolate or silently saturate density.

Continuing requires a physics decision for the FEM density state: for example,
whether `rho_cap` is a physical saturation law (and, if so, its exact evolution
at saturation) or merely a numerical validity limit requiring a newly audited
range.  Raising the cap, clipping the accepted state, or inventing a recovery
law would alter frozen physics and was not done.

Consequently the requested 1e10, 1e12, and 1e14 stable-birth survivals and the
remaining class sequence are not scientifically reportable in this milestone.
Peak stable-birth endurance remains `undetermined`; the separate cleavage-
attempt classification remains `no_endurance` from the previously accepted
attempt-floor proof.
