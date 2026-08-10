# Canonical v3 300 K diagnostic figure audit

## Scope and provenance

This figure package describes the completed canonical elastic-FEM/signed-MPZ
v3 campaign at source commit `5d6b1eef331ca88a4a55c8a34fc7d988f2882ac9`.
No physical condition was run or advanced. The accepted parameter rows,
geometry, 300 K temperature, energy-gated stable-birth endpoint, and atomic
checkpoint generations are unchanged.

The plotting preflight verifies all 48 rows of `300K_quiet_tail_state.csv`
against the active atomic checkpoint generation. For every row it requires
exact equality among:

- the CSV `H_cleave` and `cycle_hazard`;
- `summary.json` `H_cleave` and `cycle_hazard`;
- `state_metadata.json` birth-capsule cumulative cleavage hazard; and
- the condition key (class, stress, and cycle count).

The machine-readable result of that preflight is
`figures/v9_canonical_four_class_elastic_v3_300K/figure_data_audit.json`.

## Reconstruction correction

An actual summary-table defect was found. For weak-T at 310 and 329 MPa,
`300K_stable_birth_survival.csv` and the reconstructed descent table had mixed
the `qualified_stationary_marked_renewal_quadrature` representation into rows
whose keys denoted physical 10^10 and 10^12-cycle checkpoints. The underlying
atomic trajectories were correct and were not modified.

The affected values were:

| class / stress | N | incorrect mixed H | direct atomic H |
|---|---:|---:|---:|
| weak-T / 310 MPa | 10^10 | 8.488388614248607e-2 | 8.587739043261424e-6 |
| weak-T / 310 MPa | 10^12 | 8.488388614248607e-2 | 8.501674896487475e-4 |
| weak-T / 329 MPa | 10^10 | 7.645746665230013 | 7.744564855240640e-4 |
| weak-T / 329 MPa | 10^12 | 7.645746665230013 | 7.660309387636118e-2 |

At 10^14 cycles the direct and stationary-horizon values coincide. The
reconstruction now takes `H_cleave` only from direct checkpoint rows and
combines it with the evaluated checkpoint energy-gate envelope. Because all
sampled phase/Xi marks are admitted in these completed conditions,
`S_stable_birth = exp(-H_cleave)`. A partially open sampled gate fails closed
in this reconstruction and requires the checkpoint-kernel renewal solver;
it is not silently approximated.

The canonical-v3 runner no longer writes stationary future projections into
the physical checkpoint survival table. A regression test reproduces the
weak-T failure mode and proves that direct checkpoint H wins over a conflicting
stationary representation.

## Figure-by-figure data audit

1. `01_probabilistic_SN_300K`: direct N10/N50/N90 crossings are log-log
   interpolations between adjacent *actual checkpoint* H values in
   `300K_stable_birth_survival.csv`. Hollow markers at 10^14 cycles come from
   `300K_practical_survival_stress_estimates.csv`; they are bracketed
   finite-horizon stress interpolations, not endurance limits.

2. `02_cumulative_cleavage_hazard`: direct `H_cleave` at 10^8, 10^10, 10^12,
   and 10^14 cycles for both descent stresses of each class, read from
   `300K_endurance_stress_descent.csv` after exact cross-check against
   `300K_quiet_tail_state.csv` and the atomic checkpoint.

3. `03_instantaneous_cycle_hazard`: direct atomic/checkpoint
   `cycle_hazard` for the same conditions. No `H/N` surrogate is used.

4. `04_physical_MPZ_state_evolution`: direct arrays in each active
   `state_arrays.npz` plus direct scalar metadata for the lower descent stress
   (Peak 631, DBTT 648, weak-T 310, ceramic 332 MPa). It shows Taylor
   backstress, the magnitude and recorded negative sign of signed shielding,
   mesh-tip radius/blunting, both signed-channel mobile and retained contents,
   accumulated slip, wake content, and advance/advection. Positive and
   negative channel magnitudes overlap for these symmetric cyclic histories.
   Wake and advance remain exactly zero at all four checkpoints.

5. `05_phase_resolved_cleavage_response`: direct
   `kernel_cleavage_log_rate_effective_s`, `kernel_phase`, and
   `kernel_root_opening_Pa` arrays from those same 16 atomic checkpoints.
   The root opening is the cached elastic-FEM phase response; the cleavage
   rate includes the checkpoint signed-MPZ state.

6. `06_energy_gate_diagnostics`: checkpoint summaries reconstructed from
   `300K_energy_gate_envelope.csv` into
   `300K_endurance_stress_descent.csv`. Admission probability, rejected
   phase/Xi fraction, admitted-length range, and energy residual are derived
   reductions of the exact post-first-passage v10.2.30 gate evaluations. They
   are stable-birth admission diagnostics, not PD damage accumulation.

7. `07_initial_geometry_mesh`: the actual canonical-v3 mesh reconstructed
   from the accepted K360 run arguments through
   `CanonicalElasticFourClassFEMCondition`. The mesh-resolved root radius is
   45.538 micrometres, compared with the nominal 45 micrometre elliptical
   notch radius. The MPZ arrow uses the selected Peak row's audited MPZ length;
   the orange corridor is the actual local mesh-refinement corridor.

Every figure is committed as PNG and vector PDF. The reproducible source is
`scripts/plot_v9_300K_diagnostics.py`.

## Physical interpretation

At both new descent stresses for every class, cumulative H is nearly linear
in N over the decade checkpoints and the actual instantaneous cycle hazard is
nearly constant and positive. Thus these finite trajectories support a
practical stationary-orbit interpretation over 10^8--10^14 cycles, but they
do not by themselves prove a strict infinite-life classification.

- Peak: the 631 and 650 MPa paths separate by roughly 1.6 decades in H.
  Backstress, mobile content, and slip continue growing across the stored
  checkpoints, while their effect on phase-resolved cleavage rate is below
  plot resolution. At 10^14, H is 0.07468 at 631 MPa and 3.08063 at 650 MPa.
- DBTT: the 648 and 660 MPa paths similarly retain positive, almost constant
  cycle hazard. Signed shielding is the largest of the four classes in
  magnitude but remains mechanically tiny here. At 10^14, H is 0.10358 and
  1.33724 respectively.
- weak-T: this is the important transient-state case. From 10^8 to 10^12 at
  310 MPa, backstress rises from about 0.846 to 1.039 GPa, mobile content,
  slip, tip radius, and the phase-resolved rate evolve; by 10^12--10^14 the
  stored MPZ state is effectively unchanged. The *integrated* direct H still
  grows from 8.69e-8 to 8.59e-6 to 8.50e-4 to 8.49e-2. Repeating the final
  stationary value at earlier horizons was therefore physically wrong.
- ceramic: the signed-MPZ response is extremely weak, with zero wake and
  advance and negligible shielding, but the direct cleavage cycle hazard is
  positive and nearly constant. At 10^14, H is 0.07440 at 332 MPa and 0.92507
  at 342 MPa.

All sampled post-first-passage phase/Xi gate points remain admitted at these
stresses (conditional admission one, rejected fraction zero). That makes the
finite-horizon stable-birth survival equal to no-attempt survival for the
completed envelope, while strict endurance remains `undetermined`: the plots
do not establish an all-future lower bound on accepted-event measure or an
all-future closed energy gate.
