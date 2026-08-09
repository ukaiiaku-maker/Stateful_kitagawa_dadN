# Current state

## Active architecture

The current consolidated finite-feature solver is the v8.7 stateful FEM/peridynamics implementation:

```text
SN_2D_intact_FEM_stateful_local_peridynamics_v8_7_generalized_features_local_front_spacing_fixed_geometry
```

The active source on the local validated installation has been organized around:

```text
arrhenius_fracture/sn_pd2d_stateful_v8_7_generalized_features.py
arrhenius_fracture/stateful_peridynamics_v8_7_local_front_spacing.py
arrhenius_fracture/sn_feature_geometry_v8_7.py
```

v8.7 extends the validated v8.3 persistent-first-passage solver with generalized ellipse / rounded-V geometry and a local-front-spacing treatment.

## Physical architecture

- Cyclic global response: intact 2-D FEM.
- Local fracture process: stateful nonlocal/peridynamic patch.
- Candidate crack-initiation sites: persistent exponential first-passage thresholds.
- Local state includes plasticity/dislocation evolution, residual stress, candidate-site hazard, embryos, stable defects, bond damage, active-front state, and morphology metrics.
- Geometry is fixed for the production baseline; crack formation is represented by the evolving local PD network.
- PD degradation is presently one-way coupled: it redistributes the local patch response but does not change the global FEM stiffness.
- A run is a valid finite-feature failure only at a resolved `physical_handoff`; a runout-reaching no-handoff case is `right_censored`.

## Verified K360 gate

The v2.8 / v2.8.1 gate established a real end-to-end execution path for the generalized solver.

### Failure endpoint

```text
sigma_a = 735.9214951373579 MPa
status = physical_handoff
N_handoff = 3,695,962.9621467358 cycles
```

Representative event sequence:

```text
first embryo       ~1.8313e6 cycles
first stable       ~1.8313e6 cycles
first softening    ~1.9166e6 cycles
root connection    ~2.4486e6 cycles
front capture      ~2.5348e6 cycles
physical handoff   ~3.6960e6 cycles
```

Final morphology was localized and admissible: axial coverage 1.0, orientation coherence ~0.978, width/length ~0.24, 79 connected broken bonds, and no off-front broken fraction.

### Lower-stress runout endpoint

```text
sigma_a = 690.4432004940379 MPa
status = right_censored
N_end = 1.0e8 cycles
```

No embryo or broken bond was realized. The run nevertheless accumulated substantial internal-state evolution and the instantaneous birth rate fell strongly with cycling.

## What this proves

It proves that the current solver can execute both a real failure and a real runout case with valid geometry and summary/reuse behavior.

It does **not** prove that the model has a true endurance limit. The lower-stress hazard at `1e8` cycles is small but nonzero. If it approaches a positive constant, or decays too slowly, cumulative hazard still diverges and eventual failure probability remains one.
