# Corrected FEM-v3 / PD-m1 local-tip reference audit (v2)

This is a read-only reconstruction from completed checkpoints. It preserves the
original `fem_pd_local_tip_reference` audit unchanged and supersedes only its
stress-concentration labels and geometry-transfer explanation. No FEM or PD
mechanics trajectory was run, resumed, or modified.

## Stress-concentration definitions

Every row in `FEM_PD_LOCAL_TIP_BASELINE.csv` now reports three distinct factors:

- `Kt_root_opening`: maximum root opening stress divided by remote `sigma_max`;
- `Kt_root_principal`: maximum principal stress at the exact clock-driving root
  node divided by remote `sigma_max`;
- `Kt_hotspot_principal`: maximum principal stress anywhere in the accepted FEM
  field divided by remote `sigma_max`.

The legacy `FEM_Kt` column is retained for compatibility and is an alias of
`Kt_root_principal`; it must not be interpreted as the field-hotspot value.

For the blunt PD mesh the clock-driving root node is node 11 at
`(149.9162346, -10.0244931) um`. The principal hotspot is node 124 at
`(160.2, 0) um`, 14.3613 um from the root clock node (the associated Gauss-point
hotspot is in element 4659). Thus 1.707 is the exact-root principal factor and
2.034 is the nearby field-hotspot factor. The root opening factor is 1.683.

For K360, the clock node and nodal principal hotspot are both node 23. Its root
opening factor is 5.56548. The relevant opening-stress geometry transfer is

`(7804 / 1402) / (7481 / 4444) = 5.565 / 1.683 = 3.307`,

not the former generic 3.260 Kt ratio. This explains the near-exact nominal
mapping `2000 / 605 = 3.306`.

## Interpretation and acceptance gate

The approximately 9% PD-2000 extrapolation agreement remains a geometry-transfer
consistency check, not an independent constitutive validation. PD-1500 remains
outside the qualified K360 local-stress range. Constitutive acceptance is the
side-effect-free same-tensor test:

`same local tensor -> same signed-MPZ evolution -> same raw m=1 action`.

It does not require equal nominal stress, equal canonical m=3 birth time, or
equal PD front-capture time. Terminal PD states are post-seed MPZ/topology
diagnostics and are not valid pre-attempt hazard-comparison states.

## Files and provenance

- `FEM_PD_LOCAL_TIP_BASELINE.csv`: direct atomic FEM checkpoints and
  deterministic reconstruction from hash-verified PD checkpoints.
- `FEM_PD_PHASE_TENSORS.npz`: reconstructed phase tensors used by the comparison.
- `FEM_PD_CONSTITUTIVE_PARITY.json`: initial same-tensor parity qualification.
- `00_Kt_definitions.png`: explicitly separates the three Kt definitions.
- `01_effective_stress_vs_nominal.png`: labels the root-opening transfer basis.
- `source_inventory.json`: hashes of immutable source tables and generated files.

The event-boundary capsule extraction and expanded tensor/capsule parity matrix
are separately generated qualifications and use pre-attempt and exactly localized
attempt capsules rather than terminal states.
