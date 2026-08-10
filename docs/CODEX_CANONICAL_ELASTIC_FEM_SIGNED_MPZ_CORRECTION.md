# Correction: canonical four-class S-N branch must use elastic FEM bulk + signed MPZ only

## Governing correction

The Peak quiet-tail calculation stopped because the inherited v9/K360 FEM state drives a scalar dislocation density `rho_gp` toward the frozen `rho_cap = 1e17 m^-2` boundary near 4.06e9 cycles.

Do **not** resolve this by raising the cap, clipping at the cap, or inventing a saturation/recovery law.

The authoritative canonical four-class 2-D PF/FEM-CZM parity closure explicitly defines:

```text
continuum_bulk_role = elastic_fem_only
bulk_state_evolves_in_fem = false
bulk_hardening_law = none_bulk_elastic
moving_crack_tip_mpz_active = true
tip_mpz_pt_model_active = true
bulk_pt_model_active = false
bulk_scalar_rho_used_for_signed_shielding = false
```

Therefore the evolving `ep_gp/rho_gp` continuum plasticity inherited from the K360 Stateful-PD mechanics scaffold is **not part of the canonical four-class pre-birth constitutive model**.

All evolving pre-birth plasticity for the canonical four-class S-N branch must come from the audited signed persistent moving-tip MPZ:

- aggregate persistent emission;
- signed mobile populations;
- signed retained populations;
- Peierls transport;
- Taylor transport/correlation;
- Taylor backstress;
- signed retained shielding;
- accumulated slip/blunting;
- wake/escape/advection state.

The FEM bulk supplies elastic equilibrium and the resolved root stress/opening tensor only.

## Consequence for the rho-cap blocker

The `rho_cap` boundary is not a new physics decision for the canonical four-class branch. It is evidence that a noncanonical inherited bulk-plasticity state remained active.

The correct production repair is to remove that state evolution from the canonical branch, not to define behavior beyond the cap.

Preserve the current rho-cap trajectory as a diagnostic named, for example:

```text
noncanonical_K360_bulk_plasticity_control
```

Do not delete or overwrite it.

## Required implementation

Create a separately versioned canonical mechanics mode, for example:

```text
v9_canonical_four_class_elastic_fem_signed_mpz_energy_gated_birth_v3
```

In this mode:

1. FEM geometry, mesh, boundary conditions, elastic constitutive tensor, loading waveform, R ratio, frequency, and root stress extraction may reuse the qualified v9/K360 mechanics infrastructure.
2. The continuum FEM plastic constitutive update is disabled.
3. `ep_gp` must not accumulate plastic strain.
4. `rho_gp` must not evolve or contribute a bulk hardening/backstress law.
5. If placeholder arrays are required structurally by the assembler, keep them immutable at an explicitly documented neutral/reference value and prove that changing/removing the placeholder does not alter the elastic solution.
6. No FEM density cap may control accepted physical time in this canonical mode.
7. All plasticity-induced shielding/backstress/blunting must be derived only from the signed MPZ state.
8. The energy-admissible stable-birth endpoint remains the v2 definition: first canonical cleavage attempt whose existing post-first-passage v10.2.30 fixed-opening gate admits a nonzero mesh-resolved crack event.
9. Stable birth remains terminal; no post-birth crack propagation is performed.

## Parity proof before new production

Before rebuilding the S-N/endurance campaign, require a direct parity audit against the authoritative four-class 2-D tip-only source.

For at least Peak and one control class under identical prescribed loading, compare:

```text
FEM bulk role = elastic only
root-local stress/opening history
signed channel stresses
signed mobile/retained/wake state
Peierls/Taylor rates
Taylor backstress
signed shielding
blunting
cleavage log rate
cumulative cleavage action
energy-gate proposal/admitted length
```

The audit must explicitly assert:

```text
bulk_state_evolves_in_fem = false
continuum_bulk_role = elastic_fem_only
bulk_hardening_law = none_bulk_elastic
bulk_scalar_rho_used_for_signed_shielding = false
```

Then re-run >10x block-partition and atomic restart equivalence for the corrected canonical mode.

## Status of existing canonical S-N results

The 300 K Peak/DBTT/weak-T/ceramic results generated with evolving v9/K360 continuum `ep_gp/rho_gp` must not remain labeled as final canonical four-class production results.

Preserve them intact, but relabel them as:

```text
noncanonical_bulk_plasticity_controls
```

because their long-N stress history and hazard were affected by a constitutive state that is absent from the intended canonical four-class parity model.

Do not discard them; they are useful ablations showing the effect of adding bulk plasticity.

The final canonical 300 K S-N/endurance curves must be regenerated under elastic FEM + signed MPZ only.

## Efficiency implication

This correction should substantially simplify the VHCF problem.

Before stable birth and with fixed notch geometry:

- the elastic FEM phase response is deterministic and cycle-repeatable for a given applied amplitude;
- the only evolving mechanical material state is the signed MPZ;
- the FEM root phase tensor can therefore be precomputed/cached exactly, subject to any explicit effective-opening/blunting feedback required by the audited source;
- long-N acceleration should evolve the signed MPZ and cleavage/energy-gate state rather than repeatedly integrating a noncanonical continuum `rho_gp` field.

Use this structure to target 1e8, 1e10, 1e12 and 1e14 cycles without brute-force FEM cycling.

## 300 K production sequence after correction

Stay at 300 K.

Use only the final canonical rows:

```text
Peak     = v913_paper_peak01_0242980_persistent_sites
DBTT     = v913_paper_dbtt01_0202500_persistent_sites
weak-T   = v913_paper_weakT01_0129902_persistent_sites
ceramic  = v913_paper_ceramic01_0077080_persistent_sites
```

Proceed one condition at a time:

1. implement elastic-FEM/signed-MPZ v3 mode;
2. prove source parity, partition invariance, and restart equivalence;
3. run one Peak 300 K finite-life point to verify stable-birth event equivalence;
4. run one Peak low-stress VHCF tail to 1e8 and qualify the long-N accelerator;
5. extend Peak toward 1e10 -> 1e12 -> 1e14 where valid;
6. repeat sequentially for DBTT, weak-T, ceramic;
7. construct stable-birth survival using the marked renewal + post-passage energy gate;
8. classify strict endurance only from the corrected canonical model.

## Endurance semantics remain

Keep separate:

```text
cleavage_attempt_asymptotic_classification
energy_admissible_stable_crack_birth_asymptotic_classification
```

The positive intrinsic cleavage-attempt floor can still imply `cleavage_attempt = no_endurance`.

Stable-crack-birth endurance depends on the energy-admissibility marked-renewal process in the corrected elastic-FEM/signed-MPZ state.

## Stop conditions

Stop only for a genuine audited-source ambiguity, inability to reproduce the canonical tip-only source with elastic FEM bulk, or failure of partition/restart equivalence that cannot be repaired without changing physics.

Do not treat the inherited K360 `rho_cap` as a physical saturation decision for this branch.
