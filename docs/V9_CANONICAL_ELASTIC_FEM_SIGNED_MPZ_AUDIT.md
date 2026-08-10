# v9 Canonical Elastic-FEM + Signed-MPZ Audit

## Qualified model

The canonical four-class branch is now separately versioned as:

`v9_canonical_four_class_elastic_fem_signed_mpz_energy_gated_birth_v3`

Its generated contract asserts:

```text
continuum_bulk_role = elastic_fem_only
bulk_state_evolves_in_fem = false
bulk_hardening_law = none_bulk_elastic
moving_crack_tip_mpz_active = true
tip_mpz_pt_model_active = true
bulk_pt_model_active = false
bulk_scalar_rho_used_for_signed_shielding = false
```

The v9 geometry, mesh, elastic tensor, boundary conditions, waveform, R ratio,
frequency, and phase-resolved root tensor extraction are retained.  The full
elastic phase response is solved once and cached.  Structurally required
`ep_gp`, `rho_gp`, and accumulated-plastic-strain arrays are immutable
placeholders.  `ep_gp` and accumulated plastic strain remain identically zero;
`rho_gp` remains at its reference value and enters no constitutive law.  There
is no active FEM rho-cap boundary in v3.

All evolving pre-birth plastic state is the directly loaded audited v10.2.21
signed MPZ: aggregate persistent emission, signed mobile/retained/wake
populations, Peierls/Taylor kinetics, Taylor backstress, signed shielding,
blunting, escape, and advection.

## Parity and invariance

At the identical initial elastic state, v3 and the previous FEM scaffold give
bitwise-identical phase-resolved root tensors.  Peak and DBTT prescribed-
history comparisons use identical block partitions and reproduce the signed
MPZ arrays, cumulative hazard, and log cumulative hazard exactly.  Existing
all-four direct-audited-closure tests separately compare emission, Peierls,
Taylor, barrier rates, signed state, backstress, shielding, and blunting.

The zero-`ep_gp` v3 and zero-`ep_gp` control produce the same threshold-scaled
energy proposal, admitted length, and arrest reason.  Eleven-way partition
tests meet the fixed `5e-4` MPZ/hazard tolerance.  Atomic restart preserves the
complete signed state, threshold/RNG state, cached elastic request contract,
and subsequent trajectory.

## Finite-life event check

Peak, 300 K, `804.139 MPa`, seed 1720 stopped at
`N = 255.065512355417` on its first cleavage attempt.  The fixed-opening gate
admitted `2.2973400956248734e-6 m`; no post-birth state was advanced.

## Low-stress control comparison

At Peak `735.9214951373579 MPa` and 300 K:

- noncanonical evolving-bulk control: `H(1e8) = 0.051523913898147686`;
- canonical elastic-FEM v3: `H(1e8) = 8.829782039936095`.

The v3/control hazard ratio is `171.37250204615242`.  This is an ablation
result, not a fitted correction: removing noncanonical bulk plastic relaxation
leaves the fixed elastic root waveform driving the authoritative tip MPZ.

## Four-class 300 K horizons

Each class was advanced sequentially with exact atomic restarts.  The values
below are deterministic cleavage action and energy-gated marked-renewal
survival on the evaluated phase/Xi grid.

| class | stress MPa | H(1e8) | S_stable(1e8) | H(1e10) | S_stable(1e10) | H(1e12) | H(1e14) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Peak | 735.9214951373579 | 8.829782039936095 | 1.463101257737059e-4 | 882.9782039804766 | 0 | 88297.82026662635 | 8829780.713400843 |
| DBTT | 700 | 0.005913441733257106 | 0.9941040082498795 | 0.5913441733257107 | 0.5535826733401542 | 59.13441733255116 | 5913.441733070966 |
| weak-T | 350 | 0.0008723891221042562 | 0.9991279912986529 | 0.08601797107692588 | 0.917577741211807 | 8.500396125468388 | 848.0888483407676 |
| ceramic | 400 | 0.3843090994510605 | 0.6809209223644433 | 38.43090994510604 | 2.040177136680802e-17 | 3843.0909945106005 | 384309.09945106093 |

All evaluated phase/Xi marks were energy-admissible at these checkpoints, so
the tabulated stable-birth survival equals no-attempt survival to quadrature
precision.  No rejection was hidden or reclassified.

## Asymptotic classification

The separate cleavage-attempt classification remains `no_endurance` for all
four rows from the accepted strictly positive intrinsic renewal floor.

The signed MPZ pointwise convergence metric has not certified a stationary
orbit, despite the fixed elastic root waveform and nearly constant late-cycle
hazard.  Therefore the strict energy-admissible stable-birth classification is
retained as `undetermined` for all four classes.  Finite admission through
`1e14` is not used as an asymptotic proof.

## Prior results

All earlier outputs are preserved.  Rows that evolved continuum `ep_gp` and
`rho_gp` are labeled `noncanonical_bulk_plasticity_control` and are excluded
from canonical v3 curves.
