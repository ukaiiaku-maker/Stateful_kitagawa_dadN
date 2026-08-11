# FEM-v3 / PD-m1 local-tip baseline audit

## Scope and provenance

This is a read-only comparison. No canonical FEM-v3 condition was run,
resumed, extended, or modified. The source campaign is the accepted
`v9_canonical_four_class_elastic_fem_signed_mpz_energy_gated_birth_v3` tree at
commit `5d6b1eef331ca88a4a55c8a34fc7d988f2882ac9`. The material is Peak
`v913_paper_peak01_0242980_persistent_sites`, at 300 K, R=0.1 and 1000 Hz.
Registry-row and constitutive hashes are retained in every atomic capsule and
the source-file hashes are recorded in `source_inventory.json`.

The exact K360 mechanics request is
`runs/sn_v9_single_seed_skeleton/K360_FAILURE/solver_output/shielded/sigmaA_735p921MPa/run_args.json`.
It specifies an elliptical edge notch with a=360 um, b=127.279 um, nominal
rho=45 um, and a mesh-resolved root radius of 45.5379 um. The blunt PD geometry
has a=150 um, b=300 um, nominal rho=600 um, and mesh radius 598.745 um.

## Direct, reconstructed, and derived data

Direct atomic FEM-v3 quantities are the phase-resolved 2x2 root tensors, root
opening, signed channel stresses, raw and m=3 cleavage log rates, cumulative
action, instantaneous cycle hazard, signed mobile/retained/wake arrays,
accumulated slip, backstress, shielding, blunted radius, displacement phases,
and checkpoint generation. Peak checkpoints at 735.921495, 650 and 631 MPa and
N=1e8, 1e10, 1e12 and 1e14 are included.

The 804.139 MPa finite-life event is direct event-ledger evidence: stable birth
at N=255.065512355417, one admitted attempt, Xi=0.3753690264 and admitted length
2.29734 um. Its event file does not contain a signed-MPZ or phase-array
checkpoint, so checkpoint-state fields were not fabricated or inferred from a
different trajectory.

For the 2000 and 1500 MPa blunt-notch cases, the accepted mesh, neutral/evolved
FEM arrays, displacement, run arguments and shared-root capsule were read from
their hash-verified checkpoints. The exact 32-phase root tensor was
deterministically reconstructed through the same cached FEM stress-history
path used by the existing read-only stress-scale audit. This consumed no cycles
and did not write into either source result directory. Raw m=1 and diagnostic
m=3 actions were evaluated side-effect-free from copied signed-MPZ states.

Analytical Kt is the ellipse estimate `1+2a/b`: 6.657 for K360 and 2.0 for the
blunt notch. Direct FEM Kt is 5.5655 and 1.7070 respectively. Effective stress
uses the accepted signed shielding and blunted radius. The only interpolation
used in the interpretation below is log(action) versus effective stress between
existing 631 and 650 MPa FEM checkpoints. Both PD points lie outside that local
stress range; any extension is explicitly called extrapolation.

## Layer A: mechanics transfer

At N=1e14 the selected local states are:

| model | sigma_a | sigma_max | FEM Kt | max opening | max effective stress |
|---|---:|---:|---:|---:|---:|
| K360 FEM, 631 MPa | 631 | 1402 | 5.565 | 7804 MPa | 7804 MPa |
| K360 FEM, 650 MPa | 650 | 1444 | 5.565 | 8039 MPa | 8039 MPa |
| K360 FEM, 735.921 MPa | 735.921 | 1635 | 5.565 | 9102 MPa | 9102 MPa |
| blunt PD m1, 2000 MPa | 2000 | 4444 | 1.707 | 7481 MPa | 7481 MPa |
| blunt PD m1, 1500 MPa | 1500 | 3333 | 1.707 | 5610 MPa | 5610 MPa |

The FEM-Kt ratio is 3.260. Extrapolating the elastic K360 response to the 7481
MPa local stress of the PD 2000 MPa condition gives sigma_a about 605 MPa; the
nominal-amplitude ratio 2000/605=3.31 is close to the Kt ratio. Thus essentially
all of the nominal-scale difference at the upper PD point is explained by the
macronotch transfer, not a different cleavage barrier. The mesh-radius ratio
is 13.15 and is already embodied in these two resolved FEM fields; no second
Inglis or K-radius multiplier was applied.

## Layer B: constitutive parity

`FEM_PD_CONSTITUTIVE_PARITY.json` starts both implementations from the same
copied N=1e8 K360 signed-MPZ capsule and supplies the same stored phase tensors.
The canonical and shared-root m=3 action is bit-identical at
8.829782039934813e-8 per cycle. Mobile, retained and accumulated-slip arrays,
aggregate emission, shielding and tip radius are also bit-identical after the
cycle. The independent raw m=1 evaluation differs by 2.01e-14 absolute
(9.90e-14 relative), solely because one path averages linear rates and the
other uses log-sum-exp. Tensor opening and signed shear projections are
identical. This passes the 5e-13 numerical parity tolerance.

There is therefore no detected difference in root-node convention, tensor
components, crystal projection, phase ordering, initial radius, registry row,
source multiplicity, barrier, signed-state update, shielding transform, or m=3
renewal transform.

As a stress-law cross-check, a short extrapolation of the 631--650 MPa FEM raw
m1 curve predicts 7.23e-5 action/cycle at the PD-2000 effective stress; the
actual PD checkpoint gives 6.60e-5 (9% lower). This is close agreement given
that the comparison is just outside the FEM bracket and uses differently aged
MPZ states. The 1500 MPa point is far outside the FEM local-stress range; its
corresponding extrapolation is not scientifically qualified.

## Layer C: endpoint conversion

FEM-v3 integrates the canonical m=3 completed-event process, evaluates the
post-first-passage fixed-opening energy gate, and stops at the first admitted
stable birth. PD production integrates raw m=1 elementary attempts, draws a
normalized spatial mark, injects a reversible embryo, then resolves healing or
stabilization, bond softening, root connection and front capture. These are not
the same lifetime random variable. FEM stable-birth life is used only to
benchmark local hazard and the energy gate; PD determines conditional
seed-to-front behavior. `exp(-H_attempt)` is not PD stable-front survival.

The stored FEM energy-gate envelope is fully open over the evaluated phase/Xi
grid for these Peak conditions. PD does not use that gate; it uses its separate
embryo-transition and topology physics. This distinction is machine-readable
in `FEM_PD_ENDPOINT_SEMANTICS.json`.

## Long-cycle behavior

The direct K360 instantaneous m=3 cycle hazards remain positive through 1e14:
approximately 8.83e-8 at 735.921 MPa, 3.08e-14 at 650 MPa and 7.47e-16 at 631
MPa. Their changes over the final decades are tiny, while the signed state is
still formally evolving; the saved diagnostics do not claim an exact periodic
orbit. Cumulative action consequently grows approximately linearly at the
observed tail. This supports a positive practical attempt floor over the saved
horizon, not an infinite-life theorem.

Before its first attempt, the PD 1500 MPa trajectory accumulates 0.0373635
action over 1.82037e8 cycles, consistent with its reconstructed late-state raw
m1 action of 2.05e-10 per cycle. After stabilization the production policy
stops the global clock, so post-seed PD topology time cannot be interpreted as
additional attempt-hazard exposure. Neither finite censoring nor a stopped
post-seed clock proves endurance.

## Stress-bracket implication

The sharper-tip reference does not justify copying FEM nominal stresses into
the blunt geometry. It does validate the local constitutive mapping: the 2000
MPa blunt condition reaches almost the same local stress/hazard regime as a
roughly 605 MPa K360 state, primarily through the 3.26 Kt ratio. The existing
1500--2000 MPa Peak m1 bracket therefore remains the correct smooth-blunt
baseline for further numerical/probabilistic qualification. The FEM reference
does not by itself recommend a new PD stress, and no new condition was run.
