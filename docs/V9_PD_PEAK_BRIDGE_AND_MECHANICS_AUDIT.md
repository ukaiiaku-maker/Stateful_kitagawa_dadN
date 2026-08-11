# Peak PD mechanics and emission-delivery bridge audit

## Scientific classification

Neither preserved Peak calculation is a production S-N point.

- `Peak_blunt_4000MPa_K2_suppressed_runout_diagnostic` is the historical
  per-source-delivery trajectory. It omitted the audited aggregate persistent
  source multiplicity and is retained only for K=2/tail numerics.
- `Peak_blunt_12000MPa_overstress_topology_diagnostic` has remote
  `sigma_max = 26.667 GPa`, FEM `sigma1_max = 54.232 GPa`, and
  `max(sigma1/E) = 0.1323`. It is outside the conservative small-strain audit
  bound and did not reach the declared six-bond front-capture endpoint by
  `N=3e6`. Stable-site and root-connection events are not relabeled as the
  primary endpoint.

The mechanics audit uses `max(sigma/E) <= 0.05` as a documented formulation
bound: linearized kinematics neglect terms of scale `epsilon^2`, bounded here
by `0.0025`. This is not a fitted fatigue threshold. The 4 GPa diagnostic has
remote `sigma_max/E = 0.02168`, terminal FEM `sigma1_max = 18.077 GPa`, and
local `sigma1/E = 0.04409`; its mechanics pass this bound, but its delivery
bridge does not, so it remains non-production.

## Executable bridge ledger

At the controlling 4 GPa node and phase, the exact terminal reconstruction
gives the FEM tensor (Pa)

```
[[ 8.98528065e8, -5.70944157e9],
 [-5.70944157e9,  1.53184158e10]]
```

The legacy scalar-equivalent projection is `16.464 GPa`. The audited
crystallographic projection instead gives two signed drives of `7.210 GPa`.
With zero initial signed backstress, the audited per-source rates are
`1.571e-18 s^-1` per channel. The exact v10.2.21 source geometry is:

- source density: `1.410252084e16 m^-2`;
- reference source area: `2.5e-11 m^2` (25 um2);
- multiplicity per system: `M = 352563.0211`;
- aggregate signed-source rate summed over two channels: `1.108e-12 s^-1`.

The historical PD adapter instead used the scalar per-source emission rate
`8.679e-15 s^-1`. Merely multiplying that scalar rate by M would produce
`3.060e-9 s^-1`, 2,762 times the authoritative tensor-resolved aggregate.
Therefore the omitted M is real, but an `M*scalar-equivalent-rate` patch is not
constitutive parity.

PD lumped point area is `2.094e-10 m2`; candidate-site density is `5e10 m^-2`
(candidate area `2e-11 m2`). Neither is another emission multiplier: M already
contains the audited source density and source area. Candidate density controls
the separate specimen first-passage population.

At the same legacy state:

- raw cleavage rate: `1.653e10 s^-1`;
- delivery memory: `1.321e-18`;
- K=2 completion: `1.746e-36`;
- actual persistent-site birth action: `2.897e-30 cycle^-1`;
- raw-rate to gated-birth separation: 39.76 decades.

For small Lambda, K=2 is quadratic. The omitted M alone would shift the scalar
route by 11.09 birth-rate decades, but the tensor-resolved authoritative rate
is 3.44 delivery decades below that scalar aggregate, a further 6.88-decade
K=2 difference. This is why no fitted multiplier is admissible.

## Terminal tail diagnostics

The side-effect-free exact cycle leaves RNG, thresholds, cumulative actions,
and topology hashes unchanged. It repairs the missing 4 GPa FEM/Kt fields.
For the preserved legacy tail:

| horizon | max node H | 1-S_no_embryo | minimum remaining action | controlling rate /cycle |
|---:|---:|---:|---:|---:|
| 1e10 | 2.897e-20 | 2.897e-20 | 4.247e-4 | 2.897e-30 |
| 1e12 | 2.896e-18 | 2.896e-18 | 4.247e-4 | 2.902e-30 |
| 1e14 | 2.887e-16 | 2.887e-16 | 4.247e-4 | 2.052e-30 |

These are finite-horizon numerical diagnostics, not endurance limits.

## Fail-closed production status

The repository resolves per-source versus aggregate emission and rules out the
scalar equivalent-stress projection. It does not resolve how the authoritative
root-local signed MPZ state is to couple to the spatial PD patch. At least two
physically distinct implementations remain:

1. evolve one authoritative root-local signed MPZ and broadcast/distribute its
   aggregate tensor-resolved delivery to the PD process region;
2. evolve node-local signed MPZ states under each node tensor, which changes
   source geometry, backstress, shielding, and delivery spatially;
3. use canonical cleavage attempts directly as PD embryo attempts, retaining
   stabilization/healing/topology but removing K=2 delivery from birth.

The authoritative fracture source is root-local and supplies no spatial-PD
allocation rule. Alternatives 1 and 2 already produce different state and area
semantics; alternative 3 is separately represented in the repository's
canonical-cleavage branch. Consequently corrected Peak production is blocked
fail-closed pending a physics choice. DBTT, weak-T, and ceramic remain unstarted.

## Figures and files

The seven compact plots under
`runs/v9_pd_production/figures/Peak_bridge_diagnostics/` use the packaged 4 GPa
terminal exact audit, the packaged 12 GPa history, the site-matched
`1e10/1e12/1e14` table, and the accepted HC controller history. They visibly
identify both stress conditions as diagnostics rather than a physical bracket.
The complete dimensional ledger is
`runs/v9_pd_production/diagnostics/Peak_emission_delivery_bridge/Peak_4000MPa_executable_bridge_ledger.json`.
