# Physics decision: authoritative aggregate persistent emission for four-class S-N

## Decision

For the present four-class stable-crack-birth S-N campaign, use the **exact audited v10.2.21 aggregate persistent-emission closure** recovered from `persistent_site_source_v10221`.

Do **not** invent per-source stochastic identities, Exp(1) emission thresholds, emission outcome uniforms, or an emission RNG state. Those coordinates are not part of the audited executable closure and adding them would define a new constitutive model rather than transfer the existing one.

The source audit is authoritative over earlier shorthand prose that described emission as stochastic. In this exact closure:

- stochastic cleavage first passage remains stochastic;
- embryo stabilization/healing remains the persistent competing stochastic process already used by v9 through the stable-crack-birth endpoint;
- emission is an Arrhenius, stress/backstress-dependent **aggregate mean activation law** solved implicitly and deposited continuously into the signed mobile/retained state;
- `site_capacity` and `available_sites` are compatibility/diagnostic fields, not a discrete stochastic ledger.

`persistent` therefore means that the source/MPZ state and its consequences persist and evolve; it does not imply a population of individually sampled stochastic emission sites.

## Exact emission law to preserve

The audited implementation solves the implicit mean activation count

```text
N = M * lambda(sigma_drive - sigma_back(N), T) * dt
```

and deposits the resulting continuous signed line content according to the existing signed-channel state-update rules.

Port this executable law directly. Preserve its implicit backstress coupling and all active audited parameters. Do not replace it with a Poisson draw or first-passage source clock.

## Scope of stochasticity

The v9 four-class state should explicitly distinguish:

```text
aggregate_emission_mean_activation       # deterministic conditional constitutive response
signed_mobile_plus/minus
signed_retained_plus/minus
accumulated_slip / blunting state
Taylor/backstress state
shielding state

cleavage_first_passage_thresholds        # stochastic
embryo/stabilization-healing clocks       # stochastic/persistent as already qualified
stable_crack_birth                        # primary S-N endpoint
```

Do not add an `emission_rng_state` merely for symmetry with cleavage.

## Why this is the production choice

The scientific requirement is to use the **same four audited parameter rows and executable constitutive physics** as the fracture/fatigue framework without class-specific refitting. Exact executable transfer takes precedence over a prior prose expectation of stochastic emission when the audited source itself is deterministic aggregate emission.

A different audited discrete-emission implementation may be investigated later as a separately named model if provenance can be established, but its possible existence does not block the present campaign. It must not silently replace the qualified v10.2.21 closure.

## Required Codex action

Continue the signed-MPZ port using the exact aggregate persistent-emission closure.

1. Preserve the immutable audited executable package and hashes.
2. Use the exact four canonical rows:

```text
Peak     = v913_paper_peak01_0242980_persistent_sites
DBTT     = v913_paper_dbtt01_0202500_persistent_sites
weak-T   = v913_paper_weakT01_0257068_persistent_sites
ceramic  = v913_paper_ceramic01_0189364_persistent_sites
```

3. Preserve independent cleavage, emission, Peierls and Taylor parameter surfaces and the signed mobile/retained state.
4. Keep emission deterministic/aggregate exactly as the audited source defines it.
5. Keep stochastic cleavage first passage and the qualified stochastic stabilization/healing process unchanged.
6. End each S-N realization at the first stable-crack-birth event.
7. Prove prescribed-history constitutive parity and real block/restart equivalence with the aggregate-emission port.
8. Then proceed to four-class 300 K qualification one condition at a time, followed by the temperature-dependent stable-crack-birth S-N plan.

## Reporting

Campaign manifests and summaries should identify the source closure explicitly, for example:

```text
emission_model = audited_v10221_aggregate_persistent_emission
emission_stochastic = false
cleavage_first_passage_stochastic = true
fatigue_endpoint = stable_crack_birth
```

Report aggregate source activation/current signed content diagnostics rather than fictional per-source identities.

## No new physics in this decision

This decision does not change the audited material rows or equations. It removes an erroneous requirement to preserve stochastic source coordinates that never existed in the authoritative implementation.
