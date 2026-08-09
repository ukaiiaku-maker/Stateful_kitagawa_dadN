# Definition of done

The project is not complete when a single S-N plot is produced. The following gates are required.

## Gate A: preserve current baseline

- Current v8.7 source hashes and K360 regression fixtures must be recorded before physics changes.
- K360 failure remains a valid `physical_handoff` near `3.696e6` cycles for the frozen baseline.
- K360 lower-stress fixture remains a valid `right_censored` result at `1e8` cycles for the frozen baseline.

## Gate B: hazard-tail mathematics

Implement a standalone endurance-classification module and unit tests using synthetic hazards:

- constant positive hazard -> `no_endurance`;
- exponential decay -> `endurance_supported`;
- power law `p=0.5` -> `no_endurance`;
- power law `p=1.0` -> `no_endurance`;
- power law `p=1.5` -> `endurance_supported`;
- insufficient/noisy tail -> `undetermined`.

Classification must expose the fitted/bounded tail, residual-hazard integral, and reason code.

## Gate C: long-N numerical invariance

For one real stress/seed trajectory:

- run with at least two block-size policies;
- demonstrate consistent first-passage/handoff result within stated tolerance;
- interrupt and resume from a checkpoint;
- demonstrate restart equivalence;
- verify no cumulative hazard or persistent-threshold state is double-committed.

## Gate D: two mechanism classes

Implement two explicit model configurations:

1. current independent-cleavage baseline;
2. activity-conditioned initiation closure.

Each configuration must report the physical terms entering its hazard. Do not hard-code the desired endurance label.

## Gate E: small S-N pilot

Before any large sweep, run only:

```text
3 stresses x 3 seeds x 2 mechanism classes
```

with stress points chosen to include clear failure, transition, and very-long-life behavior.

Pass criteria:

- no invalid infrastructure results;
- failure probability is nondecreasing with stress within statistical uncertainty;
- all censors are true runout censors;
- tail classification is reproducible from stored diagnostics;
- outputs are restartable and hash-addressed.

## Gate F: production S-N campaign

Only after A-E pass, run a broad stress/seed campaign spanning approximately `1e3` to `1e12` cycles as needed.

Required outputs:

```text
sn_results.csv
survival_results.csv
endurance_diagnostics.json
campaign_manifest.json
```

Required plots:

1. stress amplitude vs log10(N), with failures and censors distinguished;
2. survival/failure probability vs N for selected stresses;
3. tail hazard and cumulative hazard demonstrating why each model is classified as endurance / no-endurance / undetermined;
4. state evolution (plastic activity, shielding/back stress, relevant dislocation state) for representative high- and low-stress trajectories.

## Gate G: optional da/dN / total-life continuation

Only after the initiation S-N calculation is stable, add state-preserving sharp-front continuation. Report initiation life, propagation life, and total life separately.
