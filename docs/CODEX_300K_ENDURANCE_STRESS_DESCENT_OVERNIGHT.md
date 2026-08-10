# 300 K endurance-limit stress-descent overnight plan

## Governing objective

The canonical v3 architecture is accepted:

`v9_canonical_four_class_elastic_fem_signed_mpz_energy_gated_birth_v3`

Do not change the model architecture, material rows, post-first-passage energy gate, or stable-birth semantics.

The existing four low-stress conditions were advanced through 1e14 cycles and are **not low enough to probe a practical endurance regime**. At those stresses every evaluated phase/Xi mark is energy-admissible and stable-birth survival collapses to effectively zero by sufficiently large N.

The next work package should therefore move **downward in stress**, not densify the earlier S-N knee and not spend the night trying to prove pointwise signed-MPZ stationarity at the same stresses.

300 K only.

## Existing canonical v3 anchors

Retain these results unchanged:

- Peak 735.921495 MPa
- DBTT 700 MPa
- weak-T 350 MPa
- ceramic 400 MPa

They establish an above-endurance reference regime in which the evaluated gate is fully open.

## Scientific question

For each class, determine whether lowering stress leads to either:

1. an energy-admissible stable-birth process with finite nonzero event probability through 1e14 and, asymptotically, no endurance; or
2. a regime where the post-first-passage energy gate closes for all physically possible marks after pre-birth shakedown, supporting a true stable-crack endurance limit.

Do not infer the answer from cleavage-attempt hazard alone.

## Stress-selection strategy

Do not use a uniform dense stress grid.

Use a sequential downward stress search for each class.

For a candidate stress:

1. Use the cached elastic phase response.
2. Advance the signed MPZ with the already-qualified large-N machinery.
3. Evaluate at 1e8 first, then 1e10, 1e12, and 1e14 only as needed.
4. At each checkpoint evaluate the complete phase/Xi energy-gate envelope.
5. Track both cumulative cleavage action H and energy-admission measure.

The immediate target is a stress where at 1e14 either:

- `H_cleave` is O(0.1-10), so survival is nontrivial; or
- the gate begins to reject a finite fraction of marks; or
- the gate is completely closed.

If H(1e14) is still >>10 and the gate is fully open, lower stress again.

If H(1e14) <<0.1 and the gate is fully open, a slightly higher stress may be useful only if needed to bracket the practical 1e14 median; do not densify unnecessarily.

The stress search should be driven by the model response, not by fixed arbitrary increments. A multiplicative or secant/log-hazard predictor may be used only for selecting the next stress; every reported physical result must come from the exact v3 trajectory.

## Practical endurance outputs

For every accepted candidate stress report:

- `S_stable_birth(1e8)`
- `S_stable_birth(1e10)`
- `S_stable_birth(1e12)`
- `S_stable_birth(1e14)`
- `H_cleave` at the same horizons
- conditional attempt-admission probability at each horizon
- fraction of phase/Xi envelope rejected
- minimum and maximum admitted length
- minimum energy margin across the phase/Xi envelope

Also compute, where bracketed, the stress associated with:

- `S_stable_birth(1e14)=0.9`
- `S_stable_birth(1e14)=0.5`
- `S_stable_birth(1e14)=0.1`

These are practical VHCF endurance/survival stresses, not strict infinite-life proofs.

## Strict endurance classification

Keep the strict infinite-life classification separate from the finite 1e14 survival result.

`endurance_supported` requires a model-derived argument that after some N* no physically possible cleavage attempt can produce a nonzero admitted stable crack event.

A strong route is a fully closed energy-gate envelope plus a proven monotone/conservative bound showing that later signed-MPZ evolution cannot reopen it.

`no_endurance` requires a model-derived positive lower bound on the accepted-event measure/rate, not just one admitted transient event.

If neither is established, retain `undetermined` even if S(1e14) is numerically close to zero or one.

## Overnight sequence

Proceed sequentially, one class at a time:

1. Peak
2. DBTT
3. weak-T
4. ceramic

For each class, descend stress until at least one genuinely informative 1e14 condition is obtained (nontrivial survival and/or partial/full gate closure).

Do not spend the night adding many nearby stresses around the earlier sharp S-N drop.

If a class reaches a fully closed gate before 1e14, use the remaining time to test one slightly higher stress for the boundary and one lower stress for robustness.

If a class remains fully open even after H(1e14) becomes <<1, preserve that result; it shows practical long life without a strict gate-defined endurance limit in the tested range.

## Required canonical rows

- Peak = `v913_paper_peak01_0242980_persistent_sites`
- DBTT = `v913_paper_dbtt01_0202500_persistent_sites`
- weak-T = `v913_paper_weakT01_0129902_persistent_sites`
- ceramic = `v913_paper_ceramic01_0077080_persistent_sites`

## Required outputs

Maintain atomically:

- `runs/sn_v9_canonical_four_class_elastic_v3/300K_endurance_stress_descent.csv`
- `runs/sn_v9_canonical_four_class_elastic_v3/300K_stable_birth_survival.csv`
- `runs/sn_v9_canonical_four_class_elastic_v3/300K_energy_gate_envelope.csv`
- `runs/sn_v9_canonical_four_class_elastic_v3/endurance_diagnostics.json`
- `runs/sn_v9_canonical_four_class_elastic_v3/condition_registry.json`
- `runs/sn_v9_canonical_four_class_elastic_v3/campaign_manifest.json`

Preserve all existing v3 checkpoints and all prior noncanonical bulk-plasticity controls.

## Stop conditions

Stop only for:

- a genuine constitutive/provenance ambiguity;
- failure of restart or block-partition equivalence;
- failure of the already-qualified large-N accelerator outside its validated error bounds;
- a new physical state boundary not defined by the canonical signed-MPZ model.

Do not stop merely because one class has completed a useful lower-stress point. Continue through the four-class sequence overnight when computationally practical.

Do not run temperature dependence or post-birth crack growth in this work package.
