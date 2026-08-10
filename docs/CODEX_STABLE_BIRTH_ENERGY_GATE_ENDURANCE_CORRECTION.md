# Correction: stable-crack birth must include the existing post-first-passage energy-admissibility gate

## Scope

This document supersedes the previous 300 K endurance interpretation that classified the canonical four-class branch solely from the intrinsic zero-drive cleavage first-passage hazard floor.

The positive zero-drive cleavage renewal rate proves only that **cleavage attempts** continue with nonzero probability for arbitrarily long time. It does **not** by itself prove that a **stable crack is born**.

The scientific endpoint remains:

```text
N_stable_crack_birth
```

but for the canonical four-class fatigue model this must be defined consistently with the existing v10.2.30 event physics:

```text
stable crack birth = first cleavage first-passage attempt for which the existing
post-first-passage fixed-opening elastic-energy balance admits a nonzero,
physically resolvable crack event.
```

A cleavage first-passage attempt whose existing energy gate admits zero crack extension is a **nonpropagating cleavage attempt**, not stable-crack birth.

Do not alter the cleavage hazard or add an athermal fracture trigger.

## Why the previous no-endurance proof is insufficient

The current canonical branch intentionally gives a finite thermally activated cleavage-attempt rate even at zero opening stress. Because the multihit renewal transform is strictly positive for every positive raw rate, the cumulative *attempt* hazard diverges as N->infinity.

That is mathematically correct, but if every first passage is identified with stable-crack birth then even an unloaded specimen eventually develops a stable crack by definition. This collapses the physical distinction between thermal cleavage attempts and energetically sustainable crack formation.

The existing fracture/fatigue model already preserves that distinction. First passage creates a stochastic event-length proposal; only afterward is the proposed crack event truncated/admitted by the fixed-opening elastic-energy balance. The continuum K^2/E' comparison remains diagnostic and must not gate the hazard before first passage.

Therefore retain two separate classifications:

```text
cleavage_attempt_asymptotic_classification
stable_crack_birth_asymptotic_classification
```

The existing result

```text
cleavage_attempt_asymptotic_classification = no_endurance
```

may remain as a diagnostic. It must not be reported as stable-crack-birth no-endurance until the post-first-passage energy-admissibility process is included.

## Exact physics to reuse

Reuse the existing qualified v10.2.30 event sequence without inventing new parameters:

1. Accumulate the canonical cleavage first-passage action using the exact four-class signed-MPZ state.
2. Localize the first-passage threshold crossing exactly.
3. Use the same sampled threshold to construct the existing `threshold_scaled`, mean-preserving event-length proposal.
4. Apply the existing **post-first-passage** fixed-opening elastic-energy balance at the event state/Kmax.
5. If the admitted event length is zero / no mesh-resolved admissible increment exists:
   - classify as `nonpropagating_cleavage_attempt`;
   - consume the threshold/attempt exactly as the qualified fatigue model does;
   - draw the next threshold from the continuing canonical RNG stream;
   - do not modify crack geometry;
   - continue the pre-birth cyclic state.
6. If a nonzero physically admissible event is obtained:
   - classify that instant as `stable_crack_birth`;
   - record proposal length, admitted length, energy release, resistance, event state, threshold, and cycle;
   - STOP. Do not propagate the crack further.

No post-birth da/dN, front capture, long crack growth, or handoff is part of this S-N endpoint.

## No athermal fracture criterion

Do not suppress or rescale the cleavage hazard based on K^2/E' or any deterministic toughness threshold.

The energy test acts only on the **reward/admissibility after a thermally activated first passage**. A rejected attempt remains a real consumed stochastic attempt.

Do not add `Gc0_athermal`, a fitted fatigue limit, a hard stress threshold, or a stress cutoff.

## Stable-birth survival semantics

After this correction,

```text
exp(-H_cleave)
```

is the probability of **no cleavage first-passage attempt**, not the stable-crack-birth survival probability.

The stable-birth process is a marked renewal process because:

- attempt times are controlled by successive Exp(1) first-passage thresholds;
- the same threshold also controls the threshold-scaled event-length proposal;
- the event gate may reject some attempts as nonpropagating;
- rejected attempts consume the threshold and advance the RNG stream.

Therefore do not use `S_stable = exp(-H_cleave)` after the energy-admissibility correction.

Provide separately:

```text
S_no_attempt(N)
S_stable_birth(N)
```

For stable-birth survival, either derive the exact accepted-event renewal law or use a validated threshold-sequence Monte Carlo/quadrature over a deterministic no-crack physical trajectory. Do not brute-force the full FEM separately for every threshold realization if the pre-birth physical state is identical for rejected attempts.

## Efficient 1e12-1e14 strategy

The user still wants the 300 K VHCF regime through approximately 1e12-1e14 cycles where useful.

Because nonpropagating attempts do not change geometry, exploit the structure of the problem:

1. Build/validate a deterministic pre-birth physical trajectory for each class/stress:
   - FEM cyclic state;
   - signed mobile/retained/wake MPZ state;
   - shielding/backstress/blunting;
   - canonical cleavage action rate;
   - event-energy admissibility as a function of state and sampled threshold/event proposal.
2. Qualify a quiet-tail acceleration against direct continuation over overlapping windows.
3. Use the accelerated state trajectory to 1e10, 1e12, and, when justified, 1e14 cycles.
4. Evaluate many threshold sequences cheaply against that trajectory rather than rerunning FEM for every stochastic realization.
5. Retain exact event localization and exact energy-gate evaluation at candidate accepted births.

## Endurance criterion after correction

The stable-crack-birth endurance question is now:

```text
Does the cumulative probability/intensity of ENERGY-ADMISSIBLE cleavage events
remain finite, or does an admissible stable-birth event occur with probability 1?
```

Possible asymptotic outcomes include:

- `endurance_supported`: after cyclic shakedown the post-first-passage energy gate becomes permanently closed, or the accepted-event cumulative measure is provably finite;
- `no_endurance`: the settled state retains a nonzero accepted-event probability/rate so eventual stable birth is certain;
- `undetermined`: finite-horizon evidence is insufficient to establish either.

A positive intrinsic cleavage-attempt floor alone is not sufficient for `no_endurance` of stable crack birth.

## 300 K campaign sequence

Use only the final canonical rows:

```text
Peak     = v913_paper_peak01_0242980_persistent_sites
DBTT     = v913_paper_dbtt01_0202500_persistent_sites
weak-T   = v913_paper_weakT01_0129902_persistent_sites
ceramic  = v913_paper_ceramic01_0077080_persistent_sites
```

300 K only. Do not resume temperature work.

Proceed one class/condition at a time:

1. implement and validate the stable-birth energy-admissibility endpoint on a short condition with both admitted and rejected attempts;
2. prove restart and block-partition invariance including threshold renewal after rejected attempts;
3. reclassify existing 300 K low-stress trajectories as `cleavage_attempt_tail` diagnostics, not final stable-birth endurance results;
4. select low stresses for Peak, DBTT, weak-T, ceramic;
5. continue to 1e8 -> 1e10 -> 1e12 -> 1e14 where the validated acceleration permits;
6. report stable-birth survival/endurance separately from cleavage-attempt survival.

Do not alter the four parameter rows and do not tune the energy gate to create or remove an endurance limit.

## Required outputs

Maintain:

```text
runs/sn_v9_canonical_four_class/300K_cleavage_attempt_tail_summary.csv
runs/sn_v9_canonical_four_class/300K_stable_birth_endurance_summary.csv
runs/sn_v9_canonical_four_class/stable_birth_survival_results.csv
runs/sn_v9_canonical_four_class/endurance_diagnostics.json
runs/sn_v9_canonical_four_class/condition_registry.json
runs/sn_v9_canonical_four_class/campaign_manifest.json
```

Every endurance diagnostic must state whether it concerns:

```text
cleavage_attempt
or
energy_admissible_stable_crack_birth
```

## Stop conditions

Continue autonomously unless:

- the exact v10.2.30 energy-admissibility calculation cannot be transferred without a new constitutive choice;
- stable-birth restart/partition equivalence fails for a physical reason;
- the quiet-tail acceleration cannot be bounded against direct trajectories;
- audited parameter/source provenance becomes ambiguous.

Do not stop merely because another short audit passes.
