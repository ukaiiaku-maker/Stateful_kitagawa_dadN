# Canonical four-class cleavage-to-stable-birth decision

## Decision

For the canonical four-class S-N parity branch, use the **audited sharp-front cleavage first-passage process itself as the stable-crack-birth event**.

Do not splice the older v9 spatial PD candidate-site/K=2 embryo-birth law into the audited four-class closure.

The four-class production branch therefore uses:

- exact audited signed-MPZ / aggregate persistent-emission evolution before fracture;
- one canonical front/root-local stochastic cleavage first-passage clock;
- the canonical cleavage barrier and multiplicity exactly as defined by the audited sharp-front executable;
- stable-crack-birth endpoint at the first canonical cleavage crossing.

The older v9 `completion_gated_independent_cleavage` PD candidate-site branch remains a separately named historical/legacy S-N model. It is not the canonical four-class parity branch.

## Why this is the only exact-transfer choice

The recovered audited closure already defines the complete pre-fracture state evolution and the stochastic fracture process:

1. emission is the audited deterministic aggregate persistent-source law;
2. signed mobile/retained MPZ state evolves through the audited Peierls/Taylor, correlation, escape/advection, shielding, backstress and blunting closure;
3. cleavage is a stochastic front-local first-passage process;
4. when the cleavage clock fires in the canonical sharp-front model, its reward is irreversible crack advance.

That cleavage event is therefore already an irreversible/stable fracture event. For the S-N problem whose endpoint is only the **birth of a stable crack**, the same first cleavage event is the natural stopping event.

By contrast, the old v9 PD initiation branch introduces additional constitutive structure that the canonical four-class executable does not define:

- many spatial candidate sites;
- per-site Exp(1) embryo-birth thresholds;
- an additional K=2 delivery-completion multiplier;
- reversible embryo stabilization/healing;
- a candidate-site density/multiplicity law.

Multiplying canonical cleavage by the K=2 gate, distributing one scalar MPZ shielding value over PD nodes, or assigning a scalar cleavage event to a PD site would all be new constitutive mappings. None is required when the endpoint is stable crack birth.

## Exact event definition

For this branch define

```text
fatigue_endpoint = canonical_cleavage_stable_crack_birth
```

and

```text
N_stable_crack_birth = N_first_canonical_cleavage_first_passage
```

The event is evaluated at the fixed notch/root coordinate before crack birth. The signed MPZ is root/front-local and remains centered on that fixed root until the event occurs.

No PD embryo-site identity is created. Record the deterministic root/front coordinate and the canonical front/kernel orientation/identity instead.

## Cleavage clock semantics

Use the audited canonical stochastic cleavage machinery directly or an exact state-preserving port of it.

Requirements:

- one canonical root/front-local renewal/first-passage clock as in the authoritative executable;
- persistent threshold/state sampled exactly according to the canonical implementation;
- integrated canonical cleavage hazard advanced transactionally;
- exact/conservative localization of the first crossing inside a proposed large-N block;
- event committed once at an accepted boundary;
- no candidate-site-density multiplier unless the canonical executable itself contains it;
- no K=2 completion factor in this canonical branch;
- no additional embryo stabilization/healing stage after cleavage crossing.

The stable event is irreversible by definition because the canonical cleavage event is the event that would produce crack advance in the fracture solver. Since the present S-N calculation stops at stable crack birth, no crack-advance distance needs to be sampled or executed.

## Interaction with signed MPZ state

Before the cleavage event preserve the exact audited state and coupling:

- independent cleavage barrier;
- aggregate persistent emission;
- independent emission barrier;
- Peierls kinetics;
- Taylor kinetics and correlation;
- positive/negative mobile state;
- positive/negative retained state;
- accumulated slip and blunting;
- unsigned active mobile+retained Taylor backstress;
- signed retained shielding;
- escape/advection/wake state active before fracture;
- theta=0 production signed kernel;
- exact four audited material rows.

Do not spatially distribute the scalar/root-local signed shielding over PD candidate nodes. Use the canonical front/root-local effective loading and shielding exactly as the sharp-front engine does.

## Survival/endurance consequence

For the canonical four-class branch, if the only stochastic fracture event before the endpoint is the canonical cleavage first-passage clock and the signed-MPZ evolution is deterministic conditional on no cleavage, then the stable-crack-birth survival law is directly

```text
S_stable(N) = exp[-H_cleave(N)]
```

where `H_cleave` is the cumulative canonical cleavage hazard along the no-fracture state trajectory.

This is stronger and cleaner than the legacy multistage embryo model:

- finite `H_cleave(infinity)` directly supports a nonzero stable-crack survival asymptote;
- divergent `H_cleave(infinity)` implies eventual canonical stable-crack birth probability 1 under that model;
- finite runout alone still does not establish either classification.

If additional stochastic pre-birth coordinates are found in the authoritative canonical source, include them exactly and update this statement. Do not invent them.

## Legacy v9 results

Do not delete or reinterpret the completed legacy v9 PD-candidate-site results.

Label them explicitly, for example:

```text
model_class = legacy_completion_gated_PD_candidate_sites
```

Their `N_stable_crack_birth` values and right-censored 690.443 MPa trajectory remain useful as a separate mechanistic comparison, but they are **not** parameter-parity results for the canonical four fracture classes.

The canonical four-class S-N campaign must be rerun using the root-local canonical cleavage endpoint because the governing event law is different.

## Four canonical rows

Use exactly:

```text
Peak     = v913_paper_peak01_0242980_persistent_sites
DBTT     = v913_paper_dbtt01_0202500_persistent_sites
weak-T   = v913_paper_weakT01_0257068_persistent_sites
ceramic  = v913_paper_ceramic01_0189364_persistent_sites
```

The later v10.4 weak-T/ceramic alternatives remain excluded from this campaign.

## Next autonomous implementation and run sequence

1. Implement a new versioned canonical four-class endpoint path; do not overwrite the legacy PD candidate-site branch.
2. Port/reuse the canonical cleavage clock exactly and checkpoint all its authoritative stochastic state.
3. Demonstrate prescribed-history equality of cumulative cleavage hazard and crossing cycle against the audited sharp-front executable for all four rows at representative temperatures.
4. Demonstrate block-partition and atomic restart equivalence for the full signed-MPZ + cleavage-clock trajectory.
5. Run one real 300 K condition for Peak first, stopping immediately at the first canonical cleavage crossing.
6. If valid, qualify DBTT, weak-T and ceramic sequentially at 300 K.
7. Build an adaptive single-condition-at-a-time S-N skeleton for each class using stable-crack birth as the endpoint.
8. Then begin temperature comparison with a compact subset (e.g. 300, 900, 1300 K) before filling the audited temperature grid adaptively.
9. Maintain exact class provenance, request hashes, signed-state diagnostics, cumulative cleavage hazard and endurance diagnostics.
10. Do not run post-birth crack growth, front capture, handoff or da/dN in this campaign.

## Scientific comparison to preserve

The central hypothesis is allowed to emerge from the transferred physics:

- Peak may be strongly nonmonotonic in temperature;
- DBTT may show a strong transition;
- weak-T may remain much less temperature sensitive;
- ceramic may remain a brittle/weak-plasticity control.

Do not tune these responses toward one another.
