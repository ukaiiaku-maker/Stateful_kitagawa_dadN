# Four-class parameter provenance correction

## Why this correction exists

The first full-FEM canonical S-N calculations exposed an unphysical zero-load cleavage saturation for the row

```text
v913_paper_ceramic01_0189364_persistent_sites
```

Do **not** modify the cleavage law or insert an ad hoc zero-stress suppression to repair this result.

A later and more specific project handoff resolves the provenance problem. The document **Four-Class 2-D PF/CZM Parameter Handoff**, prepared 25 July 2026, defines the final four-class paper/parity transfer after focused 2-D validation. Its authoritative PF source is:

```text
ukaiiaku-maker/PF-fracture-fatigue
branch: v10.2.22-physical-front-width-top5-dbtt-screen
canonical four-class commit: 0a340f6
```

with registry/selection files:

```text
arrhenius_fracture/data/materials/v10_2_27_v913_four_class_paper_registry.csv
arrhenius_fracture/data/materials/v10_2_27_v913_four_class_paper_selection.json
```

The final production rows in that handoff are:

```text
Peak     = v913_paper_peak01_0242980_persistent_sites
DBTT     = v913_paper_dbtt01_0202500_persistent_sites
weak-T   = v913_paper_weakT01_0129902_persistent_sites
ceramic  = v913_paper_ceramic01_0077080_persistent_sites
```

The earlier rows

```text
v913_paper_weakT01_0257068_persistent_sites
v913_paper_ceramic01_0189364_persistent_sites
```

are explicitly documented in that final handoff as **rejected earlier candidates** after long 2-D transfer behavior did not match the intended classes. They may be retained only as provenance/failed-transfer controls, not as paper-class production primaries.

This correction supersedes earlier Stateful-SN handoff text that incorrectly called 0257068/0189364 the canonical weak-T/ceramic pair.

## Immediate consequence

The current production campaign must **not** repair or retune 0189364. Preserve its zero-load result as a diagnostic/provenance artifact and mark it non-production.

Likewise, the completed 0257068 weak-T S-N point is a valid calculation of that historical candidate, but it is not part of the final canonical four-class paper/parity production set. Do not delete it or relabel it as 0129902.

Peak 0242980 and DBTT 0202500 remain unchanged and their completed canonical results may be retained.

## Required source qualification

Before running new weak-T or ceramic S-N conditions:

1. Load the final registry at `v10_2_27_v913_four_class_paper_registry.csv` from the authoritative commit/branch above.
2. Verify the registry and selection hashes against the source manifest.
3. Verify that the four rows resolve exactly to 0242980, 0202500, 0129902, and 0077080.
4. Compare the full active-parameter fingerprints for 0129902 and 0077080 against the source model.
5. Fail closed if a local registry still maps the paper-class names to 0257068 or 0189364.
6. Do not manually transcribe or partially substitute parameter fields.

The active transferred coordinates remain the full audited cleavage, emission, Peierls, Taylor, source-density/correlation, and blunting row; no refitting or parameter rescaling is authorized.

## Sanity check before production

Run a constitutive zero-load audit at 300 K for the final four rows before any S-N trajectory. At minimum report:

```text
G_cleave(sigma=0,T=300 K)
raw cleavage rate
renewal effective cleavage rate
emission rate at zero drive
signed-MPZ initial state
```

A final production row that still predicts renewal-ceiling stable-cleavage hazard at exactly zero applied load must fail closed for source/constitutive review. Do not suppress it numerically.

Also compare these diagnostics with the authoritative sharp-front executable under the same zero-load state, so a Stateful-SN adapter bug cannot be mistaken for a bad material row.

## Existing results classification

Keep the current results, but classify them explicitly:

```text
Peak_0242980_*       canonical production
DBTT_0202500_*       canonical production
weakT_0257068_*      historical/rejected-transfer control
ceramic_0189364_*    historical/rejected-transfer control; zero-load diagnostic invalid for production
```

Do not mix the historical controls into canonical four-class S-N plots or class comparisons.

## Resume sequence

After the final registry/fingerprint and zero-load audits pass:

1. Qualify `weakT 0129902` at 300 K, one condition at a time.
2. Qualify `ceramic 0077080` at 300 K, one condition at a time.
3. Build their adaptive S-N quantile skeletons using the already-qualified canonical cleavage/full-FEM path.
4. Retain existing Peak/DBTT canonical results and fill additional stresses only when needed for useful `log10(N50)` coverage.
5. Once all four final classes have useful 300 K skeletons, resume the 300/900/1300 K comparison.

## Scientific interpretation

The central comparison remains unchanged:

- Peak: potentially strong nonmonotonic temperature response.
- DBTT: strong transition response.
- weak-T 0129902: weak-temperature/FCC-like control.
- ceramic 0077080: weak-plasticity/brittle control.

Do not tune weak-T toward DBTT and do not tune ceramic simply to obtain a convenient S-N range.

The previous 0257068/0189364 runs remain scientifically useful as failed-transfer controls showing why focused 2-D validation changed the final paper-class selection.
