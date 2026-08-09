# v9 stable-crack-birth and four-class audit

## Endpoint qualification

The v9 transaction driver now supports `--fatigue-endpoint stable_crack_birth`.
It limits an accepted transaction to the exact residual competing-transition
clock, commits the stabilization/healing transition once, and stops only when
the first stabilization branch is present at the accepted physical boundary.
The persistent site clocks, outcome uniform, site identity, accepted FEM state,
and stochastic generator state are checkpointed.  The equivalence regression
compares the endpoint state with the same event in a continuing full trajectory.

All five completed 300 K conditions were reanalyzed without reruns.  The four
finite stable-birth lives are 1793124.873503177, 108922.8058062798,
5602.142229371419, and 112.15473848253177 cycles.  The 690.443 MPa trajectory
is right-censored at 100000000 cycles.  Front capture and physical handoff are
retained only in `legacy_physical_handoff_diagnostics.csv`.

For the no-birth censored trajectory, the exact aggregate active-site birth
hazard is 0.5582458473004149 and its no-embryo survival is
0.5722119311041882.  No-embryo survival is deliberately left unidentified on
realized post-embryo paths, and it is never labeled stable-crack survival.

## Canonical registry qualification

The versioned bridge calls the audited `parameter_registry_v9111.select_option`
implementation after verifying the selector, registry, and selection-manifest
SHA-256 hashes.  It uniquely resolves the mandated Peak 0242980, DBTT 0202500,
weak-T 0257068, and ceramic 0189364 rows.  The later weak-T 0129902 and ceramic
0077080 options are explicitly rejected.

The audited campaign temperature grid is:

`300, 600, 800, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300 K`.

## Fail-closed physical-adapter boundary

Class calculations have not been launched because exact row identity does not
by itself define an exact map into the current v9 pre-birth state model.  The
audited rows drive independent cleavage and emission EXP-floor surfaces,
Peierls/Taylor transport surfaces, Taylor correlation, persistent source sites,
and separate signed mobile/retained populations.  Current v9 instead drives
continuum plastic strain and scalar density through scaled copies of one common
barrier surface.  Directly assigning only the similarly named constants would
discard active audited coordinates and would not be parameter transfer.

A physics decision is required between two architectures:

1. Port the audited signed persistent-site MPZ/emission/transport state through
   stable crack birth, leaving every post-stable process absent; or
2. Define and audit a reduced v9 constitutive projection, including explicit
   rules for mobile/retained state, Taylor correlation, source-site density,
   encounter efficiency, and blunting.

Until one is selected, four-class numerical qualification is intentionally
fail-closed.  No class-specific fitting or post-stable growth work was done.
