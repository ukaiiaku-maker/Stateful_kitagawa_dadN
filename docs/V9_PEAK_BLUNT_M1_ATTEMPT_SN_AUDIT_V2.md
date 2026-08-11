# Peak blunt-notch m1 attempt S-N audit v2

This audit preserves and supersedes the campaign-status portions of V1. Commit
`edb32a2` is the approved seven-stress global-attempt S-N skeleton; commit
`b956d88` is its approved exact 1500 MPa H=ln(2) analysis boundary.

The seven actual-action conditions are 1172, 1290, 1460, 1500, 1630, 1820,
and 2000 MPa. Paths, action support, accepted numerical modes, material/source
hashes, atomic generation/checkpoint hashes, controller hashes, and mechanics
validity are recorded per condition in
`runs/sn_v9_shared_root_m1_peak/attempt_sn_v1/manifest.json`.

Two N50 rows have independent ordered-phase boundary localization:

- 1500 MPa: exact 3377055123.1199546 cycles versus interpolated
  3377055123.7406697 cycles; absolute error 0.6207151 cycle and relative error
  1.838e-10.
- 1820 MPa: exact 638586.2170300657 cycles versus interpolated
  638586.5758892479 cycles; absolute error 0.3588592 cycle and relative error
  5.620e-7.

The production V2 attempt table retains both exact and interpolated values and
uses `exact_ordered_phase_action_boundary_localization` for these rows. These
are elementary global-attempt distributions, not stable-seed or front-capture
distributions.
