# Peak blunt-notch m1 attempt S-N audit v1

This versioned audit supersedes only the campaign-status statement in
`V9_PD_M1_NUMERICAL_PROBABILISTIC_QUALIFICATION.md`; it does not rewrite that
historical numerical-gate document. Commits `f8ea4ee` and `8466833` are the
accepted 1500 and 1820 MPa actual-action milestones.

The campaign is Peak at 300 K, R=0.1, 1000 Hz, an ideal smooth 600 um nominal
blunt notch, and one 10 um correlated front segment. Conditional trajectories
are analysis-only: their threshold is not a physical draw, mark RNG is not
consumed, and topology continuation is forbidden.

Production action paths and complete provenance are recorded in
`runs/sn_v9_shared_root_m1_peak/attempt_sn_v1/manifest.json`. Cumulative action
comes from accepted exact private windows with split-window state/hazard checks;
stationary propagation is used only after its independent certificate. The
quantile CSV records the actual checkpoint bracket for every interpolation.

The 1820 MPa N50 interpolation was checked by exact ordered-phase action-boundary
localization: 638586.575889 interpolated versus 638586.217030 exact cycles, an
absolute error of 0.358859 cycles and relative error 5.62e-7. No mark, embryo,
threshold renewal, or topology continuation occurred at that analysis boundary.

The local-stress acceptance reference remains approved commit `7b5837d`:
Kt_root_opening is 5.56548 for K360 and 1.68315 for the blunt notch. Same-tensor
constitutive parity, rather than nominal-stress agreement, remains the local
physics gate.
