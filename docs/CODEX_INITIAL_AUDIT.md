# Codex initial baseline audit — 2026-08-09

## Repository identity and scope

- Reference tree: `/Volumes/Data/Data/Nanopillar_calculation/Stateful_kitagawa_dadN_canonical_reconstruction` (read only; 17 files inventoried).
- Development tree: `/Volumes/Data/Data/Nanopillar_calculation/Stateful_kitagawa_dadN`.
- Required upstream branch was clean and current at `bef62ed0d4909cbd3e3dbf45f15c97055b21b8ad` before creating `codex/sn-endurance-limit-v9`.
- GitHub Issue #1 was read on 2026-08-09. It agrees with the local handoff: audit and tail mathematics precede large-N solver changes or campaigns.
- No production S-N run was launched and no frozen baseline physics was edited.

## Provenance and dependency recovery

All 17 reference files were enumerated. The reference tree contains the three v8.7 modules, v8.3 parents, v2.8.1 orchestration correction, and checksum/provenance documents; it does not contain the common Python closure.

AST traversal found this recursive local closure:

```text
v8.7 driver
├── config
├── sn_arrhenius_chain ── fatigue_v1 ── config
├── sn_feature_geometry_v8_7 ── mesh ── config
├── sn_intact_fem ── config, mesh, sn_arrhenius_chain, sn_geometry
└── stateful_peridynamics_v8_7_local_front_spacing
    └── config, sn_geometry, sn_intact_fem
```

The missing modules were recovered byte-for-byte from the historical standalone at `Fatigue-PF/dist/stateful_pd_kitagawa_production_v2_standalone`, whose `SOURCE_MANIFEST.json`, `SHA256SUMS_STANDALONE.txt`, and `BUILD_ENVIRONMENT.json` identify their original `Fatigue-PF` sources and build environment. No compatibility substitutions were required. Exact paths, origins, hashes, and edges are in `reference/BASELINE_DEPENDENCY_MANIFEST.json`.

Frozen v8.7 identities after recovery:

```text
09e993183efff4857819b4ede0d668c909e20810936c775638edbc9b2c5809d6  sn_pd2d_stateful_v8_7_generalized_features.py
bb8341a3d8ed883f605bd99f6d74b4b011837d436d913cc1914a882a133f05b3  stateful_peridynamics_v8_7_local_front_spacing.py
8ae3e8c454297fcbc61f179b747ffbd05129aa3c8598528113aa2edf8913396e  sn_feature_geometry_v8_7.py
```

The exact v8.3 parents and their two recorded regression files were also restored as audit fixtures. The v2.8.1 orchestration files were copied from `REFERENCE_ROOT`; these are not part of the recursive solver import closure.

## Environment and validation

Validated executable and packages:

```text
/opt/homebrew/Caskroom/miniconda/base/envs/arrhenius-stateful-pd/bin/python
Python 3.12.13
numpy 2.5.1
scipy 1.18.0
matplotlib 3.11.0 (historical build manifest; imported by smoke output)
```

The base Python 3.13.2 environment is unsuitable: importing `scipy.sparse.linalg` fails because `_propack` cannot provide `_spropack`. This is an environment defect, not a recovered-source defect.

Validation results in the historical environment:

- `compileall`: passed (bytecode cache redirected to `/private/tmp`).
- Imports: v8.7 driver, `StatefulPDConfig`, `StatefulPDPatch`, and generalized `make_blunt_edge_notch_mesh`: passed.
- `python -m unittest discover -s tests -p 'test_*.py' -v`: 63 passed in 5.087 s (56 frozen v8.3 core/sweep tests plus 7 endurance tests).
- `python orchestration/v2_8_1/test_v2_8_1_summary_discovery.py -v`: 4 passed in 0.005 s.

### One-block deterministic physical smoke

One custom-resolution smoke only was run under `/private/tmp/stateful_pd_v9_smoke`: seed 42, shielded case, 700 MPa, one cycle, one accepted block, 18 x 36 base grid, zero jitter, 35 micrometre root spacing, 105 micrometre PD horizon. It constructed generalized geometry and an 806-bond PD patch, advanced FEM/plastic/PD state once, and wrote a 42-array atomic checkpoint, history, final PD state, handoff audit, plots, and `summary.json`. Result: `right_censored` at `N=1`, as expected for infrastructure smoke. The custom mesh reported `delta/h=3.00` and failed the production-only active-region-clearance preflight; the solver explicitly allowed it only as smoke/debug. No morphology or geometry invalidity occurred.

## State, stochastic, checkpoint, and outcome map

- Persistent exponential thresholds are sampled once in `StatefulPDPatch.initial_state` (`stateful_peridynamics_v8_7_local_front_spacing.py`, around lines 1430–1453) into `site_birth_threshold`; cumulative site hazards start in `birth_cumulative_hazard`.
- `_next_discrete_birth_wait` computes threshold margins; `_advance_discrete_birth_clocks` (around lines 1563–1590) finds within-block crossing fractions and commits cumulative hazard.
- `StatefulPDPatch.update` (line 2178 onward) performs the physical state transition. The driver invokes it once per accepted block and advances `cycles_total` after the update. This baseline has a within-block birth clock but has not yet demonstrated full v9 block-policy/restart invariance on a real v8.7 trajectory.
- `_save_case_checkpoint` (driver lines 284–340) writes a temporary compressed NPZ then `os.replace`s it atomically. Metadata includes model, source hashes, signature, RNG states, cycle/block state, FEM variables, mesh, and every PD dataclass array. `_load_case_checkpoint` verifies version/model/source/signature before restoration.
- Final status selection (driver lines 1288–1298) distinguishes clean `physical_handoff`, geometry-saturated handoff, clean `right_censored`, geometry-saturated censor, and invalid/max-block modes. Summary discovery/reuse is corrected by v2.8.1 metadata-based lookup rather than formatted stress-directory identity.

## Large-N caps, floors, and asymptotic risks

- The independent cleavage barrier contains `crack_floor_frac` (default 0.010). This is an energy-barrier floor, not a hazard floor, but it can leave a positive residual nucleation hazard and therefore materially favors `H(infinity)=infinity`; it must be included in tail diagnostics/ablations.
- Arrhenius exponent arguments are clipped to `[-700, 0]`. The lower clip prevents underflow and leaves a tiny positive rate; it does not create zero hazard, but can impose an artificial positive asymptote at extreme states.
- Birth probability uses `1-exp(-clip(H_block, 0, 700))`; discrete first passage uses the unclipped nonnegative hazard increment. The probability clip saturates at one and does not create endurance.
- Negative computed rates/hazard increments are projected to zero with `maximum(..., 0)`. These are domain guards; any physical positive rate that becomes negative numerically would be hidden, so v9 must audit their activation.
- Crossing division uses `1e-300`; tolerances use machine epsilon. Neither is a hazard cutoff, but both need block-invariance tests.
- `rho_floor=1e8 m^-2`, `rho_cap=1e17 m^-2`, back stress maximum 1 GPa, delivery-rate cap (default infinity), PD amplification caps, and front state-shift clipping can change the long-time state and hence hazard asymptote.
- The baseline hardening/shielding is strongly one-sided in the observed K360 censor. Dislocation evolution includes bounded state and dynamic terms, while the fracture pathway lacks a demonstrated competing slow degradation/recovery mechanism. Whether the actual tail is integrable remains undetermined because no stored long-tail history was present in the reconstruction.

## Tail analyzer

`arrhenius_fracture.endurance` fits constant, exponential, and power-law positive tails over nested windows. A result is supported only when model, integrability class, and fit quality are window-stable. It reports model/parameters, window sensitivity, remaining integrated hazard (conservative maximum across finite fits), and survival asymptote. Short, noisy, clipped/nonpositive, poor-fit, model-sensitive, or near-`p=1` tails return `undetermined`. Exact `p=1` is correctly classified as divergent.

The existing K360 `1e8` censor cannot presently be classified from the reconstructed evidence: summary statements say its hazard declined and remained positive, but no sufficiently resolved hazard history was supplied. Continuing a merely positive finite endpoint is mathematically insufficient; a fitted tail and window sensitivity are required.

## Proposed v9 large-N architecture (not implemented in this deliverable)

1. Freeze v8.7 as an adapter/reference and introduce a versioned v9 state capsule containing cycle, physical fields, cumulative per-site hazards, persistent thresholds, every RNG bit-generator state, accepted-step counters, and request/source hashes.
2. Use a propose/evaluate/accept loop with logarithmically growing blocks. Bound state/hazard interpolation error; reject/subdivide a block when the bound fails.
3. Locate the earliest persistent-threshold crossing inside a proposed block, truncate to it, and commit exactly once. Make the accepted interval the sole unit of state mutation and checkpoint generation.
4. Write checkpoints and summaries atomically with generation IDs and content hashes. Resume only exact request/source identities and preserve complete stochastic state.
5. Prove block-partition invariance and interrupted/resumed equivalence before adding two explicitly named closures: frozen independent cleavage and activity-conditioned initiation based on dimensional irreversible activity (for example plastic dissipation/event delivery), with every hazard factor emitted.
6. Feed no-failure hazard/state histories to the standalone tail analyzer. Keep `undetermined` when extrapolation is not defensible. Only after these gates run the approved 3 x 3 x 2 pilot.

## Unresolved issues

- No raw K360 long-tail history or accepted K360 checkpoint is included in `REFERENCE_ROOT`; the two scientific anchors therefore remain documented but were not rerun (doing so would violate the short-run scope).
- v8.7-specific regression tests beyond import and one-block execution were not present; current focused physics tests are the exact v8.3 parent suite.
- The smoke's intentionally coarse custom mesh is not production-qualified. A short production-resolution restart/invariance fixture should be designed for v9 without rerunning the million-cycle gate.
- Summary JSON is written directly rather than with the checkpoint's temporary-file/replace protocol; v9 should make all terminal identity artifacts atomic.
