# Canonical source reconstruction — 2026-08-09

## Decision
The active development baseline is **Stateful-PD v8.7 generalized features + local-front spacing**, with the **v2.8.1 gated-pilot/summary-discovery orchestration**.

## Evidence
The user-provided local archive `Archive(20260809-160146).zip` contains three top-level packages:

- `stateful_pd_v8_5_standalone_builder_v1_1/`
- `stateful_pd_kitagawa_v2_8_gated_pilot/`
- `stateful_pd_kitagawa_v2_8_1_summary_discovery_hotfix/`

It does **not** contain the v8.7 core solver source files. However, both v2.8 runners explicitly invoke `arrhenius_fracture.sn_pd2d_stateful_v8_7_generalized_features`, confirming that v8.7 was the later solver used by the latest gated-pilot work.

The exact v8.7 source is reproducibly reconstructed from the exact validated v8.3 baseline and the v2.6 add-only installer transformation preserved in the project artifacts. Reconstructed hashes match the July runtime source hashes exactly:

```text
09e993183efff4857819b4ede0d668c909e20810936c775638edbc9b2c5809d6  arrhenius_fracture/sn_pd2d_stateful_v8_7_generalized_features.py
bb8341a3d8ed883f605bd99f6d74b4b011837d436d913cc1914a882a133f05b3  arrhenius_fracture/stateful_peridynamics_v8_7_local_front_spacing.py
8ae3e8c454297fcbc61f179b747ffbd05129aa3c8598528113aa2edf8913396e  arrhenius_fracture/sn_feature_geometry_v8_7.py
```

Exact v8.3 parent hashes:

```text
10286f2d659c2253df33fd84c7514e22bfb4b651dcbe8719ce2d93a76ae1fbf8  sn_pd2d_stateful_v8_3.py
c08e79213266bdb1d975f0e21435b547be6dd603fd412b3633e1c49c3d8b08b3  stateful_peridynamics_v8_3.py
```

## Important limitation
The local archive is not a complete standalone solver checkout. Before the repository is declared runnable, restore and verify the exact common dependency closure used by the v8.7 standalone (`config.py`, `mesh.py`, `sn_arrhenius_chain.py`, `sn_geometry.py`, `sn_intact_fem.py`, and any transitive imports). Do not silently substitute unrelated newer implementations.

## Historical v8.5
The local v8.5 standalone builder predates the v8.7 generalized-feature/local-front-spacing solver and should be retained only as historical/reference material, not as the active solver baseline.
