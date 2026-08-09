# K360 authoritative checkpoint recovery before v9 physical adapter

## Why this note exists

The committed `reference/k360_v2_8_tail/` directory intentionally contains only derived CSV/JSON evidence. That directory is insufficient for exact continuation.

However, the original July v2.8 gated run was performed in the persistent local registry. Before declaring K360 continuation impossible, check whether that original run directory is still present on the user's machine.

## Exact local path to check

```text
/Volumes/Data/Data/Nanopillar_calculation/stateful_pd_kitagawa_production_v2_persistent_registry/runs/stateful_pd_kitagawa_v2_8_gated_pilot/K360_CENSOR/solver_output/shielded/sigmaA_690p443MPa/
```

Expected files include:

```text
checkpoint_latest.npz
pd_state_final.npz
run_args.json
sn_stateful_pd_history.csv
summary.json
```

The enclosing gate directory should also contain:

```text
K360_CENSOR/REQUEST_V2_8.json
K360_CENSOR/EXIT_V2_8.json
K360_CENSOR/command.txt
K360_CENSOR/run.log
```

Do not recreate any of these if absent. Missing authoritative state must fail closed.

## Expected checkpoint identity from the preserved July artifact

A preserved copy of the original censor checkpoint has these properties:

```text
cycles = 100000000.0
next_block = 190
```

Its metadata contains:

```text
checkpoint_version
model_id
source_sha256
signature
next_block
cycles
Wp_total
pd_scalars
candidate_rng_state
event_rng_state
rows
```

The recorded v8.7 source hashes are:

```text
driver:
09e993183efff4857819b4ede0d668c909e20810936c775638edbc9b2c5809d6

pd_module:
bb8341a3d8ed883f605bd99f6d74b4b011837d436d913cc1914a882a133f05b3
```

The checkpoint contains 42 NPZ entries, including the evolving FEM and PD state required by the v8.7 resume implementation. Important entries include:

```text
mesh_nodes
root_xy
ep_gp
rho_gp
epsp_acc_gp
u
last_residual
pd__site_node_index
pd__site_status
pd__site_birth_threshold
pd__site_birth_cycle
pd__site_stable_cycle
pd__birth_cumulative_hazard
pd__delivery_memory
pd__completion
pd__growth
pd__bond_damage
pd__active_front_*
pd__front_*
metadata_json
```

The metadata preserves both candidate and event RNG bit-generator states.

## Important distinction

The original v8.7 checkpoint was designed by its own save/load implementation as the authoritative state for exact **v8.7** continuation. That does not automatically make it a complete `V9StateCapsule`.

The next gate is therefore not to run K360 immediately. It is to build and test a fail-closed `v8.7 checkpoint -> v9 capsule` adapter.

The adapter must:

1. validate the checkpoint version, model ID, source hashes and complete v8.7 run signature;
2. load only authoritative fields that are actually present;
3. reconstruct only deterministic geometry/connectivity that the original v8.7 loader itself reconstructs from the same run request; never reconstruct stochastic state;
4. preserve candidate/event RNG state exactly;
5. preserve persistent site thresholds, statuses and cumulative hazards exactly;
6. map every evolving FEM/plastic/PD field into the v9 capsule;
7. identify any genuinely new v9 controller fields as newly initialized **adapter/controller state**, not historical physical state;
8. prove adapter round-trip / v8.7 resume equivalence on a shortened checkpoint before using the K360 checkpoint;
9. fail closed if any field required for physical continuation cannot be justified from the checkpoint plus immutable run identity.

## Local preflight for Codex

Run this before changing the conclusion in `V9_NUMERICAL_ARCHITECTURE_AUDIT.md`:

```bash
K360=/Volumes/Data/Data/Nanopillar_calculation/stateful_pd_kitagawa_production_v2_persistent_registry/runs/stateful_pd_kitagawa_v2_8_gated_pilot/K360_CENSOR/solver_output/shielded/sigmaA_690p443MPa

for f in \
  checkpoint_latest.npz \
  pd_state_final.npz \
  run_args.json \
  sn_stateful_pd_history.csv \
  summary.json; do
  if test -f "$K360/$f"; then
    shasum -a 256 "$K360/$f"
  else
    echo "MISSING $K360/$f"
  fi
done
```

If `checkpoint_latest.npz` exists, inspect its NPZ keys and `metadata_json` before asserting insufficiency.

Do not launch continuation merely because the checkpoint exists. First complete the physical adapter, array codec, block-partition proof and restart-equivalence proof required by the v9 numerical architecture gate.
