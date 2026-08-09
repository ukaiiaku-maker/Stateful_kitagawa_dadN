# V9 physical-adapter audit — fail-closed gate result

## Branch and scope

- Branch: `codex/sn-endurance-limit-v9`
- Implementation base: merge `b93a016` (the final audit/implementation commits are reported in the handoff message).
- Frozen v8.7 physics files were not edited.
- No K360 continuation, K360 rerun, S-N pilot, optimization, sweep, `1e12` run, or da/dN work was performed.

This gate produced a strict checkpoint validator/converter and production array codec. It also found a real-physics block-partition failure that prevents claiming a completed v9 physical integrator. The gate therefore stops **failed closed**, as required.

## K360 authoritative checkpoint provenance

The exact documented path was checked:

```text
/Volumes/Data/Data/Nanopillar_calculation/stateful_pd_kitagawa_production_v2_persistent_registry/runs/stateful_pd_kitagawa_v2_8_gated_pilot/K360_CENSOR/solver_output/shielded/sigmaA_690p443MPa/
```

All expected solver files were absent:

```text
checkpoint_latest.npz
pd_state_final.npz
run_args.json
sn_stateful_pd_history.csv
summary.json
```

The enclosing `REQUEST_V2_8.json`, `EXIT_V2_8.json`, `command.txt`, and `run.log` were also absent. Nothing was recreated. Consequently the requested K360 hash/key/metadata checks (`version=4`, model/source/signature, RNG, `cycles=1e8`, `next_block=190`) cannot be performed locally.

**K360 eligibility decision: not eligible for exact v9 continuation.** The authoritative physical/stochastic checkpoint is missing.

## Native v8.7 save/load audit

`_save_case_checkpoint` atomically writes seven evolving FEM/geometry arrays, every NumPy array in `StatefulPDState`, and scalar/RNG/request/history metadata. `_load_case_checkpoint` validates checkpoint version, model, source hashes and complete signature; it regenerates mesh/PD structure from the immutable request, overwrites evolving nodes/root, rebuilds geometry caches, initializes a state object only as a container, replaces every stored field, and restores both RNG states.

The new `v9_v87_adapter.py` accepts the same exact schema only. It rejects missing/unknown NPZ or metadata fields, wrong version/model/source/signature, missing thresholds/RNG/physical arrays, invalid cycle/block coordinates, geometry/shape inconsistencies, invalid statuses and inconsistent site indices. It opens NPZ with `allow_pickle=False` and makes imported arrays read-only.

### Field mapping

| v8.7 checkpoint field | Physical meaning | v9 capsule/array field | Category | Proof required |
|---|---|---|---|---|
| `mesh_nodes` | evolving FEM geometry | array + geometry identity | A exact | bitwise loader/adapter comparison; finite `(nn,2)` |
| `root_xy` | current feature root | array + geometry identity | A exact | bitwise; finite `(2,)` |
| `ep_gp` | FEM plastic-strain tensor components | FEM/plastic arrays | A exact | bitwise, `(3,ne)` |
| `rho_gp` | integration-point density | rho/back-stress state | A exact | bitwise, `(ne,)` |
| `epsp_acc_gp` | accumulated equivalent plastic strain | FEM/plastic arrays | A exact | bitwise, `(ne,)` |
| `u` | accepted displacement | FEM state | A exact | bitwise, `(2*nn,)` |
| `last_residual` | accepted residual-stress projection | FEM state | A exact | bitwise, `(nn,)` |
| `pd__available` | mean-field available fraction | embryo/stable state array | A exact | bitwise, point shape |
| `pd__embryo` | mean-field embryo fraction | embryo/stable state array | A exact | bitwise, point shape |
| `pd__stable` | mean-field stable fraction | embryo/stable state array | A exact | bitwise, point shape |
| `pd__inactive` | mean-field inactive fraction | embryo/stable state array | A exact | bitwise, point shape |
| `pd__candidate_sites` | realized candidate count/node | candidate state array | A exact | integer point shape/conservation |
| `pd__available_sites` | available realized count/node | candidate state array | A exact | integer point shape/conservation |
| `pd__embryo_sites` | embryo realized count/node | candidate state array | A exact | integer point shape/conservation |
| `pd__stable_sites` | stable realized count/node | candidate state array | A exact | integer point shape/conservation |
| `pd__inactive_sites` | inactive realized count/node | candidate state array | A exact | integer point shape/conservation |
| `pd__born_sites_cumulative` | realized births/node | candidate state array | A exact | bitwise point shape |
| `pd__healed_sites_cumulative` | realized healing/node | transition state array | A exact | bitwise point shape |
| `pd__site_node_index` | persistent site identity→node | candidate identity and v9 clocks | A exact | bounds and site-shape checks |
| `pd__site_status` | available/embryo/stable/inactive code | candidate status and clock active flag | A exact | codes restricted to 0–3 |
| `pd__site_birth_threshold` | persistent Exp(1) threshold | `PersistentClock.threshold` | A exact | positive finite; never reconstructed |
| `pd__site_birth_cycle` | realized birth event time | candidate event state | A exact | bitwise including NaNs |
| `pd__site_stable_cycle` | realized stabilization time | candidate event state | A exact | bitwise including NaNs |
| `pd__birth_cumulative_hazard` | cumulative birth hazard/node | node array and per-site clock lookup | A exact | finite nonnegative, node shape |
| `pd__delivery_memory` | Lambda/node | delivery-memory array | A exact | bitwise point shape |
| `pd__completion` | Q(K,Lambda)/node | completion array | A exact | bitwise point shape |
| `pd__growth` | stable-defect growth state | progression array | A exact | bitwise point shape |
| `pd__crack_normal_c2` | crack-normal director cosine | PD state array | A exact | bitwise point shape |
| `pd__crack_normal_s2` | crack-normal director sine | PD state array | A exact | bitwise point shape |
| `pd__crack_orientation_weight` | director confidence | PD state array | A exact | bitwise point shape |
| `pd__bond_damage` | cohesive/PD damage | PD damage array | A exact | bitwise bond shape |
| `pd__primary_seed_rejected` | rejected primary sites | front state array | A exact | bitwise point shape |
| `pd__active_front_contact_xy` | front/surface contact | front state array | A exact | bitwise `(2,)` |
| `pd__active_front_tip_xy` | active tip | front state array | A exact | bitwise `(2,)` |
| `pd__active_front_bonds` | selected crack mask | front state array | A exact | bitwise bond shape |
| `pd__front_backbone_bonds` | backbone mask | front state array | A exact | bitwise bond shape |
| `pd__front_wake_bonds` | wake mask | front state array | A exact | bitwise bond shape |
| `pd__front_process_bonds` | process-zone mask | front state array | A exact | bitwise bond shape |
| `pd__active_front_path_xy` | ordered active path | front state array | A exact | bitwise `(n,2)` |
| `pd__healed_cumulative` | mean-field healed amount | transition array | A exact | bitwise point shape |
| `pd__born_cumulative` | mean-field births | transition array | A exact | bitwise point shape |
| `metadata.checkpoint_version` | schema identity | adapter validation | A exact | equals 4 |
| `metadata.model_id` | frozen model identity | solver identity | A exact | equals frozen `MODEL_ID` |
| `metadata.source_sha256` | driver/PD identities | source hashes | A exact | exact two v8.7 hashes |
| `metadata.signature` | immutable physical request | request hash/source identity | A exact | complete dictionary equality |
| `metadata.next_block` | native controller boundary | provenance; not reused as v9 policy | A exact | integer expected boundary |
| `metadata.cycles` | accepted cycle coordinate | `capsule.cycle` | A exact | finite expected value |
| `metadata.Wp_total` | accumulated plastic work | FEM/plastic metadata | A exact | exact scalar |
| every `metadata.pd_scalars` field | event times/front counters/scalars | embryo/front JSON metadata | A exact | exact schema/value comparison |
| `candidate_rng_state` | candidate RNG stream | RNG state/digest | A exact | valid bit-generator state; never reconstructed |
| `event_rng_state` | event RNG stream | RNG state/digest | A exact | valid bit-generator state; never reconstructed |
| `metadata.rows` | accepted diagnostic history | historical metadata | A exact | exact JSON-equivalent sequence |
| mesh elements/connectivity/boundaries | immutable structural topology, not stored | regenerated structural cache | B deterministic | reproduce native loader construction from exact signature |
| feature indices, PD neighborhoods/bonds/weights | deterministic derived structure | regenerated adapter cache | B deterministic | exact request/source and geometry-cache equivalence |
| `controller_next_block_cycles` | new adaptive proposal | v9 controller | C new | labeled new; no physical/RNG mutation |
| commit/rejection counters | v9 transaction bookkeeping | v9 controller | C new | initialized at boundary; accounting invariants |
| generation ID/manifest hashes | atomic storage identity | v9 generation metadata | C new | component hash verification |

All other `StatefulPDState` scalar fields are required by exact name in `metadata.pd_scalars`; missing or unknown names are rejected.

## Production array codec

`v9_array_codec.py` implements generation directories containing:

```text
state_arrays.npz
state_metadata.json
summary.json
manifest.json
```

Arrays are contiguous, non-object, and written with `allow_pickle=False` on load. Metadata and manifest record stable dtype/shape schema version `V9_ARRAY_GENERATION_1`. Each component has SHA-256 and byte length. Files are completed and fsynced before the temporary directory is atomically renamed; `ACTIVE.json` is atomically published last. Loading verifies every component hash/size, schema, key, dtype and shape. Missing-manifest generations are rejected/ignored and component corruption is detected.

## Shortened real checkpoint and round-trip equivalence

The earlier accepted one-cycle physical smoke checkpoint was used read-only:

```text
/private/tmp/stateful_pd_v9_smoke/shielded/sigmaA_700MPa/checkpoint_latest.npz
SHA-256 d7001eb0d87525ba70a3d38fb17329698fa2da8121a76bd35ed1dd25fc900ace
cycles=1, next_block=1, 42 NPZ entries, 3328 persistent sites
```

Native `_load_case_checkpoint` and the adapter were applied independently. All 41 stored array fields matched bitwise (`equal_nan=True` for never-realized event-time arrays); candidate and event RNG states matched exactly. The adapter round-trip writer preserved every array and metadata value.

Path A resumed the original checkpoint natively for one controlled cycle. Path B converted it to a v9 capsule, round-tripped only historical fields to a v8.7-compatible checkpoint, then used the same native continuation. At `cycles=2`, `next_block=2`:

- all 41 arrays matched bitwise including FEM, density, displacement, residual, Lambda, Q, statuses, thresholds, cumulative hazards, embryo/stable state, damage/front state;
- `Wp_total`, all PD scalars, both RNG states, rows, signature and source hashes matched exactly;
- event sequences were empty and identical.

This proves checkpoint mapping/round-trip equivalence. It does **not** prove a v9 log-domain physical continuation, because both continuation paths still use frozen native v8.7 physics.

## Real-physics block-partition diagnostic — failed gate

The shortened real case was run from the same initial request to `N=2` with fixed candidate blocks 1.0 and 0.05 cycles (>10x separation). Both reached the same endpoint, and thresholds, site statuses, growth, damage, candidate RNG and event RNG were exact. Continuous results were:

| State | Max absolute difference | Relative difference |
|---|---:|---:|
| `ep_gp` | `1.473e-21` | `2.322e-10` |
| `rho_gp` | `2.441e-3 m^-2` | `2.441e-15` |
| `epsp_acc_gp` | `1.533e-21` | `2.267e-10` |
| displacement | `1.271e-21 m` | `1.791e-16` |
| residual state | `7.805e-11 Pa` | `3.683e-10` |
| delivery memory Lambda | `3.133e-7` | **`0.587`** |
| completion Q | `1.175e-13` | **`0.825`** |
| cumulative birth hazard | `1.135e-14` | **`0.965`** |

The large relative memory/completion/hazard differences occur at tiny absolute values but directly control endurance asymptotics. They cannot be waived as harmless. Native v8.7 fixed-block integration is therefore not an acceptable v9 physical trajectory or partition-invariance reference.

No v9 real-physics partition result is claimed. The physical adapter must replace this integration with log-domain hazard/memory quadrature plus embedded error bounds.

## Real restart result

The A/B round-trip native continuation proves exact adapter state preservation across terminate/reload. The earlier synthetic transactional engine proves v9 controller restart. A combined real-physics v9 restart proof is **not available**, because the real v9 physical integrator failed the prerequisite partition gate and has not been implemented. This remains a blocker.

## Log-rate and embedded error status

The existing v9 log-rate module traces completed delivery, cleavage, completion-gated birth, stabilization, healing, growth and linkage without exponent clipping or `1e-300` physical floors. The new checkpoint adapter preserves their required physical state. However, the actual shortened continuation still invokes native v8.7 rate/integration code.

Therefore these requested items remain incomplete and are not misrepresented as passing:

- log-domain physical integration of all seven rate paths;
- stable accumulation of ultra-small hazards without conversion to zero/floor;
- embedded errors for FEM/plastic, density/back stress, Lambda, birth and progression hazards;
- localization of density/cap/governing-law boundaries before commit;
- real v9 event-driven partition and restart equivalence.

The generic v9 engine already rejects/subdivides based on model errors and reports cap activations, but a real physical model has not yet supplied those estimates/crossing functions.

## Final gate decision

| Requirement | Result |
|---|---|
| Exact K360 checkpoint found/validated | **FAIL — path absent** |
| Fail-closed v8.7 checkpoint validator/converter | PASS |
| Historical/deterministic/new-controller categories | PASS |
| Production array generation codec | PASS |
| Short real checkpoint mapping/round-trip | PASS, bitwise |
| Native A/B resume after adapter round-trip | PASS, bitwise |
| Real v9 log-domain physical continuation | **INCOMPLETE** |
| Real >10x block-partition invariance | **FAIL — native diagnostic exposes Lambda/Q/hazard sensitivity** |
| Real v9 restart equivalence | **BLOCKED by physical integrator** |
| K360 eligible for continuation | **NO** |
| S-N pilot eligible | **NO** |

## Remaining barriers to an S-N pilot

1. Implement a new physical `TransactionModel` rather than delegating accepted evolution to v8.7.
2. Integrate delivery memory and all hazards in log space with embedded quadrature/state error estimates.
3. Localize all state/cap law changes and persistent stochastic events.
4. Pass real >10x partition invariance at tolerances justified for tiny hazards, not only bulk FEM fields.
5. Pass real atomic restart with identical RNG/event sequence.
6. Recover an authoritative K360 checkpoint or accept that K360 cannot be continued exactly.

No further scientific run is authorized by this audit.
