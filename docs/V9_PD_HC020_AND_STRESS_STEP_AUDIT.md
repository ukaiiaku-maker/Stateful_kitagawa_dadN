# V9 PD HC-020 and K360 stress-step audit

## Scope

This is an intermediate numerical/physics qualification. It preserves the
K360 calculations and does not classify them as the final `a=150 um`, `b=300
um`, `rho=600 um` four-class campaign.

## Packaged virgin protocol

`virgin_840MPa` is packaged at generation
`generation_744b3f0fe9614b6db832b215aa22d990`. The package manifest hashes the
atomic arrays/metadata/summary, legacy checkpoint, run arguments, event
history, final state, crack-handoff audit, and controller/mode histories. It
contains no image sequence.

Its event cycles are 36,019.351986 to first embryo, 36,021.574947 to stable
site, 45,746.295543 to first softening, and 867,148.341194 to both root
connection and front capture.

## Aged stress-step protocol

`aged_690MPa_then_840MPa` starts from the hash-verified generation
`generation_c18519735d8640a8b6cad153c26158ff` at
`N=168634947.0289538`. Loading validates the complete original 690 MPa run
signature and accepted source fingerprint. Only `sigma_a` changes to
840.4432004940379 MPa; checkpoint arrays, site actions and thresholds, both RNG
states, site identities, delivery memory, FEM/PD state, geometry, and topology
are restored rather than reconstructed.

At `N=180000000` the aged protocol remains dormant. Additional cycles to
embryo, stabilization, and front capture are all right-censored above
11,365,052.9710462 cycles. No stress above 840.443 MPa was used.

## HC-020 operator

The production operator now trains on exact private macro-windows. Each
candidate window is compared with an independently evaluated two-part exact
partition. Adoption requires bounded active-state, cumulative log-action, and
monotone-ledger errors; unchanged transition signatures; a conservative
first-passage action guard; and the configured minimum accepted cycles per
exact map. The right partition begins at the exact physical coordinate
`N + dN_left`, rather than merely restoring the left terminal state at `N`;
a cycle-coordinate-dependent regression prevents this distinction from being
lost in nominally stationary tests.

At the current 690 MPa terminal head, three accepted 64-cycle windows advance
192 cycles with nine exact maps (21.333 cycles/map). Across the three windows:

- maximum active-state partition error: `7.795486238666265e-12`;
- maximum cumulative log-action error: `5.134393596151199e-5`;
- protected threshold/RNG/topology state: unchanged;
- direct versus accelerated maximum cumulative action:
  `0.14025654034944365` versus `0.14025654034944357`;
- continuous accelerated versus persisted 64+128 restart: every stored array
  identical, maximum absolute and relative difference zero.

The earlier trial with `accepted_projected_cycles=0` is retained as a rejected
qualification and is not called acceleration.

## Persistent-site audit

The matched table is reconstructed by loading each atomic state and evaluating
one authoritative private v9 cycle. It uses PD-local coordinates, correcting
the ambiguity of indexing a PD-local site node directly into the global FEM
node array. It reports all 4,349 persistent sites at both checkpoints and a
separate extrema JSON. The four requested extrema are not assumed to be one
site: site 683 has minimum remaining action, site 1189 minimum instantaneous
wait, site 13 maximum cumulative action, and the maximum-rate identity changes
from site 786 at `N=1e8` to site 614 at the terminal checkpoint.
