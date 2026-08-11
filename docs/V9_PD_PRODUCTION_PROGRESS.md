# V9 PD production progress

## Current capability

- Synthetic milestone `f7c8fb1`: protected-state dormant engine, periodic and
  projective modes, event guard, atomic mode diagnostics.
- Historical finite-life/front qualification is complete. The production
  repair now includes monotone-ledger closure, log-domain event calculations,
  field-specific active coordinates, and an adaptive exact private-window
  training operator.
- Focused kernel qualification includes formal `1e300` log-domain skipping,
  exact-window partition guarding, and persisted restart equivalence.
- A real cached-FEM/spatial-PD exact-cycle callback now executes through the
  production phase ordering. A 10-cycle direct/callback comparison agrees in
  FEM state to roundoff, PD memory/population state to about `1e-12` relative,
  and cumulative log action to about `2e-13` relative. No threshold, RNG,
  site-status, bond, topology, geometry, or physical-cycle state changes during
  private evaluation.
- The real high-cycle orchestration path is opt-in via `--pd-high-cycle`, writes
  mode/controller histories, returns to the direct loop at an event guard, and
  uses no PD images.

## Real production status

The historical K360 finite-life/front control, 690 MPa terminal generation,
virgin 840 MPa finite-life condition, aged-state stress-step comparison, and
HC-020 current-head overlap are complete development qualifications. They are
not the final blunt-notch four-class campaign. The production target is now the
300 K `a=150 um`, `b=300 um`, `rho=600 um` half-elliptical notch, beginning
with Peak. Existing v3 and historical PD trajectories remain unchanged.

## Historical K360 real replay

The exact archived 735.921495 MPa request was replayed with `pd_image_policy =
none`. The integrated path used high-cycle qualification attempts with
transactional fallback through curved transients and reproduced the ordered
sequence at node 23:

| event | preserved control N | replay N | relative shift |
|---|---:|---:|---:|
| first embryo | 1,793,122.651 | 1,805,851.926 | 0.71% |
| stable site | 1,793,124.874 | 1,805,854.149 | 0.71% |
| first softening | 1,817,260.678 | 1,821,610.036 | 0.24% |
| root connection | 2,300,120.914 | 2,367,711.399 | 2.94% |
| front capture | 2,351,351.153 | <=2,400,000 | <=2.07% |

The replay initially exposed a sub-ULP event-time stall and loss of controller
scale across localization. Both are repaired. The corrected trajectory retains
primary seed 23 with zero reselections, reaches maximum bond damage 0.974, 18
broken bonds, root connection, and active-front capture.

## Real low-stress high-cycle qualification

The archived shielded K360 690.443200494 MPa trajectory has been preserved and
continued atomically beyond `N=1.60e7` without a realized embryo, stable site,
or topology change. Real private-map qualification exposed and repaired three
numerical-coordinate defects: spatial ledger variation was confused with time
variation; structural-zero expected populations were extrapolated as an
exponential source; and the quasistatic FEM displacement warm start was treated
as constitutive memory. Expected birth/healing ledgers are now closed by the
validated population conservation identities rather than an inaccurate
independent extrapolation.

The real engine has accepted and persisted projective segments (38 cycles in
the first fully conservation-closed qualification, with state error at most
`3.53e-6` and log-hazard validation error `1.73e-4`). This proves adoption of a
real projected physical state, but one-cycle training still limits useful
segment growth. The next numerical task is a multi-cycle exact training burst;
the direct embedded macro-stepper remains the accepted fallback and no
accelerated low-stress material result is yet claimed.

The same archived 690.443200494 MPa trajectory subsequently completed an
atomic `N=1e8` boundary (generation
`generation_3655f6ccf3044ebf8a826261d9987b23`). It remains dormant: 4,349
available realized sites, zero realized births/stable sites, zero bond damage,
and no seed/front/topology transition. The maximum nodal cumulative birth
action is `0.140255573333055`, the maximum final per-cycle birth intensity is
`6.978364374432825e-13`, and the minimum remaining persistent-site threshold
action is `9.647106699392403e-6`. This is a finite runout/checkpoint, not an
endurance classification. Expected-population diagnostics reach about 0.627
births/stable sites and remain distinct from the realized persistent clocks.

Production fallback was bounded with a private-map efficiency budget. Once
the projective fit ceased providing useful cycles per exact map, control
returned to the authoritative embedded macro-stepper and deferred another
qualification until the next logarithmic checkpoint. Mode-history writes are
append-only across subsequent atomic resumes.

An extension of the same condition reached atomic generation
`generation_c18519735d8640a8b6cad153c26158ff` at
`N=168634947.0289538`, again with zero realized births, stable sites, damage,
or topology change. The maximum cumulative action increased only from
`0.140255573333055` to `0.140256540348356` over the additional 6.86e7 cycles;
the nearest remaining threshold changed from `9.647106699392403e-6` to
`9.645864666437046e-6`. The actual final maximum birth rate is
`1.4023253737703178e-13` per cycle and the direct controller estimates a
`2.00129794714e10`-cycle next wait. This strong rate decay is scientifically
important but is not by itself a proof of endurance.

The former single-cycle projective limitation is superseded by HC-020. The
adaptive operator evaluates an exact macro-window and an independent two-part
exact partition before adoption. On the current 690 MPa head it accepted three
64-cycle windows: 192 physical cycles with nine exact maps, or 21.33 accepted
cycles per map. Maximum state partition error was `7.80e-12`; maximum
log-action error was `5.14e-5`. Direct and accelerated terminal maximum action
agree within `8e-17`, and a persisted 64+128-cycle restart is array-identical to
the continuous accelerated result. Acceleration is reported only for this
accepted result; the earlier zero-accepted-cycle trial remains a rejection.

## Current-head A/B and requested stress increase

A fresh current-head 735.921495137 MPa direct/high-cycle A/B exposed that 196
locally accepted but inefficient projected cycles shifted the first embryo by
86.85 cycles and front capture by 4,330.74 cycles. Production admission now
requires at least 16 projected cycles per exact map; sub-64-cycle proposals
fail closed pending the multi-cycle operator. The repaired high-cycle leg then
matched the direct reference exactly at every reported event and terminal
quantity: embryo `1805740.0338697084`, stable site `1805742.256830253`, first
softening `1824173.1357521766`, root connection `2322912.5057363925`, front
capture `2382710.901963831`, seed node 23, zero reselections, 16 broken bonds,
and maximum damage `0.9701612494297833`.

The accepted dormant 690.443 MPa generation was packaged with all restart
arrays, context, hashes, controller/mode history, and all 4,349 available site
clocks under `runs/v9_pd_production/packages/`. Following the requested rule,
the stress was increased by 150 MPa to `840.4432004940379 MPa` with unchanged
K360 geometry, temperature, seeds, and physics. This condition formed a stable
site: first embryo `36019.35198636835`, stable site `36021.57494691292`, first
softening `45746.295542562904`, and root/front capture
`867148.3411936606` cycles. The terminal calculation contains two realized
births/stable sites, 33 broken bonds, and maximum damage `0.956203983658257`.

The virgin condition is packaged independently as protocol `virgin_840MPa`.
The exact 690 MPa terminal generation was also restored and continued with only
the applied stress changed, protocol `aged_690MPa_then_840MPa`. The aged state
has reached `N=180000000`, or 11,365,052.971 additional cycles at 840 MPa,
without an embryo, stable site, damage, or front capture. Thus the aged
additional event times are right-censored beyond 11.365 million cycles, versus
36,019.352 cycles to embryo, 36,021.575 to stabilization, and 867,148.341 to
front capture in virgin material. Stress was not raised above 840 MPa.

Matched persistent-site tables at `N=1e8` and `N=168634947.0289538` retain all
site IDs and separately identify the minimum-action (site 683), maximum-rate
(786 then 614), minimum-wait (site 1189), and maximum-cumulative-action (site
13) clocks. The table reports the exact one-cycle log/linear birth action,
remaining action, wait, delivery rate/memory, K=2 completion, raw cleavage
rate, effective opening stress, equivalent-stress emission drive, backstress,
and state shift.

## Blunt-notch Peak diagnostic start

The target geometry is active at 300 K with `a=150 um`, `b=300 um`,
`rho=600 um`, the audited Peak row, and no PD images. The 12,000 MPa overstress
topology diagnostic records first embryo at `81050.8859089888` cycles,
stabilization at `81053.2753569623`, first softening at
`92053.27535696232`, and root connection at `1528000`. At the preserved
`N=3e6` checkpoint the original seed (node 297) has zero reselections and three
broken bonds. Six-bond front capture remains right-censored above `3e6`; this
checkpoint is not represented as a completed front-capture endpoint.

This run exposed an update-partition defect in primary-seed stall detection:
small, monotone damage increments could be rejected before accumulating the
coarse diagnostic increment. Numerically resolved positive progress now resets
the consecutive stall counter, while the original coarse milestone remains the
reported progress reference. A regression exercises the slowly advancing seed
under small controller blocks.

The independent lower-stress Peak condition at 4,000 MPa was continued
logarithmically from `N=1e6` through `1e8` to `1e10`. It has no embryo, stable
site, bond damage, root connection, or front capture. On the final continuation
the HC-020 operator accepted 2,688 exact-window cycles with 130 exact-map
evaluations, qualified the periodic state, and advanced the remaining
`8,999,997,312` cycles by the stationary map. This is a preserved finite
runout, not an endurance classification.

Neither calculation is a qualified physical S-N bracket. The 4 GPa package is
now labeled `Peak_blunt_4000MPa_K2_suppressed_runout_diagnostic`; the 12 GPa
package is `Peak_blunt_12000MPa_overstress_topology_diagnostic`. The executable
bridge audit found omitted aggregate source multiplicity and, separately, a
non-parity scalar equivalent-stress projection. Corrected Peak production is
fail-closed pending selection of a root-local versus node-local signed-MPZ/PD
delivery allocation (or direct canonical-cleavage embryo attempts).

Both Peak checkpoints have hash-verified packages under
`runs/v9_pd_production/packages/`. The direct post-event 12,000 MPa package
explicitly records the absence of high-cycle controller/mode files rather than
inventing a controller history. The complete regression suite passes: 178
tests and 7 subtests.
