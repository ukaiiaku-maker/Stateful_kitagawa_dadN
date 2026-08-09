# Endurance asymptotic hardening audit — 2026-08-09

## Scope and conclusions

This gate changed no v8.7 physics and ran no production solver, optimizer, or S-N pilot. It hardened finite-tail interpretation, independently reproduced the committed K360 tail, derived the completion-gated asymptotics, separated birth from physical handoff, and audited large-N floors/caps.

The main conclusions are:

1. The K360 censor provides strong **empirical** evidence for an integrable aggregate birth-hazard segment through `1e8` cycles. It is not a physics proof of endurance.
2. The K=2 completion-gated architecture is structurally capable of finite cumulative birth hazard when cyclic delivery decays sufficiently quickly.
3. The frozen numerical implementation gives Arrhenius mechanisms strictly positive rates through exponent clipping and finite state caps. Taken literally to infinite cycles, those lower tails can leave positive delivery memory/completion and hence a divergent birth hazard. This is a numerical/constitutive asymptote that must not be mistaken for physical no-endurance without an underflow-independent rate treatment.
4. Birth endurance and physical-handoff endurance are distinct. Finite total birth hazard is sufficient for a nonzero no-handoff population; divergent birth hazard alone does not prove eventual handoff.

## Hardened classification contract

`arrhenius_fracture.endurance` now reports separately:

```text
empirical_tail_evidence:
  integrable_tail_observed | divergent_tail_observed | ambiguous_tail_observed

physics_asymptotic_classification:
  endurance_supported | no_endurance | undetermined
```

Finite-window constant/exponential/power-law fits use nested 50/65/80 percent windows, fit-quality and dynamic-range checks. Their finite extrapolation is named `empirical_extrapolated_remaining_integrated_hazard`; the corresponding survival value is explicitly `empirical_extrapolated_survival_asymptote`.

An empirical fit alone always leaves the physics classification `undetermined`. `endurance_supported` requires a supplied model-derived finite remaining-hazard upper bound. `no_endurance` requires a supplied divergence proof. The sole exception is a test whose analytical generating law is explicitly marked exact and closed. Zero/nonpositive histories remain ambiguous because clipping or underflow cannot manufacture endurance.

`EndpointEnduranceDiagnostics` contains both `birth_endurance_diagnostic` and `physical_handoff_endurance_classification`. It propagates finite birth hazard to handoff endurance, but does not propagate divergent birth hazard to no-handoff endurance without proof for the full stabilization/growth/linkage process.

Adversarial coverage includes an exponential plus hidden positive floor, `p=1.5` plus hidden floor, represented late crossover from `p>1` to `p=0.8`, insufficient dynamic range, and a zero/clipped tail.

## Independently reproduced K360 metrics

The repeatable implementation is `scripts/analyze_k360_endurance_tail.py`; its regression compares computed values against the committed JSON audit. It reads 190 interval rows and reproduces:

```text
cumulative expected births at runout = 0.6223159419460373

tail fraction  p                   R^2                 extrapolated remaining integral
0.50           2.174184685808068   0.999066470078393   0.004121774905290739
0.65           2.249218577541490   0.999855752947341   0.003641359511422709
0.80           2.296280533132353   0.999947775777379   0.003429157811243042

naive conservative fit-based S(infinity) = 0.5344924270789695
```

Thus:

```text
empirical_tail_evidence = integrable_tail_observed
physics_asymptotic_classification = undetermined
physical_handoff_endurance_classification = undetermined
```

The last label is additionally unresolved because this history measures expected births, not the complete probability of stabilization, persistent growth, connected damage, and physical handoff.

## K=2 completion-gated asymptotics

At a site, v8.7 implements

```text
dLambda/dt = r_delivery(t) - Lambda/tau
h_birth(t) = h_cleave(t) Q(2, Lambda(t))
Q(2, Lambda) = gammainc(2, Lambda) ~ Lambda^2 / 2  as Lambda -> 0.
```

With fixed cycle frequency, time and cycle number differ only by a constant, so integrability exponents are unchanged. The stable memory solution is the exponential convolution

```text
Lambda(t) = exp(-t/tau) Lambda(0)
          + integral_0^t exp(-(t-s)/tau) r_delivery(s) ds.
```

For nonnegative delivery, this gives the following sufficient conditions.

### Positive delivery asymptote

If `r_delivery -> r_inf > 0`, then `Lambda -> tau r_inf` (or a positive periodic orbit for phase-resolved cyclic delivery). If also `h_cleave -> c > 0`, `Q(2,Lambda)` has a positive asymptote and `integral h_birth dN` diverges. There is no birth endurance.

### Delivery shuts off or decays exponentially

If delivery becomes zero after finite time, memory subsequently decays as `exp(-t/tau)` and `Q ~ exp(-2t/tau)`; with bounded cleavage hazard the remaining birth hazard is integrable. The same is true when delivery itself has a sufficiently fast exponential bound. This supports a nonzero no-birth—and therefore no-handoff—population.

### Algebraic delivery

For delivery varying slowly compared with finite `tau`, `Lambda(N) ~ tau r_delivery(N)`. If

```text
r_delivery(N) = O(N^-alpha)
h_cleave(N)   = O(N^beta)
```

then

```text
h_birth(N) = O(N^(-2 alpha + beta)).
```

A sufficient finite-hazard condition is `2 alpha - beta > 1`; divergence is expected when a matching positive lower asymptote has `2 alpha - beta <= 1`. For bounded nonvanishing cleavage (`beta=0`), K=2 therefore changes the activity threshold to:

```text
alpha > 1/2  -> finite cumulative birth hazard (sufficient)
alpha <= 1/2 -> divergent when corresponding lower bounds hold.
```

The K360 observed aggregate exponent near `2.17–2.30` is consistent with an effective delivery exponent near `1.09–1.15` if cleavage is asymptotically bounded and completion dominates. That inference is empirical and is not a bound on unobserved cycles.

### Code mapping

- Phase-resolved delivery source: driver `phase_resolved_delivery_rate`, around lines 174–205.
- Cleavage barrier: driver `ScratchExpFloorBarrier`, around lines 83–119.
- Exact finite-memory phase update and `gammainc(K,Lambda)`: PD module around lines 1774–1905.
- Long-block periodic-memory acceleration: PD module around lines 1907–1980.
- Cleavage rate and completion-gated birth rate: PD module around lines 2000–2065.
- Persistent exponential thresholds/cumulative hazard: PD module around lines 1430–1590.

## Structural capability versus frozen implementation

The differential-equation **form** is capable of true birth endurance: it has no required direct birth channel after delivery vanishes, and K=2 quadratically suppresses small delivery memory.

The exact frozen implementation, however, evaluates plastic delivery mechanisms using EXP-floor barriers and `exp(clip(exponent, -700, 0))`. Series-rate denominators are guarded at `1e-300`, and density/stress-related states are finitely capped/floored. With positive prefactors and delivery scale, this produces a strictly positive, though potentially fantastically small, completed-flow delivery rate at every finite capped state. Cleavage uses the same `-700` exponent lower clip and is also strictly positive. If these numerical tails are interpreted literally for infinite N, memory approaches a positive orbit and cumulative birth hazard diverges.

Therefore:

- structural completion-gated physics: **capable of true endurance** under a physically justified vanishing-delivery bound;
- K360 observed birth tail: **strong empirical integrable-tail evidence**;
- frozen implementation taken literally beyond all resolved scales: **numerically biased toward no birth endurance** by positive Arrhenius tails;
- physically justified asymptotic classification of K360: **undetermined**, pending a dimensionally and numerically defensible treatment of rates far below resolvable physical activity.

Removing the exponent clip or allowing floating-point underflow to zero would not solve this: numerical zero would manufacture endurance. v9 must carry log-rates/asymptotic bounds and distinguish physical process extinction from representational underflow.

## Birth versus physical-handoff endurance

The endpoint chain is:

```text
delivery + cleavage -> birth -> embryo stabilization/healing
                    -> stable growth -> bond linkage/front capture
                    -> physical_handoff
```

If total birth hazard is finite, exponential first-passage thresholds give a nonzero probability of no birth; no birth implies no handoff. Thus finite birth hazard is sufficient for an endurance population with respect to handoff.

If birth hazard diverges, births eventually occur under the site-clock assumptions, but embryos can heal and be returned/inactivated; stabilization, stable growth, link activity, morphology acceptance, and front capture all evolve. Divergent births alone are insufficient to prove eventual physical handoff. v9 must record separate instantaneous and cumulative hazards/activities for birth, stabilization/healing, growth, linkage, and final topology.

## Floors and caps affecting the asymptote

| Item | Class | Code behavior | Can change finite/divergent class? |
|---|---|---|---|
| `crack_floor_frac`, `max(1e-4, ...)`, 0.95 fraction ceiling | Physical constitutive assumption with fixed numerical minimum | Bounds cleavage activation energy at high stress; finite zero-stress `G0` still gives positive thermal cleavage | Yes indirectly; guarantees cleavage does not vanish through barrier growth |
| Arrhenius exponent clip `[-700,0]` in cleavage and plastic mechanisms | Numerical safety bound | Prevents underflow and imposes a strictly positive represented rate | **Yes**; can convert a physically vanishing delivery/rate into a positive asymptote and divergent hazard |
| `G >= 1e-12 eV`, `G0 >= 1e-12`, `sigc >= 1e6 Pa`, shear-modulus ratio >= 0.05 | Numerical/constitutive safety bounds | Prevent invalid barriers/scales | Can bound rates away from intended asymptotic behavior; audit in log-rate form |
| Plastic EXP-floor parameters (`Gfloor_fraction`, `Gfloor_min_eV`, `Gfloor_max_fraction`) | Physical constitutive assumptions | Positive activation rates remain at all finite stress | Yes through delivery asymptote |
| Series-rate `maximum(rate,1e-300)` and inverse denominator guard | Numerical safety bound | Prevents divide-by-zero | **Yes** if treated as a physical lower rate; must remain algebraic/log-domain only |
| `rho_floor`, `rho_cap`; local `rho >= 1e6` | Physical state-domain / numerical bounds | Clips dislocation density after every accepted block | Yes; forces a bounded state and can prevent continued decay/growth that would change delivery exponent |
| `sigma_back_max_GPa * P`, with bounded `P` | Physical constitutive saturation | Caps shielding/back stress | Yes; may leave positive driving stress/rates or prevent unbounded shielding |
| `pd_amplification_cap`, `front_amplification_cap`, amplification floor 0.25 | Finite-resolution/local-coupling bound | Bounds PD/FEM traction amplification/unloading | Yes through bounded cleavage/link driving; the 0.25 floor can prevent complete local unloading in the amplification ratio |
| `front_link_state_shift_z_clip` | Physical/numerical constitutive cap | Bounds barrier shift applied to linkage | Yes for post-birth progression/handoff, not initial birth directly |
| `max_dep_phase`, `max_rho_rel_phase`, `target_dep_eq_block`, `target_rho_rel_block` | Numerical integration controls | Limit proposed block/phase increments, not the constitutive per-cycle solution | Should not change the class if convergence/block invariance holds; can do so if clipping rather than subdivision occurs, hence mandatory invariance tests |
| `rho` clipping after `drho_cycle*dN` | Numerical state projection | Clips accepted state rather than rejecting/subdividing | Yes at cap/floor activation; v9 must surface activation and use event/subdivision semantics |
| `delivery_rate_cap_s` | Optional physical/numerical cap (default infinity) | Upper-bounds delivery | Cannot create endurance by itself; can alter transient/tail scale and exponent |
| `delivery_scale >= 0`, negative delivery/rates projected to zero | Constitutive switch/domain guard | Exact zero scale disables activity; negative artifacts are zeroed | Yes if activated unintentionally; exact zero must be identified as an explicit mechanism choice, never an inferred floor |
| `birth_scale`, pre/post-capture birth scales clipped/selected; post-capture default zero | Model competition assumption | Suppresses further births after selected crack states | Yes for post-capture birth counting, but handoff may already be governed by progression; must not be used as endurance evidence |
| Hazard increments/rates projected with `maximum(...,0)` | Numerical domain guard | Removes negative numerical values | Could manufacture zero hazard if a positive physical value is lost; activation counters required |
| Probability exponent clip `[0,700]` | Numerical safety bound | Saturates transition probability near one | Does not create endurance; may affect crossing accuracy unless persistent clocks/subdivision dominate |
| `1e-300` crossing/transition denominators and `1e-30` time/frequency/strain guards | Numerical safety bounds | Avoid division by zero | Not physical rates, but can alter asymptotics if allowed to enter model quantities |
| 512-cycle delivery transient cutoff and `1e-9/1e-11` periodic-memory tolerances | Numerical acceleration | Switches to periodic/Simpson long-block approximation | Can bias tiny tails; requires error bounds and partition invariance before v9 large-N use |
| Geometry radius/area limits, ALE taper/freeze | Finite-domain/model-resolution bounds | Invalidates, freezes, or terminates unresolved geometry | Cannot establish endurance; geometry-limited results must remain invalid/provisional |
| Front/growth/link rate and state clipping to `[0,1]` | Probability/state-domain bounds | Bounds damage and progression state | Can affect handoff asymptote; must be tested separately from birth endurance |

No finite runout, rate underflow, plotting `log10(max(rate,1e-300))`, or output precision is admissible as an endurance criterion.

## Required v9 observables and architecture

Each accepted interval must preserve and output:

```text
cleavage log-rate/hazard
delivery log-rate and physical activity units
Lambda and Q(K,Lambda)
aggregate/per-site birth hazard, cumulative hazard, thresholds and margins
stabilization and healing hazards
stable-growth and linkage/front activities
birth_endurance_diagnostic
physical_handoff_endurance_classification
all floor/cap activation counters
```

The proposed solver remains transactional:

1. immutable v8.7 adapter plus versioned state capsule;
2. log-domain rates and asymptotic bounds, never underflow-defined zero;
3. propose/evaluate/error-bound/accept loop with logarithmically growing blocks;
4. reject/subdivide on state, memory, hazard, or cap-crossing error;
5. locate earliest persistent-threshold crossing inside the block and truncate;
6. commit physical and stochastic state exactly once per accepted interval;
7. atomic checkpoint, summary, and manifest generations with request/source hashes;
8. prove block-partition and interruption/restart equivalence before mechanisms or pilot.

The baseline mechanism should be named `completion_gated_independent_cleavage`. A separately selectable no-endurance comparison may be `persistent_direct_cleavage`: a physically justified nonzero direct thermally activated initiation channel that survives delivery shakedown. It must expose units, rate law and hazard decomposition; it must not be an arbitrary stress cutoff or damage-per-cycle constant.

## Gate status and unresolved proof obligations

This asymptotic-analysis gate is complete when the full regression suite passes. Before any pilot, v9 still needs:

- a physical/log-domain decision about ultra-small delivery rates and what constitutes genuine activity extinction;
- a model-derived bound on the unobserved K360 delivery/cleavage tail;
- real block-partition and restart-equivalence tests;
- separate progression/handoff asymptotic diagnostics;
- independent evidence for any persistent-direct-cleavage comparison channel.
