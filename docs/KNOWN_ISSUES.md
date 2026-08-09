# Known issues and lessons from prior campaign attempts

## 1. Finite runout is not an endurance-limit criterion

The K360 lower-stress gate reached `1e8` cycles without a realized birth, but its instantaneous hazard remained nonzero. The central unresolved scientific issue is the asymptotic hazard tail.

## 2. Strong hardening/shielding can create apparent endurance

The lower-stress gate evolved toward a strongly stabilized state and the instantaneous birth rate decreased by orders of magnitude. This may be physical fatigue training/shakedown, or it may be too strong because the model lacks a competing slow degradation/recovery/softening channel. Long-N diagnostics must separate these possibilities.

## 3. Embryo stabilization is very fast in the current baseline

In the K360 failure gate, embryo-to-stable evolution occurred in only a few cycles compared with a million-cycle first-birth time. The stabilization stage is therefore effectively instantaneous for this condition and should not be over-interpreted as independently calibrated physics.

## 4. Current global/local coupling is one-way

PD damage redistributes the local patch response but is not fed back into the global FEM stiffness. This is acceptable for initiation-to-handoff studies but must be documented if longer crack growth is attempted.

## 5. v8.7 changed the geometry/resolution implementation

The original v8.3 production audit used a global patch median spacing. Deep features could be rejected even when the root/front was locally resolved. v8.7 uses local front/root spacing and a refined crack corridor. This change is intentional, but regression against the K360 gate is required.

## 6. Previous optimizer wrappers failed before doing physical work

Two important orchestration failures occurred:

- an early optimizer advanced/ranked candidates even though solver jobs had failed;
- a later wrapper pre-created output directories and used `--skip-existing`, causing solver jobs to exit without producing summaries.

Campaign infrastructure must fail closed. A missing summary, model mismatch, corrupt checkpoint, geometry-invalid result, or nonzero return code is never candidate performance.

## 7. Summary discovery cannot depend on display formatting

A completed v2.8 result used a six-digit stress directory tag while the wrapper predicted a higher-precision tag. v2.8.1 fixed this by discovering summaries from their stored physical metadata. New campaign code should use request hashes / manifests, not floating-point directory names, as identity.

## 8. Long jobs require restart equivalence before production use

Any campaign expected to run for hours or days must demonstrate interruption/restart equivalence on a real short physical case before large sweeps are authorized.
