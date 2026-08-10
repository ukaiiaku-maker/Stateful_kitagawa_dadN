# v9 canonical four-class cleavage/stable-birth implementation audit

## Implemented boundary

`v9_canonical_four_class_birth.py` defines the separately versioned model
`v9_canonical_four_class_cleavage_stable_birth_v1`.  Its endpoint is
`canonical_cleavage_stable_crack_birth`: the first canonical sharp-front
cleavage first passage.  It neither imports nor calls the legacy PD
candidate-site/K=2/stabilization path, and it executes no crack-advance reward.

The state owns the exact audited signed-MPZ adapter, one root-local cleavage
threshold, cumulative physical cleavage hazard, and the canonical RNG state.
The threshold is the authoritative unit-exponential draw from
`SeedSequence([hazard_seed, engine_id])`.  No new threshold is drawn at birth.

## Constitutive and hazard parity

The raw cleavage rate comes directly from the exact registry-selected audited
manifest.  The canonical multiplicity is the sharp-front renewal law

```text
lambda_eff = gammainc(m, lambda_raw*tau_c) / tau_c
m = 3, tau_c = 1e-6 s
```

with the same `1e12` argument bound as `UnifiedMPZFrontEngine.lambda_cleave`.
Tests cover all four canonical rows at 300, 900, and 1300 K.  The signed state,
aggregate emission, Peierls/Taylor transport, backstress, shielding, blunting,
and wake coordinates remain owned by the directly loaded immutable executable.

The root-local stress transformation reproduces the canonical engine ordering:
nominal root stress is converted to nominal K using the initial root radius,
signed active shielding is subtracted in K space, and the result is divided by
the current blunted-radius factor.

## Low-rate hardening

The representable linear renewal rate is bit-for-bit the SciPy `gammainc`
expression used by the source.  An algebraically identical log-domain small-x
form is retained in parallel.  Thus an underflowed linear diagnostic is never
interpreted as a physical zero.  Cumulative hazard is accumulated with
`logaddexp`; `S_stable=exp(-H_cleave)` has the canonical stable-birth survival
semantics because there is no other stochastic pre-endpoint event.

Finite runout is still not endurance evidence.  Endurance classification must
use the asymptotic behavior of the cumulative canonical hazard.

## Transaction and restart evidence

An interval uses the canonical plastic-half-step / midpoint-cleavage /
plastic-half-step ordering.  A proposed crossing is localized by replaying
private signed-state copies; only the localized accepted state is committed.
The event consumes no post-birth RNG and leaves unused proposed time uncommitted.

The versioned capsule contains every signed active/wake array and scalar, the
threshold, cumulative linear and log hazard, time/event coordinates, and the
complete NumPy bit-generator state.  Atomic restart matches an uninterrupted
multi-load history exactly.  A 10-way partition comparison passes the signed
state and hazard convergence gate.  The complete repository suite passes 134
tests.

## Scope still to execute

The canonical endpoint state is implemented and qualified under prescribed
root histories.  It is not yet wired into the full v9 FEM cyclic driver, and no
canonical Peak/DBTT/weak-T/ceramic production condition has therefore been
claimed.  The next implementation seam is the existing accepted FEM phase
history: supply its root-local opening and signed shear history to this state,
persist the canonical capsule in the atomic generation, and stop immediately
at the localized crossing.  Post-birth PD, front capture, handoff, and da/dN
must remain disabled for that driver.
