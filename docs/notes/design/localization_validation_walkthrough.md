# The θ = 0 validation: input drains, epochs, and the localization trigger

*Companion explainer, non-normative. The spec ([§8.4][s8-4], [§8.6][s8-6]; context in [§2.1][s2-1], [§9.4][s9-4]) wins
wherever the two disagree. Written against the design as of rows 179–182: detection
policy declared by the guard's return type, honest priors under budgeted re-firing,
`localization_budget`/`firing_budget` as the two event budgets.*

This document explains one small mechanism — a single extra guard evaluation at the
start of every localization — and why it is load-bearing: it closes a correctness hole
in the localization trigger, supplies a value the root-finders always needed but the
spec never sourced, and draws the jurisdictional line between the two kinds of edges a
guard can exhibit. The mechanism is easy to state and easy to implement; what is worth
preserving is the *reasoning*, because every piece of it follows from timing rules
settled elsewhere in the spec, and a future change to any of those rules should be
checked against this chain.

## Section 1 — The problem: the prior and the probes straddle an input epoch

The ordering of operations at a grid instant tₙ is:

> boundary sequence of the completed frame → input drain → integrate the next frame

Priors — each event's record of its predicate at the last boundary — are sampled at
the boundary's *quiescence*, before the drain. Localization probes during the next
frame evaluate guards by running sweeps against the *post-drain* inputs. So the prior
and the probes sit on opposite sides of an input step:

```
      boundary tₙ                        frame [tₙ, tₙ₊₁]                boundary tₙ₊₁
──────┬────────────────────┬──────────────────────────────────────┬──────
      quiescence:          drain:          integrate under u_new  guards see u_new
      priors sampled       u_old → u_new   probes: σ(x̂(θ), u_new)
      under u_old
```

The localization trigger is "prior not-holding, holding now" ([§8.4][s8-4]). Its original
endpoint argument claimed the bracket's left end stays strictly not-holding because
the interpolant reproduces the endpoint exactly (x̂(0) = xₙ) and probe sweeps are
deterministic. That argument silently assumes the θ = 0 probe reproduces the prior's
*sample*. Determinism only guarantees it reproduces the prior's *state*. The two
coincide only if everything else the guard reads is also unchanged — and one input is
not.

## Section 2 — Only `u` can differ: the availability audit

Enumerate every channel a guard can read, and compare its value in the prior's
evaluation context (tₙ's quiescence) against its value at a θ = 0 probe:

| guard input | at θ = 0, comes from | consistent with the prior's context? |
|---|---|---|
| `x` | xₙ, written into the state buffer | yes — x̂(0) = xₙ exactly (Hermite endpoint) |
| `y`, `w`, `ws` | recomputed by the interior sweep from xₙ | yes, *if* the sweep's other inputs match |
| `m` | current mode stores | yes — handlers only run at boundaries, and priors are sampled at quiescence, after the last one |
| discrete cells | the signal table, ZOH-held | yes — interior sweeps exclude discrete entries, so the cells still hold the frame's values |
| `t` | computed: t = tₙ (indexed grid time, never accumulated) | yes |
| `u` | the current root input slots — **post-drain** | **no** — the prior was sampled pre-drain |

Sweeps are deterministic, so the audit is exhaustive: if the θ = 0 evaluation
disagrees with the prior, the drain is the only possible culprit. (Under the retired
deferral design there was a second source — the manufactured not-holding prior. With
row 181's honest priors, the prior is unconditionally the quiescent sample, and the
drain is the *sole* source of disagreement. This is what makes the validation probe a
conclusive discriminator rather than a heuristic.)

One ordering rule makes the "still consistent" rows true, and the spec makes it
explicit ([§8.4][s8-4]): the localized-event trigger checks run against the **arrival sweep**
at tₙ₊₁ — *before* the boundary's due-gated sweep refreshes any discrete cells. So
probes see the same held cell values the frame actually evolved under, and t\*
firings strictly precede tₙ₊₁'s boundary sequence.

## Section 3 — Input epochs, and the two kinds of edges

The drain runs once per frame, at the frame top, so `u` is piecewise constant in
model time. Call each maximal span of constant `u` an **input epoch**; the drain at a
frame top is an epoch seam. A guard σ(x, u, …), viewed as a function of model time
along the run, can change for exactly two reasons:

- **Within an epoch** — only the trajectory moves. A not-holding → holding transition
  inside a frame is a **trajectory-caused edge**: there is a genuine crossing instant
  t\* ∈ (tₙ, tₙ₊₁] at which σ passed through zero. Root-finding it is meaningful.
- **At an epoch seam** — `u` steps discontinuously, and σ can jump from negative to
  positive *at the seam itself*, with no state motion and no zero crossing anywhere.
  It teleported across zero rather than passing through it. This is an
  **epoch-caused edge**.

A concrete pair, for a localized event with guard σ = x.ω − u.ω_ref (state crossing a
commanded threshold):

- *Trajectory-caused:* ω_ref stays 100; ω climbs 95 → 103 during the frame. Prior:
  σ = −5. Arrival: σ = +3. The state moved past the threshold; ω(t\*) = 100 exists
  and is localizable.
- *Epoch-caused:* ω holds at ~95; the drain lands a new command ω_ref = 90. Prior
  (pre-drain): σ = 95 − 100 = −5. Arrival: σ = 95 − 90 = +5. An edge is observed —
  but the *threshold moved past the state*. Under u_new, σ was ≈ +5 for the entire
  frame. There is no root; a bracketing root-finder cannot even initialize (both ends
  positive).

The gate idiom for localizing mixed predicates — `(gate) ? σ : -one(σ)`, with the
gate reading `u`/`m` — is blessed by row 179, so `u`-reading localized guards are a
supported class, not a corner. Epoch-caused edges are therefore ordinary events, and
the trigger must handle them by design, not by accident.

## Section 4 — The validation sequence

Nothing per-guard is stored across the boundary — the detection registers hold
predicate samples and a counter, never values (detection bookkeeping, not model
memory). The validation is a *reconstruction*, using the standard probe mechanics
(write a state into the buffer, run the interior sweep, call the guard):

1. Integration lands at tₙ₊₁; arrival sweep; each localized event's guard is
   evaluated. Event E: prior not-holding, σ₁ = σ(tₙ₊₁) holding → trigger.
2. **θ = 0 probe, before anything else is built**: write xₙ into the state buffer —
   xₙ is already retained by the stepper, because the Hermite interpolant needs
   (xₙ, ẋₙ, xₙ₊₁, ẋₙ₊₁), so the probe requires no new retention — run one interior
   sweep, evaluate E's guard → σ₀. No interpolant is needed for this probe:
   x̂(0) = xₙ *identically*, so it runs before the interpolant exists.
3. **σ₀ < 0 (not-holding)** — trajectory-caused. Genuine bracket. Now pay the lazy
   costs — one sweep for ẋₙ₊₁, build the interpolant — and hand ITP/Brent the
   bracket **with both endpoint values**, (tₙ, σ₀) and (tₙ₊₁, σ₁), σ₁ retained
   transiently from step 1's evaluation. This is where the left bracket value comes
   from: ITP and Brent are value-based, the prior stores only a predicate, and the
   spec previously never said who supplies σ(tₙ). The validation probe *is* its
   source.
4. **σ₀ ≥ 0 (holding)** — epoch-caused. The drain flipped the guard at the frame
   top; no in-frame crossing exists. Discard the localization: E fires inside
   tₙ₊₁'s ordinary iteration. Mechanically this is almost nothing — *not localizing
   is the action*: skip and fall through, and the boundary's event phase detects the
   edge (prior not-holding, sample holding) and fires E exactly as a
   boundary-detected event.
5. The buffer is overwritten by the next probe or restored as usual — probes own the
   buffer during localization.

Cost accounting for the degenerate path: one interior sweep total. ẋₙ₊₁'s sweep is
never paid, the interpolant is never built, `localization_budget` is untouched (it
counts localizations performed; none was), and no warning is emitted.

## Section 5 — Why boundary firing is the *correct* outcome, not a concession

The cause of an epoch-caused edge is the input step, and input timing is a frame fact
by construction: the input arrived asynchronously during the previous frame's wall
time and was quantized to the frame top by the drain. This is the same doctrine that
forbids draining at t\* boundaries — replay determinism must not depend on
localization arithmetic. An event whose cause is only defined to frame resolution has
a firing time that is only defined to frame resolution. **There is no finer truth for
localization to recover — t\* doesn't exist**, not merely "isn't worth computing."

This is the discrete-driven exactness doctrine ([§2.1][s2-1], row 179) resurfacing: for edges
caused by `u`/`m` steps, boundary detection is not a cheap approximation, it *is* the
semantics. Localization's jurisdiction is trajectory facts; the θ = 0 probe is the
jurisdictional check at its door. That is also why the degenerate path warns nothing:
an input-gated event firing at boundary granularity is correct sampled-input timing,
and a warning would fire on every such event — pure noise.

Note the symmetry with the *right*-end degeneracy already in [§8.4][s8-4]: t\* = tₙ₊₁
exactly → localization result discarded, event fires inside tₙ₊₁'s ordinary
iteration. The θ = 0 validation is its left-end mirror, discharging into the same
mechanism: when localization has nothing to say, the event falls back into the
ordinary boundary iteration — one boundary, one snapshot, no special path.

## Section 6 — What the validation repairs in the spec's own argument

Before: "t\* = tₙ is structurally impossible" rested on the prior's *testimony*
(not-holding at quiescence) plus determinism — and the testimony is about the wrong
epoch for `u`-reading guards. After: the argument rests on the *observed* left end.
The validation probe either establishes σ₀ strictly not-holding — from which
everything downstream follows unconditionally: the returned endpoint (the smallest
probed point where the predicate holds) is strictly later than tₙ, the guard
observably holds at t\*, and the post-fire prior records an actual observation — or
it routes the event to boundary firing, where none of those claims are needed.
Observation replaces testimony at exactly one point, and the whole endpoint-policy
chain becomes unconditional.

A closing remark on design method: the probe never asks *why* the prior disagrees
with the left end — drain-staled, or (in the retired deferral design)
manufactured. It replaces the register's account of the left end with an observation
of the left end under the current epoch, and lets the ordinary machinery act on it:
bracket → localize, no bracket → boundary firing. One observation, no case analysis
over the register's history. Mechanisms with that shape tend to survive redesigns of
the bookkeeping around them — this one survived the replacement of deferral by
budgeted re-firing (row 181) without changing a single step.

<!-- citation link definitions — generated by tools/linkify.jl; do not edit -->
[s2-1]: framework_spec.md#21-events-two-detection-policies
[s8-4]: framework_spec.md#84-localization-mechanics
[s8-6]: framework_spec.md#86-event-iteration-at-boundaries-to-quiescence-budgeted
[s9-4]: framework_spec.md#94-inbound-per-device-staging-representation-and-the-drain
