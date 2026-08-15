# Event visibility at a boundary, from the ground up

*A companion explainer, not normative text. The ground truth is
`framework_spec.md` [§5.3][s5-3] (step-boundary semantics), [§8.6][s8-6] (event iteration)
and decision D-154, which settled the visibility rule below (the table is
written only by sweeps; a component fires at most one event per round) on
2026-08-07, superseding the round-3 rule of D-016/D-100/D-152. If this document
and the spec ever disagree, the spec wins.*

Everything here answers one question: **when several events fire at the same
boundary, what world does each handler see?** The spec's one-sentence answer
([§5.3][s5-3]) is: *a handler executes against exactly the world its guard fired
on; one component's transition reaches everyone — itself included — through
the next sweep.* This document builds up to that sentence and shows why it
needs no mechanism at all.

## 1. The primitive fact: cells hold references to immutable values

A signal-table cell does not hold bytes that get mutated in place — it holds
a *reference* to an immutable value ([§4.1][s4-1], [§7][s7]). "Writing a cell" means
replacing that reference with a reference to a *new* value:

```julia
table[battery_bus_voltage] = 24.0    # cell → ref to the value 24.0
# ... a later sweep, after the breaker tripped:
table[battery_bus_voltage] = 0.0     # cell now → ref to 0.0
```

The old value is untouched by the overwrite. Anyone who already loaded the
old reference keeps a perfectly coherent view of `24.0` for as long as they
hold it, and the GC keeps it alive. This is exactly the mechanism [§9.2][s9-2] uses
to publish snapshots across tasks; here the same property is what lets a
handler's bundle stay meaningful no matter what later sweeps do.

## 2. Bundles and gathers

Every user function is called as `fn(comp, args)`, where `args` is a
NamedTuple bundle of views ([§5.2][s5-2]). Two of its fields are **gathers** — small
compiled loops that load cell references into a NamedTuple — and they
partition the table by ownership:

- **`y` gathers the component's *own* cells** (its published ports and
  locals);
- **`u` gathers *foreign* cells** (other components' ports, routed through
  the wiring).

No cell is in both sets. The partition is still worth knowing — it is why a
component can never see another's state — but it no longer carries a
freshness distinction: within a round, `y` and `u` alike read the table that
round's sweep left, because that is the only table there is.

## 3. The round, and who writes what

At a step boundary the event phase iterates ([§8.6][s8-6]): rounds of

> re-run the boundary sweep → evaluate **all** guards once → fire the
> eligible events, **at most one per component** → repeat until a round fires
> nothing,

with each declared event firing at most `firing_budget` times per boundary
([§8.6][s8-6], default 4, eligibility read against its last-observed sample), and declaration
order ([§11.2][s11-2]) picking among a component's simultaneously-eligible events.
Firing an event means `handler → project`. That is all.

Now the structural fact the whole design leans on: **the signal table has a
single writer, and it is the sweep.** A handler is pure like every other user
function: it returns `(; x, m)` as immutable values, and the *framework*
latches them into the component's state stores immediately after the return,
before `project`. Latching writes *stores*, not cells; `project` normalizes
those stores; auto-publication is a stage-1 sweep act like any other
([§12.5][s12-5]). Nothing — no user code, no framework step — writes the table
between the sweep that opened a round and the sweep that opens the next.

So the set of table cells that can change between a round's start and its
end is *empty*.

## 4. Why cross-component order dissolves

The classic hazard: components A and B both fire in the same round, and B's
handler reads a port of A through its `u`. If A's transition could reach the
table before B's handler runs, then whether B sees pre- or post-transition A
would depend on *execution order* — and "declaration order" orders events
only within one component. Cross-component order would fall to the build's
executor order, a schedule artifact that rewiring silently permutes: model
semantics leaking from a build detail, the same disease [§8.6][s8-6] diagnosed when
it rejected the single-pass cascade.

Under [§3](#3-the-round-and-who-writes-what) the hazard has no way to arise. A's transition reaches nothing
until the next sweep; B's bundle reads a table that has not moved. There is
no mechanism here — no freezing, no pre-materialization, no copy — because
there is nothing mid-round for order to observe. Swapping the execution order
of A and B changes nothing observable, so **cross-component handler order
stops being a semantic decision at all**. The framework still fixes a
canonical order (the executor's component order, declaration order within a
component) — but only so the [§13.4][s13-4] execution cursor ("event round *r*,
component *c*") and the diagnostics stream are deterministic and nameable,
not because trajectories depend on it.

## 5. One boundary, end to end

Two components, wired together:

```julia
# Battery: publishes the bus voltage; trips its breaker on overload
h_x(b::Battery, (; m))          = (; bus_voltage = m.tripped ? 0.0 : 24.0)
guard_ovld(b::Battery, (; u))   = u.i_load > b.i_max          # wired load current
handler_ovld(b::Battery, (; m)) = (; m = (; m..., tripped = true))

# Avionics: on reaching target altitude, latches what it saw at capture
guard_cap(a::Avionics, (; x))   = x.alt - a.alt_target        # boundary-detected: σ ≥ 0
handler_cap(a::Avionics, (; m, u)) =
    (; m = (; m..., captured = true, v_at_capture = u.bus_voltage))
```

`avionics.u.bus_voltage` is wired from the battery's port. At the boundary
`t = 12.34`, the integration step has *both* pushed the load current over
`i_max` and carried `alt` across `alt_target` — two independent events on one
boundary.

**Round 1:**

1. *Sweep.* The table is boundary-consistent: `battery/bus_voltage → 24.0`,
   the load current above `i_max`, `alt` past target.
2. *Guards, all at once.* Both hold, both priors were not-holding → fired
   set `{Battery.ovld, Avionics.cap}`. They belong to different components,
   so both fire this round.
3. *Battery fires* (canonical order — see [§4](#4-why-cross-component-order-dissolves); it doesn't matter): the
   handler returns `(; m = (; m..., tripped = true))`, the framework latches
   it into the mode store, `project` normalizes. The table is untouched:
   `battery/bus_voltage` still holds `24.0`.
4. *Avionics fires*: `handler_cap` reads `u.bus_voltage` from that same
   round-1 table and gets `24.0`. So `v_at_capture = 24.0`: the avionics
   latched the world *as it stood at the boundary*, which is what "two
   simultaneous events" ought to mean.

**Round 2:** full re-sweep against post-transition state — `h_x` now decodes
`tripped = true`, so `bus_voltage = 0.0` propagates everywhere. Guards
re-evaluate; if the avionics also declares a low-voltage guard
(`u.bus_voltage < 18`), it is *newly* enabled and fires now, seeing `0.0`.
Cross-component causality is quantized into rounds, and a genuine cascade
takes one round per causal link.

**Round 3:** re-sweep, nothing new fires → quiescence → due `g` updates →
publication.

Same trajectory as the old design produced, with strictly simpler mechanics.

**The same-component twin.** Now suppose one component declares `A1` and `A2`
(in that order) and both guards hold in round 1. Only `A1` fires: a component
fires at most one event per round. `A2` is not lost — it is *blocked*, and
round 2 re-evaluates its guard against the post-`A1` sweep. If the premise
still holds, `A2` fires in round 2, reading a table that already reflects
`A1`'s transition — sequential composition, exactly as intended, one round
later. If `A1`'s transition falsified `A2`'s premise, `A2` does not fire at
all, which is the honest answer; a within-round sequence would have fired it
on a premise that no longer held. Declaration order is thus a **priority with
re-decision**, not a simultaneity.

## 6. The epoch rule

Everything above is one sentence: **a handler executes against exactly the
world its guard fired on.** Own `y`, foreign `u`, own `x` and `m` are all the
firing round's sweep — a single epoch, no bundle straddling two of them. In
particular `y = h(x)` holds at every handler entry, the coherence a
component author is entitled to assume.

Serialization is what makes this free. A component's state stores are written
only by its own handlers, and it fires at most one event per round, so no
same-round writer precedes any handler's entry: there is no stale-decode
corner to patch, because there is no second writer to be stale with respect
to. Sequential composition still exists in full — it simply happens across
rounds ([§5](#5-one-boundary-end-to-end)), each later event re-decided against a settled table.

## 7. Costs, and the rejected shapes

- **Mechanism: none.** Nothing is materialized early, staged, carried or
  copied; the executor builds each handler's bundle at dispatch, from the
  live table, exactly as it builds every other bundle. "No shadow table, no
  allocation" is true by construction.
- **One extra sweep per serialized same-component event.** A component with
  two genuinely simultaneous events takes one more intra-boundary round than
  it used to. Sweeps are microseconds and the case is rare.
- **Per-event re-decode with a frozen round-start `u`** (rejected; the prior
  design of record, D-016/D-100/D-152): two freshness classes inside one
  bundle, delivered by pre-materializing the fired components' gathers —
  standing staging machinery whose only purchase was same-round multi-handler
  own-`y` freshness, which serialization provides for free.
- **Live table + canonical order** (rejected): deterministic, free, but only
  under a design where transitions *do* reach the table mid-round — which
  makes cross-component handler semantics depend on the executor's schedule
  order.
- **Copy the table at round start** (rejected): identical semantics, pays an
  allocation for nothing.
- **Handlers stripped of their own `y`** (rejected): makes the incoherence
  unwritable rather than absent, and keeps both the two-epoch bundle and the
  fire-on-a-falsified-premise hole.

## 8. The doctrine limit, stated openly

A handler *cannot* opt into seeing a same-round foreign transition. If two
components genuinely need same-instant sequential coupling, the design's
answer is the existing one: that coupling is a cascade — one round per link,
deterministic. If even one round is too slow, the two FSMs belong in one
component, where declaration order now gives priority with re-decision across
rounds: still exact, still deterministic, and each step taken against a
coherent world rather than a mid-transition one. This is the
synchronous-languages position (transitions in a micro-step see the
pre-state; effects appear at the next micro-step), and it is the trade this
rule makes.

<!-- citation link definitions — generated by tools/linkify.jl; do not edit -->
[s11-2]: framework_spec.md#112-the-declaration-inventory
[s12-5]: framework_spec.md#125-the-always-on-conformance-check
[s13-4]: framework_spec.md#134-runtime-failures-one-catch-site-an-execution-cursor
[s4-1]: framework_spec.md#41-immutable-value-semantics
[s5-2]: framework_spec.md#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws
[s5-3]: framework_spec.md#53-structural-feedthrough-stage-roles-schedule-and-step-boundaries
[s7]: framework_spec.md#7-state-and-data-representation
[s8-6]: framework_spec.md#86-event-iteration-at-boundaries-to-quiescence-budgeted
[s9-2]: framework_spec.md#92-outbound-snapshot-publication
