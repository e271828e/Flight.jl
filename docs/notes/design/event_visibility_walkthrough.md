# Event visibility at a boundary, from the ground up

*A companion explainer, not normative text. The ground truth is
`framework_spec.md` [§5.3](framework_spec.md#53-structural-feedthrough-stage-roles-schedule-and-step-boundaries) (step-boundary semantics), [§8.6](framework_spec.md#86-event-iteration-at-boundaries-to-quiescence-once-per-event) (event iteration)
and decision row 100, which settled the visibility rule below (round-start
`u`, live `y`, pre-materialized bundles) on 2026-07-31 during the round-3
cluster-4 adjudication. If this document and the spec ever disagree, the
spec wins.*

Everything here answers one question: **when several events fire at the same
boundary, what world does each handler see?** The spec's one-sentence answer
([§5.3](framework_spec.md#53-structural-feedthrough-stage-roles-schedule-and-step-boundaries)) is: *all guards and handlers read the same boundary snapshot, plus
their own component's refreshed ports; one component's transition reaches
others through the next sweep.* This document builds the mechanism that makes
that sentence true — and shows why it costs nothing.

## 1. The primitive fact: cells hold references to immutable values

A signal-table cell does not hold bytes that get mutated in place — it holds
a *reference* to an immutable value ([§4.1](framework_spec.md#41-immutable-value-semantics), [§7](framework_spec.md#7-state-and-data-representation)). "Writing a cell" means
replacing that reference with a reference to a *new* value:

```julia
table[battery_bus_voltage] = 24.0    # cell → ref to the value 24.0
# ... later, a re-decode:
table[battery_bus_voltage] = 0.0     # cell now → ref to 0.0
```

The old value is untouched by the overwrite. Anyone who already loaded the
old reference keeps a perfectly coherent view of `24.0` for as long as they
hold it, and the GC keeps it alive. This is exactly the mechanism [§9.2](framework_spec.md#92-outbound-snapshot-publication) uses
to publish snapshots across tasks; here the same trick is applied *within*
one boundary, between the event phase's writers and readers.

## 2. Bundles and gathers: the ownership partition

Every user function is called as `fn(comp, args)`, where `args` is a
NamedTuple bundle of views ([§5.2](framework_spec.md#52-two-stage-outputs-signatures-bundles-and-the-hand-off-laws)). Two of its fields are **gathers** — small
compiled loops that load cell references into a NamedTuple — and they
partition the table by ownership:

- **`y` gathers the component's *own* cells** (its published ports and
  locals);
- **`u` gathers *foreign* cells** (other components' ports, routed through
  the wiring).

No cell is in both sets. That disjointness is what makes the visibility rule
below expressible without policing anything.

## 3. The round, and who writes what

At a step boundary the event phase iterates ([§8.6](framework_spec.md#86-event-iteration-at-boundaries-to-quiescence-once-per-event)): rounds of

> re-run the boundary sweep → evaluate **all** guards once → fire the
> newly-fired events → repeat until a round fires nothing,

with each declared event firing at most once per boundary. Firing an event
means `handler → project → auto-publish → re-run the component's own output
stages` (the per-event **re-decode**; auto-published cells are rewritten from
the just-latched stores at their stage-1 position, [§5.3](framework_spec.md#53-structural-feedthrough-stage-roles-schedule-and-step-boundaries)).

Now the structural fact the whole design leans on: **no user code writes
anything, anywhere.** A handler is pure like every other user function: it
returns `(; x, m)` as immutable values, and the *framework* latches them into
the component's state stores — immediately after the return, before `project`
and before the next fired event runs (that immediacy is what makes
same-component sequential composition real, [§6.1](framework_spec.md#61-connections-and-hierarchy)). The only mid-round *table*
writes are the framework's re-decodes, confined to the *firing component's
own* cells. So the set of cells that can change between a round's start and
its end is exactly "the cells of components that fired this round" — small,
and known the moment the guards have been evaluated.

## 4. The hazard the rule exists for

Suppose components A and B both fire in the same round, and B's handler reads
a port of A through its `u`. Executed naively against the live table:

1. A's handler runs, A's re-decode refreshes A's cells;
2. B's handler runs, its `u` gather loads A's cells — and sees
   **post-transition A**.

Whether B observes pre- or post-transition A now depends on *execution
order* — and "declaration order" ([§11.2](framework_spec.md#112-the-declaration-inventory)) orders events only within one
component. Cross-component order would fall to the build's executor order, a
schedule artifact that rewiring silently permutes: model semantics leaking
from a build detail, the same disease [§8.6](framework_spec.md#86-event-iteration-at-boundaries-to-quiescence-once-per-event) diagnosed when it rejected the
single-pass cascade.

## 5. The rule, and its nearly-free mechanism

**`u` binds at round start; `y` binds live.** Concretely, the round becomes:

> sweep → evaluate all guards → determine the fired set →
> **materialize the `u` gathers of every fired component now**, before any
> handler runs → execute the fired events, each with its per-event re-decode
> reading that same round-start `u`.

"Materialize" means: build the args NamedTuple early. The executor already
builds one per call; this moves the `u` reference-loads a few microseconds
earlier, once per fired event. Because the values are immutable ([§1](#1-the-primitive-fact-cells-hold-references-to-immutable-values)), loading
the references early *is* a capture of the round-start world — no copy, no
shadow table, no allocation. Nothing that happens to the cells afterwards can
change what the bundle already holds.

`y` stays bound to the live table, which is what delivers the "plus their own
component's refreshed ports" clause — see [§6.1](framework_spec.md#61-connections-and-hierarchy) below.

## 6. One boundary, end to end

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
   set `{Battery.ovld, Avionics.cap}`. (The fired set is snapshot-consistent
   under any option — guards run before any handler. The question is only
   what the handlers see.)
3. *Pre-materialization.* Battery's `args.u` loads the load-current
   reference; Avionics' `args.u` loads the reference to `24.0`. These
   NamedTuples now own those references.
4. *Battery fires* (canonical order — see [§4.4](framework_spec.md#44-function-valued-signals-environment-access); it doesn't matter): the handler
   returns `(; m = (; m..., tripped = true))`, which the framework latches
   into the mode store; project; the re-decode runs `h_x` against the new
   mode → the cell `battery/bus_voltage` is overwritten with a reference to
   `0.0`.
5. *Avionics fires*: `handler_cap` runs with the bundle built in step 3 —
   `u.bus_voltage` dereferences the `24.0` it already holds. The cell
   overwrite in step 4 replaced the *cell's* reference, not the value the
   bundle captured. So `v_at_capture = 24.0`: the avionics latched the world
   *as it stood at the boundary*, which is what "two simultaneous events"
   ought to mean.

**Round 2:** full re-sweep against post-transition state — `bus_voltage =
0.0` now propagates everywhere. Guards re-evaluate; if the avionics also
declares a low-voltage guard (`u.bus_voltage < 18`), it is *newly* enabled
and fires now, with round-2 bundles, seeing `0.0`. This is "one component's
transition reaches others through the next sweep", working as written:
cross-component causality is quantized into rounds, and a genuine cascade
takes one round per causal link.

**Round 3:** re-sweep, nothing new fires → quiescence → due `g` updates →
publication.

## 7. Why cross-component order dissolves

Ask what Avionics' handler could possibly observe about Battery having fired
earlier in the same round:

- Battery's latched `x`/`m`? Invisible — bundles carry only *own* state; no
  cross-component state view exists anywhere in the design.
- Battery's refreshed cells? Invisible — Avionics' `u` was materialized
  before any re-decode.

Swapping the execution order of Battery and Avionics changes nothing
observable, so **cross-component handler order stops being a semantic
decision at all**. The framework still fixes a canonical order (the
executor's component order, declaration order within a component) — but only
so the [§13.4](framework_spec.md#134-runtime-failures-one-catch-site-an-execution-cursor) execution cursor ("event round *r*, component *c*") and the
diagnostics stream are deterministic and nameable, not because trajectories
depend on it.

## 8. The own-component exception: why `y` binds live

[§5.3](framework_spec.md#53-structural-feedthrough-stage-roles-schedule-and-step-boundaries) promises handlers "their own component's refreshed ports", and it must.
If one component declares events `A1`, `A2` (declaration order) and both
fire, `A2`'s handler sees post-`A1` state stores — `x`/`m` are live; that is
the sequential-composition guarantee — so its own `y` must be the decode of
*that* state, not a stale pre-`A1` decode; otherwise `A2` receives an
incoherent bundle where `y ≠ h(x)`.

The ownership partition ([§2](#2-bundles-and-gathers-the-ownership-partition)) makes the split exact: `y` = own cells (live),
`u` = foreign cells (round-start). One corner, stated for completeness: a
component whose input port is wired, through some contortion, from its own
stage-1 output still gets round-start `u` — the binding rule is uniform, and
the freshness rule is about `y` only.

The re-decode gets the same treatment: the fired component's output stages
also read `u`, and they reuse the *same* pre-materialized field. Its
refreshed ports are therefore the decode of (new `x`/`m`, round-start `u`) —
coherent with the round's semantics — and the next round's full re-sweep
re-decodes everything against settled values anyway, so nothing stale
survives the round.

## 9. Costs, and the rejected shapes

- **Allocations: zero new.** A NamedTuple of references to immutable values
  is stack-shaped; the executor builds it anyway, just later. Rounds where
  nothing fires — the overwhelming majority of boundaries — never reach the
  materialization step.
- **Live table + canonical order** (rejected): deterministic, free, but makes
  cross-component handler semantics depend on the executor's schedule order —
  a rewiring that permutes it silently changes trajectories.
- **Copy the table at round start** (rejected): identical semantics, pays an
  allocation per firing round; unnecessary because only fired components'
  cells can change mid-round, and early loading sidesteps even those.

## 10. The doctrine limit, stated openly

Round-start visibility means a handler *cannot* opt into seeing a same-round
foreign transition. If two components genuinely need same-instant sequential
coupling, the design's answer is the existing one: that coupling is a
cascade — it takes one round, and it is deterministic. If even one round is
too slow, the two FSMs belong in one component, where declaration order gives
exact sequencing. This is the synchronous-languages position (transitions in
a micro-step see the pre-state; effects appear at the next micro-step), and
it is the trade this rule makes.
