# Sample-Time Extensions: A Worked Proposal

**Status: non-normative companion; adopted into the spec.** This document is the
worked sample-time proposal for the framework: the declaration surface, compilation
pipeline, runtime representation and diagnostics for **phase offsets** and **absolute
rate declarations**. It supersedes section 2 of the extensions charter
(`framework_extensions.md`), which sketched the same territory before the pipeline was
worked out. **Adopted 2026-08-12 as decision [D-185][d-185]–[D-187][d-187]**, with the section 9
amendments applied ([§10.5][s10-5], [§8.7][s8-7], [§9.1][s9-1], [§9.2][s9-2], [§9.7][s9-7], [§14.5][s14-5], [§14.8][s14-8], Appendices B/C —
the spec wins on any conflict) and two deltas from the proposal as written: the
`ratespec` normalization of section 2 was **dropped** — the wrappers are the whole
value vocabulary, a bare integer or bare quantity is a declaration error, and an
unlisted discrete child defaults to `Relative(1)` (bare-value entries surviving in the
examples below read as `Relative(K)`/`Absolute(q)` under the adopted surface) — and
the declaring method was renamed `rates` → `sample_times` (this document is swept to
match; quoted pre-adoption spec text keeps its original wording). What remains
non-normative here is the exposition: the worked examples, timelines and derivations
live only in this document. Everything below stays on the **static lattice**: tick
instants remain build-time constants, and nothing from the closed tick-queue axis
(`framework_extensions.md` section 1.4) is dragged in.

Provenance: the gap-analysis discussion of 2026-08-08/09 and the sample-time pipeline
sessions of 2026-08-12.

## Contents

- [1. Where the parent-relative paradigm hurts](#1-where-the-parent-relative-paradigm-hurts)
- [2. The declaration surface: two registers, one concept](#2-the-declaration-surface-two-registers-one-concept)
- [3. The composition law](#3-the-composition-law)
- [4. Anchors: absolute declarations anywhere in the tree](#4-anchors-absolute-declarations-anywhere-in-the-tree)
- [5. The pipeline, stage by stage](#5-the-pipeline-stage-by-stage)
- [6. Worked example: the staggered sensor suite](#6-worked-example-the-staggered-sensor-suite)
- [7. Worked example: nested anchors and the refusal path](#7-worked-example-nested-anchors-and-the-refusal-path)
- [8. Grid diagnostics: suggestion, attribution, warnings](#8-grid-diagnostics-suggestion-attribution-warnings)
- [9. Spec amendments this proposal requires](#9-spec-amendments-this-proposal-requires)

---

## 1. Where the parent-relative paradigm hurts

Inside an authored assembly, relative integers are the right language: an inner/outer
cascade at 1:5 *is* a 1:5 design, and whoever writes the assembly owns that
arithmetic. The convention that makes `K ≥ 1` livable should be stated explicitly:
**a scope's base rate is its fastest member** (that member gets `K = 1`). With that
discipline, refactoring never fights the tree.

The degenerate case is the top. Suppose avionics must run at 50 Hz (the actuator bus
says so) and the nav filter at 30 Hz (the GPS says so). Today the user does the
arithmetic the framework should do:

```julia
# by hand: the common grid is 150 Hz, so...
sample_times(::Aircraft) = (avionics = 3, nav = 5)
Simulation(world; Δt_base = 1//150, h = 1//150)
```

Those `3` and `5` are not design ratios — nothing in the aircraft is "3:5"; they are
residue of a GCD computation. And they are *unstable* residue: change nav to 40 Hz
and the grid becomes 200 Hz, `avionics = 4`, `nav = 5`, `Δt_base = 1//200` — three
coordinated edits at two sites to express one changed fact. The pain sits at the
boundary where relative structure meets absolute facts, and the spec's own doctrine
already names that boundary: "absolute rates are deployment decisions made once at
the root" ([§10.5][s10-5]).

There is a second, subtler pain the root-only framing creates. Some absolute rates
are not deployment decisions at all — they are **facts about the modeled system**. A
GPS receiver emits at 1 Hz because the device does; a 1553 bus has its schedule; an
ADC pipeline has a fixed conversion offset. Forcing those declarations to the root
*breaks* encapsulation: the root must know device internals in order to re-declare
them. Section 4 develops this into the case for absolute declarations anywhere in
the tree.

## 2. The declaration surface: two registers, one concept

A rate entry can be declared in one of two **registers** — relative or absolute —
and the proposal makes each an explicit wrapper type, with `Period`/`Hz` as pure
quantity types underneath:

```julia
struct Period; T::Rational{Int}; end          # Period(1//50)
struct Hz;     f::Rational{Int}; end          # Hz(50), Hz(1//2)
period(q::Period) = q.T
period(q::Hz)     = 1 // q.f

struct Relative
    K::Int                                    # fire every K-th scope tick...
    Φ::Int                                    # ...starting from the Φ-th
    Relative(K::Integer, Φ::Integer = 0) = new(K, Φ)
end

struct Absolute
    T::Rational{Int}                          # normalized: always a period
    τ::Rational{Int}                          # offset, in exact rational seconds
    Absolute(q::Union{Period, Hz}, τ = 0//1) = new(period(q), τ)
end
```

The symmetry is the point, and worth stating as the mental model: **both registers
are the same pair — (period, phase) — expressed in different unit systems, and the
type name declares the unit system.** `Relative(K, Φ)` is period and phase in *scope
ticks*; `Absolute(q, τ)` is period and phase in *seconds*. Even the validity
constraint has the same shape in both: `0 ≤ phase < period` (`0 ≤ Φ < K` and
`0 ≤ τ < T` respectively).

Design notes, each with its reason:

- **`Absolute` normalizes at construction.** `Hz` vs `Period` is a spelling choice at
  the declaration site, not a fact the pipeline needs to remember; the struct stores
  the rational period, and diagnostics can always print both renderings
  (`1//50 s (50 Hz)`). Keeping `Period`/`Hz` as real types still pays off elsewhere:
  deployment's `Δt_base` should accept them too — `Simulation(world; Δt_base =
  Hz(1000))` — so one quantity vocabulary serves both sites.
- **Terse sugar for the zero-phase cases.** The overwhelmingly common declaration is
  a bare ratio, and `Relative(1)` everywhere would be ceremony. One normalization
  function at the Stratum A intake keeps today's ergonomics:

  ```julia
  ratespec(x::Integer)           = Relative(x)
  ratespec(q::Union{Period, Hz}) = Absolute(q)      # bare quantity = zero offset
  ratespec(r::Union{Relative, Absolute}) = r
  ```

  so `imu = 1` and `b = Hz(50)` stay legal, and the explicit wrappers appear exactly
  when a phase does — which is also when the reader most needs the field names.
- **Exact rationals, never floats.** GCD over floats is ill-defined; `Rational`
  makes every derivation below exact and every diagnostic printable. `Rational{Int}`
  fields enforce it structurally — but someone *will* type `Period(0.02)`, and a
  missing-method error is a poor teacher. A dedicated float method that throws
  "periods are exact rationals; write `Period(1//50)`" is one line and worth it.
  Note `Hz(::Rational)` earns its keep here: `Hz(1//2)` is the honest spelling of
  0.5 Hz.
- **Range validation belongs to Stratum A, not the constructors.** Constructor-side
  checks give instant REPL feedback but fail the *evaluation* of a `sample_times` body with
  a raw exception, whereas Stratum A's charter ([§9.1][s9-1], [§13.1][s13-1]) is to *collect*
  declaration defects with path attribution ("`gnss` in `sample_times(::Sensors)`: phase 20
  not in `0 ≤ Φ < K = 20`"). Today's `K ≥ 1` check already lives there. One
  validation site, collected diagnostics; the constructors are plain data carriers.
- **Naming caveat, flagged and dismissed**: `Relative` and `Absolute` are generic
  adjectives to claim as exported names; `RelativeRate`/`AbsoluteRate` would be
  clash-safer. But the short forms read beautifully inside a `sample_times` NamedTuple,
  where their context is unambiguous. Cosmetic; decide at build time.

## 3. The composition law

**Relative form.** `Relative(K, φ)` with `0 ≤ φ < K` means "fire on every K-th tick
of my scope, starting from its φ-th" — scope-relative, exactly like the bare
multiplier, so a library assembly can stagger its internal children without knowing
any absolute rate.

The reason this is clean and not merely plausible is that the composition law works
out affinely. Write a scope's compiled schedule as an absolute divisor and phase
`(D_s, Φ_s)`, both in base ticks: it ticks at base indices `n = Φ_s + m·D_s`. A
child declared `Relative(K_c, φ_c)` inside it fires on scope ticks with
`m ≡ φ_c (mod K_c)`, so its base indices are

```
n = Φ_s + (φ_c + l·K_c)·D_s  =  (Φ_s + φ_c·D_s) + l·(K_c·D_s)
```

That is: `D_child = K_c·D_s` (the existing multiplicative rule, untouched) and
`Φ_child = Φ_s + φ_c·D_s`. Multipliers compose multiplicatively, phases affinely,
and both collapse at build into one `(D, Φ)` pair per discrete component — the same
"all scoping compiles away" shape [§10.5][s10-5] already has. The boundary sweep's gating test
changes from `idx % D == 0` to `(idx − Φ) % D == 0`. One subtraction; the lattice
stays static.

One property of the fold deserves proving once, because the runtime leans on it in
section 5.5: **composition preserves `0 ≤ Φ < D`**. If `Φ_s ≤ D_s − 1` and
`φ ≤ K − 1`, then `Φ_child ≤ (D_s − 1) + (K − 1)·D_s = K·D_s − 1 = D_child − 1`. So
every compiled phase is a *canonical* residue of its divisor, by construction, with
no normalization pass.

**Absolute form.** `Absolute(q, τ)` declares tick instants `t = τ + k·T`, with
`0 ≤ τ < T` and τ an exact rational time. Since ticks may only land on base ticks,
**absolute periods and offsets jointly constrain the base grid** — they join the GCD
pool of section 5.3. This is the one subtlety with teeth: an offset of `T/2` can
cost a 2× finer grid; an offset of `T/1000` can cost 1000×. Offsets should be simple
fractions of their period — though section 8 shows why that authoring rule cannot be
the engine's actual test.

Three properties distinguish the registers, and they are what make the split
principled rather than cosmetic:

1. **Relative offsets never refine the base grid.** `φ` is an integer count of scope
   ticks — it selects among instants that already exist. Grid cost is confined
   entirely to the absolute register.
2. **Relative offsets cannot leave the scope grid.** A child can fire on the 3rd or
   19th scope tick of a cycle, never *between* scope ticks. Staggering off-grid
   requires expressing the offset where the grid is finer: in absolute time (paying
   refinement there), or by declaring the scope base finer than its fastest member
   so unused slots exist.
3. **`K = 1` admits no stagger** (`0 ≤ φ < K` forces `φ = 0`): a component occupying
   its scope's full grid cannot be phase-shifted against it. Two same-rate siblings
   are staggered by giving the scope twice their rate and declaring `Relative(2, 0)`
   and `Relative(2, 1)` — the same "stagger costs resolution" rule one level down.

## 4. Anchors: absolute declarations anywhere in the tree

An absolute entry may appear in any scope's `sample_times`, not only the root's:

```julia
sample_times(::Avionics) = (computing = Relative(2), sensors = Absolute(Hz(500)))
```

Call the `(T, τ)` pair such an entry establishes an **anchor**. The mechanics cost
nothing — section 5 shows the composition fold handles anchors at any depth with one
re-seed rule — but the semantics deserve their own section, because three
consequences follow, one of them with teeth.

**Severing.** An absolute entry **severs the child from the enclosing scope's
grid**. The child is not at any multiple or phase of the scope's ticks; if the
assembly is instantiated under a `K = 5` scope somewhere, its relative children
scale with it and the anchored child does not. Three corollaries the spec text would
have to state:

- The `K ≥ 1` constraint vanishes for anchored children — an anchored child may tick
  faster than its scope, because there is no relation left to constrain. "A child
  cannot tick faster than its scope" becomes "…than the scope it is *relative* to."
- The section 1 convention "a scope's base rate is its fastest member" excludes
  anchored children from "member."
- Phase relationships between an anchored child and its relative siblings are
  **deployment-emergent**: whether their ticks ever coincide depends on how the grid
  derivation works out. Not wrong — it is what "absolute" means — but it is exactly
  the kind of fact the printable schedule of section 5.3 exists to answer (section 7
  shows a coincidence structure silently rewired by a 1/300 s offset edit).

**The doctrinal line.** [§10.5][s10-5] rejected absolute-first declaration because it welds
deployment rates into reusable definitions. Mid-tree anchors are coherent only if
the line moves to something sharper: **an absolute declaration inside a library type
is legitimate when the rate is a fact about the modeled system, not a preference
about the simulation.** A GPS receiver emitting at 1 Hz, a bus schedule, an ADC
pipeline with a fixed conversion offset — those rates are as intrinsic to the
assembly as its wiring, and forcing them to the root breaks encapsulation. Meanwhile
"run the controller at 400 Hz in this study" remains a deployment choice, and the
existing idiom — the assembly exposes `K` as a constructor parameter — remains the
answer for it. Absolute pinning *from outside* a subtree's contract stays rejected
as action at a distance. The framework cannot police the distinction; it becomes
authoring doctrine, one paragraph next to the declaration forms.

**What this does not damage.** [§10.5][s10-5]'s structural argument for never caching
`Δt`-derived coefficients survives fully intact, and in a satisfying way: the
pinning happens in the *enclosing assembly's* `sample_times` — the same site where `K`
lives today. The component type itself remains rate-agnostic; its author still
cannot know the rate, still consumes the bundle's `Δt`. What changes is only which
register the parent's declaration uses. Nesting is likewise unproblematic — an
anchor inside an anchored subtree just seeds again; anchors need no relation to each
other beyond sharing the base lattice, which the grid derivation guarantees.

## 5. The pipeline, stage by stage

### 5.1 The inversion: gate data, not due sets

A run has millions of boundaries; no table of due sets is ever built. [§10.5][s10-5]'s "the
due set is a pure function of the frame index" means the due set exists only
*intensionally*: compilation produces, per discrete component, a small amount of
integer **gate data** — the `(D, Φ)` pair — and the due set at any boundary is "the
entries whose gate passes at this boundary's tick index," evaluated inline where the
entry sits in the compiled sweep. So "how does the compiler determine the due set
for each boundary" decomposes into: how the pipeline derives `(D, Φ)` from the
declarations (5.2–5.3), where it lives at each lifecycle stage (5.3–5.4), and how
the runtime gate resolves every boundary flavor (5.5).

### 5.2 Stratum A: validation and the composition fold

Stratum A already does "rates validation and compilation of relative multipliers
into absolute divisors — everything except binding `Δt_base`" ([§9.1][s9-1]). The proposal
keeps that shape and enriches the fold. During the tree walk, each scope's `sample_times`
NamedTuple is read, normalized through `ratespec`, and validated per entry — `K ≥ 1`,
`0 ≤ Φ < K`, `T > 0`, `0 ≤ τ < T`, keys matching the discrete/scope children — all
collected with path attribution, per [§9.1][s9-1]'s usual register.

The fold itself: each node carries a triple `(anchor, m, c)` — which anchor it hangs
from, and its divisor and phase *in that anchor's tick units*. Treat the base grid
itself as **anchor 0**, with `(T, τ) = (Δt_base, 0)`, symbolic until deployment.
Then:

- **Seed**: the root scope is `(A0, 1, 0)`.
- **Relative child** `Relative(K, φ)` under scope `(a, m_s, c_s)`:
  `(a, K·m_s, c_s + φ·m_s)` — the affine law of section 3, applied in anchor-tick
  units.
- **Absolute child** `Absolute(q, τ)`: **sever and re-seed** — allocate a fresh
  anchor `A_k = (period(q), τ)` and continue below with `(A_k, 1, 0)`.

The canonical-residue invariant of section 3 holds within each anchor's subtree by
the same induction (`c < m` throughout, from the seed `(1, 0)`). The output, stored
in the `Build` artifact as plain printable data ([§9.2][s9-2]): an **anchor table** — each
anchor's exact `(T, τ)` rationals plus the declaring scope's path — and a
**component table** of `(anchor, m, c)` triples with their declaration provenance
(the `Relative`/`Absolute` chain down the tree). For a fully relative model the only
anchor is `A0` and the triples *are* the final base-tick `(D, Φ)` pairs, exactly as
today; nothing about the current design's output changes until an absolute entry
exists.

Note what the `Build` genuinely cannot carry when anchors exist: final divisors —
they don't exist until `Δt_base` binds. That is consistent with the `Build`'s
charter, since the same `Build` already backs many `Simulation`s with different
deployment parameters. The structural shift is small and localized: divisor
resolution for anchored entries moves from Stratum A to deployment binding, and it
is one multiply-add per component, not a second composition walk.

### 5.3 Deployment binding: the pool, the rule, the resolution

At `Simulation` construction, where `Δt_base`, `h`, `n` already bind and validate
([§9.1][s9-1]), three new steps slot in.

**The pool.** Collect every anchor's period and every nonzero offset — the
**constraint pool**. A `Δt_base` is admissible **iff it divides every pool entry,
i.e. iff it divides their GCD** (`gcd(a//b, c//d) = gcd(a, c)//lcm(b, d)` — exact
over rationals, which is why section 2 insists on them). The admissible set is
exactly `gcd(pool)/k` for integer `k ≥ 1`; the *maximum* (coarsest) admissible
`Δt_base` is the GCD itself.

**The derive-vs-declare rule.** May the engine set `Δt_base` itself? Only sometimes,
and the boundary is principled. An **unanchored** component (one hanging from
anchor 0) has period `m·Δt_base` — proportional to the deployment knob. If
`Δt_base` were silently derived from the pool, then editing any anchor anywhere
would rescale every unanchored component in the model:

```
root: logger = 10 (relative);  deep inside: b = Absolute(Hz(50))
→ pool {1//50}, Δt_base = 1//50, logger at 5 Hz
now add offset 1//100 to b
→ pool {1//50, 1//100}, Δt_base = 1//100, logger silently at 10 Hz
```

Action at a distance of exactly the kind the design refuses elsewhere. The rule:

> **Derive `Δt_base` only when every discrete component is anchored** (anchor 0 is
> unpopulated) — then `Δt_base` is pure bookkeeping and no component's period
> depends on it. **If any unanchored component exists, deployment must declare
> `Δt_base` explicitly**, and the anchors validate against it.

The refusal is not a dead end: the engine has everything needed to make it
constructive — section 8 specifies the suggestion. Either way, validation failures
are collected `DeploymentInvalid`s in [§9.1][s9-1]'s style, and with anchors the diagnostic
gains error locality: "period `1//500` does not divide declared `Δt_base` — declared
`Absolute(Hz(500))` at `Avionics`, for `sensors`," straight from the provenance
column.

**Resolution.** With `Δt_base` in hand, one division pair per anchor and one
multiply-add per component:

```julia
# per anchor k:      D_k = T_k / Δt_base       Φ_k = τ_k / Δt_base   (both exact Ints)
# per component:     D   = m · D_k             Φ   = Φ_k + c · D_k
#                    Δt  = D · Δt_base
```

The residue invariant survives this last step too (`Φ_k < D_k` from `τ < T`, `c < m`
from the fold, so `Φ ≤ (D_k − 1) + (m − 1)·D_k = D − 1`). The output is the **bound
schedule**: per discrete component, `(D, Φ, Δt)` plus the anchor and provenance
columns carried through. This deserves to be a **named, printable artifact on the
`Simulation`** — it is the single source of truth for `Δt` ([§10.5][s10-5]), the substrate of
every diagnostic in section 8, and the table a user reads to answer "when does what
run, and what coincides with what."

### 5.4 The execution form: entry fields and the gate

At activation compile, the bound schedule bakes into the compiled executor, and
[§9.7][s9-7] has already made the load-bearing choice: *"an entry carries what selects
code — component type, stage — in type parameters, and what is plain data — tick
divisor, layout offsets — in fields."* `D`, `Φ` and `Δt` are exactly the second
kind. Two instances of one controller type at different rates and phases differ only
in field values, share an entry type, and compile to **one** body — the property the
cell-store benchmark ([D-162][d-162]) was run to protect. Schematically:

```julia
struct HxuBoundaryEntry{C, ...}   # C, stage → select code
    D::Int; Φ::Int                # gate data
    Δt::Float64                   # injected into the bundle it constructs
    offsets::...                  # cell addresses, as today
end
```

Inside the compile-time-unrolled boundary body, the gate is one test wrapped around
the entry's evaluation:

```julia
(idx - e.Φ) % e.D == 0 && eval_stage(e, table, stores)
```

The same gate wraps the entry's appearances in all three tick-sensitive blocks —
boundary `sweep_hx(idx)`, boundary `sweep_hxu(idx)`, and `ticks(idx)` — since
due-ness is per component, per boundary, not per stage. With the schedule tuple
concretely typed, those field loads constant-fold in practice; either way it is the
"one subtraction" over today's `idx % D == 0`. The interior (zero-arg) arities are
untouched: discrete entries are *absent* from them at compile time ([§10.5][s10-5]'s static
split), so the phase machinery adds literally nothing to the RK/localization hot
path. `Δt` semantics are likewise unchanged: the bundle's `Δt` is still
`D·Δt_base` — the offset shifts firing instants, never the period — so [§15.2][s15-2]'s
discretized laws and the never-cache-`Δt` rule are unaffected.

### 5.5 Boundary flavors: one gate, no special cases

The loop owns one `Int` tick counter, advanced at frame tops. All four boundary
flavors fall out of the same gate — the main evidence the representation is right:

- **Frame top**: the loop passes the frame's tick index; the gates evaluate; the due
  set is whatever passes.
- **Quiescence re-sweeps** ([§10.6][s10-6]): every re-sweep of the boundary passes the *same*
  index, so the gates return the same answers. [§10.5][s10-5]'s "the due set is computed once
  for the boundary" is a semantic fact — a pure function of an index that does not
  change within the boundary — not a memoized bitmask. Nothing to store, nothing to
  invalidate.
- **Boundary zero**: `idx = 0` makes the gate `(0 − Φ) % D == 0`, which by the
  canonical-residue invariant (`0 ≤ Φ < D`) holds **iff `Φ = 0`**. That is exactly
  the amended boundary-zero rule — "everything with `Φ = 0`" — falling out of the
  ordinary gate rather than being implemented anywhere. In today's no-offset design
  every `Φ` is 0, so "at boundary zero everything is due" is the degenerate case of
  the same identity. This is what the residue invariant of section 3 buys: it makes
  index 0 honest. An offset component's first tick is at `Φ·Δt_base`; until then its
  cells hold the values the build probe populated from `init_s` ([§9.3][s9-3]) — a coherent
  ZOH story, since those are exactly the values a tick at `t₀⁻` would have produced.
- **`t*` boundaries**: the empty due set is **arity selection, not an index trick**.
  There is no sentinel index that fails all gates — a `D = 1, Φ = 0` component
  passes at *every* index — so the `t*` event iteration runs the zero-arg interior
  sweep arities (which contain no discrete entries) rather than the boundary arities
  with a magic index. This is presumably *why* [§10.5][s10-5] phrases the rule as "the due set
  is a property of the boundary" rather than of the index: at `t*` the counter has
  not advanced, and the modulo image of the unadvanced index would wrongly re-admit
  the previous frame's due set. The two-arity structure [§9.7][s9-7] already commits to is
  the mechanism that implements the emptiness.

### 5.6 Why the due set is never materialized

For completeness, the two designs this displaces. A **precomputed due-set table over
the hyperperiod** (the schedule repeats with period `lcm(Dᵢ)`) sounds attractive
until you notice the lcm grows multiplicatively where the gcd-derived base grid
grows slowly — staggered prime-ish divisors blow it up fast — and it converts two
arithmetic ops into a dynamically-indexed load of runtime data. A **per-frame-top
recomputed mask** threads mutable state through the sweep for the privilege of
replacing a subtraction and a remainder with a load. Both lose to the pure gate.

The one place a materialized "due set as data" genuinely earns its way in is the
extensions charter's triggered-updates item (`framework_extensions.md` section 1.2,
item 7): a triggered entry's due-ness is *not* a function of the index — its gate
would read a staged trigger flag instead of doing modulo. The gate-per-entry
structure generalizes to that; a `(D, Φ)` table does not, and does not need to.

## 6. Worked example: the staggered sensor suite

The sensor suite runs at a 0.002 s period with a 0.001 s offset; the control laws at
0.02 s, offset 0:

```julia
sample_times(::Vehicle) = (
    sensors = Absolute(Period(1//500), 1//1000),   # 0.002 s, offset 0.001 s
    ctrl    = Absolute(Period(1//50)),             # 0.02 s
)
sample_times(::Sensors) = (imu = 1, gnss = Relative(20, 3))
```

Every discrete component here hangs from an anchor (`imu` and `gnss` under the
`sensors` anchor, `ctrl` its own), so this model takes the **derivation path**: the
engine may set `Δt_base` itself, printing the derivation (section 8).

**Grid derivation.** The pool is `{1//500, 1//1000, 1//50}`, giving
`Δt_base = 1//1000`. Note that neither period alone needs a 1 kHz grid —
`gcd(1//500, 1//50) = 1//500` — it is the offset, half the sensor period, that
refines it (leave-one-out factor 2, in section 8's terms). Deployment then needs
`h ≤ 1 ms`; take `h = 1 ms`, `n = 1`.

**Compiled pairs.** The `sensors` anchor resolves to `(D = 2, Φ = 1)`; `ctrl` to
`(D = 20, Φ = 0)`. Cascading into the scope with the composition law — the scope
arrives as `(D_s, Φ_s) = (2, 1)`, its grid `t = 1, 3, 5, … ms`, and that pair is
what the children inherit; a child never sees seconds, only scope ticks:

- **`imu`**: `Relative(1)` → `D = 2, Φ = 1`. Every scope tick: 1, 3, 5, … ms — the
  scope's declared rate *is* its fastest member, per the section 1 convention.
- **`gnss`**: `Relative(20, 3)` → `D = 40, Φ = 1 + 3·2 = 7`. Ticks at 7, 47, 87, …
  ms — 25 Hz, on the 3rd scope tick of each cycle of 20.

The `Sensors` author declared ratios and a tick-count stagger; they did not know —
could not know, preserving the [§10.5][s10-5] property — that a scope tick would turn out to
be 2 ms wide or displaced by 1 ms.

**Timeline** (base index in milliseconds):

```
t (ms):    0   1   2   3   4  ...  19  20  21  ...
ctrl:      C                           C
sensors:       S       S      ...  S       S
```

The suite's and ctrl's due sets are *disjoint at every boundary* — the designed
effect of a half-period stagger. Boundary by boundary:

- **t = 0**: due = {ctrl} (the boundary-zero rule of section 5.5). The sensors have
  never ticked; ctrl reads their build-probed `init_s` cells for the first
  millisecond.
- **Odd ms**: due = {sensors}; ctrl's cells hold.
- **Even ms, not multiples of 20**: due set empty — the boundary exists for events
  and frame structure only. These 1 kHz boundaries in a model whose fastest task is
  500 Hz are the overhead the offset-refined grid charges.
- **Multiples of 20 ms**: due = {ctrl}, whose `h_su` reads sensor cells written at
  the previous odd boundary — data exactly 1 ms old, every time, deterministically.

**What the offset bought.** With offset 0, sensors and ctrl coincide at every 20 ms
boundary and [§10.5][s10-5]'s simultaneous-tick rule gives ctrl same-instant sensor data via
topological order — the idealized synchronous-sampling picture. With the offset,
that becomes "pipelined, deterministically 1 ms stale" — which is precisely the
knob: a real ADC-plus-bus pipeline delivers data some fixed fraction of a sample
period old, and the offset expresses that latency *structurally*, in the schedule,
rather than through delay blocks. The second classic use is load shaping under
real-time pacing ([§10.7][s10-7]): with the stagger, no single 1 ms frame executes both the
sensor stack and the control stack, so worst-case frame cost is the `max` of the two
rather than their sum.

**The merged schedule.** At 7, 47, 87 … ms both `imu` and `gnss` are due — gnss's
instants are a subset of the scope grid, which imu fully occupies, so they always
coincide, and the simultaneous-tick rule orders them topologically within one sweep
(if imu feeds gnss, gnss reads same-boundary fresh imu output: a clean fusion
cascade with zero intra-suite latency). Meanwhile `ctrl`, reading at 20, 40, 60, 80
… ms, sees gnss data aged alternately 13 ms and 33 ms (writes at 7, 47, …) — not an
artifact but the honest schedule of a multi-rate system, and *fully determined at
build* from the `(D, Φ)` table. That determinism suggests a diagnostic: the [§5.6][s5-6]
feedthrough tracer could print per-wire data ages ("this connection delivers data
aged 13–33 ms"), turning the staleness question from a debugging session into a
table lookup.

## 7. Worked example: nested anchors and the refusal path

A root `Aircraft` with a continuous `Vehicle` and a discrete `Avionics` — a
root-relative chain, a mid-tree anchor, a phased relative *under* an anchor, and a
nested anchor with an offset, all in one tree:

```julia
sample_times(::Aircraft) = (avionics = Relative(1),)
sample_times(::Avionics) = (computing = Relative(2), sensors = Absolute(Hz(500)))
sample_times(::Sensors)  = (imu = Relative(2, 1), gnss = Absolute(Hz(10), 1//100))
```

(`vehicle` is continuous, appears in no `sample_times`, and none of this concerns it — it
lives in the interior sweep and integrates with `h`.)

### 7.1 Stratum A

Validation passes throughout (`0 ≤ 1 < 2` for imu; `0 ≤ 1//100 < 1//10` for gnss;
mid-tree `Absolute` legal per section 4). The fold:

| step | declaration | action | result |
|---|---|---|---|
| root scope | — | seed | `(A0, 1, 0)` |
| `avionics` | `Relative(1)` | `m = 1·1, c = 0 + 0·1` | `(A0, 1, 0)` |
| `computing` | `Relative(2)` | `m = 2·1, c = 0 + 0·1` | `(A0, 2, 0)` |
| `sensors` | `Absolute(Hz(500))` | **sever**, anchor `A1 = (1//500, 0)` | `(A1, 1, 0)` |
| `imu` | `Relative(2, 1)` | `m = 2·1, c = 0 + 1·1` | `(A1, 2, 1)` |
| `gnss` | `Absolute(Hz(10), 1//100)` | **sever**, anchor `A2 = (1//10, 1//100)` | `(A2, 1, 0)` |

Stored in the `Build`:

```
anchors:     A0 = base grid (deployment's)
             A1 = (T = 1//500, τ = 0)       — Absolute(Hz(500)) at Avionics, for sensors
             A2 = (T = 1//10, τ = 1//100)   — Absolute(Hz(10), 1//100) at Sensors, for gnss

component            anchor   m    c    provenance
avionics/computing   A0       2    0    Relative(1)@Aircraft · Relative(2)@Avionics
…/sensors/imu        A1       2    1    Relative(2,1)@Sensors
…/sensors/gnss       A2       1    0    (anchor is its own rate)
```

Note what the severing did: `computing`'s two-level relative chain composed through,
while `imu`'s pair is relative to `A1` and *no trace* of the `avionics` scope
survives in it — and `gnss`, a leaf declared absolute, is the degenerate anchor
whose subtree is itself.

### 7.2 Binding: refusal, suggestion, resolution

`computing` hangs from anchor 0, so the derive-vs-declare rule refuses to pick
`Δt_base` — but constructively:

```
pool = {1//500, 1//10, 1//100}       # A1.T, A2.T, A2.τ (zero offsets contribute nothing)
gcd  = 1//500                        # lcm(500, 10) = 500; lcm(500, 100) = 500
```

> `Δt_base` must be declared: unanchored discrete components exist
> (`avionics/computing`). The 2 anchor declarations admit `Δt_base = 1//500 s`
> (500 Hz) at coarsest, or any integer division of it.

The offset `1//100` happened to be grid-friendly (`5 × 1//500`) — it joined the pool
but did not refine the grid. Say the user declares `Δt_base = 1//500`, `h = 1//500`
(`n = 1`). Anchor resolution: `A1 → (D = 1, Φ = 0)`; `A2 → (D = 50, Φ = 5)`. (Had
the user declared `Δt_base = 1//300`: `(1//500)/(1//300) = 3//5`, not an integer →
`DeploymentInvalid` naming `A1`, *declared at `Avionics` for `sensors`*.) Then the
bound schedule:

```
component            D    Φ    Δt        rate     first tick
avionics/computing   2    0    1//250 s  250 Hz   t = 0
…/sensors/imu        2    1    1//250 s  250 Hz   t = 2 ms
…/sensors/gnss       50   5    1//10 s   10 Hz    t = 10 ms
```

`imu`'s `Δt` is 4 ms even though its first tick is at 2 ms — the phase shifts
instants, never the period. And `computing`'s 250 Hz is a **deployment artifact**:
declare `Δt_base = 1//1000` instead (also admissible) and `computing` silently
becomes 500 Hz while `imu` and `gnss` do not move — which is exactly why the engine
refuses to pick the value itself and only suggests.

### 7.3 The grid, drawn

With `n = 1`, every base tick is a frame top. Base index `k`, `t = 2k` ms;
hyperperiod `lcm(2, 2, 50) = 50` base ticks = 100 ms:

```
k:          0   1   2   3   4   5   6   7   8   9  10  ...  54  55  56
t (ms):     0   2   4   6   8  10  12  14  16  18  20  ...  108 110 112
computing:  C       C       C       C       C       C   ...  C       C
imu:            I       I       I       I       I       ...      I
gnss:                           G                        ...      G
```

- **k = 0** (boundary zero): due = {computing} — the only `Φ = 0` row, straight from
  the gate at index 0. `imu` and `gnss` hold their build-probed values until their
  first ticks.
- **Even k / odd k**: the two 250 Hz components are perfectly staggered, so no
  boundary executes both — the load-shaping effect of section 6, obtained here
  *accidentally*, purely from how the phases resolved. No boundary is ever empty.
- **k = 5, 55, 105, …**: due = {imu, gnss}. `Φ_gnss = 5` is odd and `D_gnss = 50` is
  even, so **every gnss tick lands on an imu tick and never on a computing tick** —
  parity is preserved forever. Nobody declared that; it fell out of GCD arithmetic
  and an offset. If gnss consumes imu, it reads same-boundary-fresh imu data at
  every one of its ticks — a happy outcome you would want to *verify*, not assume,
  which is precisely the printable table's job.

### 7.4 Variant: an awkward offset

Change one declaration: `gnss = Absolute(Hz(10), 1//150)`. Stratum A is untouched —
only `A2 = (1//10, 1//150)` changes. The offset now carries a prime factor (3) the
rest of the pool does not have:

```
pool = {1//500, 1//10, 1//150}
gcd(1//500, 1//150) = 1 // lcm(500, 150) = 1//1500      # 500 = 2²·5³, 150 = 2·3·5²
```

Coarsest admissible `Δt_base = 1//1500 s ≈ 0.667 ms` — a 3× refinement. Binding at
the suggestion (`h = 1//1500`, `n = 1`): `A1 → (3, 0)`, `A2 → (150, 10)`, and:

```
component            D    Φ    Δt        rate     first tick
avionics/computing   2    0    1//750 s  750 Hz   t = 0
…/sensors/imu        6    3    1//250 s  250 Hz   t = 2 ms
…/sensors/gnss       150  10   1//10 s   10 Hz    t = 1//150 s ≈ 6.67 ms
```

The comparison with 7.2 exposes four things:

- **Anchored components' physical schedules are `Δt_base`-invariant.** `imu` still
  ticks at 2, 6, 10 … ms, `gnss` still every 100 ms — their `(D, Φ)` rescaled in
  *tick units* but denote the same instants. The anchor contract working: seconds
  in, seconds out.
- **The unanchored `computing` rescaled again**: 250 → 750 Hz, purely because the
  suggested grid tripled and its period is `2·Δt_base`. A `1/300 s` edit to a GNSS
  *offset* three scopes away is what tripled the flight computer's rate — the
  refusal rule forces the user to look at this number, and if 250 Hz was a
  requirement, `computing` should be anchored too.
- **The coincidence structure silently rewired.** gnss's base indices `10 + 150k`
  are always even and never `≡ 3 (mod 6)`, so gnss now *always* coincides with
  `computing` and *never* with `imu` — the exact opposite of 7.3. If gnss consumes
  imu, its data went from same-boundary-fresh to deterministically 2/3 ms stale: an
  observable behavioral change caused by an offset edit, invisible without the
  table.
- **Empty boundaries appeared** (odd indices not `≡ 3 (mod 6)` — a third of all
  boundaries), and since `h ≤ Δt_base` the `Vehicle` integrator takes 3× the RHS
  evaluations per simulated second. That is the actual bill: **3× stepping
  throughput plus a rewired staleness pattern.**

The grid (base tick = 2/3 ms; hyperperiod `lcm(2, 6, 150) = 150` ticks = 100 ms):

```
k:          0    1    2    3    4    5    6    7    8    9   10   11   12
t (ms):     0   2/3  4/3   2   8/3 10/3   4  14/3 16/3   6  20/3 22/3  8
computing:  C         C         C         C         C         C         C
imu:                       I                             I
gnss:                                                         G
empty:           ·              ·              ·                   ·
```

## 8. Grid diagnostics: suggestion, attribution, warnings

### 8.1 Attribution must be global, and it can be exact

The instinctive rule — "warn when an offset is not a simple fraction of its
period" — is authoring guidance (section 3), but it cannot be the *engine's* test:
in section 7, `τ/T = 1/10` cost nothing and `τ/T = 1/15` cost 3×, though they are
comparably "non-simple." Demand is **relational**: `1//100`'s denominator was
already inside the pool's prime content (`100` divides `500`); `1//150` brought a
new prime. A declaration is expensive only relative to everything else declared, so
blame must be computed against the actual pool — which the engine can do *exactly*,
two ways:

- **Leave-one-out refinement factors.** For each pool entry `p`,
  `r_p = gcd(pool ∖ p) / gcd(pool)` — always an integer ≥ 1 (the full GCD divides
  every partial one, since `gcd(pool) = gcd(gcd(pool ∖ p), p)`), read as "how much
  coarser the grid would be without this entry." For section 7.4:
  `r(offset 1//150) = 3`, `r(period 1//500) = 10`, `r(period 1//10) = 1`. Note the
  honest surprise: leave-one-out fingers the *sensors anchor* as the bigger marginal
  driver (without it, `1//150` alone would suffice). Both entries are jointly
  responsible — which is the truth, and why the report lists all `r_p > 1` rather
  than crowning one culprit.
- **Prime attribution**, the sharper cut: `1/Δt_base = 1500 = 2²·3·5³`, and each
  prime power traces to the entry (or entries) whose denominator supplies it — `2²,
  5³` from the 500 Hz anchor, **the factor 3 from the gnss offset alone**. This
  pinpoints exactly what an edit changed: versus the `1//100` configuration, the
  regression is the single prime 3, owned by `1//150`. (For unit-fraction periods
  `1/Δt_base` is the lcm of denominators; general rationals work the same through
  `gcd(numerators)/lcm(denominators)`.) Both computations are one-liners over a pool
  of a handful of rationals.

### 8.2 Where it surfaces

Resist warning thresholds almost everywhere — the design already creates the right
moments:

1. **The suggestion message** (refusal path — guaranteed an audience, since the
   build will not proceed without a declared `Δt_base`): the coarsest admissible
   value, the admissible set (`gcd/k`), and the attribution table. "Coarsest
   `Δt_base = 1//1500`; drivers: sensors period `1//500` (×10), gnss offset
   `1//150` (×3)."
2. **The derivation info line** (fully-anchored path — the one place refinement can
   happen *silently*): always print the derived value with its drivers. This is
   where an actual advisory earns a threshold, and the principled metric is **grid
   utilization**: `min_i D_i`, the number of base ticks between the *fastest*
   component's ticks. It directly counts what a fine grid costs (empty boundaries,
   forced-small `h`); the charter's `T/1000`-offset horror case shows up as
   utilization ~1000 where section 7.4 is a mild 2. One caveat built in: a scope
   deliberately declared faster than its fastest member to buy stagger room
   (section 3, property 3 — legal and sometimes intended) inflates the metric, so
   the advisory should say "grid is N× finer than the fastest declared work" and
   name the drivers, not scold.
3. **The printable bound schedule** carries the attribution and utilization columns
   permanently, for post-hoc reading.
4. **A constructive companion**: when an offset is a driver, compute the nearest
   offsets that *would not* refine — the multiples of what the rest of the pool
   supports. For section 7.4: "on the `1//500` grid the admissible offsets near
   `1//150 ≈ 6.67 ms` are `3//500` (6 ms) and `4//500` (8 ms); declaring `3//500`
   keeps `Δt_base = 1//500`." If the 6.67 ms was physically meaningful, the user
   keeps it and pays knowingly; if it was "about 6–7 ms," the snap just saved 3×
   throughput. The diagnostic becomes a repair, not a complaint.

### 8.3 The hyperperiod chart

"Show me the first ticks" generalizes into a cheap, well-defined print utility on
the bound schedule, because the whole pattern is **periodic with the hyperperiod
`lcm(Dᵢ)`** — the gate is pure modulo arithmetic, so one hyperperiod is the complete
truth, not a sample. A `show`-style method rendering exactly the charts of sections
6–7 over `k = 0 … lcm(D) − 1` (with a guard for silly hyperperiods) answers "when
does what run, and what coincides with what" as a table lookup. It shares its
substrate with the per-wire data-age idea of section 6: both read the `(D, Φ)`
table, and together they turn the two classic multi-rate debugging sessions — "why
didn't it fire" and "why is my data stale" — into lookups.

## 9. Spec amendments this proposal requires

All small, all honest; collected so a future increment can check them off:

- **[§10.5][s10-5], rate declaration.** The relative-declaration paragraph gains the two
  registers and the severing rule (section 4); `K ≥ 1` is scoped to the relative
  register; the fastest-member convention is stated; the recorded limitation "no
  phase offsets in the first cut (no demonstrated use)" is replaced by a pointer
  here.
- **[§10.5][s10-5], boundary zero.** "At boundary zero the due set is everything" refines to
  "everything with `Φ = 0`" — implemented by nothing, per section 5.5. [§14.5][s14-5]'s
  boundary-zero transitions and the trim read-back need a sentence each
  acknowledging that offset components are not due there.
- **[§9.1][s9-1].** Stratum A's "compilation of relative multipliers into absolute
  divisors" becomes compilation into `(anchor, m, c)` triples; deployment binding
  gains the pool, the derive-vs-declare rule, and anchor resolution.
- **[§9.2][s9-2].** The `Build` gains the anchor table and the rate columns of section
  5.2; the bound schedule becomes a named printable artifact on the `Simulation`
  (section 5.3).
- **[§9.7][s9-7].** Entry fields gain `Φ`; the gating sentence becomes
  `(idx − Φ) % D == 0`.
- **`Δt` semantics**: explicitly unchanged — the bundle's `Δt` is still
  `D·Δt_base`; offsets shift firing instants, not the period.
- **Documentation next to [§10.5][s10-5]'s simultaneous-tick rule.** The choice between
  coincident ticks (fresh same-instant reads via topological order) and staggered
  ticks (pipelined, deterministically aged reads) is a modeling decision with
  observable consequences — sections 6 and 7.4 — and should be documented where a
  user debugging "why does my controller see stale data" will find it.
- **[Appendix C][sC]**: no new error kinds — `DeploymentInvalid` covers the new
  validation failures; one new *advisory* kind if the derivation path adopts the
  grid-utilization warning of section 8.2.

<!-- citation link definitions — generated by tools/linkify.jl; do not edit -->
[d-162]: framework_decisions.md#d-162--adopt-per-eltype-homogeneous-cell-stores-over-per-instance
[d-185]: framework_decisions.md#d-185--adopt-the-phased-two-register-sample-time-declaration
[d-187]: framework_decisions.md#d-187--make-the-bound-schedule-a-named-artifact-with-exact-grid-diagnostics
[s10-5]: framework_spec.md#105-multi-rate-tick-scheduling
[s10-6]: framework_spec.md#106-event-iteration-at-boundaries-to-quiescence-budgeted
[s10-7]: framework_spec.md#107-real-time-pacing
[s13-1]: framework_spec.md#131-reporting-policy-collect-the-checks-fail-the-evaluations-fast
[s14-5]: framework_spec.md#145-boundary-zero-an-ordinary-boundary-with-authored-incoming-transitions
[s14-8]: framework_spec.md#148-the-trim-service-solver-seam-scratch-stores-commit-and-report
[s15-2]: framework_spec.md#152-torture-tests-for-the-52-interfaces-pistonengine-and-the-fcs-pid-cascade
[s5-6]: framework_spec.md#56-diagnostics-feedthrough-tracing
[s8-7]: framework_spec.md#87-rate-scopes
[s9-1]: framework_spec.md#91-three-strata
[s9-2]: framework_spec.md#92-the-build-artifact
[s9-3]: framework_spec.md#93-probing-and-input-synthesis
[s9-7]: framework_spec.md#97-the-compiled-executor
[sC]: framework_spec.md#appendix-c-the-diagnostic-kind-set
