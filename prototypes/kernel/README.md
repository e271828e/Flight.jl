# Kernel prototype

The walking skeleton for the framework in `docs/notes/design/framework_spec.md`,
built to keepable standards and grown one increment at a time. Increment 1 (the
cell-store representation bench) lives in `../cellstore_bench` and stays frozen
there — D-162 cites its numbers.

**Increment 2 — the continuous tier walks.** A model builds, schedules itself
from its own feedthrough structure, integrates, and does it without allocating.

**Increment 3 — the discrete tier, at one rate.** Both tiers now share the
model: tier is read off the declaration shape, discrete state and modes live in
their own stores, `g` runs at boundaries after the output stages, and the ZOH
holds mid-step because the interior sweep has no discrete entries to gate.
Every discrete component sits at `D = 1, Φ = 0`; the grid is increment 4.

    julia --project=. check.jl

## What is real here

| piece | spec | file |
| --- | --- | --- |
| leaf walk: flatten / reconstruct / the activation retype | §7.1, §7.2 | `src/leaves.jl` |
| declaration layer, both tiers' arities, the bundle law, `probe_value` | §5.2, §8.2, §9.3 | `src/declare.jl` |
| per-eltype cell stores and the store bundle | §9.7, D-162 | `src/store.jl` |
| entries, chunked unrolled walk, the interior/boundary split | §9.7, §10.5 | `src/executor.jl` |
| tier classification, probe, feedthrough graph, layout, embed-accept | §8.2, §9.3, §5.3, D-166 | `src/build.jl` |
| `Simulation`, RK4, the boundary macro-sequence, `phase_bodies` | §10.2, §10.3, §9.7 | `src/sim.jl` |
| the coverage set: `Plant`, `Gain`, `Sum`, `DiscreteIntegrator`, `TickCounter`, `Smoother`, `WorkGain`, `ModedSource` | §5.2, §6.2, §7.3 | `src/library.jl` |

The properties the tests pin down, each of which is a spec claim rather than a
programming convenience:

- **The schedule is derived, not authored.** Stage-1 (`h_x`) ports carry no
  input dependence, so consuming one adds no edge to the feedthrough graph. The
  reference model's feedback path closes through the plant's `h_x` port and
  schedules as `sum → ctl → plant`; rewiring the same loop through the plant's
  stage-2 port makes it a genuine algebraic loop and the build rejects it by
  name.
- **Derivative completeness is structural.** `f` returns a value shaped like
  `init_x`; the probe checks the shape once and the runtime scatter into the
  flat `ẋ` block is then safe by construction. A forgotten field is impossible,
  not merely detectable.
- **The phase-body roster is fixed and total.** `phase_bodies(sim)` always
  returns all four bodies. A model with no discrete components still gets
  `ticks`, empty — it compiles to a no-op, and its `@ballocated` assertion
  passes vacuously, which is the point: consumers iterate the roster with no
  per-model branching.
- **Tier is read off the declaration shape, never announced.** For a stateful
  leaf the update law carries it (`f` continuous, `g` discrete); for a
  stateless one the contract arity does. Every other tier-implying declaration
  must agree, and disagreement names the offending one — `g` beside a
  two-argument `output_types`, `init_m` on a discrete leaf, both arities of one
  contract.
- **The ZOH is not implemented.** It is the absence of any way to change a
  discrete cell mid-step: the interior sweep is compiled from continuous
  entries alone, so the hot path carries no gating test, and a discrete cell
  cannot move across a step because nothing in that walk writes it.

**The store bundle is in, and with it the *plural* in "per-eltype stores".**
The bench that settled the representation (D-162) measured exactly one buffer;
`StoreBundle` now holds one `CellStore` per leaf element type, keyed by the
eltype's name. Selection is static — an address carries its port type, whose
leaf eltype names the buffer at compile time — so a deliberately pinned
`Float64` leaf (D-166) keeps a buffer of its own beside the `Dual` one instead
of being flattened into it as a zero-partial, which is what increment 2 had to
refuse by name. The bundle is one concrete type per model, keeping `Chunk`'s
store parameter single: chunk-type count, not model size, is what bounds the
compile curve D-162 blessed.

Correctness is checked against analytically integrated references — the
continuous loop by matrix exponential, the sampled loop by its exact ZOH
discretization `q[k+1] = Ad q[k] + Bd s[k]` — with a tolerance, never `==`
(D-163). The sampled-data reference is the sharp one: it only matches if the
hold is real *and* output stages run before updates within a boundary.

## What is deliberately absent

**Increment 4 — multi-rate:** the harmonic grid, the two-register `sample_times`
declaration (`Relative` composing affinely down the tree, `Absolute` severing
and re-seeding), its compilation to one `(D, Φ)` pair per component, the
`Relative(1)` default, rate scopes (§8.7) and the boundary sweep's
`(idx - Φ) % D` gate (D-185; the one-arg phase-body arity exists and is
exercised, but gates nothing yet). At `D = 1, Φ = 0` the gate is identically
true, so increment 4 must leave increment 3's tests passing unchanged.

Beyond those: events (guards, handlers, `project`, localization), hierarchy and
assemblies, computed connections, auto-published ports, §9.5's always-on
conformance check, §13.2's diagnostic framing (build errors here are a plain
`BuildError` with a good message, not the structured carrier), and the entire
runtime periphery.

**Two seams are collapsed rather than missing.** `build(spec; Δt)` stands in for
the deployment-binding stratum: `Δt` is entry-field data, so it must be fixed
before the executor exists, but in the design it arrives at binding time and the
`Build` artifact itself is deployment-free. And a frozen discrete component's
inputs are *synthesized* here, where the design carries the nominal activation's
cell contents across to the non-nominal one (§9.4). Both matter only in that a
reader should not mistake the prototype's shape for the spec's.

**One absence here is a refusal rather than a silence.** A cell must be
leaf-homogeneous. Mixed-leaf cells are legal in the design — a pinned `Float64`
beside `T` fields inside one declared struct — but their addresses need a
cursor per eltype, and this increment lays out one offset per cell. The layout
rejects them by name rather than mislaying them; multi-cursor addressing is a
scope cut, not doctrine.

## Authoring caveat found while building this

Declarations written in a **local scope** never reach the framework. Inside a
`let`, a function body or a `@testset`, `h_x(::MyComp, (; x)) = …` does not add
a method to the global `h_x`: it binds a *new local function* of that name.
Calls inside the block see it and look fine; the global generic function the
build dispatches on is untouched, so `build` sees a component that declares
nothing at all and proceeds silently. Ordinary top-level definitions — the
script, module and notebook-cell cases — are unaffected. `check.jl` defines its
malformed test components at top level for exactly this reason.

This is the local-scope sibling of the `using Flight` trap §8.1 already
documents, and §8.1's shadowing check does not reach it: there is no
parent-module binding to compare against, because the shadow is a local binding
that disappears with its block. Ratified as D-164: a component that declares
nothing and defines no stage is a build error — an inert component cannot be
intentional — spec'd as the inert-component check in §8.1's stage register,
raised as `DeadStage` at the probe (§9.3). Not built here.
