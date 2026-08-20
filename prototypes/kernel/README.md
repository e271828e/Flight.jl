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
Every discrete component sits at `D = 1, Φ = 0`; the grid arrives in
increment 5, after assemblies.

**Increment 4 — hierarchy and assemblies.** `build` takes the root component
instance. Class is read off declaration shape, children are struct fields and
container elements, and the three connection declarations resolve against
slash-separated paths into flat primitives with one producer per input.
Assemblies are virtual for execution: what reaches the executor is exactly what
increment 3 already ran. The whole-tree obligation model replaces the silent
unwired face — an input fed by nothing is a build error, and the root
assembly's own input faces are the slots.

**Increment 5 — the multi-rate grid.** The two-register `sample_times`
declaration folds to one `(anchor, m, c)` triple per component in Stratum A;
deployment binding — now a real seam, `Simulation(build(root); h, n, Δt_base)`
— resolves the triples to `(D, Φ, Δt)` by exact GCD arithmetic and compiles the
executor behind it, one deployment-free `Build` backing many `Simulation`s. The
boundary sweep gates discrete entries by `(idx − Φ) % D`, and a step boundary
that is not a base tick runs the zero-arg bodies: an empty due set is arity
selection, never a sentinel index failing every gate (D-185).

    julia --project=. check.jl

## What is real here

| piece | spec | file |
| --- | --- | --- |
| leaf walk: flatten / reconstruct / the activation retype | §7.1, §7.2 | `src/leaves.jl` |
| declaration layer, both tiers' arities, the bundle law, `probe_value`, `AbstractComponent`, the three connection declarations, and the two rate registers (`Relative`/`Absolute` over `Period`/`Hz`, plain data carriers) with `sample_times` | §5.2, §8.2, §8.6, §8.7, §9.3, D-185 | `src/declare.jl` |
| class by declaration shape, children and containers, paths and the reach rule, endpoint and face resolution, the flatten pass with the obligation model and the sample-time fold to `(anchor, m, c)` triples | §8.5, §8.6, §6.1, §9.1 | `src/assembly.jl` |
| per-eltype cell stores and the store bundle | §9.7, D-162 | `src/store.jl` |
| entries, chunked unrolled walk, the interior/boundary split, the `(idx − Φ) % D` gate on discrete boundary entries | §9.7, §10.5 | `src/executor.jl` |
| tier classification, probe, feedthrough graph, layout, embed-accept, the `Build` artifact, the activation seam (nominal at build, `activation(b, T)` as a cached Stratum-C re-run, frozen products carried across, the eager `activations` keyword), deployment binding (three cross-validated `Δt_base` sources, the GCD constraint pool, per-anchor division pairs) and per-deployment entry compilation | §8.2, §9.3, §5.3, §9.1, §9.2, §9.4, D-166 | `src/build.jl` |
| `Simulation` with its deployment constructor and bound schedule, RK4, the boundary macro-sequence at ticks and off ticks, `phase_bodies`, the path- and face-addressed table accessors | §10.2, §10.3, §10.5, §9.7, §11.3 | `src/sim.jl` |
| the coverage set: `Plant`, `Gain`, `Sum`, `DiscreteIntegrator`, `TickCounter`, `Smoother`, `WorkGain`, `ModedSource`, `ZOH`, `Ramp`, plus `Group` and the named `SampledLoop`/`Vehicle`/`MultiRate` assemblies | §5.2, §6.2, §7.3, §8.5, §10.5 | `src/library.jl` |

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
- **A path reaches exactly as far as the declaration knows.** The reach rule is
  about the declaring type's knowledge, not the instance's: the same
  `"inner/sum/a"` against the same `SampledLoop` value resolves when the field
  is declared `::SampledLoop` and is a build error when it is declared `::L`,
  because deep-routing *past* a generically held child hard-codes one
  implementation. Resolving *to* one is face-level access and legal (§13.3): a
  route through concretely-declared structure may end at a generic child's
  face, and a single hop always may — which is what lets a container's
  elements be addressed by their parent's own declarations.
- **Being fed is a whole-tree obligation, not a per-declaration one.** An unfed
  child input inside one assembly is merely awaiting a claim from above — a
  sibling wire, an ancestor's deep route, or an `input_connections` entry
  handing it up a level. The error fires at the root, for the chain that never
  terminates; the one legitimate terminus is a root input face, which is a
  slot. The one-producer rule spans levels the same way, so an ancestor's route
  onto an input a sub-assembly already wires is caught where the two claims
  meet, with both entries named.
- **A face is derived, never declared.** It owns no cell: it resolves — through
  as many levels of re-export as there are — to its ultimate internal endpoint,
  and takes that endpoint's type and tier. `Vehicle`'s `y` and `cmd` are one
  assembly's faces sourced from the two tiers, and at a `Dual` activation `y`
  walks while `cmd` stays pinned, with nothing anywhere declaring either.
- **An anchored divisor does not exist until `Δt_base` binds.** The `Build`
  carries `(anchor, m, c)` triples against exact rational anchors, and the same
  `Build` lands `MultiRate`'s gnss at `D = 10` or `D = 5` depending on the
  deployment keywords — with nothing writable shared between the two
  `Simulation`s, because each materializes its own stores and buffers. The
  worked example compiles to exactly §9.2's table.
- **Due-ness is one subtraction and one remainder, at boundaries only.** The
  interior walk is untouched by increment 5 — no index reaches it — and a step
  boundary that is not a base tick runs the zero-arg bodies, which *are* the
  boundary walk with every discrete entry gated out. The multi-rate sampled
  loop only matches its exact discretization at `Δt_ctl = D·Δt_base = 4h` if
  the gate admits the controller at exactly its own ticks, the hold spans the
  sub-ticks, and the bundle's `Δt` is the compiled schedule's — one test, three
  ways to fail it.
- **A cell holds what the build probe populated until a sweep first writes it.**
  Entry compilation seeds every cell from the probe products, so an offset
  component's pre-first-tick reads and a frozen component's pinned cells are
  the same story (§10.5, §9.3) rather than a lucky zero.
- **The nominal activation runs at build; any other is a Stratum-C re-run.**
  `activation(b, T)` materializes at first request and caches on the `Build`
  (§9.4); structure and schedule are computed once and shared. A frozen
  component sits outside the non-nominal executable set: its stages are never
  probed there, and its complete products — cells included — are carried
  across from the nominal activation, holding what a tick at `t₀⁻` computed
  from real upstream values rather than anything synthesized. That is what
  makes a discrete stage reading `t` lawful (it is `Float64` whenever the
  stage actually runs), and why a pinned-leaf lurk detonates at the `Dual`
  activation — or eagerly, via `build(root; activations = (Float64, D8))`,
  the CI idiom — rather than at `build`.

**The store bundle is in, and with it the *plural* in "per-eltype stores".**
The bench that settled the representation (D-162) measured exactly one buffer;
`StoreBundle` now holds one `CellStore` per leaf element type, keyed by the
eltype's name. Selection is static — an address carries its port type, whose
leaf eltypes name the buffers at compile time, with one cursor per eltype in
address *fields* — so a deliberately pinned `Float64` leaf (D-166) keeps a
buffer of its own beside the `Dual` one instead of being flattened into it as
a zero-partial, which is what increment 2 had to refuse by name, and a
mixed-leaf cell (an `Int` tag beside `T` leaves, or a pin inside a declared
struct) simply spans several buffers, its leaves each living with their own
kind. The `K = 1` homogeneous cell is a case, not a special case; the C2M
point in `../cellstore_bench` (2026-08-20) re-confirmed D-162's flat compile
curve for exactly this address shape. The bundle is one concrete type per model, keeping `Chunk`'s
store parameter single: chunk-type count, not model size, is what bounds the
compile curve D-162 blessed.

Correctness is checked against analytically integrated references — the
continuous loop by matrix exponential, the sampled loop by its exact ZOH
discretization `q[k+1] = Ad q[k] + Bd s[k]` — with a tolerance, never `==`
(D-163). The sampled-data reference is the sharp one: it only matches if the
hold is real *and* output stages run before updates within a boundary.

## What is deliberately absent

**D-187's grid diagnostics:** the bound schedule here is plain data on the
`Simulation` (path, `D`, `Φ`, `Δt`), not the named printable artifact — no
hyperperiod-chart `show`-form, no anchor/provenance columns carried through, no
leave-one-out refinement factors or prime attribution, no nearest-non-refining
offset suggestions, no `GridUtilization` advisory. The refusal messages name
the offending anchor and the pool's GCD, and stop there.

Beyond that: events (guards, handlers, `project`, localization), computed
connections and the §8.8 helpers (`input_passthrough` and the generic-holding
sugar), visibility enforcement (§8.3), did-you-mean suggestion lists (an error
names the offender plainly), auto-published ports, §9.5's always-on conformance
check, §13.2's diagnostic framing (build errors here are a plain `BuildError`
with a good message, not the structured carrier), partial advance and the
per-run overrides (§12.6 — `run!(sim, t_end)` takes `t_end` to the nearest step
boundary), and the entire runtime periphery.

## Stand-ins: where the prototype's shape is not the spec's

**Rule: nothing deviates silently.** Every construct a reader could mistake
for the design's is accounted for in exactly one of three places: the "what is
real" table, the absence list above, or a row here naming the spec shape it
replaces. The table is transactional — an increment that introduces a stand-in
adds its row in the same commit, and one that retires a stand-in deletes it. A
deviation in none of the three places is a defect in this file, not a liberty
the prototype has. `prototypes/` sits outside the design tools' checked
rosters, so no battery enforces this section: the diff review is the
enforcement.

The table is currently empty: every stand-in introduced through increment 5
was retired on 2026-08-20 (the `Tuple` container, mixed-leaf cells, the
whole-rebuild activation and the synthesized frozen inputs). The rule stands
for whatever the next increment brings — a new stand-in re-creates the table
with its row, in the introducing commit.

## Authoring caveat found while building this

Declarations written in a **local scope** never reach the framework. Inside a
`let`, a function body or a `@testset`, `h_x(::MyComp, (; x)) = …` does not add
a method to the global `h_x`: it binds a *new local function* of that name.
Calls inside the block see it and look fine; the global generic function the
build dispatches on is untouched, so `build` sees a component that declares
nothing at all. Ordinary top-level definitions — the
script, module and notebook-cell cases — are unaffected. `check.jl` defines its
malformed test components at top level for exactly this reason.

This is the local-scope sibling of the `using Flight` trap §8.1 already
documents, and §8.1's shadowing check does not reach it: there is no
parent-module binding to compare against, because the shadow is a local binding
that disappears with its block. Ratified as D-164: a component that declares
nothing and defines no stage is a build error — an inert component cannot be
intentional — spec'd as the inert-component check in §8.1's stage register,
raised as `DeadStage` at the probe (§9.3). Increment 4 catches the case where
the trap lands, one stratum earlier and under the other rule: a component that
declares nothing has no *class* to read either (§8.5), so the build names both
families rather than reaching the probe at all. `DeadStage` itself is not built.
