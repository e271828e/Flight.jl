# Kernel prototype

The walking skeleton for the framework in `docs/notes/design/framework_spec.md`,
built to keepable standards and grown one increment at a time. Increment 1 (the
cell-store representation bench) lives in `../cellstore_bench` and stays frozen
there — decision row 162 cites its numbers.

**Increment 2 — the continuous tier walks.** A model builds, schedules itself
from its own feedthrough structure, integrates, and does it without allocating.

    julia --project=. check.jl

## What is real here

| piece | spec | file |
| --- | --- | --- |
| leaf walk: flatten / reconstruct / the activation retype | §7.1, §7.2 | `src/leaves.jl` |
| declaration layer, the bundle law, `probe_value` | §5.2, §11.2, §12.3 | `src/declare.jl` |
| per-eltype cell store with offsets in entry fields | §12.7, row 162 | `src/store.jl` |
| entries, chunked unrolled walk, phase bodies in both arities | §12.7 | `src/executor.jl` |
| probe, feedthrough graph, algebraic-loop rejection, layout | §12.3, §5.3, §5.5 | `src/build.jl` |
| `Simulation`, RHS evaluation, RK4, `phase_bodies` | §8.2, §8.3, §12.7 | `src/sim.jl` |
| the coverage set: `Plant`, `Gain`, `Sum` | §5.2, §6.2 | `src/library.jl` |

Three properties the tests pin down, each of which is a spec claim rather than a
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
  returns all four bodies. This model has no discrete components, so `ticks` is
  empty — it compiles to a no-op, and its `@ballocated` assertion passes
  vacuously, which is the point: consumers iterate the roster with no per-model
  branching.

Correctness is checked against an analytically integrated closed loop (matrix
exponential), with a tolerance — never `==` (row 163).

## What is deliberately absent

Increment 3 and beyond: the discrete tier (`z`, `h_z`/`h_zu`/`g`, modes,
workspace), multi-rate tick scheduling and the boundary sweep's gating (the
one-arg phase-body arity exists and is exercised, but gates nothing yet), events
(guards, handlers, `project`, localization), hierarchy and assemblies, computed
connections, auto-published ports, §12.5's always-on conformance check, §13.2's
diagnostic framing (build errors here are a plain `BuildError` with a good
message, not the structured carrier), and the entire runtime periphery.

## Authoring caveat found while building this

Declarations written in a **local scope** never reach the framework. Inside a
`let`, a function body or a `@testset`, `h_x(::MyComp, (; x)) = …` does not add
a method to the global `h_x`: it binds a *new local function* of that name.
Calls inside the block see it and look fine; the global generic function the
build dispatches on is untouched, so `build` sees a component that declares
nothing at all and proceeds silently. Ordinary top-level definitions — the
script, module and notebook-cell cases — are unaffected. `check.jl` defines its
malformed test components at top level for exactly this reason.

This is the local-scope sibling of the `using Flight` trap §11.1 already
documents, and §11.1's shadowing check does not reach it: there is no
parent-module binding to compare against, because the shadow is a local binding
that disappears with its block. Recorded as row 164 with a candidate mitigation
(reject a component that declares nothing and defines no stage — an inert
component cannot be intentional), not built here.
