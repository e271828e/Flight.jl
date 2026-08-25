# The simulation object with its deployment binding (§9.1, §9.2), the stepper
# seam (§10.2) and the phase-body accessor (§9.7). The framework owns the loop;
# the one delegated operation is "advance the continuous state from `t` by `h`".

struct Simulation{T,S,B,CL,EV,M}
    store::S
    xbuf::Vector{T}
    ẋbuf::Vector{T}
    clock::CL
    bodies::B
    events::EV                    # the compiled event set with its §10.6 registers
    layout::Layout
    flat::Flat
    h::Float64                    # the continuous step, bound at deployment
    n::Int                        # steps per base tick: Δt_base = n·h (§10.5)
    Δt_base::Float64
    firing_budget::Int            # per-event firings per boundary (§10.6)
    localization_tol::Float64     # relative bracket-width stop (§10.4)
    localization_budget::Int      # t* boundaries permitted per frame (§10.4)
    has_localized::Bool           # any localized event compiled in: the frame loop's fast-path key
    sched::Vector{@NamedTuple{path::String, D::Int, Φ::Int, Δt::Float64}}   # the bound schedule (§9.2)
    sstores::Vector{Any}          # discrete state stores, by component index
    mstores::Vector{Any}          # mode stores, by component index
    stepper::M                    # the seam's backend (§10.2), owning its own scratch
    xnext::Vector{T}              # the retained arrival pair (§10.4): xₙ₊₁ saved before trials clobber
    ẋnext::Vector{T}              # the buffer, ẋₙ₊₁ paid only past a validated trigger
end

"""
    Simulation(build::Build, T = Float64; h, n = 1, Δt_base = nothing,
               method = RK4, firing_budget = 4, localization_tol = 1e-6,
               localization_budget = 8, chunk_size = 16)
    Simulation(root, T = Float64; …)

Deployment binding at construction (§9.1, §9.2): `Δt_base` from one of three
cross-validated sources — explicit (a `Rational` or `Period`/`Hz` value), the
`n·h` product (the default path), or GCD derivation over the anchors' constraint
pool, requested as `Δt_base = :derive` and permitted only with every discrete
component anchored. The scalar picks the activation the entries compile over —
the nominal one directly, any other via `activation(b, T)`'s cached Stratum-C
re-run (§9.4). The convenience form is *defined as* `Simulation(build(root), T;
…)`; entry compilation lives behind the binding because `Δt`, `D` and `Φ` are
entry data, and one `Build` backs many `Simulation`s.

`method` selects the integration backend across the stepper seam (§10.2): a
stepper type — `RK4`, the default, or `Heun` — materialized against the flat
buffer at binding like every other deployment product. The method is
trajectory-determining and grid-independent, exactly like the keywords below;
nothing outside the backend's own struct knows which one ran.

`firing_budget` is §10.6's deployment keyword: how many times each declared
event may fire at one boundary, an integer ≥ 1 defaulting to 4 — a legitimate
re-enable is one or two firings deep, a toggling FSM pair chatters without
bound, and 4 separates them without ever binding on a healthy model.

`localization_tol` and `localization_budget` are §10.4's: the relative bracket
width at which root-finding stops (positive, defaulting to 1e-6 — the event
time can only ever be as accurate as the `O(h⁴)` interpolant, so anything
tighter buys nothing while every trial evaluation costs a full sweep), and how
many localizations one frame admits (an integer ≥ 1 defaulting to 8 — a
legitimate multi-event frame needs three or four, chattering needs tens). All
three are trajectory-determining and grid-independent: they stand beside `h`
and `n`, validated here with their siblings, and enter none of the grid
arithmetic.
"""
function Simulation(b::Build, ::Type{T} = Float64; h = nothing, n = nothing,
                    Δt_base = nothing, method = RK4, firing_budget = 4,
                    localization_tol = 1e-6, localization_budget = 8,
                    chunk_size::Int = 16) where {T}
    method isa Type && method <: AbstractStepper ||
        throw(BuildError("method must be a stepper type — RK4 or Heun — got $method (§10.2)"))
    firing_budget isa Integer && firing_budget ≥ 1 ||
        throw(BuildError("firing_budget must be an integer ≥ 1, got $firing_budget (§10.6)"))
    localization_tol isa Real && localization_tol > 0 ||
        throw(BuildError("localization_tol must be a positive real, got $localization_tol (§10.4)"))
    localization_budget isa Integer && localization_budget ≥ 1 ||
        throw(BuildError("localization_budget must be an integer ≥ 1, got $localization_budget (§10.4)"))
    act = activation(b, T)
    bound = bind_schedule(b, h, n, Δt_base)
    c = compile(b, act, bound.D, bound.Φ, bound.Δt; chunk_size)
    stepper = method(T, length(c.xbuf))
    Simulation{T,typeof(c.store),typeof(c.bodies),typeof(c.clock),typeof(c.events),
               typeof(stepper)}(
        c.store, c.xbuf, c.ẋbuf, c.clock, c.bodies, c.events, act.layout, b.flat,
        bound.h, bound.n, bound.Δt_base, Int(firing_budget), Float64(localization_tol),
        Int(localization_budget), any(c.events.localized), bound.sched,
        c.sstores, c.mstores, stepper, zeros(T, length(c.xbuf)), zeros(T, length(c.xbuf)))
end

Simulation(root::AbstractComponent, ::Type{T} = Float64; kw...) where {T} =
    Simulation(build(root), T; kw...)

"""
    phase_bodies(sim)

The compiled bodies of the nominal activation, bound over this simulation's own
buffers — **these are the bodies the loop runs**, not re-derivations, which is
what makes the §7.5 measurement honest. The roster is fixed and total: a model
with no discrete components still gets `ticks`, empty, compiling to a no-op
whose `@ballocated` assertion passes vacuously, so consumers iterate uniformly
with no per-model branching.
"""
phase_bodies(sim::Simulation) = sim.bodies

# --- evaluation ---------------------------------------------------------------

"""
One RHS evaluation: *evaluating the RHS means running the sweep* (§5.3). The
interior variant of each sweep block, then the `f` block against the complete
fresh table. Leaves `ẋbuf` holding the derivative of whatever `xbuf` holds.
"""
@inline function evaluate!(sim::Simulation)
    sim.bodies.sweep_1()
    sim.bodies.sweep_2()
    sim.bodies.rhs()
    nothing
end

"""
The boundary macro-sequence at a base tick, final form (§5.3, §10.6):

> integrate → project → [sweep → guards → handlers] iterated to quiescence
> (under the firing budget) → all due `g` updates

Integration has just written the state, so projection runs first — between the
write and its decode; the event phase then iterates with the due set fixed for
the whole boundary, and the due updates run last, after quiescence, reading
post-transition values off the settled table. Output stages before updates, so
a discrete component's cells carry `y[k]` computed from `s[k]` while `g`
produces `s[k+1]` — the sampled-data recursion, ordered by construction rather
than by convention.
"""
@inline function boundary!(sim::Simulation, tick::Int)
    _projects!(sim.events)
    event_phase!(sim, tick)
    sim.bodies.ticks(tick)
    nothing
end

"""
The macro-sequence at a step boundary that is *not* a base tick: the tick
counter has not advanced, so the due set is empty — and emptiness is arity
selection, never a sentinel index failing every gate (§10.5, D-185). The
zero-arg interior bodies are exactly the boundary walk with every discrete
entry gated out, so consistency is restored and nothing discrete can move;
projection and the event phase run in full — every step boundary is a boundary
(§10.4).
"""
@inline function offtick_boundary!(sim::Simulation)
    _projects!(sim.events)
    event_phase!(sim, nothing)
    nothing
end

# One iteration round's sweep: the whole gated schedule (§10.6), in the due-set
# arity the boundary fixed — never the update laws, which wait for quiescence.
@inline _round!(sim::Simulation, tick::Int) =
    (sim.bodies.sweep_1(tick); sim.bodies.sweep_2(tick); nothing)
@inline _round!(sim::Simulation, ::Nothing) =
    (sim.bodies.sweep_1(); sim.bodies.sweep_2(); nothing)

"""
The event phase (§10.6): [sweep → guards → handlers] iterated to quiescence,
the fixed point where a round of handlers fires nothing. The first round always
runs — it is the boundary sweep that restores table consistency — so a model
with no events degenerates to exactly that sweep, with no register work at all.

An event fires in a round iff its predicate is observed holding, the sample
observed before it was not-holding, and its firing count for this boundary is
below the budget; at most one event fires per component per round, the first
eligible in declaration order — priority with re-decision. The one register
subtlety is D-191's: an eligible-but-blocked event's last-observed sample is
*not* overwritten, so the edge it presented stands into the next round. Budget
exhaustion degrades rather than throwing: the event's further edges at this
boundary are lost under a `FiringBudget` warning, at most once per event per
boundary, while every other event iterates untouched. At quiescence the prior
is updated unconditionally from the final samples — every prior an honest
observation of a settled boundary.
"""
function event_phase!(sim::Simulation, tick)
    es = sim.events
    _round!(sim, tick)
    n = length(es.prior)
    n == 0 && return nothing
    copyto!(es.last, es.prior)
    fill!(es.count, 0)
    fill!(es.warned, false)
    budget = sim.firing_budget
    while true
        _guards!(es)
        fill!(es.comp_fired, false)
        any_fired = false
        for i in 1:n
            edge = !es.last[i] && es.now[i]
            eligible = edge && es.count[i] < budget
            if edge && !eligible && !es.warned[i]
                es.warned[i] = true
                (path, name) = es.names[i]
                @warn "FiringBudget: event `$name` at `$path` has fired $budget times at " *
                      "the boundary t = $(sim.clock.t); its further edges here are lost (§10.6)"
            end
            firing = eligible && !es.comp_fired[es.owner[i]]
            es.fire[i] = firing
            if firing
                es.comp_fired[es.owner[i]] = true
                es.count[i] += 1
                any_fired = true
            end
            # A blocked edge stays unconsumed (D-191): only the eligible-but-
            # blocked sample stands; every other event takes this round's.
            (eligible && !firing) || (es.last[i] = es.now[i])
        end
        any_fired || break
        _fire!(es)
        _round!(sim, tick)
    end
    copyto!(es.prior, es.last)
    nothing
end

# --- the stepper seam, framework side (§10.2) ----------------------------------

"""
    step!(sim, h)

Advance the continuous state from `t` by `h`, delegated across the stepper
seam to whichever backend the deployment bound (stepper.jl). The seam is never
entered empty: with no continuous state the step degenerates to advancing `t`,
the backend is simply not called, and no backend contract has to say what it
would do at N = 0.
"""
@inline function step!(sim::Simulation, h)
    isempty(sim.xbuf) ? (sim.clock.t += h) : step!(sim.stepper, sim, h)
    nothing
end

"""
Initialize: state at its declared values, table consistent, clock at `t₀`.

Boundary zero is an ordinary boundary with an empty integrate (§10.5, §14.5):
the gate at index 0 admits exactly the components with `Φ = 0`, implemented by
nothing. An offset component holds its probe-populated cells until its first
tick at `Φ·Δt_base`.

Boundary zero also establishes every event prior as not-holding (§10.6), so a
predicate already holding in the authored state fires at `t₀` — derived, not
asserted — and a warm restart (`init!` again) resets all three registers from
scratch: such predicates fire again at the new `t₀`.
"""
function init!(sim::Simulation{T}; t₀::T = zero(T)) where {T}
    sim.clock.t = t₀
    sim.clock.t₀ = t₀
    sim.clock.step = 0
    fill!(sim.events.prior, false)
    boundary!(sim, 0)
    nothing
end

"""
Advance to `t_end` one frame at a time — `frame!` carries each grid step
`[tₖ₋₁, tₖ]` through the §10.4 localization loop, firing any `t*` boundaries it
brackets on the way — then the frame-top boundary (§10.3): a base tick every
`n` frames, where the gate reads the tick index, and the empty-due-set boundary
in between. The grid is driven by the step counter, so `t_end` is taken to the
nearest frame top; partial advance is the periphery's (§12.6) and absent here.
"""
function run!(sim::Simulation{T}, t_end::Real) where {T}
    target = round(Int, Float64(t_end) / sim.h)
    while sim.clock.step < target
        k = (sim.clock.step += 1)
        frame!(sim, k)
        k % sim.n == 0 ? boundary!(sim, k ÷ sim.n) : offtick_boundary!(sim)
    end
    nothing
end

# --- reading and writing the table outside the loop ---------------------------
# Path-addressed, dictionary-driven, deliberately off the measured path: these
# are boundary-time operations, never called from inside a phase body. Paths are
# §8.6's canonical strings, and an assembly's own path addresses its faces —
# which resolve to the cells they derive from.

port(sim::Simulation, path::String, name::Symbol) = gather(sim.store, sim.layout.addr[(path, name)])

"""
Write a root slot, by the root input face's name (§11.3): the write surface is
the root's own faces, the one terminal fed by no component.
"""
function set_slot!(sim::Simulation, face::AbstractString, v)
    scatter!(sim.store, sim.layout.addr[("", Symbol(face))], v)
    nothing
end

"""
State at `path`, from whichever home owns it: `x` in the flat buffer on the
continuous tier, `s` in the component's own store on the discrete one (§7.3).
"""
function state(sim::Simulation{T}, path::String) where {T}
    ci = index_of(sim.flat, path)
    sim.sstores[ci] === nothing || return sim.sstores[ci][]
    _tier(i) = classify_tier(sim.flat.paths[i], sim.flat.comps[i])
    _decls(i) = declarations(sim.flat.comps[i], _tier(i), T)
    off = 0
    for i in 1:(ci-1)
        _tier(i) === CONTINUOUS && (off += nleaves(typeof(_decls(i).x)))
    end
    reconstruct(typeof(_decls(ci).x), sim.xbuf, off)
end

"""Modes at `path` (§7.3). Read-only here: modes are written by handlers alone."""
modes(sim::Simulation, path::String) = sim.mstores[index_of(sim.flat, path)][]
