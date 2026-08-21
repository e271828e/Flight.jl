# The simulation object with its deployment binding (§9.1, §9.2), the stepper
# seam (§10.2) and the phase-body accessor (§9.7). The framework owns the loop;
# the one delegated operation is "advance the continuous state from `t` by `h`".

struct Simulation{T,S,B,CL}
    store::S
    xbuf::Vector{T}
    ẋbuf::Vector{T}
    clock::CL
    bodies::B
    layout::Layout
    flat::Flat
    h::Float64                    # the continuous step, bound at deployment
    n::Int                        # steps per base tick: Δt_base = n·h (§10.5)
    Δt_base::Float64
    sched::Vector{@NamedTuple{path::String, D::Int, Φ::Int, Δt::Float64}}   # the bound schedule (§9.2)
    sstores::Vector{Any}          # discrete state stores, by component index
    mstores::Vector{Any}          # mode stores, by component index
    work::NTuple{5,Vector{T}}     # RK4 scratch: x₀, k₁..k₄
end

"""
    Simulation(build::Build, T = Float64; h, n = 1, Δt_base = nothing, chunk_size = 16)
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
"""
function Simulation(b::Build, ::Type{T} = Float64; h = nothing, n = nothing,
                    Δt_base = nothing, chunk_size::Int = 16) where {T}
    act = activation(b, T)
    bound = bind_schedule(b, h, n, Δt_base)
    c = compile(b, act, bound.D, bound.Φ, bound.Δt; chunk_size)
    work = ntuple(_ -> zeros(T, length(c.xbuf)), 5)
    Simulation{T,typeof(c.store),typeof(c.bodies),typeof(c.clock)}(
        c.store, c.xbuf, c.ẋbuf, c.clock, c.bodies, act.layout, b.flat,
        bound.h, bound.n, bound.Δt_base, bound.sched, c.sstores, c.mstores, work)
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
    sim.bodies.sweep_hx()
    sim.bodies.sweep_hxu()
    sim.bodies.rhs()
    nothing
end

"""
The boundary macro-sequence at a base tick (§10.3, §10.6 in miniature): the
boundary sweep restores signal-table consistency with the due set gated off
`tick`, then the due updates run. Output stages before updates, so a discrete
component's cells carry `y[k]` computed from `s[k]` while `g` produces `s[k+1]`
— the sampled-data recursion, ordered by construction rather than by convention.
"""
@inline function boundary!(sim::Simulation, tick::Int)
    sim.bodies.sweep_hx(tick)
    sim.bodies.sweep_hxu(tick)
    sim.bodies.ticks(tick)
    nothing
end

"""
The macro-sequence at a step boundary that is *not* a base tick: the tick
counter has not advanced, so the due set is empty — and emptiness is arity
selection, never a sentinel index failing every gate (§10.5, D-185). The
zero-arg interior bodies are exactly the boundary walk with every discrete
entry gated out, so consistency is restored and nothing discrete can move.
"""
@inline function offtick_boundary!(sim::Simulation)
    sim.bodies.sweep_hx()
    sim.bodies.sweep_hxu()
    nothing
end

# --- the stepper (§10.2): in-house fixed-step RK4 over the flat buffer ---------
# One-step method, arbitrary `h`, zero allocation. The seam is never entered
# empty: with no continuous state the step degenerates to advancing `t`.

function step!(sim::Simulation{T}, h::T) where {T}
    x, ẋ = sim.xbuf, sim.ẋbuf
    n = length(x)
    if n == 0
        sim.clock.t += h
        return nothing
    end
    x₀, k₁, k₂, k₃, k₄ = sim.work
    t₀ = sim.clock.t
    copyto!(x₀, x)

    evaluate!(sim); copyto!(k₁, ẋ)
    _advance!(x, x₀, k₁, h / 2); sim.clock.t = t₀ + h / 2
    evaluate!(sim); copyto!(k₂, ẋ)
    _advance!(x, x₀, k₂, h / 2)
    evaluate!(sim); copyto!(k₃, ẋ)
    _advance!(x, x₀, k₃, h); sim.clock.t = t₀ + h
    evaluate!(sim); copyto!(k₄, ẋ)

    @inbounds for i in 1:n
        x[i] = x₀[i] + (h / 6) * (k₁[i] + 2k₂[i] + 2k₃[i] + k₄[i])
    end
    nothing
end

@inline function _advance!(x, x₀, k, dt)
    @inbounds for i in eachindex(x)
        x[i] = x₀[i] + dt * k[i]
    end
end

"""
Initialize: state at its declared values, table consistent, clock at `t₀`.

Boundary zero is an ordinary boundary with an empty integrate (§10.5, §14.5):
the gate at index 0 admits exactly the components with `Φ = 0`, implemented by
nothing. An offset component holds its probe-populated cells until its first
tick at `Φ·Δt_base`.
"""
function init!(sim::Simulation{T}; t₀::T = zero(T)) where {T}
    sim.clock.t = t₀
    sim.clock.step = 0
    boundary!(sim, 0)
    nothing
end

"""
Advance to `t_end` in steps of the deployment's `h`, one boundary per step
(§10.3): a base tick every `n` steps, where the gate reads the tick index, and
the empty-due-set boundary in between. The grid is driven by the step counter,
so `t_end` is taken to the nearest step boundary; partial advance is the
periphery's (§12.6) and absent here.
"""
function run!(sim::Simulation{T}, t_end::Real) where {T}
    h = T(sim.h)
    target = round(Int, Float64(t_end) / sim.h)
    while sim.clock.step < target
        step!(sim, h)
        k = (sim.clock.step += 1)
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

"""Modes at `path` (§7.3). Read-only here: handlers are the event system."""
modes(sim::Simulation, path::String) = sim.mstores[index_of(sim.flat, path)][]
