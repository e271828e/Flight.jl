# The simulation object, the stepper seam (§10.2) and the phase-body accessor
# (§9.7). The framework owns the loop; the one delegated operation is "advance
# the continuous state from `t` by `h`".

struct Simulation{T,S,B,CL}
    store::S
    xbuf::Vector{T}
    ẋbuf::Vector{T}
    clock::CL
    bodies::B
    layout::Layout
    flat::Flat
    Δt::Float64                   # the base tick period; 0.0 with no discrete tier
    xstores::Vector{Any}          # discrete state stores, by component index
    mstores::Vector{Any}          # mode stores, by component index
    work::NTuple{5,Vector{T}}     # RK4 scratch: x₀, k₁..k₄
end

function Simulation(store, xbuf, ẋbuf, clock, bodies, layout, flat, Δt, xstores, mstores,
                    ::Type{T}) where {T}
    work = ntuple(_ -> zeros(T, length(xbuf)), 5)
    Simulation{T,typeof(store),typeof(bodies),typeof(clock)}(
        store, xbuf, ẋbuf, clock, bodies, layout, flat, Δt, xstores, mstores, work)
end

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
The boundary macro-sequence (§10.3, §10.6 in miniature): the boundary sweep
restores signal-table consistency, then the due updates run. Output stages
before updates, so a discrete component's cells carry `y[k]` computed from
`x[k]` while `g` produces `x[k+1]` — the sampled-data recursion, ordered by
construction rather than by convention.
"""
@inline function boundary!(sim::Simulation, tick::Int)
    sim.bodies.sweep_hx(tick)
    sim.bodies.sweep_hxu(tick)
    sim.bodies.ticks(tick)
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
everything at `Φ = 0` is due, which at one rate is everything, so the discrete
tier takes its first tick here.
"""
function init!(sim::Simulation{T}; t₀::T = zero(T)) where {T}
    sim.clock.t = t₀
    boundary!(sim, 0)
    nothing
end

"""
Advance to `t_end` in steps of `h`, running the macro-sequence at each boundary.

With a discrete tier every step boundary is a tick, so `h` must be the base tick
period: the harmonic grid at one rate, `Δt_base = h` (§10.5). Increment 5 lifts
this to `Δt_base = n·h` with the gate.
"""
function run!(sim::Simulation{T}, t_end::T, h::T) where {T}
    sim.Δt == 0.0 || h ≈ sim.Δt ||
        throw(ArgumentError("this model ticks at Δt = $(sim.Δt); step it at that period, " *
                            "not $h (§10.5)"))
    tick = 0
    while sim.clock.t < t_end - eps(t_end)
        step!(sim, min(h, t_end - sim.clock.t))
        tick += 1
        boundary!(sim, tick)
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
State at `path`, from whichever home owns it: the flat buffer on the continuous
tier, the component's own store on the discrete one (§7.3).
"""
function state(sim::Simulation{T}, path::String) where {T}
    ci = index_of(sim.flat, path)
    sim.xstores[ci] === nothing || return sim.xstores[ci][]
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
