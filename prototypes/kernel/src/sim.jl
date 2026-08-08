# The simulation object, the stepper seam (§8.2) and the phase-body accessor
# (§12.7). The framework owns the loop; the one delegated operation is "advance
# the continuous state from `t` by `h`".

struct Simulation{T,S,B,CL}
    store::S
    xbuf::Vector{T}
    ẋbuf::Vector{T}
    clock::CL
    bodies::B
    layout::Layout
    spec::ModelSpec
    work::NTuple{5,Vector{T}}     # RK4 scratch: x₀, k₁..k₄
end

function Simulation(store, xbuf, ẋbuf, clock, bodies, layout, spec, ::Type{T}) where {T}
    work = ntuple(_ -> zeros(T, length(xbuf)), 5)
    Simulation{T,typeof(store),typeof(bodies),typeof(clock)}(
        store, xbuf, ẋbuf, clock, bodies, layout, spec, work)
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

"""The boundary sweep: restores signal-table consistency at an accepted step (§8.3)."""
@inline function boundary_sweep!(sim::Simulation, tick::Int)
    sim.bodies.sweep_hx(tick)
    sim.bodies.sweep_hxu(tick)
    nothing
end

# --- the stepper (§8.2): in-house fixed-step RK4 over the flat buffer ---------
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

"""Initialize: state at its declared values, table consistent, clock at `t₀`."""
function init!(sim::Simulation{T}; t₀::T = zero(T)) where {T}
    sim.clock.t = t₀
    boundary_sweep!(sim, 0)
    nothing
end

"""Advance to `t_end` in steps of `h`, restoring table consistency at each boundary."""
function run!(sim::Simulation{T}, t_end::T, h::T) where {T}
    tick = 0
    while sim.clock.t < t_end - eps(t_end)
        step!(sim, min(h, t_end - sim.clock.t))
        tick += 1
        boundary_sweep!(sim, tick)
    end
    nothing
end

# --- reading and writing the table outside the loop ---------------------------
# Path-addressed, dictionary-driven, deliberately off the measured path: these
# are boundary-time operations, never called from inside a phase body.

port(sim::Simulation, path::Symbol, name::Symbol) = gather(sim.store, sim.layout.addr[(path, name)])

function set_slot!(sim::Simulation, path::Symbol, face::Symbol, v)
    scatter!(sim.store, sim.layout.addr[(path, face)], v)
    nothing
end

function state(sim::Simulation{T}, path::Symbol) where {T}
    ci = index_of(sim.spec, path)
    d = declarations(sim.spec.comps[ci], T)
    off = 0
    for i in 1:(ci-1)
        off += nleaves(typeof(declarations(sim.spec.comps[i], T).x))
    end
    reconstruct(typeof(d.x), sim.xbuf, off)
end
