# The synthetic component set and the model-spec generator.
#
# Deliberately not the real declaration layer (§8.2): declarations here take
# the activation scalar `T` explicitly instead of being walked out of nominal
# `Float64` declarations, because the walk is not what this bench measures.
# Everything else is faithful: stage functions take state and input views and
# return named tuples of port values (§4.3, §5.2), and the port vocabulary
# covers a scalar, a static vector, and a custom isbits struct — the three
# shapes a cell store has to carry.

using StaticArrays
using LinearAlgebra

# --- a custom struct port value ----------------------------------------------

struct Pose3{T}
    p::SVector{3,T}
    ψ::T
end

# --- component types ----------------------------------------------------------
# Declarations: `input_types(c, T)`, `output_types(c, T)`, `init_x(c, T)`.
# Stage function: `h_xu(c, x, u)` — the interior continuous stage, the body the
# `sweep_2()` measurement (§9.7) runs.

"""Root source: no inputs, drives the chain. One instance per model."""
struct Source
    ω::Float64
end

input_types(::Source, ::Type{T}) where {T} = NamedTuple()
output_types(::Source, ::Type{T}) where {T} =
    (cmd = T, wind_a = SVector{3,T}, wind_b = SVector{3,T}, pose = Pose3{T})
init_x(::Source, ::Type{T}) where {T} = (φ = T(0.1),)

@inline function h_xu(c::Source, x, u)
    s, cs = sincos(c.ω * x.φ)
    cmd = 0.5 * s + 0.25 * cs
    wind_a = SVector(s, cs, s * cs)
    wind_b = SVector(cs - s, s + cs, 0.5 * s)
    pose = Pose3(SVector(x.φ, s, cs), c.ω * cs)
    (; cmd, wind_a, wind_b, pose)
end

"""Scalar lag: the one-in one-out scalar shape. One instance per model."""
struct Filter
    k::Float64
end

input_types(::Filter, ::Type{T}) where {T} = (in = T,)
output_types(::Filter, ::Type{T}) where {T} = (out = T,)
init_x(::Filter, ::Type{T}) where {T} = (z = T(0.0),)

@inline h_xu(c::Filter, x, u) = (out = c.k * (u.in - x.z) + x.z,)

"""Summing junction (§6.2): stateless, vector in, vector out."""
struct Junction end

input_types(::Junction, ::Type{T}) where {T} = (a = SVector{3,T}, b = SVector{3,T})
output_types(::Junction, ::Type{T}) where {T} = (sum = SVector{3,T},)
init_x(::Junction, ::Type{T}) where {T} = NamedTuple()

@inline h_xu(::Junction, x, u) = (sum = u.a + u.b,)

"""
The workhorse: the component type replicated `N` times to scale the model.
Mixed port shapes on both faces and a ~20-op body, in the size class of the
2026-07 compile-cost anchors (§9.7).
"""
struct Workhorse
    k::Float64
    τ::Float64
end

input_types(::Workhorse, ::Type{T}) where {T} =
    (cmd = T, wind = SVector{3,T}, pose_in = Pose3{T})
output_types(::Workhorse, ::Type{T}) where {T} =
    (pose = Pose3{T}, rate = SVector{3,T}, load = T)
init_x(::Workhorse, ::Type{T}) where {T} =
    (v = SVector{3,T}(1.0, 0.0, -0.5), ω = SVector{3,T}(0.0, 0.1, 0.0), e = T(0.0))

@inline function h_xu(c::Workhorse, x, u)
    vr = x.v - u.wind
    q = 0.5 * dot(vr, vr)
    ω = x.ω + c.k * cross(vr, u.pose_in.p)
    ψ = u.pose_in.ψ + c.τ * x.e
    p = u.pose_in.p + vr * c.τ + c.k * ω
    load = q * c.k - u.cmd * x.e
    (pose = Pose3(p, ψ), rate = ω, load = load)
end

# --- model spec ---------------------------------------------------------------
# Plain printable data, in the spirit of §9.2's `Build`: what the store layout
# and the schedule are compiled *from*. Components are listed in topological
# order; an input names its producer by component position and port name.

struct CompSpec
    path::Symbol
    comp::Any
    inputs::Vector{Pair{Symbol,Tuple{Int,Symbol}}}
end

struct ModelSpec
    comps::Vector{CompSpec}
end

Base.length(spec::ModelSpec) = length(spec.comps)

"""
    chain_spec(N)

Source → Filter, Junction → a chain of `N` `Workhorse`s, each taking the
previous one's `pose`. The chain is what forces a real topological order and
keeps the sweep from folding away; the fixed head keeps the model genuinely
heterogeneous at every `N`.
"""
function chain_spec(N::Int)
    comps = CompSpec[
        CompSpec(:src, Source(0.7), []),
        CompSpec(:fil, Filter(0.25), [:in => (1, :cmd)]),
        CompSpec(:jun, Junction(), [:a => (1, :wind_a), :b => (1, :wind_b)]),
    ]
    for i in 1:N
        pose_src = i == 1 ? (1, :pose) : (3 + i - 1, :pose)
        push!(comps, CompSpec(Symbol(:wrk, i), Workhorse(0.1 + 0.01i, 0.02),
            [:cmd => (2, :out), :wind => (3, :sum), :pose_in => pose_src]))
    end
    ModelSpec(comps)
end

# --- derived layout facts (build time) ---------------------------------------

out_types(cs::CompSpec, ::Type{T}) where {T} = output_types(cs.comp, T)
in_types(cs::CompSpec, ::Type{T}) where {T} = input_types(cs.comp, T)
x_value(cs::CompSpec, ::Type{T}) where {T} = init_x(cs.comp, T)

"""Every cell the model needs, in schedule order: `(component index, port, type)`."""
function cell_inventory(spec::ModelSpec, ::Type{T}) where {T}
    inv = Tuple{Int,Symbol,Type}[]
    for (ci, cs) in enumerate(spec.comps)
        for (port, P) in pairs(out_types(cs, T))
            push!(inv, (ci, port, P))
        end
    end
    inv
end

"""Per-component offsets into the shared flat state buffer, and the buffer."""
function state_layout(spec::ModelSpec, ::Type{T}) where {T}
    offsets = Int[]
    buf = T[]
    for cs in spec.comps
        push!(offsets, length(buf))
        append!(buf, collect(T, _leaf_values(x_value(cs, T))))
    end
    (offsets = offsets, buf = buf)
end

_leaf_values(v::Real) = (v,)
_leaf_values(v::StaticArray) = Iterators.flatten(map(_leaf_values, Tuple(v)))
_leaf_values(v::NamedTuple) = Iterators.flatten(map(_leaf_values, values(v)))
_leaf_values(v::Tuple) = Iterators.flatten(map(_leaf_values, v))
_leaf_values(v) = Iterators.flatten(map(_leaf_values, ntuple(i -> getfield(v, i), fieldcount(typeof(v)))))
