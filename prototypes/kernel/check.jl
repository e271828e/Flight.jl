# Acceptance tests for increment 2: the continuous-tier walking skeleton.
#
#   julia --project=. check.jl

using Test, StaticArrays, LinearAlgebra, BenchmarkTools, ForwardDiff

include("src/leaves.jl")
include("src/declare.jl")
include("src/store.jl")
include("src/executor.jl")
include("src/build.jl")
include("src/sim.jl")
include("src/library.jl")

const D8 = ForwardDiff.Dual{Nothing,Float64,8}

@testset "the activation walk (§7.2)" begin
    @test retype(D8, Float64) === D8
    @test retype(D8, SVector{3,Float64}) === SVector{3,D8}
    @test retype(D8, Int) === Int
    @test retype(Float64, SVector{2,Float64}) === SVector{2,Float64}
    @test retype_value(D8, (q = SVector(1.0, 2.0),)).q isa SVector{2,D8}
end

@testset "the bundle law (§5.2)" begin
    # A name appears iff the store or fact exists: the stateless gain sees no
    # `x`, the no-feedthrough stage sees no `u`, `t` is always there.
    @test bundle_names(h_x, Plant(), ()) === (:x, :t)
    @test bundle_names(h_xu, Plant(), (:y,)) === (:x, :u, :y_x, :t)
    @test bundle_names(f, Plant(), (:y,)) === (:x, :u, :y, :t)
    @test bundle_names(h_xu, Gain(1.0), ()) === (:u, :t)
end

@testset "the schedule follows the feedthrough graph (§5.3)" begin
    sim = build(feedback_model())
    # sum first (both its inputs are loop-breaking), then ctl, then plant
    paths = [e.comp isa Sum ? :sum : e.comp isa Gain ? :ctl : :plant
             for c in sim.bodies.sweep_hxu.chunks for e in c.entries]
    @test paths == [:sum, :ctl, :plant]
    @test length([e for c in sim.bodies.sweep_hx.chunks for e in c.entries]) == 1
    @test length([e for c in sim.bodies.rhs.chunks for e in c.entries]) == 1
end

@testset "an algebraic loop is a build error (§5.5)" begin
    err = try
        build(feedback_model(feedback_port = :power))
        nothing
    catch e
        e
    end
    @test err isa BuildError
    @test occursin("algebraic loop", err.msg)
    @test occursin("plant", err.msg) && occursin("ctl", err.msg) && occursin("sum", err.msg)
end

@testset "the loop integrates the right trajectory" begin
    # Closed-loop reference: ẋ = (A - B k C) x + B k r, integrated exactly.
    ω, ζ, k, r = 2.0, 0.1, 4.0, 0.7
    A = SMatrix{2,2}(0.0, -ω^2, 1.0, -2ζ * ω)
    B = SVector(0.0, 1.0)
    Acl = A - B * k * SVector(1.0, 0.0)'
    exact(t) = exp(Acl * t) * (Acl \ (B * k * r)) - Acl \ (B * k * r)

    sim = build(feedback_model(; k, ω, ζ))
    set_slot!(sim, :sum, :a, r)
    init!(sim)
    run!(sim, 2.0, 1e-3)

    # Tolerance, never `==` (D-163): RK4 truncation dominates at ~1e-12 here.
    @test state(sim, :plant).q ≈ exact(2.0) rtol = 1e-8
    @test port(sim, :plant, :y) ≈ exact(2.0)[1] rtol = 1e-8
    # the table is consistent at the boundary: `power` is a fresh decode
    @test port(sim, :plant, :power) ≈ port(sim, :ctl, :out) * exact(2.0)[2] rtol = 1e-8
end

@testset "the phase-body roster is fixed and total (§9.7)" begin
    sim = build(feedback_model())
    b = phase_bodies(sim)
    @test keys(b) === (:sweep_hx, :sweep_hxu, :rhs, :ticks)
    @test b.ticks() === nothing            # empty body: legal, a no-op
    @test b.ticks(3) === nothing           # and total in both arities
    for name in keys(b)
        body = b[name]
        body(); body(0)
        @test @ballocated($body()) == 0
        @test @ballocated($body(1)) == 0
    end
end

@testset "gate 1: stepping does not allocate (§7.5)" begin
    sim = build(feedback_model())
    init!(sim)
    step!(sim, 1e-3)
    @test @ballocated(step!($sim, 1e-3)) == 0
    @test @ballocated(evaluate!($sim)) == 0
end

@testset "the whole continuous path is generic over the scalar (§7.2)" begin
    sim = build(feedback_model(), D8)
    set_slot!(sim, :sum, :a, D8(0.7))
    init!(sim)
    run!(sim, D8(0.05), D8(1e-3))
    @test state(sim, :plant).q isa SVector{2,D8}
    @test ForwardDiff.value(port(sim, :plant, :y)) != 0.0
end

@testset "instances of one component type share one compiled body (D-162)" begin
    # Two independent loops: eight components, still one entry type per stage
    # per component type — the store's addressing keeps offsets in fields.
    two = ModelSpec(vcat(feedback_model().comps,
                         [ComponentSpec(:plant2, Plant(), [:u => (:ctl2, :out)]),
                          ComponentSpec(:ctl2, Gain(3.0), [:e => (:sum2, :e)]),
                          ComponentSpec(:sum2, Sum(), [:b => (:plant2, :y)])]))
    sim = build(two)
    types(body) = unique(typeof(e) for c in body.chunks for e in c.entries)
    @test length(types(sim.bodies.sweep_hx)) == 1     # two Plants, one h_x body
    @test length(types(sim.bodies.sweep_hxu)) == 3    # Plant, Gain, Sum
    @test length(types(sim.bodies.rhs)) == 1
end

# Malformed components for the probe tests. Defined at top level, not inside the
# testset: a declaration written in a local scope binds a *new local function*
# of that name rather than adding a method to the global one (D-164), so `build`
# would dispatch on the untouched global and silently see the fallback
# declarations instead.
struct Undeclared end
output_types(::Undeclared, ::Type{T}) where {T <: Real} = (a = T,)
h_x(::Undeclared, (; t)) = (a = 1.0, b = 2.0)

struct Unproduced end
output_types(::Unproduced, ::Type{T}) where {T <: Real} = (a = T, b = T)
h_x(::Unproduced, (; t)) = (a = 1.0,)

struct BadDerivative end
init_x(::BadDerivative) = (q = SVector(0.0, 0.0),)
f(::BadDerivative, (; x)) = (q = 0.0,)

struct NoFlow end
init_x(::NoFlow) = (q = 1.0,)

# The whole-signature forgotten-`T`: the retired one-argument form, which the
# `::Any` fallback would otherwise swallow into a component declaring nothing.
struct OneArgDecl end
output_types(::OneArgDecl) = (a = Float64,)
h_x(::OneArgDecl, (; t)) = (a = 1.0,)

@testset "the probe rejects malformed components (§9.3)" begin
    one(c) = ModelSpec([ComponentSpec(:c, c)])
    @test_throws BuildError build(one(Undeclared()))
    @test_throws BuildError build(one(Unproduced()))
    @test_throws BuildError build(one(BadDerivative()))
    @test_throws BuildError build(one(NoFlow()))
    @test_throws BuildError build(one(OneArgDecl()))
end

# The constant-branch idiom (D-166): a literal `Float64` returned into a
# declared-`T` port is a lawful arrival, embedded as a zero-partial.
struct ConstantBranch end
input_types(::ConstantBranch, ::Type{T}) where {T <: Real} = (in = T,)
output_types(::ConstantBranch, ::Type{T}) where {T <: Real} = (out = T, vec = SVector{2,T})
h_xu(::ConstantBranch, (; u)) = (out = u.in > 0 ? u.in : 0.0, vec = SVector(0.0, 1.0))

# A `Dual` arriving at a deliberately pinned leaf: the one honest cause, and the
# one that earns the didactic hint.
struct PinnedGetsDual end
output_types(::PinnedGetsDual, ::Type{T}) where {T <: Real} = (frozen = Float64,)
h_x(::PinnedGetsDual, (; t)) = (frozen = t,)

@testset "embed-accept keeps the constant branch legal (D-166)" begin
    one(c) = ModelSpec([ComponentSpec(:c, c)])
    # Both ports return literal `Float64`s at a `Dual` activation — the scalar
    # through a branch not taken, the `SVector` wholesale.
    sim = build(one(ConstantBranch()), D8)
    init!(sim)
    # What the table holds is the cell's type, the constant embedded into it.
    @test port(sim, :c, :out) isa D8
    @test port(sim, :c, :vec) isa SVector{2,D8}
    @test ForwardDiff.value(port(sim, :c, :vec)[2]) == 1.0

    # The converse is not accepted: a `Dual` at a pinned leaf is an error, with
    # the hint that names the one honest cause.
    err = try
        build(one(PinnedGetsDual()), D8)
        nothing
    catch e
        e
    end
    @test err isa BuildError
    @test occursin("participates in differentiation", err.msg)
end

# A deliberately pinned leaf (D-166): `frozen` is declared `Float64` rather than
# `T`, so it must not follow the activation scalar.
struct PinnedLeaf end
output_types(::PinnedLeaf, ::Type{T}) where {T <: Real} = (a = T, frozen = Float64)
h_x(::PinnedLeaf, (; t)) = (a = t, frozen = 2.0)

@testset "a pinned leaf is refused, not silently promoted (D-166)" begin
    one(c) = ModelSpec([ComponentSpec(:c, c)])
    # Nominally the pin and the activation scalar coincide, so it builds.
    @test build(one(PinnedLeaf())) isa Simulation
    # Off nominal it needs its own store, which this increment does not build.
    # The point of the check is that this is an error rather than a `Float64`
    # quietly flattened into the `Dual` buffer as a zero-partial.
    err = try
        build(one(PinnedLeaf()), D8)
        nothing
    catch e
        e
    end
    @test err isa BuildError
    @test occursin("frozen", err.msg) && occursin("own store", err.msg)
end
