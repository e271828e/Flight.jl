# --- the continuous walking skeleton (increment 2) ------------------------------

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
    @test bundle_names(h_x, Plant(), CONTINUOUS, ()) === (:x, :t)
    @test bundle_names(h_xu, Plant(), CONTINUOUS, (:y,)) === (:x, :u, :y_x, :t)
    @test bundle_names(f, Plant(), CONTINUOUS, (:y,)) === (:x, :u, :y, :t)
    @test bundle_names(h_xu, Gain(1.0), CONTINUOUS, ()) === (:u, :t)
end

@testset "the schedule follows the feedthrough graph (§5.3)" begin
    sim = Simulation(feedback_model(); h = 1//1000)
    # sum first (both its inputs are loop-breaking), then ctl, then plant
    paths = [e.comp isa Sum ? :sum : e.comp isa Gain ? :ctl : :plant
             for e in walked(sim.bodies.sweep_2)]
    @test paths == [:sum, :ctl, :plant]
    @test length(walked(sim.bodies.sweep_1)) == 1
    @test length(walked(sim.bodies.rhs)) == 1
end

@testset "an algebraic loop is a build error (§5.5)" begin
    # `build` alone: rejection needs no deployment, which is the strata split.
    err = failure(() -> build(feedback_model(feedback_port = "power")))
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

    sim = Simulation(feedback_model(; k, ω, ζ); h = 1//1000)
    init!(sim, fragment(inputs = (ref = r,)))
    run!(sim; t_end = 2.0)

    # Tolerance, never `==` (D-163): RK4 truncation dominates at ~1e-12 here.
    @test state(sim, "children/plant").q ≈ exact(2.0) rtol = 1e-8
    @test port(sim, "children/plant", :y) ≈ exact(2.0)[1] rtol = 1e-8
    # the table is consistent at the boundary: `power` is a fresh decode
    @test port(sim, "children/plant", :power) ≈
          port(sim, "children/ctl", :out) * exact(2.0)[2] rtol = 1e-8
end

@testset "the phase-body roster is fixed and total (§9.7)" begin
    sim = Simulation(feedback_model(); h = 1//100)
    b = phase_bodies(sim)
    @test keys(b) === (:sweep_1, :sweep_2, :rhs, :ticks)
    @test b.ticks() === nothing            # empty body: legal, a no-op
    @test b.ticks(3) === nothing           # and total in both arities
    for name in keys(b)
        body = b[name]
        body(); body(0)
        @test @ballocated($body()) == 0
        @test @ballocated($body(1)) == 0
    end
end

@testset "the chunk walk is allocation-free at any width (§9.7)" begin
    # Six independent loops at chunk_size = 1: sweep_2 walks 18 one-entry
    # chunks — a chunk count no other fixture approaches — so the §7.5 canary
    # covers the outer walk over the chunk tuple itself, not just the entry
    # walks within one chunk.
    six = Group(NamedTuple{ntuple(i -> Symbol(:m, i), 6)}(ntuple(_ -> feedback_model(), 6)),
                (), ("ref" => ntuple(i -> "children/m$(i)/ref", 6),), ())
    sim = Simulation(six; h = 1//100, chunk_size = 1)
    @test length(sim.bodies.sweep_2.interior) > 16
    for name in keys(phase_bodies(sim))
        body = phase_bodies(sim)[name]
        body(); body(0)
        @test @ballocated($body()) == 0
        @test @ballocated($body(1)) == 0
    end
end

@testset "gate 1: stepping does not allocate (§7.5)" begin
    sim = Simulation(feedback_model(); h = 1//1000)
    init!(sim, fragment(inputs = (ref = 0.0,)))
    step!(sim, 1e-3)
    @test @ballocated(step!($sim, 1e-3)) == 0
    @test @ballocated(evaluate!($sim)) == 0
end

@testset "the whole continuous path is generic over the scalar (§7.2)" begin
    sim = Simulation(feedback_model(), D8; h = 1//1000)
    init!(sim, fragment(inputs = (ref = D8(0.7),)))
    run!(sim; t_end = 0.05)
    @test state(sim, "children/plant").q isa SVector{2,D8}
    @test ForwardDiff.value(port(sim, "children/plant", :y)) != 0.0
end

@testset "instances of one component type share one compiled body (D-162)" begin
    # Two independent loops, each a sub-assembly of one root: eight components,
    # still one entry type per stage per component type — the store's addressing
    # keeps offsets in fields. The root's one face fans out to both.
    two = Group((; a = feedback_model(), b = feedback_model(k = 3.0)), (),
                ("ref" => ("children/a/ref", "children/b/ref"),), ())
    sim = Simulation(two; h = 1//100)
    types(body) = unique(typeof(e) for e in walked(body))
    @test length(types(sim.bodies.sweep_1)) == 1     # two Plants, one h_x body
    @test length(types(sim.bodies.sweep_2)) == 3    # Plant, Gain, Sum
    @test length(types(sim.bodies.rhs)) == 1

    # The discrete tier keeps the property: a state store is a `Ref` whose
    # *type* every instance of a component type shares, so the store lives in a
    # field and two counters still compile to one `g` body.
    counters = Simulation(Group((; c1 = TickCounter(), c2 = TickCounter()), (), (), ());
                          h = 1//10)
    @test length(walked(counters.bodies.ticks)) == 2
    @test length(types(counters.bodies.ticks)) == 1
    @test length(types(counters.bodies.sweep_1)) == 1
    # And one bundle type per model, whatever the eltype count (D-162).
    @test counters.store isa StoreBundle
end

# Malformed components for the probe tests. Defined at top level, not inside the
# testset: a declaration written in a local scope binds a *new local function*
# of that name rather than adding a method to the global one (D-164), so `build`
# would dispatch on the untouched global and silently see the fallback
# declarations instead.
struct Undeclared <: AbstractComponent end
output_types(::Undeclared, ::Type{T}) where {T <: Real} = (a = T,)
h_x(::Undeclared, (; t)) = (a = 1.0, b = 2.0)

struct Unproduced <: AbstractComponent end
output_types(::Unproduced, ::Type{T}) where {T <: Real} = (a = T, b = T)
h_x(::Unproduced, (; t)) = (a = 1.0,)

struct BadDerivative <: AbstractComponent end
init_x(::BadDerivative) = (q = SVector(0.0, 0.0),)
f(::BadDerivative, (; x)) = (q = 0.0,)

struct NoFlow <: AbstractComponent end
init_x(::NoFlow) = (q = 1.0,)

struct Inert <: AbstractComponent end

@testset "the probe rejects malformed components (§9.3)" begin
    @test_throws BuildError build(single(Undeclared()))
    @test_throws BuildError build(single(Unproduced()))
    @test_throws BuildError build(single(BadDerivative()))
    @test_throws BuildError build(single(NoFlow()))
end
