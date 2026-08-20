# Acceptance tests for the kernel prototype: the continuous-tier walking
# skeleton (increment 2), the discrete tier at one rate (increment 3), the
# hierarchy (increment 4) and the multi-rate grid (increment 5).
#
#   julia --project=. check.jl

using Test, StaticArrays, LinearAlgebra, BenchmarkTools, ForwardDiff

include("src/leaves.jl")
include("src/declare.jl")
include("src/assembly.jl")
include("src/store.jl")
include("src/executor.jl")
include("src/build.jl")
include("src/sim.jl")
include("src/library.jl")

const D8 = ForwardDiff.Dual{Nothing,Float64,8}

# The entries a phase body walks, per variant (§10.5): the boundary variant is
# the full list — its discrete entries wearing their `Gated` wrapper, unwrapped
# here — the interior one the continuous entries alone.
walked(body, variant = :boundary) =
    [e isa Gated ? e.e : e for c in getfield(body, variant) for e in c.entries]

# The gated entries of a body's boundary variant.
gated(body) = count(e isa Gated for c in body.boundary for e in c.entries)

# A model of one component, wrapped in the assembly every model's root must be.
single(c) = Group((; c = c), (), (), ())

# The error a build raises, for the tests that read the message.
failure(f) =
    try
        f()
        nothing
    catch e
        e
    end

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
             for e in walked(sim.bodies.sweep_hxu)]
    @test paths == [:sum, :ctl, :plant]
    @test length(walked(sim.bodies.sweep_hx)) == 1
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
    set_slot!(sim, "ref", r)
    init!(sim)
    run!(sim, 2.0)

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
    sim = Simulation(feedback_model(); h = 1//1000)
    init!(sim)
    step!(sim, 1e-3)
    @test @ballocated(step!($sim, 1e-3)) == 0
    @test @ballocated(evaluate!($sim)) == 0
end

@testset "the whole continuous path is generic over the scalar (§7.2)" begin
    sim = Simulation(feedback_model(), D8; h = 1//1000)
    set_slot!(sim, "ref", D8(0.7))
    init!(sim)
    run!(sim, 0.05)
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
    @test length(types(sim.bodies.sweep_hx)) == 1     # two Plants, one h_x body
    @test length(types(sim.bodies.sweep_hxu)) == 3    # Plant, Gain, Sum
    @test length(types(sim.bodies.rhs)) == 1

    # The discrete tier keeps the property: a state store is a `Ref` whose
    # *type* every instance of a component type shares, so the store lives in a
    # field and two counters still compile to one `g` body.
    counters = Simulation(Group((; c1 = TickCounter(), c2 = TickCounter()), (), (), ());
                          h = 1//10)
    @test length(walked(counters.bodies.ticks)) == 2
    @test length(types(counters.bodies.ticks)) == 1
    @test length(types(counters.bodies.sweep_hx)) == 1
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

# --- class (§8.5) -------------------------------------------------------------
# Class is read off *which* well-known declarations a type defines:
# `child_connections` the assembly marker, any leaf declaration a primitive's.

struct HoldsComponents <: AbstractComponent      # components, but no class to read
    inner::Gain
end

struct BothFamilies <: AbstractComponent         # assembly marker beside a contract
    inner::Gain
end
child_connections(::BothFamilies) = ()
output_types(::BothFamilies, ::Type{T}) where {T <: Real} = (a = T,)
h_x(::BothFamilies, (; t)) = (a = 1.0,)

@testset "class is read off the declaration shape (§8.5)" begin
    @test classify("c", Gain(1.0)) === PRIMITIVE
    @test classify("c", feedback_model()) === ASSEMBLY
    @test classify("c", Vehicle()) === ASSEMBLY

    # A component that declares nothing and defines no stage cannot be
    # intentional (D-164) — and now says so as a missing class, naming both
    # families rather than failing later and elsewhere.
    err = failure(() -> classify("c", Inert()))
    @test err isa BuildError
    @test occursin("child_connections", err.msg) && occursin("output_types", err.msg)

    # Sharpened when the type holds components: the likely omission, named.
    err = failure(() -> classify("c", HoldsComponents(Gain(1.0))))
    @test occursin("holds components but declares no `child_connections`", err.msg)

    # Both families on one type: an assembly owns no state and no contract of
    # its own, so this is a build error too.
    err = failure(() -> classify("c", BothFamilies(Gain(1.0))))
    @test err isa BuildError && occursin("output_types", err.msg)

    # A primitive at the root is a Stratum A error (§9.1): a model's root is an
    # assembly, whose input faces are the root slots.
    err = failure(() -> build(Plant()))
    @test err isa BuildError && occursin("root", err.msg)
end

# --- container children (§8.5) ------------------------------------------------

struct MixedContainer <: AbstractComponent
    kids::NamedTuple
end
child_connections(::MixedContainer) = ()

struct EmptyRoster <: AbstractComponent          # parametric code needs no special case
    roster::NamedTuple
    src::ModedSource
end
child_connections(::EmptyRoster) = ()

struct TupleRoster{U <: Tuple} <: AbstractComponent
    units::U
end
child_connections(::TupleRoster)  = ("units/1/out" => "units/2/e",)
input_connections(::TupleRoster)  = ("in" => "units/1/e",)
output_connections(::TupleRoster) = ("units/2/out" => "y",)

@testset "container children are path-named `field/key` and `field/1` (§8.5)" begin
    sim = Simulation(Group((; c1 = TickCounter(), c2 = TickCounter()), (), (), ());
                     h = 1//10)
    @test sim.flat.paths == ["children/c1", "children/c2"]
    @test state(sim, "children/c2") === (n = 0,)

    # An empty container contributes zero children, and is not an error.
    b = build(EmptyRoster((;), ModedSource()))
    @test b.flat.paths == ["src"]

    # A container mixing components with anything else is one, by name.
    err = failure(() -> build(single(MixedContainer((a = Gain(1.0), b = 2.0)))))
    @test err isa BuildError
    @test occursin("mixes components", err.msg) && occursin("kids", err.msg)

    # The `Tuple` form: the same rule with index segments, `"field/1"…"field/N"`
    # (§8.5), addressable by the parent's declarations like any child name.
    tsim = Simulation(TupleRoster((Gain(2.0), Gain(3.0))); h = 1//10)
    @test tsim.flat.paths == ["units/1", "units/2"]
    set_slot!(tsim, "in", 1.0)
    init!(tsim)
    @test port(tsim, "units/2", :out) === 6.0
    @test port(tsim, "", :y) === port(tsim, "units/2", :out)

    # The mixing rule is form-blind.
    err = failure(() -> build(TupleRoster((Gain(1.0), 2.0))))
    @test err isa BuildError
    @test occursin("mixes components", err.msg) && occursin("units", err.msg)
end

# --- paths and the reach rule (§6.1, §8.6) ------------------------------------
# The same instance, held two ways: the deep route is legal through concretely
# declared fields and a build error past a generically held child, which is the
# rule being about the declaration's knowledge rather than the build's.

struct ConcreteHold <: AbstractComponent
    inner::SampledLoop
end
child_connections(::ConcreteHold) = ()
input_connections(::ConcreteHold) = ("ref" => "inner/sum/a",)
output_connections(::ConcreteHold) = ("inner/plant/y" => "y",)

struct GenericHold{L <: AbstractComponent} <: AbstractComponent
    inner::L
end
child_connections(::GenericHold) = ()
input_connections(::GenericHold) = ("ref" => "inner/sum/a",)
output_connections(::GenericHold) = ("inner/plant/y" => "y",)

@testset "a path stops at the first generically held child (§6.1)" begin
    # Two levels down into declared structure, bypassing the sub-assembly's own
    # face: legal, and the routed input is fed exactly once.
    sim = Simulation(ConcreteHold(SampledLoop()); h = 1//50)
    set_slot!(sim, "ref", 1.0)
    init!(sim)
    @test port(sim, "", :y) === port(sim, "inner/plant", :y)

    # The identical paths against the identical instance, held generically: the
    # concrete child would resolve them, and the declaration may not.
    err = failure(() -> build(GenericHold(SampledLoop())))
    @test err isa BuildError
    @test occursin("generically", err.msg) && occursin("inner", err.msg)
end

# The rule's other half (§13.3): resolving *to* a generic child is face-level
# access. `mid` is concretely declared and `inner` generically, so a route from
# the outer type may end at `inner`'s faces and may not go past them.

struct MidHold{L <: AbstractComponent} <: AbstractComponent
    inner::L
end
child_connections(::MidHold) = ()
input_connections(::MidHold) = ("ref" => "inner/ref",)
output_connections(::MidHold) = ("inner/y" => "y",)

struct FaceReach <: AbstractComponent            # ends at the generic child's face
    mid::MidHold{SampledLoop}
end
child_connections(::FaceReach) = ()
input_connections(::FaceReach) = ("ref" => "mid/inner/ref",)
output_connections(::FaceReach) = ("mid/inner/y" => "y",)

struct PastReach <: AbstractComponent            # one segment further: past it
    mid::MidHold{SampledLoop}
end
child_connections(::PastReach) = ()
input_connections(::PastReach) = ("ref" => "mid/inner/sum/a",)
output_connections(::PastReach) = ("mid/inner/plant/y" => "y",)

@testset "a route may end at a generic child's face, never go past it (§6.1, §13.3)" begin
    sim = Simulation(FaceReach(MidHold(SampledLoop())); h = 1//50)
    set_slot!(sim, "ref", 1.0)
    init!(sim)
    @test port(sim, "", :y) === port(sim, "mid/inner/plant", :y)

    err = failure(() -> build(PastReach(MidHold(SampledLoop()))))
    @test err isa BuildError
    @test occursin("generically", err.msg) && occursin("mid/inner", err.msg)
end

# --- the three connection declarations (§8.6) ---------------------------------

struct BackwardsWire <: AbstractComponent        # a consumer endpoint on a port
    a::ModedSource
    b::Gain
end
child_connections(::BackwardsWire) = ("a/out" => "b/out",)

struct BackwardsFace <: AbstractComponent        # a producer endpoint on a face
    a::ModedSource
    b::Gain
end
child_connections(::BackwardsFace) = ("a/out" => "b/e",)
output_connections(::BackwardsFace) = ("b/e" => "y",)

struct SlashedFace <: AbstractComponent
    a::ModedSource
end
child_connections(::SlashedFace) = ()
output_connections(::SlashedFace) = ("a/out" => "sensors/out",)

struct CollidingFaces <: AbstractComponent
    a::ModedSource
    b::Gain
end
child_connections(::CollidingFaces) = ("a/out" => "b/e",)
input_connections(::CollidingFaces) = ("y" => "b/e",)
output_connections(::CollidingFaces) = ("b/out" => "y",)

@testset "direction is declared by the method, endpoints cross-check it (§8.6)" begin
    err = failure(() -> build(BackwardsWire(ModedSource(), Gain(1.0))))
    @test err isa BuildError
    @test occursin("child_connections", err.msg) && occursin("consumer", err.msg)

    err = failure(() -> build(BackwardsFace(ModedSource(), Gain(1.0))))
    @test err isa BuildError
    @test occursin("output_connections", err.msg) && occursin("producer", err.msg)
end

@testset "face names carry two invariants (§8.6)" begin
    err = failure(() -> build(SlashedFace(ModedSource())))
    @test err isa BuildError && occursin("reserved for structural paths", err.msg)

    err = failure(() -> build(CollidingFaces(ModedSource(), Gain(1.0))))
    @test err isa BuildError && occursin("unique", err.msg)
end

# --- the whole-tree obligation model (§6.1) -----------------------------------

struct Starved <: AbstractComponent              # an obligation chain that never ends
    g::Gain
end
child_connections(::Starved) = ()

struct DoubleFed <: AbstractComponent            # handed up *and* wired
    loop::SampledLoop
    src::ModedSource
end
child_connections(::DoubleFed) = ("src/out" => "loop/ref", "src/out" => "loop/sum/a")

struct DoubleFedSibling <: AbstractComponent     # an ancestor's route onto a wired input
    loop::SampledLoop
    src::ModedSource
end
child_connections(::DoubleFedSibling) = ("src/out" => "loop/sum/b",)

@testset "every input is fed exactly once, across levels (§6.1)" begin
    err = failure(() -> build(Starved(Gain(1.0))))
    @test err isa BuildError
    @test occursin("fed by nothing", err.msg) && occursin("g", err.msg)

    err = failure(() -> build(DoubleFed(SampledLoop(), ModedSource())))
    @test err isa BuildError
    @test occursin("fed twice", err.msg) && occursin("loop/sum", err.msg)

    # The same rule one level down: the sub-assembly's own wire against the
    # ancestor's deep route, and the message names both entries.
    err = failure(() -> build(DoubleFedSibling(SampledLoop(), ModedSource())))
    @test err isa BuildError
    @test occursin("child_connections at `loop`", err.msg) &&
          occursin("child_connections at the root component", err.msg)

    # The one legitimate terminus: the root's own input faces are the slots,
    # synthesized by `probe_value` and written by face name (§11.3).
    sim = Simulation(Group((; c = Gain(2.0)), (), ("in" => "children/c/e",), ());
                     h = 1//100)
    @test sim.flat.slots == [:in]
    init!(sim)
    @test port(sim, "children/c", :out) == 0.0        # the synthesized value
    set_slot!(sim, "in", 3.0)
    init!(sim)
    @test port(sim, "children/c", :out) == 6.0
end

# --- tier classification (§8.2) -----------------------------------------------
# Tier is read off the declaration shape. These components are the shapes the
# classifier has to separate, plus the four ways a declaration set can disagree.

struct DiscreteCounter <: AbstractComponent   # stateful discrete: `g` decides
end
init_x(::DiscreteCounter) = (n = 0,)
output_types(::DiscreteCounter) = (n = Int,)
h_x(::DiscreteCounter, (; x)) = (n = x.n,)
g(::DiscreteCounter, (; x)) = (n = x.n + 1,)

struct DiscreteMap <: AbstractComponent       # stateless discrete: the arity decides
end
input_types(::DiscreteMap) = (a = Int,)
output_types(::DiscreteMap) = (b = Int,)
h_xu(::DiscreteMap, (; u)) = (b = 2u.a,)

struct BothUpdates <: AbstractComponent       # `f` and `g` on one component
end
init_x(::BothUpdates) = (q = 1.0,)
output_types(::BothUpdates, ::Type{T}) where {T <: Real} = (a = T,)
h_x(::BothUpdates, (; x)) = (a = x.q,)
f(::BothUpdates, (; x)) = (q = 0.0,)
g(::BothUpdates, (; x)) = (q = x.q,)

struct WrongArity <: AbstractComponent        # `g` beside a two-argument contract
end
init_x(::WrongArity) = (n = 0,)
output_types(::WrongArity, ::Type{T}) where {T <: Real} = (a = T,)
h_x(::WrongArity, (; x)) = (a = 1.0,)
g(::WrongArity, (; x)) = (n = x.n,)

struct ModesOnDiscrete <: AbstractComponent   # `init_m` is continuous-only
end
init_x(::ModesOnDiscrete) = (n = 0,)
init_m(::ModesOnDiscrete) = (phase = :idle,)
output_types(::ModesOnDiscrete) = (a = Int,)
h_x(::ModesOnDiscrete, (; x)) = (a = x.n,)
g(::ModesOnDiscrete, (; x)) = (n = x.n,)

struct BothArities <: AbstractComponent       # a member of both contract families
end
output_types(::BothArities, ::Type{T}) where {T <: Real} = (a = T,)
output_types(::BothArities) = (a = Float64,)
h_x(::BothArities, (; t)) = (a = 1.0,)

@testset "tier is read off the declaration shape (§8.2)" begin
    # The two deciders: the update law for a stateful leaf, the contract arity
    # for a stateless one.
    @test classify_tier("c", Plant()) === CONTINUOUS
    @test classify_tier("c", Gain(1.0)) === CONTINUOUS
    @test classify_tier("c", DiscreteCounter()) === DISCRETE
    @test classify_tier("c", DiscreteMap()) === DISCRETE

    # The discrete bundle sets: `Δt` is a discrete-tier fact, `m` a continuous
    # one, both under the shared state letter (D-173).
    @test bundle_names(h_x, DiscreteCounter(), DISCRETE, ()) === (:x, :t, :Δt)
    @test bundle_names(g, DiscreteCounter(), DISCRETE, (:n,)) === (:x, :y, :t, :Δt)
    @test bundle_names(h_xu, DiscreteMap(), DISCRETE, ()) === (:u, :t, :Δt)

    # Disagreement names the offending declaration and the tier the rest announce.
    for (c, offender) in ((BothUpdates(), "g"), (WrongArity(), "output_types"),
                          (ModesOnDiscrete(), "init_m"), (BothArities(), "output_types"))
        err = failure(() -> classify_tier("c", c))
        @test err isa BuildError
        @test occursin(offender, err.msg)
    end

    # A store with no update law is §8.2's sibling of the classless component.
    @test_throws BuildError classify_tier("c", NoFlow())

    # The base tick period is deployment's, not the build's: the same `Build`
    # deploys at any admissible grid, and the executor cannot exist before one
    # binds because `Δt`, `D` and `Φ` are entry-field data (§9.1, §9.7).
    b = build(single(DiscreteCounter()))
    @test b isa Build
    err = failure(() -> Simulation(b))
    @test err isa BuildError && occursin("h", err.msg)
    @test Simulation(b; h = 1//10) isa Simulation
end

# --- the discrete tier at one rate (§10.5) ------------------------------------

@testset "the sampled loop matches the exact ZOH discretization" begin
    # The reference is the hybrid system solved exactly: over each interval the
    # plant is linear with a *constant* input, so its transition is
    # `q[k+1] = Ad q[k] + Bd s[k]` with `Ad = exp(A Δt)`, and the controller
    # advances by its own recursion. Nothing here matches unless the hold is
    # real — a mid-step re-run of `ctl` would move `u` inside the interval and
    # break the closed form.
    kI, ω, ζ, Δt, r, N = 3.0, 2.0, 0.1, 0.02, 0.7, 100
    A = SMatrix{2,2}(0.0, -ω^2, 1.0, -2ζ * ω)
    B = SVector(0.0, 1.0)
    Ad = exp(A * Δt)
    Bd = A \ ((Ad - I) * B)

    q, s = SVector(0.0, 0.0), 0.0
    for _ in 1:N
        q, s = Ad * q + Bd * s, s + kI * Δt * (r - q[1])
    end
    s_next = s + kI * Δt * (r - q[1])   # the update boundary N itself performs

    sim = Simulation(sampled_loop(; kI, ω, ζ); h = 1//50)   # h = Δt: one rate, n = 1
    set_slot!(sim, "ref", r)
    init!(sim)
    run!(sim, N * Δt)

    # The tolerance is RK4's, not the semantics': the reference integrates
    # exactly where the loop takes one RK4 step per sample, which at `ωΔt` this
    # size leaves ~1e-8 relative. A violated hold moves `u` mid-interval and
    # costs ~1e-3, so the margin still separates the two by orders of magnitude.
    @test state(sim, "children/plant").q ≈ q rtol = 1e-6
    # The cell holds `s[N]` — the output stage ran at this boundary — while the
    # store already holds `s[N+1]`. That gap *is* the sampled-data ordering:
    # output stages before updates, within one boundary (§10.5).
    @test port(sim, "children/ctl", :u) ≈ s rtol = 1e-6
    @test state(sim, "children/ctl").s ≈ s_next rtol = 1e-6
end

@testset "the ZOH holds by compile-time absence (§10.5)" begin
    sim = Simulation(sampled_loop(); h = 1//50)
    set_slot!(sim, "ref", 1.0)         # excite the loop, or nothing moves at all
    init!(sim)

    # Structural: the interior variants carry continuous entries only, so there
    # is no gating test on the hot path — the hold is not implemented, it is
    # the absence of any way to change the cell.
    @test length(walked(sim.bodies.sweep_hx, :interior)) == 1        # plant only
    @test length(walked(sim.bodies.sweep_hx)) == 2                   # plus ctl
    @test isempty(walked(sim.bodies.ticks, :interior))
    @test length(walked(sim.bodies.ticks)) == 1

    # Semantic: a step is made of interior evaluations, so the discrete cell
    # cannot move across one, while the continuous table does. Run a few
    # boundaries first — the loop starts at rest with `u` held at zero, so the
    # very first interval moves nothing at all.
    run!(sim, 0.1)
    u₀, y₀ = port(sim, "children/ctl", :u), port(sim, "children/plant", :y)
    @test u₀ != 0.0
    step!(sim, 0.02)
    @test port(sim, "children/ctl", :u) == u₀   # untouched: never gathered, never written
    @test port(sim, "children/plant", :y) != y₀
end

@testset "discrete state and modes live outside the buffer (§7.3)" begin
    sim = Simulation(Group((; counter = TickCounter(), moded = ModedSource()), (), (), ());
                     h = 1//10)
    # The flat buffer is continuous state only; the counter's `Int` is in its
    # own store, and no store mirrors another.
    @test isempty(sim.xbuf)
    @test state(sim, "children/counter") === (n = 0,)
    @test modes(sim, "children/moded") === (phase = :idle,)

    # `Int` and `Bool` cells force their own buffers — the plural in
    # "per-eltype stores", first exercised here.
    @test Set(keys(sim.store.stores)) == Set([Symbol(Int), Symbol(Bool), Symbol(Float64)])

    init!(sim)
    @test state(sim, "children/counter") === (n = 1,)   # boundary zero is a tick
    @test port(sim, "children/counter", :n) === 0       # the cell holds what it published
    @test port(sim, "children/counter", :even) === true
    run!(sim, 0.5)
    # Six boundaries: zero, then one per step. The store leads the cell by one,
    # the update having run after the output stage at each of them.
    @test state(sim, "children/counter") === (n = 6,)
    @test port(sim, "children/counter", :n) === 5
end

@testset "the workspace is scratch, on both tiers (§7.3)" begin
    sim = Simulation(Group((; sm = Smoother(0.5), src = ModedSource(), wg = WorkGain(2.0)),
                           ("children/src/out" => "children/sm/a",
                            "children/src/out" => "children/sm/b",
                            "children/src/out" => "children/wg/in"), (), ());
                     h = 1//10)
    init!(sim)
    @test port(sim, "children/wg", :out) == 0.0      # 2 × the idle-phase constant
    run!(sim, 0.3)
    @test state(sim, "children/sm").v isa SVector{2,Float64}

    # Allocation is what the idiom is for: in-place math on scratch, an isbits
    # snapshot into the store, and nothing on the measured path.
    b = phase_bodies(sim)
    for name in keys(b)
        body = b[name]
        body(); body(1)
        @test @ballocated($body()) == 0
        @test @ballocated($body(1)) == 0
    end
end

# A `g` that widens its own store: the discrete world is pinned, so this is an
# error rather than a silent conversion at the store assignment.
struct WidenedUpdate <: AbstractComponent end
init_x(::WidenedUpdate) = (n = 0,)
output_types(::WidenedUpdate) = (n = Int,)
h_x(::WidenedUpdate, (; x)) = (n = x.n,)
g(::WidenedUpdate, (; x)) = (n = x.n + 0.5,)

@testset "a discrete successor is the store's own type (§7.3)" begin
    @test_throws BuildError build(single(WidenedUpdate()))
end

@testset "the discrete tier is frozen at a non-nominal activation (§7.2)" begin
    # The plant's input is declared `T` and wired to a discrete `Float64` cell:
    # a lawful arrival, embedded as a zero-partial. `ctl`'s stages do not run
    # at all, and its cell holds what the probe wrote — which is exactly what a
    # tick at `t₀⁻` would have produced.
    simd = Simulation(sampled_loop(), D8; h = 1//50)
    @test isempty(walked(simd.bodies.ticks))
    @test length(walked(simd.bodies.sweep_hx)) == 1          # plant only; ctl frozen
    @test port(simd, "children/ctl", :u) isa Float64

    init!(simd)
    run!(simd, 0.04)
    @test state(simd, "children/plant").q isa SVector{2,D8}
    @test port(simd, "children/ctl", :u) == 0.0              # held, never recomputed
end

# --- hierarchy end to end (§8.5, §8.6) ----------------------------------------

@testset "a two-level assembly runs the sampled loop through its faces" begin
    kI, ω, ζ, Δt, r, k, N = 3.0, 2.0, 0.1, 0.02, 0.7, 2.0, 50
    A = SMatrix{2,2}(0.0, -ω^2, 1.0, -2ζ * ω)
    B = SVector(0.0, 1.0)
    Ad = exp(A * Δt)
    Bd = A \ ((Ad - I) * B)

    # The vehicle's gain scales the reference before the loop sees it.
    q, s = SVector(0.0, 0.0), 0.0
    for _ in 1:N
        q, s = Ad * q + Bd * s, s + kI * Δt * (k * r - q[1])
    end

    sim = Simulation(Vehicle(; k, kI, ω, ζ); h = 1//50)
    @test sim.flat.paths == ["loop/plant", "loop/ctl", "loop/sum", "trim"]
    set_slot!(sim, "ref", r)
    init!(sim)
    run!(sim, N * Δt)
    @test state(sim, "loop/plant").q ≈ q rtol = 1e-6
    @test port(sim, "loop", :cmd) ≈ s rtol = 1e-6
end

@testset "a face's type and tier are its internal endpoint's (§8.6)" begin
    sim = Simulation(Vehicle(); h = 1//50)
    set_slot!(sim, "ref", 1.0)
    init!(sim)
    run!(sim, 0.1)
    # A face is its endpoint: no cell of its own, at any level of re-export.
    @test port(sim, "loop", :y) === port(sim, "loop/plant", :y)
    @test port(sim, "", :y) === port(sim, "loop/plant", :y)
    @test port(sim, "", :cmd) === port(sim, "loop/ctl", :u)
    # And a deep-routed face reaches a grandchild's port directly.
    @test port(sim, "", :power) === port(sim, "loop/plant", :power)

    # Tier-neutral, and the tiers are *derived*: at a non-nominal activation the
    # continuous-sourced face walks while the discrete-sourced one stays pinned.
    simd = Simulation(Vehicle(), D8; h = 1//50)
    @test port(simd, "", :y) isa D8
    @test port(simd, "", :cmd) isa Float64
end

# The constant-branch idiom (D-166): a literal `Float64` returned into a
# declared-`T` port is a lawful arrival, embedded as a zero-partial.
struct ConstantBranch <: AbstractComponent end
input_types(::ConstantBranch, ::Type{T}) where {T <: Real} = (in = T,)
output_types(::ConstantBranch, ::Type{T}) where {T <: Real} = (out = T, vec = SVector{2,T})
h_xu(::ConstantBranch, (; u)) = (out = u.in > 0 ? u.in : 0.0, vec = SVector(0.0, 1.0))

# A `Dual` arriving at a deliberately pinned leaf: the one honest cause, and the
# one that earns the didactic hint.
struct PinnedGetsDual <: AbstractComponent end
output_types(::PinnedGetsDual, ::Type{T}) where {T <: Real} = (frozen = Float64,)
h_x(::PinnedGetsDual, (; t)) = (frozen = t,)

@testset "embed-accept keeps the constant branch legal (D-166)" begin
    # Both ports return literal `Float64`s at a `Dual` activation — the scalar
    # through a branch not taken, the `SVector` wholesale.
    sim = Simulation(Group((; c = ConstantBranch()), (), ("in" => "children/c/in",), ()),
                     D8; h = 1//100)
    init!(sim)
    # What the table holds is the cell's type, the constant embedded into it.
    @test port(sim, "children/c", :out) isa D8
    @test port(sim, "children/c", :vec) isa SVector{2,D8}
    @test ForwardDiff.value(port(sim, "children/c", :vec)[2]) == 1.0

    # The converse is not accepted: a `Dual` at a pinned leaf is an error, with
    # the hint that names the one honest cause.
    err = failure(() -> build(single(PinnedGetsDual()), D8))
    @test err isa BuildError
    @test occursin("participates in differentiation", err.msg)
end

# A deliberately pinned leaf (D-166): `frozen` is declared `Float64` rather than
# `T`, so it must not follow the activation scalar.
struct PinnedLeaf <: AbstractComponent end
output_types(::PinnedLeaf, ::Type{T}) where {T <: Real} = (a = T, frozen = Float64)
h_x(::PinnedLeaf, (; t)) = (a = t, frozen = 2.0)

@testset "a pinned leaf lives in its own store (D-166, D-162)" begin
    # Nominally the pin and the activation scalar coincide: one buffer.
    sim = Simulation(single(PinnedLeaf()); h = 1//100)
    @test keys(sim.store.stores) === (Symbol(Float64),)
    # Off nominal the pin keeps a `Float64` buffer of its own beside the `Dual`
    # one, rather than being flattened into it as a zero-partial.
    sim = Simulation(single(PinnedLeaf()), D8; h = 1//100)
    @test Set(keys(sim.store.stores)) == Set([Symbol(D8), Symbol(Float64)])
    init!(sim)
    @test port(sim, "children/c", :a) isa D8
    @test port(sim, "children/c", :frozen) isa Float64
    @test port(sim, "children/c", :frozen) == 2.0  # a stored constant, not a computed product
end

# A mixed-leaf cell: legal in the design (a pinned leaf inside a declared
# struct), but its layout needs multi-cursor addressing this increment does not
# build — a named refusal, not a mislayout.
struct TaggedValue{T}
    v::T
    n::Int
end
struct MixedCell <: AbstractComponent end
output_types(::MixedCell, ::Type{T}) where {T <: Real} = (out = TaggedValue{T},)
h_x(::MixedCell, (; t)) = (out = TaggedValue(t, 1),)

@testset "a mixed-leaf cell is refused by name" begin
    err = failure(() -> build(single(MixedCell())))
    @test err isa BuildError
    @test occursin("mixed-leaf", err.msg) && occursin("out", err.msg)
end

# --- multi-rate: the two registers (§8.7, §10.5, D-185) -------------------------

@testset "the wrappers are plain data over exact rationals (D-185)" begin
    # `Period` and `Hz` are one quantity, two spellings, normalized at
    # construction; floats are refused with the exact spelling named.
    @test period(Hz(50)) === 1//50
    @test period(Period(1//50)) === 1//50
    @test period(Hz(1//2)) === 2//1
    @test occursin("Period(1//50)", failure(() -> Period(0.02)).msg)
    @test occursin("Hz(1//2)", failure(() -> Hz(0.5)).msg)
    @test occursin("Rational", failure(() -> Absolute(Hz(50), 0.001)).msg)
    @test occursin("quantity", failure(() -> Absolute(1//50)).msg)

    # Plain data carriers: no range checks of their own — those are Stratum A's,
    # with path attribution, at the fold.
    @test Relative(5) === Relative(5, 0)
    @test Relative(0).K == 0
    @test Absolute(Period(0)).T == 0
end

@testset "the fold validates with path attribution (§8.7, §9.1)" begin
    rated(rates) = Group((; c = TickCounter()), (), (), (), rates)

    err = failure(() -> build(rated((; var"children/c" = 2))))
    @test err isa BuildError
    @test occursin("whole value vocabulary", err.msg) && occursin("children/c", err.msg)

    for (rates, fragment) in (((; var"children/c" = Relative(0)),      "K ≥ 1"),
                              ((; var"children/c" = Relative(2, 2)),   "0 ≤ φ < K"),
                              ((; var"children/c" = Absolute(Period(0))),        "T > 0"),
                              ((; var"children/c" = Absolute(Hz(50), 1//40)),    "0 ≤ τ < T"))
        err = failure(() -> build(rated(rates)))
        @test err isa BuildError && occursin(fragment, err.msg)
    end

    # A key names an immediate child, and nothing else: a stray name and a deep
    # key meet the same rule.
    err = failure(() -> build(rated((; nav = Relative(2)))))
    @test err isa BuildError && occursin("immediate child", err.msg)
    err = failure(() -> build(rated((; var"children/c/x" = Relative(2)))))
    @test err isa BuildError && occursin("immediate child", err.msg)

    # A key on a continuous child is the Δt-on-continuous error at declaration
    # time: keys name discrete or scope children (§8.7).
    err = failure(() -> build(Group((; c = Gain(1.0)), (), ("in" => "children/c/e",), (),
                                    (; var"children/c" = Relative(2)))))
    @test err isa BuildError && occursin("continuous", err.msg)

    # A bare container field name applies one declaration to every element.
    sim = Simulation(Group((; c1 = TickCounter(), c2 = TickCounter()), (), (), (),
                           (; children = Relative(2, 1))); h = 1//10)
    @test length(sim.sched) == 2
    @test all(e.D == 2 && e.Φ == 1 for e in sim.sched)
end

# --- multi-rate: the fold and the bound schedule (§9.1, §9.2, §10.5) ------------

@testset "the worked example compiles to the spec's pairs (§9.2)" begin
    # Three discrete components under two scopes, deployed at Δt_base = 2 ms:
    # inner (1, 0), outer (5, 2), gnss (10, 0) — §9.2's table, exactly.
    sim = Simulation(MultiRate(); h = 1//500)
    @test sim.n == 1 && sim.Δt_base == 0.002
    @test [(e.path, e.D, e.Φ) for e in sim.sched] ==
          [("fcs/inner", 1, 0), ("fcs/outer", 5, 2), ("gnss", 10, 0)]
    @test [e.Δt for e in sim.sched] ≈ [0.002, 0.01, 0.02]

    # The gate is structural: the interior variants carry no discrete entry, the
    # boundary variants gate every one of them, and nothing else.
    @test isempty(walked(sim.bodies.sweep_hxu, :interior))
    @test gated(sim.bodies.sweep_hxu) == 3
    @test gated(sim.bodies.sweep_hx) == 0            # the ramp is continuous
end

@testset "the hyperperiod chart is readable off the cells (§10.5)" begin
    # The ramp makes every sample carry its acquisition time (1 + t), so each
    # cell says when its owner last ticked — the chart's dots, and the
    # deterministic aging of the stagger, in the spec's own numbers.
    sim = Simulation(MultiRate(); h = 1//500)
    init!(sim)                                       # boundary zero: Φ = 0 is due
    @test port(sim, "fcs/inner", :out) == 1.0
    @test port(sim, "gnss", :out) == 1.0

    run!(sim, 2 * 0.002)                             # base tick 2: outer's first tick
    @test port(sim, "fcs/inner", :out) ≈ 1.004       # fresh at every tick
    @test port(sim, "fcs/outer", :out) == 1.0        # gnss's k = 0 sample: two ticks old
    run!(sim, 7 * 0.002)                             # base tick 7
    @test port(sim, "fcs/outer", :out) == 1.0        # the same sample, seven ticks old
    @test port(sim, "gnss", :out) == 1.0             # gnss itself holds until k = 10
    run!(sim, 10 * 0.002)
    @test port(sim, "gnss", :out) ≈ 1.02             # its second tick
    run!(sim, 12 * 0.002)
    @test port(sim, "fcs/outer", :out) ≈ 1.02        # re-aged two ticks at k = 12
end

@testset "relative scopes compose affinely; anchors sever (§10.5, §9.1)" begin
    # Under a scope at Relative(2, 1): D = K·Dₛ, Φ = Φₛ + φ·Dₛ.
    inner_g = Group((; a = TickCounter(), b = TickCounter()), (), (), (),
                    (; var"children/a" = Relative(1), var"children/b" = Relative(5, 2)))
    sim = Simulation(Group((; f = inner_g), (), (), (),
                           (; var"children/f" = Relative(2, 1))); h = 1//100)
    @test [(e.D, e.Φ) for e in sim.sched] == [(2, 1), (10, 5)]

    # Neither has Φ = 0, so boundary zero admits neither; over base ticks 1…10,
    # `a` ticks at the odd indices and `b` at 5 alone.
    init!(sim)
    @test state(sim, "children/f/children/a") === (n = 0,)
    run!(sim, 10 * 0.01)
    @test state(sim, "children/f/children/a") === (n = 5,)
    @test state(sim, "children/f/children/b") === (n = 1,)

    # An anchor severs: a relative child of an anchored subtree composes against
    # the anchor, not the enclosing grid — D = 3·D₁ with D₁ = (1//50)/(1//500).
    gps = Group((; rx = TickCounter()), (), (), (),
                (; var"children/rx" = Relative(3)))
    sim = Simulation(Group((; gps = gps), (), (), (),
                           (; var"children/gps" = Absolute(Hz(50)))); h = 1//500)
    @test [(e.D, e.Φ) for e in sim.sched] == [(30, 0)]
end

# --- multi-rate: deployment binding (§9.1, §9.2) --------------------------------

@testset "one Build backs many Simulations; Δt_base has three sources (§9.1)" begin
    b = build(MultiRate())

    # The n·h product (the default path), an explicit n, the explicit keyword
    # (Rational or quantity): the anchored divisor is deployment's, not the
    # build's — the same Build lands gnss at D = 10 or D = 5.
    s1 = Simulation(b; h = 1//500)
    s2 = Simulation(b; h = 1//500, n = 2)
    s3 = Simulation(b; h = 1//500, Δt_base = 1//250)
    s4 = Simulation(b; h = 1//500, Δt_base = Period(1//250))
    @test [e.D for e in s1.sched] == [1, 5, 10]
    @test [e.D for e in s2.sched] == [1, 5, 5]
    @test s3.n == 2 && s3.sched == s2.sched == s4.sched

    # Nothing writable is shared: each Simulation materializes its own buffers.
    init!(s1); run!(s1, 0.02)
    init!(s2)
    @test port(s1, "fcs/inner", :out) ≈ 1.02
    @test port(s2, "fcs/inner", :out) == 1.0

    # Cross-validation, collected as plain BuildErrors here (§13.2 is absent):
    @test occursin("h", failure(() -> Simulation(b)).msg)
    @test occursin("float", failure(() -> Simulation(b; h = 1e-3)).msg)
    @test occursin("disagrees", failure(() -> Simulation(b; h = 1//500, Δt_base = 1//250,
                                                         n = 3)).msg)
    @test occursin("harmonic", failure(() -> Simulation(b; h = 1//300,
                                                        Δt_base = 1//500)).msg)
    @test occursin("n must be", failure(() -> Simulation(b; h = 1//500, n = 0)).msg)

    # A non-dividing anchor is refused with its declaring scope and key, and the
    # admissible set is named off the pool.
    err = failure(() -> Simulation(b; h = 1//500, Δt_base = 3//250))
    @test err isa BuildError
    @test occursin("key `gnss`", err.msg) && occursin("gcd(pool)", err.msg)
end

@testset "Δt_base derivation demands an all-anchored model (§9.1)" begin
    # Derivation with an unanchored component present is action at a distance:
    # refused constructively, naming the components whose periods would rescale.
    err = failure(() -> Simulation(build(MultiRate()); h = 1//500, Δt_base = :derive))
    @test err isa BuildError
    @test occursin("fcs/inner", err.msg) && occursin("Δt_base", err.msg)

    # All anchored: the pool is every period and every nonzero offset, and the
    # derived value is its GCD — the offset drives the grid 2× finer here.
    anchored = Group((; c = TickCounter()), (), (), (),
                     (; var"children/c" = Absolute(Hz(50), 1//100)))
    sim = Simulation(anchored; h = 1//500, Δt_base = :derive)
    @test sim.Δt_base == 0.01 && sim.n == 5
    @test [(e.D, e.Φ) for e in sim.sched] == [(2, 1)]
end

# --- multi-rate: the gate at run time (§10.5) -----------------------------------

@testset "boundary zero admits Φ = 0; an offset component holds its probe cells" begin
    # `z` at Relative(2, 1) is not due at boundary zero: until its first tick at
    # Φ·Δt_base its cell holds what the build probe populated — the ramp's value
    # at the probe, a coherent ZOH story (§10.5, §9.3).
    late = Group((; src = Ramp(5.0), z = ZOH()),
                 ("children/src/out" => "children/z/in",), (),
                 ("children/z/out" => "y",),
                 (; var"children/z" = Relative(2, 1)))
    sim = Simulation(late; h = 1//100)
    init!(sim)
    @test port(sim, "", :y) == 5.0                   # probe-populated, not swept
    run!(sim, 0.01)                                  # base tick 1: the first tick
    @test port(sim, "", :y) ≈ 5.01
    run!(sim, 0.02)                                  # base tick 2: not due, holds
    @test port(sim, "", :y) ≈ 5.01
end

@testset "a multi-rate sampled loop matches its exact discretization" begin
    # The increment-3 reference, one level harder: the controller ticks every
    # second base tick and every base tick is two continuous steps, so the exact
    # recursion runs at Δt_ctl = D·Δt_base = 4h — and only matches if the gate
    # admits `ctl` at exactly its own ticks, the hold spans the sub-ticks and
    # off-tick boundaries, and the bundle's Δt is the compiled schedule's
    # (§10.5's single source of truth), not h or Δt_base.
    kI, ω, ζ, r, N = 3.0, 2.0, 0.1, 0.7, 50
    Δt_ctl = 0.02
    A = SMatrix{2,2}(0.0, -ω^2, 1.0, -2ζ * ω)
    B = SVector(0.0, 1.0)
    Ad = exp(A * Δt_ctl)
    Bd = A \ ((Ad - I) * B)
    q, s = SVector(0.0, 0.0), 0.0
    for _ in 1:N
        q, s = Ad * q + Bd * s, s + kI * Δt_ctl * (r - q[1])
    end

    # §10.5's exposed-multiplier idiom: the deployment preference arrives as a
    # constructor parameter, and the declaration stays the assembly's.
    sim = Simulation(SampledLoop(; kI, ω, ζ, ctl_rate = Relative(2)); h = 1//200, n = 2)
    @test [(e.path, e.D, e.Φ) for e in sim.sched] == [("ctl", 2, 0)]
    @test sim.sched[1].Δt ≈ Δt_ctl
    set_slot!(sim, "ref", r)
    init!(sim)
    run!(sim, N * Δt_ctl)
    @test state(sim, "plant").q ≈ q rtol = 1e-6
    @test port(sim, "ctl", :u) ≈ s rtol = 1e-6
end

@testset "the gated boundary walk does not allocate (§7.5)" begin
    sim = Simulation(MultiRate(); h = 1//500)
    init!(sim)
    bods = phase_bodies(sim)
    for name in keys(bods)
        body = bods[name]
        body(); body(1); body(2)
        @test @ballocated($body()) == 0
        @test @ballocated($body(2)) == 0             # a boundary where gates split
    end
    sim2 = Simulation(SampledLoop(; ctl_rate = Relative(2)); h = 1//200, n = 2)
    init!(sim2)
    @test @ballocated(step!($sim2, 0.005)) == 0
    @test @ballocated(offtick_boundary!($sim2)) == 0
    @test @ballocated(boundary!($sim2, 3)) == 0
end
