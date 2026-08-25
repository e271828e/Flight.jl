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
init_s(::DiscreteCounter) = (n = 0,)
output_types(::DiscreteCounter) = (n = Int,)
h_s(::DiscreteCounter, (; s)) = (n = s.n,)
g(::DiscreteCounter, (; s)) = (n = s.n + 1,)

struct DiscreteMap <: AbstractComponent       # stateless discrete: the arity decides
end
input_types(::DiscreteMap) = (a = Int,)
output_types(::DiscreteMap) = (b = Int,)
h_su(::DiscreteMap, (; u)) = (b = 2u.a,)

struct BothUpdates <: AbstractComponent       # `f` and `g` on one component
end
init_x(::BothUpdates) = (q = 1.0,)
output_types(::BothUpdates, ::Type{T}) where {T <: Real} = (a = T,)
h_x(::BothUpdates, (; x)) = (a = x.q,)
f(::BothUpdates, (; x)) = (q = 0.0,)
g(::BothUpdates, (; x)) = (q = x.q,)

struct WrongArity <: AbstractComponent        # `g` beside a two-argument contract
end
init_s(::WrongArity) = (n = 0,)
output_types(::WrongArity, ::Type{T}) where {T <: Real} = (a = T,)
h_s(::WrongArity, (; s)) = (a = 1.0,)
g(::WrongArity, (; s)) = (n = s.n,)

struct WrongLetter <: AbstractComponent       # a continuous stage name on a `g` leaf
end
init_s(::WrongLetter) = (n = 0,)
output_types(::WrongLetter) = (a = Int,)
h_x(::WrongLetter, (; s)) = (a = s.n,)
g(::WrongLetter, (; s)) = (n = s.n,)

struct ModesOnDiscrete <: AbstractComponent   # `init_m` is continuous-only
end
init_s(::ModesOnDiscrete) = (n = 0,)
init_m(::ModesOnDiscrete) = (phase = :idle,)
output_types(::ModesOnDiscrete) = (a = Int,)
h_s(::ModesOnDiscrete, (; s)) = (a = s.n,)
g(::ModesOnDiscrete, (; s)) = (n = s.n,)

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
    # one, and the state letters are the tiers' own — `s`/`y_s` here against
    # `x`/`y_x` above, disjoint by construction (D-195).
    @test bundle_names(h_s, DiscreteCounter(), DISCRETE, ()) === (:s, :t, :Δt)
    @test bundle_names(g, DiscreteCounter(), DISCRETE, (:n,)) === (:s, :y, :t, :Δt)
    @test bundle_names(h_su, DiscreteMap(), DISCRETE, ()) === (:u, :t, :Δt)
    @test bundle_names(h_su, DiscreteCounter(), DISCRETE, (:n,)) === (:s, :y_s, :t, :Δt)

    # Disagreement names the offending declaration and the tier the rest
    # announce — including the wrong-letter case the split families restore, a
    # continuous stage name on a leaf whose update law is `g` (D-195).
    for (c, offender) in ((BothUpdates(), "g"), (WrongArity(), "output_types"),
                          (WrongLetter(), "h_x"),
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
