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

    # Any component may be the root (D-208): a primitive one flattens to the
    # single leaf at the root path, its `input_types` keys the root inputs.
    b = build(Plant())
    @test b.flat.paths == [""]
    @test b.flat.root_inputs == [:u]
end

@testset "a primitive root is the whole model (§8.2, §9.1, D-208)" begin
    # The bare leaf, deployed and driven: the same second-order plant the
    # feedback model wraps, now with `u` a root input rather than a wired port.
    # Under a step `u` held from `t₀`, ẋ = A x + B u integrates exactly.
    ω, ζ, u = 2.0, 0.1, 0.7
    A = SMatrix{2,2}(0.0, -ω^2, 1.0, -2ζ * ω)
    B = SVector(0.0, 1.0)
    exact(t) = exp(A * t) * (A \ (B * u)) - A \ (B * u)

    sim = Simulation(Plant(; ω, ζ); h = 1//1000)
    init!(sim, fragment(inputs = (u = u,)))
    @test port(sim, "", :u) === u              # the root input's own cell
    run!(sim; t_end = 2.0)

    @test state(sim, "").q ≈ exact(2.0) rtol = 1e-8
    @test port(sim, "", :y) ≈ exact(2.0)[1] rtol = 1e-8
    @test port(sim, "", :power) ≈ u * exact(2.0)[2] rtol = 1e-8

    # Totality reaches it like any other root input, and the condition's own
    # vocabulary composes at the root with no `at` prefix in sight.
    @test occursin("UninitializedInputs", failure(() -> init!(sim)).msg)
    sim2 = Simulation(Plant(; ω, ζ); h = 1//1000)
    init!(sim2, combine(condition(Plant(); y = 1.0), fragment(inputs = (u = 0.0,))))
    @test state(sim2, "").q === SVector(1.0, 0.0)
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
    init!(tsim, fragment(inputs = (in = 1.0,)))
    @test port(tsim, "units/2", :out) === 6.0
    @test port(tsim, "", :y) === port(tsim, "units/2", :out)

    # The mixing rule is form-blind.
    err = failure(() -> build(TupleRoster((Gain(1.0), 2.0))))
    @test err isa BuildError
    @test occursin("mixes components", err.msg) && occursin("units", err.msg)
end

# --- paths and the one-level rule (§6.1, §8.6) --------------------------------
# The same instance, held two ways. Under D-207 a wiring endpoint names an
# immediate child and one of its faces, so the route is declared level by level
# and the declaration's knowledge of its own field is no longer the question:
# the concrete and the generic holder wire identically, and neither may reach
# past the child it names.

struct ConcreteHold <: AbstractComponent
    inner::SampledLoop
end
child_connections(::ConcreteHold) = ()
input_connections(::ConcreteHold) = ("ref" => "inner/ref",)
output_connections(::ConcreteHold) = ("inner/y" => "y",)

struct GenericHold{L <: AbstractComponent} <: AbstractComponent
    inner::L
end
child_connections(::GenericHold) = ()
input_connections(::GenericHold) = ("ref" => "inner/ref",)
output_connections(::GenericHold) = ("inner/y" => "y",)

struct PastReach <: AbstractComponent            # one segment further: past it
    inner::SampledLoop
end
child_connections(::PastReach) = ()
input_connections(::PastReach) = ("ref" => "inner/sum/a",)
output_connections(::PastReach) = ("inner/y" => "y",)

struct PastGenericReach{L <: AbstractComponent} <: AbstractComponent   # the same, generic
    inner::L
end
child_connections(::PastGenericReach) = ()
input_connections(::PastGenericReach) = ("ref" => "inner/sum/a",)
output_connections(::PastGenericReach) = ("inner/y" => "y",)

@testset "a wiring endpoint names one child and one of its faces (§6.1, D-207)" begin
    # Routed through the sub-assembly's own face: legal, the routed input is fed
    # exactly once, and the re-exported face aliases the port behind it.
    sim = Simulation(ConcreteHold(SampledLoop()); h = 1//50)
    init!(sim, fragment(inputs = (ref = 1.0,)))
    @test port(sim, "", :y) === port(sim, "inner/plant", :y)

    # The identical declarations against the identical instance, held
    # generically: substitutability now holds at *every* boundary, so the
    # generic holder builds too, and the concrete/generic distinction has left
    # this register entirely.
    gsim = Simulation(GenericHold(SampledLoop()); h = 1//50)
    init!(gsim, fragment(inputs = (ref = 1.0,)))
    @test gsim.flat.paths == sim.flat.paths && gsim.flat.conns == sim.flat.conns
    run!(sim; t_end = 0.2)                       # equal wiring, and equal trajectories:
    run!(gsim; t_end = 0.2)                      # the t₀ table alone would prove nothing
    @test state(gsim, "inner/plant").q === state(sim, "inner/plant").q
    @test port(gsim, "", :y) === port(sim, "", :y)

    # One segment further — the grandchild's own port, bypassing `inner`'s face
    # — is the build error, whatever the field's declared type.
    for bad in (PastReach(SampledLoop()), PastGenericReach(SampledLoop()))
        err = failure(() -> build(bad))
        @test err isa BuildError
        @test occursin("reaches past `inner`", err.msg) && occursin("§6.1", err.msg)
    end
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

struct DoubleFed <: AbstractComponent            # two producers onto one face
    loop::SampledLoop
    src::ModedSource
    src2::ModedSource
end
child_connections(::DoubleFed) = ("src/out" => "loop/ref", "src2/out" => "loop/ref")

struct Doubler <: AbstractComponent              # its own wire onto the input its face routes to
    s::ModedSource
    g::Gain
end
child_connections(::Doubler) = ("s/out" => "g/e",)
input_connections(::Doubler) = ("in" => "g/e",)
output_connections(::Doubler) = ("g/out" => "y",)

struct DoubleFedSibling <: AbstractComponent     # an ancestor's route onto a wired input
    loop::Doubler
    src::ModedSource
end
child_connections(::DoubleFedSibling) = ("src/out" => "loop/in",)

@testset "every input is fed exactly once, across levels (§6.1)" begin
    err = failure(() -> build(Starved(Gain(1.0))))
    @test err isa BuildError
    @test occursin("fed by nothing", err.msg) && occursin("g", err.msg)

    err = failure(() -> build(DoubleFed(SampledLoop(), ModedSource(), ModedSource())))
    @test err isa BuildError
    @test occursin("fed twice", err.msg) && occursin("loop/sum", err.msg)

    # The same rule one level down: the sub-assembly's own wire against the
    # ancestor's route through the face, and the message names both entries.
    err = failure(() -> build(DoubleFedSibling(Doubler(ModedSource(), Gain(1.0)),
                                               ModedSource())))
    @test err isa BuildError
    @test occursin("child_connections at `loop`", err.msg) &&
          occursin("child_connections at the root component", err.msg)

    # The one legitimate terminus: the root's own input faces are the root inputs,
    # authored by the init service's condition (§11.3, §14.6).
    sim = Simulation(Group((; c = Gain(2.0)), (), ("in" => "children/c/e",), ());
                     h = 1//100)
    @test sim.flat.root_inputs == [:in]
    init!(sim, fragment(inputs = (in = 0.0,)))
    @test port(sim, "children/c", :out) == 0.0
    init!(sim, fragment(inputs = (in = 3.0,)))
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
