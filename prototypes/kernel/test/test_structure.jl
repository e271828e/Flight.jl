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
    # A container's elements are children *of the parent*, in declaration order.
    # `Group` names them bare, its `children` being transparent (D-211); the
    # naming rule itself is the undeclared containers' below.
    sim = Simulation(Group((; c1 = TickCounter(), c2 = TickCounter()));
                     h = 1//10)
    @test sim.flat.paths == ["c1", "c2"]
    @test state(sim, "c2") === (n = 0,)

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

# --- name-transparent containers (§8.5, D-211) --------------------------------
# One container field per type may be declared name-transparent, and its
# elements then go by bare key everywhere a child name appears. `Group` is the
# library case, exercised throughout these files; below is the rule itself, on
# named types, beside `TupleRoster` — the undeclared container it leaves alone.

struct TransparentRoster{U <: NamedTuple} <: AbstractComponent
    units::U
end
child_connections(::TransparentRoster)     = ("a/out" => "b/e",)
input_connections(::TransparentRoster)     = ("in" => "a/e",)
output_connections(::TransparentRoster)    = ("b/out" => "y",)
transparent_container(::TransparentRoster) = :units

struct Colliding <: AbstractComponent            # a bare key against a sibling field
    kids::NamedTuple
    c1::TickCounter
end
child_connections(::Colliding) = ()
transparent_container(::Colliding) = :kids

struct Pathological <: AbstractComponent         # ...and against another container's element
    kids::NamedTuple
    units::Tuple
end
child_connections(::Pathological) = ()
transparent_container(::Pathological) = :kids

struct SelfNamed <: AbstractComponent            # a bare key equal to its own field's name
    kids::NamedTuple
end
child_connections(::SelfNamed) = ()
transparent_container(::SelfNamed) = :kids

struct Shadowed{K <: NamedTuple, U <: Tuple} <: AbstractComponent
    kids::K                                      # ...and equal to a sibling *container's*
    units::U                                     # — whose own children it would hide
    trim::Gain
end
transparent_container(::Shadowed) = :kids
child_connections(::Shadowed)  = ("trim/out" => "units/e",)
input_connections(::Shadowed)  = ("in" => "trim/e",)
output_connections(::Shadowed) = ("units/out" => "y",)

struct OpaqueDeclared <: AbstractComponent       # the declaration names a component field
    c::TickCounter
end
child_connections(::OpaqueDeclared) = ()
transparent_container(::OpaqueDeclared) = :c

struct AbsentDeclared <: AbstractComponent       # ...and here, no field of the type at all
    kids::NamedTuple
end
child_connections(::AbsentDeclared) = ()
transparent_container(::AbsentDeclared) = :nope

@testset "a name-transparent container contributes bare keys (§8.5, D-211)" begin
    # Naming is the only thing the declaration changes: the same two children in
    # the same declaration order, addressed without the field segment — wiring
    # endpoints, the flat list and the read path alike.
    sim = Simulation(TransparentRoster((a = Gain(2.0), b = Gain(3.0))); h = 1//10)
    @test sim.flat.paths == ["a", "b"]
    init!(sim, fragment(inputs = (in = 1.0,)))
    @test port(sim, "b", :out) === 6.0
    @test port(sim, "", :y) === port(sim, "b", :out)

    # The undeclared container keeps its key segment: `TupleRoster` above wires
    # and reads the same topology as `"units/1"`, and the default is `nothing`.
    @test build(TupleRoster((Gain(2.0), Gain(3.0)))).flat.paths == ["units/1", "units/2"]
    @test transparent_container(TupleRoster((Gain(1.0),))) === nothing
    @test transparent_container(Group((;))) === :children

    # Two children may not share a name, whatever produced them. The check is
    # general — bare keys only make the case reachable — and it names both
    # parties: here a bare key against a sibling field, and against another
    # container's composite name.
    err = failure(() -> build(Colliding((c1 = TickCounter(),), TickCounter())))
    @test err isa BuildError && occursin("two children are named `c1`", err.msg)
    @test occursin("field `c1`", err.msg) &&
          occursin("name-transparent container field `kids`, element `c1`", err.msg)
    err = failure(() -> build(Pathological((var"units/1" = TickCounter(),),
                                           (TickCounter(),))))
    @test err isa BuildError && occursin("two children are named `units/1`", err.msg)

    # The family's other two arms, both reachable only by a bare key and neither
    # of them a duplicate *child* name, so the check above can see neither. An
    # element keyed with its own field's name is indistinguishable from
    # `sample_times`' field-name sugar (§8.7)...
    err = failure(() -> build(SelfNamed((kids = TickCounter(),))))
    @test err isa BuildError && occursin("the bare key `kids`", err.msg)
    @test occursin("name-transparent container field `kids`, element `kids`", err.msg) &&
          occursin("`sample_times`' field-name sugar", err.msg)

    # ...and one equal to a sibling container field's name shadows the
    # `"field/key"` grammar that reaches *its* children: no child bears the bare
    # name, so nothing collides, yet `"units/1/e"` would resolve to the bare
    # child and the one-level rejection would blame the wrong party (§6.1).
    err = failure(() -> build(Shadowed((units = Gain(2.0),), (Gain(3.0),), Gain(3.0))))
    @test err isa BuildError && occursin("the bare key `units`", err.msg)
    @test occursin("name-transparent container field `kids`, element `units`", err.msg) &&
          occursin("collides with container field `units`", err.msg)

    # The exemption, on the same type one instantiation away: an *empty* sibling
    # container reaches no children, so there is no grammar to shadow and the
    # bare key stands. It has to be that way round — an empty `Tuple` is an empty
    # container and empty inert data at once, so reserving its name would refuse
    # over inert data, a false positive where no shadow exists. Legality is
    # per-instantiation, as every wiring judgment already is.
    sim = Simulation(Shadowed((units = Gain(2.0),), (), Gain(3.0)); h = 1//10)
    @test sim.flat.paths == ["units", "trim"]          # the bare child, and nothing under it
    init!(sim, fragment(inputs = (in = 1.0,)))
    @test port(sim, "units", :out) === 6.0             # reads resolve `units` bare
    @test port(sim, "", :y) === 6.0                    # and so does the wiring register

    # And the declaration must name a container field of the type — a component
    # field and an absent name are refused alike.
    for bad in (OpaqueDeclared(TickCounter()), AbsentDeclared((; c = TickCounter())))
        err = failure(() -> build(bad))
        @test err isa BuildError && occursin("names no container field", err.msg)
    end
end

@testset "`Group`'s keyword form normalizes a bare `Pair` (§8.5, D-211)" begin
    g = Group((; c = Gain(2.0)); inputs = "in" => "c/e", outputs = "c/out" => "y")
    @test input_connections(g) == ("in" => "c/e",)
    @test output_connections(g) == ("c/out" => "y",)
    @test child_connections(g) == () && sample_times(g) == (;)

    # A tuple passes through as written, and every unnamed keyword is empty.
    w = Group((; a = Gain(1.0), b = Gain(2.0)); wires = ("a/out" => "b/e",))
    @test child_connections(w) == ("a/out" => "b/e",)
    @test input_connections(w) == () && output_connections(w) == ()

    # The normalized declarations are the ones that build.
    sim = Simulation(g; h = 1//10)
    init!(sim, fragment(inputs = (in = 2.0,)))
    @test port(sim, "", :y) === 4.0
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

struct RootCollision <: AbstractComponent        # one key in both contracts
end
input_types(::RootCollision, ::Type{T}) where {T <: Real} = (u = T,)
output_types(::RootCollision, ::Type{T}) where {T <: Real} = (u = T, v = T)
h_xu(::RootCollision, (; u)) = (u = 2u.u, v = 1.0)

struct DeadFace <: AbstractComponent             # a face routed to nothing at all
    g::Gain
end
child_connections(::DeadFace) = ()
input_connections(::DeadFace) = ("in" => "g/e", "dead" => ())
output_connections(::DeadFace) = ("g/out" => "y",)

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

    # Uniqueness follows the root's class, not its family (D-210): a primitive
    # root's faces are its two contract declarations' keys together, and a
    # shared key would place the root input's cell over the output port's — the
    # authored value read back as the stage's, silently.
    err = failure(() -> build(RootCollision()))
    @test err isa BuildError && occursin("appear twice", err.msg) && occursin("§8.6", err.msg)

    # Below the root the same leaf is untouched: its input face aliases its
    # producer's cell and places nothing, so there is no collision to forbid.
    @test build(fed(RootCollision(), "u")) isa Build
end

@testset "an `input_connections` entry routes to at least one endpoint (§8.6, D-210)" begin
    # At the root and one level down alike, and at declaration: the entry is
    # named, and the refusal is a declaration error rather than the starvation
    # further down the tree that the shape used to surface as.
    err = failure(() -> build(DeadFace(Gain(1.0))))
    @test err isa BuildError && occursin("routes to no internal endpoint", err.msg)
    @test occursin("input_connections at the root component", err.msg)

    nested = Group((; sub = DeadFace(Gain(2.0))); inputs = ("in" => "sub/in",))
    err = failure(() -> build(nested))
    @test err isa BuildError && occursin("routes to no internal endpoint", err.msg)
    @test occursin("input_connections at `sub`", err.msg)

    # The misdiagnosis this closes: the dead face used to build, and a condition
    # addressing it read "declares no input face" — byte-identical to the message
    # a bare typo earns. Authoring one can no longer get that far.
    err = failure(() -> resolve_condition(at("sub", fragment(inputs = (dead = 1.0,))),
                                build(nested)))
    @test occursin("routes to no internal endpoint", err.msg) &&
          !occursin("declares no input face", err.msg)
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
    sim = Simulation(Group((; c = Gain(2.0)); inputs = ("in" => "c/e",));
                     h = 1//100)
    @test sim.flat.root_inputs == [:in]
    init!(sim, fragment(inputs = (in = 0.0,)))
    @test port(sim, "c", :out) == 0.0
    init!(sim, fragment(inputs = (in = 3.0,)))
    @test port(sim, "c", :out) == 6.0
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

# --- §13.3's primitives and §8.8's passthrough pair -----------------------------
# The declaration surface computing its own boundary: `input_connections` and
# `output_connections` are ordinary functions of the instance, so an assembly may
# derive its entries from a child's contract. The helpers are the pass-through
# case, and everything they return is an ordinary pair.

# The child under test: two input faces and two output faces, so both filters
# have something to bite on either side.
faced() = Group((; s = Sum(), g = Gain(2.0));
                wires = "s/e" => "g/e",
                inputs = ("a" => "s/a", "b" => "s/b"),
                outputs = ("s/e" => "sum", "g/out" => "scaled"))

# `a` is fed here, so it leaves the input face surface; `b` is passed through,
# and `scaled` is re-exported one level up. Nothing else about the two
# declarations is authored.
struct Passed{C <: AbstractComponent} <: AbstractComponent
    inner::C
    trim::Gain
end
child_connections(::Passed) = ("trim/out" => "inner/a",)
input_connections(p::Passed) = (input_passthrough(p, "inner"; except = ("a",))...,
                                "e" => "trim/e")
output_connections(p::Passed) = output_passthrough(p, "inner"; only = ("scaled",))

# The same assembly with both boundaries written out by hand: the twin the
# computed one must match entry for entry.
struct HandWired{C <: AbstractComponent} <: AbstractComponent
    inner::C
    trim::Gain
end
child_connections(::HandWired) = ("trim/out" => "inner/a",)
input_connections(::HandWired) = ("inner.b" => "inner/b", "e" => "trim/e")
output_connections(::HandWired) = ("inner/scaled" => "inner.scaled",)

# `prefix = ""` drops the prefixing, and the bare `b` then collides with the
# hand-written face beside it — the build's own uniqueness check, not the
# helper's business (§8.8).
struct Unprefixed{C <: AbstractComponent} <: AbstractComponent
    inner::C
    trim::Gain
end
child_connections(::Unprefixed) = ("trim/out" => "inner/a",)
input_connections(u::Unprefixed) =
    (input_passthrough(u, "inner"; prefix = "", except = ("a",))..., "b" => "trim/e")

# A face both wired and passed through: `except` is missing `a`, so the wire and
# the route both claim it.
struct DoubleClaimed{C <: AbstractComponent} <: AbstractComponent
    inner::C
    trim::Gain
end
child_connections(::DoubleClaimed) = ("trim/out" => "inner/a",)
input_connections(d::DoubleClaimed) = (input_passthrough(d, "inner")..., "e" => "trim/e")

# The same computed boundary over a name-transparent container: the helper's
# `child_path` is a bare key, and the `Group` it names is a child like any other
# (D-211).
struct PassedGroup{U <: NamedTuple} <: AbstractComponent
    units::U
end
transparent_container(::PassedGroup) = :units
child_connections(::PassedGroup) = ("trim/out" => "inner/a",)
input_connections(p::PassedGroup) = (input_passthrough(p, "inner"; except = ("a",))...,
                                     "e" => "trim/e")
output_connections(p::PassedGroup) = output_passthrough(p, "inner"; only = ("scaled",))

@testset "the §13.3 primitives resolve one level and list faces in order" begin
    m = feedback_model()
    @test resolve(m, "sum") === m.children.sum
    @test resolve_terminal(m, "sum/a") === (m.children.sum, "a")

    # Declaration order, both classes, and the `T`-independent key set read at
    # the nominal activation.
    @test input_faces(resolve(m, "sum")) == ["a", "b"]
    @test output_faces(resolve(m, "plant")) == ["y", "power"]
    @test input_faces(SampledLoop()) == ["ref"]
    @test output_faces(SampledLoop()) == ["y", "cmd", "power"]
    @test input_faces(faced()) == ["a", "b"] && output_faces(faced()) == ["sum", "scaled"]

    # One level, the same rule wiring resolution runs: a deeper path is a build
    # error naming the child it reaches past, and an unknown segment comes with
    # the sibling list in hand.
    err = failure(() -> resolve(m, "sum/a"))
    @test err isa BuildError && occursin("reaches past `sum`", err.msg)
    err = failure(() -> resolve(m, "nope"))
    @test err isa BuildError && occursin("its children are plant, ctl, sum", err.msg)
end

@testset "the §8.8 passthrough pair computes what a hand-wired twin declares" begin
    p, w = Passed(faced(), Gain(3.0)), HandWired(faced(), Gain(3.0))

    # The computed entries are the authored ones, pair for pair.
    @test input_connections(p) == input_connections(w)
    @test output_connections(p) == output_connections(w)
    @test input_connections(p) == ("inner.b" => "inner/b", "e" => "trim/e")
    @test output_connections(p) == ("inner/scaled" => "inner.scaled",)

    # And the two builds are the same model: same root inputs, same exported
    # faces, same schedule, same trajectory.
    bp, bw = build(p), build(w)
    @test bp.flat.root_inputs == bw.flat.root_inputs == [:var"inner.b", :e]
    @test bp.flat.out_faces == bw.flat.out_faces
    @test bp.flat.paths == bw.flat.paths
    sp, sw = Simulation(p; h = 1//10), Simulation(w; h = 1//10)
    cond = fragment(inputs = (var"inner.b" = 1.0, e = 2.0))
    init!(sp, cond); init!(sw, cond)
    @test port(sp, "", :var"inner.scaled") === port(sw, "", :var"inner.scaled")
    @test port(sp, "", :var"inner.scaled") === 2.0 * (3.0 * 2.0 - 1.0)
end

@testset "the passthrough filters, and refuses what it cannot mean (§8.8)" begin
    m = faced()

    # `except` drops, `only` keeps — in the author's order — and the labelling
    # keywords are independent of both.
    @test input_passthrough(m, "s") == ("s.a" => "s/a", "s.b" => "s/b")
    @test input_passthrough(m, "s"; except = ("a",)) == ("s.b" => "s/b",)
    @test input_passthrough(m, "s"; only = ("b", "a")) == ("s.b" => "s/b", "s.a" => "s/a")
    @test input_passthrough(m, "s"; prefix = "env", sep = "_") ==
          ("env_a" => "s/a", "env_b" => "s/b")
    @test input_passthrough(m, "s"; prefix = "") == ("a" => "s/a", "b" => "s/b")

    # The output side is the mirror, its pairs reading along the flow.
    @test output_passthrough(m, "g") == ("g/out" => "g.out",)
    @test output_passthrough(m, "s"; only = ("e",), prefix = "") == ("s/e" => "e",)

    # Exclusivity is enforced, not documented.
    err = failure(() -> input_passthrough(m, "s"; except = ("a",), only = ("b",)))
    @test err isa BuildError && occursin("mutually exclusive", err.msg)

    # A filter naming a face the child does not have errors with the list in
    # hand, on either side.
    err = failure(() -> input_passthrough(m, "s"; only = ("z",)))
    @test err isa BuildError && occursin("`z` names no face", err.msg) &&
          occursin("its faces are `a`, `b`", err.msg)
    err = failure(() -> output_passthrough(m, "g"; except = ("z",)))
    @test err isa BuildError && occursin("its faces are `out`", err.msg)

    # A deeper `child_path` meets the one-level rejection like any endpoint.
    err = failure(() -> input_passthrough(m, "s/a"))
    @test err isa BuildError && occursin("reaches past `s`", err.msg)
end

@testset "every computed entry meets the build's own checks (§8.8)" begin
    # `prefix = ""` collides with a hand-written face beside it, and the
    # uniqueness check does not care which of the two was computed.
    err = failure(() -> build(Unprefixed(faced(), Gain(3.0))))
    @test err isa BuildError && occursin("face name(s) b", err.msg) &&
          occursin("appear twice", err.msg)

    # A face both wired and passed through is a two-producers error, named at
    # both claimants.
    err = failure(() -> build(DoubleClaimed(faced(), Gain(3.0))))
    @test err isa BuildError && occursin("is fed twice", err.msg)
    @test occursin("child_connections", err.msg) && occursin("input_connections", err.msg)
end

@testset "the helpers address a transparent container's child by bare key (D-211)" begin
    g = PassedGroup((inner = faced(), trim = Gain(3.0)))
    @test input_connections(g) == ("inner.b" => "inner/b", "e" => "trim/e")
    @test output_connections(g) == ("inner/scaled" => "inner.scaled",)

    sim = Simulation(g; h = 1//10)
    @test sim.flat.paths == ["inner/s", "inner/g", "trim"]
    init!(sim, fragment(inputs = (var"inner.b" = 1.0, e = 2.0)))
    @test port(sim, "", :var"inner.scaled") === 2.0 * (3.0 * 2.0 - 1.0)
end
