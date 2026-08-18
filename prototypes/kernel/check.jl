# Acceptance tests for the kernel prototype: the continuous-tier walking
# skeleton (increment 2) and the pieces of increment 3 landed so far.
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

# The entries a phase body walks, per variant (§10.5): the boundary variant is
# the full list, the interior one the continuous entries alone.
walked(body, variant = :boundary) =
    [e for c in getfield(body, variant) for e in c.entries]

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
    sim = build(feedback_model())
    # sum first (both its inputs are loop-breaking), then ctl, then plant
    paths = [e.comp isa Sum ? :sum : e.comp isa Gain ? :ctl : :plant
             for e in walked(sim.bodies.sweep_hxu)]
    @test paths == [:sum, :ctl, :plant]
    @test length(walked(sim.bodies.sweep_hx)) == 1
    @test length(walked(sim.bodies.rhs)) == 1
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
    types(body) = unique(typeof(e) for e in walked(body))
    @test length(types(sim.bodies.sweep_hx)) == 1     # two Plants, one h_x body
    @test length(types(sim.bodies.sweep_hxu)) == 3    # Plant, Gain, Sum
    @test length(types(sim.bodies.rhs)) == 1

    # The discrete tier keeps the property: a state store is a `Ref` whose
    # *type* every instance of a component type shares, so the store lives in a
    # field and two counters still compile to one `g` body.
    counters = build(ModelSpec([ComponentSpec(:c1, TickCounter()),
                                ComponentSpec(:c2, TickCounter())]); Δt = 0.1)
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

struct Inert end

@testset "the probe rejects malformed components (§9.3)" begin
    one(c) = ModelSpec([ComponentSpec(:c, c)])
    @test_throws BuildError build(one(Undeclared()))
    @test_throws BuildError build(one(Unproduced()))
    @test_throws BuildError build(one(BadDerivative()))
    @test_throws BuildError build(one(NoFlow()))
end

# --- tier classification (§8.2) -----------------------------------------------
# Tier is read off the declaration shape. These components are the shapes the
# classifier has to separate, plus the four ways a declaration set can disagree.

struct DiscreteCounter end                 # stateful discrete: `g` decides
init_x(::DiscreteCounter) = (n = 0,)
output_types(::DiscreteCounter) = (n = Int,)
h_x(::DiscreteCounter, (; x)) = (n = x.n,)
g(::DiscreteCounter, (; x)) = (n = x.n + 1,)

struct DiscreteMap end                     # stateless discrete: the arity decides
input_types(::DiscreteMap) = (a = Int,)
output_types(::DiscreteMap) = (b = Int,)
h_xu(::DiscreteMap, (; u)) = (b = 2u.a,)

struct BothUpdates end                     # `f` and `g` on one component
init_x(::BothUpdates) = (q = 1.0,)
output_types(::BothUpdates, ::Type{T}) where {T <: Real} = (a = T,)
h_x(::BothUpdates, (; x)) = (a = x.q,)
f(::BothUpdates, (; x)) = (q = 0.0,)
g(::BothUpdates, (; x)) = (q = x.q,)

struct WrongArity end                      # `g` beside a two-argument contract
init_x(::WrongArity) = (n = 0,)
output_types(::WrongArity, ::Type{T}) where {T <: Real} = (a = T,)
h_x(::WrongArity, (; x)) = (a = 1.0,)
g(::WrongArity, (; x)) = (n = x.n,)

struct ModesOnDiscrete end                 # `init_m` is continuous-only
init_x(::ModesOnDiscrete) = (n = 0,)
init_m(::ModesOnDiscrete) = (phase = :idle,)
output_types(::ModesOnDiscrete) = (a = Int,)
h_x(::ModesOnDiscrete, (; x)) = (a = x.n,)
g(::ModesOnDiscrete, (; x)) = (n = x.n,)

struct BothArities end                     # a member of both contract families
output_types(::BothArities, ::Type{T}) where {T <: Real} = (a = T,)
output_types(::BothArities) = (a = Float64,)
h_x(::BothArities, (; t)) = (a = 1.0,)

@testset "tier is read off the declaration shape (§8.2)" begin
    one(c) = ComponentSpec(:c, c)
    # The two deciders: the update law for a stateful leaf, the contract arity
    # for a stateless one.
    @test classify_tier(one(Plant())) === CONTINUOUS
    @test classify_tier(one(Gain(1.0))) === CONTINUOUS
    @test classify_tier(one(DiscreteCounter())) === DISCRETE
    @test classify_tier(one(DiscreteMap())) === DISCRETE

    # The discrete bundle sets: `Δt` is a discrete-tier fact, `m` a continuous
    # one, both under the shared state letter (D-173).
    @test bundle_names(h_x, DiscreteCounter(), DISCRETE, ()) === (:x, :t, :Δt)
    @test bundle_names(g, DiscreteCounter(), DISCRETE, (:n,)) === (:x, :y, :t, :Δt)
    @test bundle_names(h_xu, DiscreteMap(), DISCRETE, ()) === (:u, :t, :Δt)

    # Disagreement names the offending declaration and the tier the rest announce.
    for (c, offender) in ((BothUpdates(), "g"), (WrongArity(), "output_types"),
                          (ModesOnDiscrete(), "init_m"), (BothArities(), "output_types"))
        err = try
            classify_tier(one(c))
            nothing
        catch e
            e
        end
        @test err isa BuildError
        @test occursin(offender, err.msg)
    end

    # A component that declares nothing and defines no stage cannot be
    # intentional (D-164); a store with no update law is the §8.2 sibling.
    @test_throws BuildError classify_tier(one(Inert()))
    @test_throws BuildError classify_tier(one(NoFlow()))

    # A discrete tier needs its base tick period before the executor exists:
    # `Δt` is entry-field data, not a run-time argument (§9.7).
    err = try
        build(ModelSpec([one(DiscreteCounter())]))
        nothing
    catch e
        e
    end
    @test err isa BuildError && occursin("Δt", err.msg)
    @test build(ModelSpec([one(DiscreteCounter())]); Δt = 0.1) isa Simulation
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

    sim = build(sampled_loop(; kI, ω, ζ); Δt)
    set_slot!(sim, :sum, :a, r)
    init!(sim)
    run!(sim, N * Δt, Δt)

    # The tolerance is RK4's, not the semantics': the reference integrates
    # exactly where the loop takes one RK4 step per sample, which at `ωΔt` this
    # size leaves ~1e-8 relative. A violated hold moves `u` mid-interval and
    # costs ~1e-3, so the margin still separates the two by orders of magnitude.
    @test state(sim, :plant).q ≈ q rtol = 1e-6
    # The cell holds `s[N]` — the output stage ran at this boundary — while the
    # store already holds `s[N+1]`. That gap *is* the sampled-data ordering:
    # output stages before updates, within one boundary (§10.5).
    @test port(sim, :ctl, :u) ≈ s rtol = 1e-6
    @test state(sim, :ctl).s ≈ s_next rtol = 1e-6
end

@testset "the ZOH holds by compile-time absence (§10.5)" begin
    sim = build(sampled_loop(); Δt = 0.02)
    set_slot!(sim, :sum, :a, 1.0)      # excite the loop, or nothing moves at all
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
    run!(sim, 0.1, 0.02)
    u₀, y₀ = port(sim, :ctl, :u), port(sim, :plant, :y)
    @test u₀ != 0.0
    step!(sim, 0.02)
    @test port(sim, :ctl, :u) == u₀        # untouched: never gathered, never written
    @test port(sim, :plant, :y) != y₀
end

@testset "discrete state and modes live outside the buffer (§7.3)" begin
    sim = build(ModelSpec([ComponentSpec(:counter, TickCounter()),
                           ComponentSpec(:moded, ModedSource())]); Δt = 0.1)
    # The flat buffer is continuous state only; the counter's `Int` is in its
    # own store, and no store mirrors another.
    @test isempty(sim.xbuf)
    @test state(sim, :counter) === (n = 0,)
    @test modes(sim, :moded) === (phase = :idle,)

    # `Int` and `Bool` cells force their own buffers — the plural in
    # "per-eltype stores", first exercised here.
    @test Set(keys(sim.store.stores)) == Set([Symbol(Int), Symbol(Bool), Symbol(Float64)])

    init!(sim)
    @test state(sim, :counter) === (n = 1,)      # boundary zero is a tick
    @test port(sim, :counter, :n) === 0          # the cell holds what it published
    @test port(sim, :counter, :even) === true
    run!(sim, 0.5, 0.1)
    # Six boundaries: zero, then one per step. The store leads the cell by one,
    # the update having run after the output stage at each of them.
    @test state(sim, :counter) === (n = 6,)
    @test port(sim, :counter, :n) === 5
end

@testset "the workspace is scratch, on both tiers (§7.3)" begin
    sim = build(ModelSpec([ComponentSpec(:sm, Smoother(0.5),
                                         [:a => (:src, :out), :b => (:src, :out)]),
                           ComponentSpec(:src, ModedSource()),
                           ComponentSpec(:wg, WorkGain(2.0), [:in => (:src, :out)])]);
                Δt = 0.1)
    init!(sim)
    @test port(sim, :wg, :out) == 0.0            # 2 × the idle-phase constant
    run!(sim, 0.3, 0.1)
    @test state(sim, :sm).v isa SVector{2,Float64}

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
struct WidenedUpdate end
init_x(::WidenedUpdate) = (n = 0,)
output_types(::WidenedUpdate) = (n = Int,)
h_x(::WidenedUpdate, (; x)) = (n = x.n,)
g(::WidenedUpdate, (; x)) = (n = x.n + 0.5,)

@testset "a discrete successor is the store's own type (§7.3)" begin
    one(c) = ModelSpec([ComponentSpec(:c, c)])
    @test_throws BuildError build(one(WidenedUpdate()); Δt = 0.1)
end

@testset "the discrete tier is frozen at a non-nominal activation (§7.2)" begin
    # The plant's input is declared `T` and wired to a discrete `Float64` cell:
    # a lawful arrival, embedded as a zero-partial. `ctl`'s stages do not run
    # at all, and its cell holds what the probe wrote — which is exactly what a
    # tick at `t₀⁻` would have produced.
    sim = build(sampled_loop(); Δt = 0.02)
    simd = build(sampled_loop(), D8; Δt = 0.02)
    @test isempty(walked(simd.bodies.ticks))
    @test length(walked(simd.bodies.sweep_hx)) == 1          # plant only; ctl frozen
    @test port(simd, :ctl, :u) isa Float64

    init!(simd)
    run!(simd, D8(0.04), D8(0.02))
    @test state(simd, :plant).q isa SVector{2,D8}
    @test port(simd, :ctl, :u) == 0.0                        # held, never recomputed
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

@testset "a pinned leaf lives in its own store (D-166, D-162)" begin
    one(c) = ModelSpec([ComponentSpec(:c, c)])
    # Nominally the pin and the activation scalar coincide: one buffer.
    sim = build(one(PinnedLeaf()))
    @test keys(sim.store.stores) === (Symbol(Float64),)
    # Off nominal the pin keeps a `Float64` buffer of its own beside the `Dual`
    # one, rather than being flattened into it as a zero-partial.
    sim = build(one(PinnedLeaf()), D8)
    @test Set(keys(sim.store.stores)) == Set([Symbol(D8), Symbol(Float64)])
    init!(sim)
    @test port(sim, :c, :a) isa D8
    @test port(sim, :c, :frozen) isa Float64
    @test port(sim, :c, :frozen) == 2.0    # a stored constant, not a computed product
end

# A mixed-leaf cell: legal in the design (a pinned leaf inside a declared
# struct), but its layout needs multi-cursor addressing this increment does not
# build — a named refusal, not a mislayout.
struct TaggedValue{T}
    v::T
    n::Int
end
struct MixedCell end
output_types(::MixedCell, ::Type{T}) where {T <: Real} = (out = TaggedValue{T},)
h_x(::MixedCell, (; t)) = (out = TaggedValue(t, 1),)

@testset "a mixed-leaf cell is refused by name" begin
    one(c) = ModelSpec([ComponentSpec(:c, c)])
    err = try
        build(one(MixedCell()))
        nothing
    catch e
        e
    end
    @test err isa BuildError
    @test occursin("mixed-leaf", err.msg) && occursin("out", err.msg)
end
