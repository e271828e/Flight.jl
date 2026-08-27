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
    init!(sim, fragment(inputs = (ref = r,)))
    run!(sim; t_end = N * Δt)
    @test state(sim, "loop/plant").q ≈ q rtol = 1e-6
    @test port(sim, "loop", :cmd) ≈ s rtol = 1e-6
end

@testset "a face's type and tier are its internal endpoint's (§8.6)" begin
    sim = Simulation(Vehicle(); h = 1//50)
    init!(sim, fragment(inputs = (ref = 1.0,)))
    run!(sim; t_end = 0.1)
    # A face is its endpoint: no cell of its own, at any level of re-export.
    @test port(sim, "loop", :y) === port(sim, "loop/plant", :y)
    @test port(sim, "", :y) === port(sim, "loop/plant", :y)
    @test port(sim, "", :cmd) === port(sim, "loop/ctl", :u)
    # And a two-level chain aliases the one cell all the way down: the vehicle's
    # `power` re-exports the loop's own `power` face, which re-exports the
    # plant's port — one alias per level, no cell of its own at either (D-207).
    @test port(sim, "", :power) === port(sim, "loop", :power) ===
          port(sim, "loop/plant", :power)

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
    sim = Simulation(Group((; c = ConstantBranch()); inputs = ("in" => "c/in",)),
                     D8; h = 1//100)
    init!(sim, fragment(inputs = (in = 0.0,)))
    # What the table holds is the cell's type, the constant embedded into it.
    @test port(sim, "c", :out) isa D8
    @test port(sim, "c", :vec) isa SVector{2,D8}
    @test ForwardDiff.value(port(sim, "c", :vec)[2]) == 1.0

    # The converse is not accepted: a `Dual` at a pinned leaf is an error, with
    # the hint that names the one honest cause. It fails at the `Dual`
    # activation's own Stratum-C re-run, not at `build` (§9.4's lazy lurk).
    b = build(single(PinnedGetsDual()))
    err = failure(() -> activation(b, D8))
    @test err isa BuildError
    @test occursin("participates in differentiation", err.msg)
end

# The activation seam (§9.1, §9.4): the nominal activation runs at build, any
# other is a cached Stratum-C re-run, and a frozen component's products are
# carried across from the nominal activation rather than probed or synthesized.
struct NomSource <: AbstractComponent end
output_types(::NomSource, ::Type{T}) where {T <: Real} = (val = T,)
h_x(::NomSource, (; t)) = (val = 3.0 + t,)

struct FrozenReader <: AbstractComponent end
input_types(::FrozenReader) = (in = Float64,)
output_types(::FrozenReader) = (out = Float64,)
h_su(::FrozenReader, (; u)) = (out = 2.0 * u.in,)

struct ClockStamp <: AbstractComponent end
output_types(::ClockStamp) = (stamp = Float64,)
h_s(::ClockStamp, (; t)) = (stamp = t,)

@testset "a non-nominal activation re-runs Stratum C; frozen products carry (§9.4)" begin
    pair() = Group((; src = NomSource(), rd = FrozenReader());
                   wires = ("src/val" => "rd/in",))

    # The nominal activation runs at build and *is* the Float64 activation; a
    # non-nominal one materializes at first request and is cached on the Build.
    b = build(pair())
    @test activation(b, Float64) === b.nominal
    @test activation(b, D8) === activation(b, D8)

    # The frozen reader's cells hold what the *nominal* probe computed from its
    # real upstream value — 2·(3.0 + 0.0) — not a value synthesized off its
    # declaration. Its cell pins while its producer's walks.
    simd = Simulation(b, D8; h = 1//100)
    @test port(simd, "src", :val) isa D8
    @test port(simd, "rd", :out) === 6.0

    # A discrete stage is never probed at a non-nominal activation (§9.4's
    # executable set): `t` in a discrete bundle is lawful, because it is a
    # `Float64` whenever the stage actually runs — so this must not detonate
    # as a `Dual` arriving at a pinned declaration.
    sims = Simulation(single(ClockStamp()), D8; h = 1//100)
    @test port(sims, "c", :stamp) === 0.0

    # §9.4's opt-in exhaustive mode: the listed activations materialize at
    # build time, which is where CI catches a lurking pinned leaf.
    err = failure(() -> build(single(PinnedGetsDual()); activations = (Float64, D8)))
    @test err isa BuildError && occursin("participates in differentiation", err.msg)
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
    @test port(sim, "c", :a) isa D8
    @test port(sim, "c", :frozen) isa Float64
    @test port(sim, "c", :frozen) == 2.0  # a stored constant, not a computed product
end

# Mixed-leaf cells (§7.2's per-leaf table): the ordinary route, an `Int` leaf
# beside `T` leaves, and the D-166 route, a pinned `Float64` inside a declared
# struct. The cell spans one buffer per leaf eltype; its address carries one
# cursor per eltype as an `NTuple` field (D-162's C2M point).
struct TaggedValue{T}
    v::T
    n::Int
end
struct MixedCell <: AbstractComponent end
output_types(::MixedCell, ::Type{T}) where {T <: Real} = (out = TaggedValue{T},)
h_x(::MixedCell, (; t)) = (out = TaggedValue(t, 1),)

struct PinnedPair{T}
    a::T
    ref::Float64
end
struct PinnedInside <: AbstractComponent end
output_types(::PinnedInside, ::Type{T}) where {T <: Real} = (out = PinnedPair{T},)
h_x(::PinnedInside, (; t)) = (out = PinnedPair(t, 2.0),)

@testset "a mixed-leaf cell lays out across its eltypes' buffers (§7.2, D-162)" begin
    # The Int leaf beside T: mixed at every activation. The tag must come back
    # as a stored `Int`, not a converted double in the `T` buffer.
    sim = Simulation(single(MixedCell()); h = 1//100)
    @test Set(keys(sim.store.stores)) == Set([Symbol(Float64), Symbol(Int)])
    init!(sim)
    out = port(sim, "c", :out)
    @test out isa TaggedValue{Float64} && out.n === 1
    step!(sim, 1e-2)
    @test @ballocated(step!($sim, 1e-2)) == 0

    # The pinned leaf inside a declared struct (D-166): homogeneous at nominal
    # (K = 1), mixed off it — same declaration, and at `Dual` the `T` half
    # walks while `ref` stays a pinned `Float64` in its own buffer.
    sim = Simulation(single(PinnedInside()); h = 1//100)
    @test keys(sim.store.stores) === (Symbol(Float64),)
    sim = Simulation(single(PinnedInside()), D8; h = 1//100)
    @test Set(keys(sim.store.stores)) == Set([Symbol(D8), Symbol(Float64)])
    init!(sim)
    out = port(sim, "c", :out)
    @test out isa PinnedPair{D8}
    @test out.a isa D8 && out.ref === 2.0
end
