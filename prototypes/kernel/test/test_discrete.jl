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
    init!(sim, fragment(inputs = (ref = r,)))
    run!(sim; t_end = N * Δt)

    # The tolerance is RK4's, not the semantics': the reference integrates
    # exactly where the loop takes one RK4 step per sample, which at `ωΔt` this
    # size leaves ~1e-8 relative. A violated hold moves `u` mid-interval and
    # costs ~1e-3, so the margin still separates the two by orders of magnitude.
    @test state(sim, "plant").q ≈ q rtol = 1e-6
    # The cell holds `s[N]` — the output stage ran at this boundary — while the
    # store already holds `s[N+1]`. That gap *is* the sampled-data ordering:
    # output stages before updates, within one boundary (§10.5).
    @test port(sim, "ctl", :u) ≈ s rtol = 1e-6
    @test state(sim, "ctl").acc ≈ s_next rtol = 1e-6
end

@testset "the ZOH holds by compile-time absence (§10.5)" begin
    sim = Simulation(sampled_loop(); h = 1//50)
    init!(sim, fragment(inputs = (ref = 1.0,)))    # excite the loop, or nothing moves at all

    # Structural: the interior variants carry continuous entries only, so there
    # is no gating test on the hot path — the hold is not implemented, it is
    # the absence of any way to change the cell.
    @test length(walked(sim.exec.bodies.sweep_1, :interior)) == 1        # plant only
    @test length(walked(sim.exec.bodies.sweep_1)) == 2                   # plus ctl
    @test isempty(walked(sim.exec.bodies.ticks, :interior))
    @test length(walked(sim.exec.bodies.ticks)) == 1

    # Semantic: a step is made of interior evaluations, so the discrete cell
    # cannot move across one, while the continuous table does. Run a few
    # boundaries first — the loop starts at rest with `u` held at zero, so the
    # very first interval moves nothing at all.
    run!(sim; t_end = 0.1)
    u₀, y₀ = port(sim, "ctl", :u), port(sim, "plant", :y)
    @test u₀ != 0.0
    step!(sim, 0.02)
    @test port(sim, "ctl", :u) == u₀   # untouched: never gathered, never written
    @test port(sim, "plant", :y) != y₀
end

@testset "discrete state and modes live outside the buffer (§7.3)" begin
    sim = Simulation(Group((; counter = TickCounter(), moded = ModedSource()));
                     h = 1//10)
    # The flat buffer is continuous state only; the counter's `Int` is in its
    # own store, and no store mirrors another.
    @test isempty(sim.exec.xbuf)
    @test state(sim, "counter") === (n = 0,)
    @test modes(sim, "moded") === (phase = :idle,)

    # `Int` and `Bool` cells force their own buffers — the plural in
    # "per-eltype stores", first exercised here.
    @test Set(keys(sim.exec.store.stores)) == Set([Symbol(Int), Symbol(Bool), Symbol(Float64)])

    init!(sim)
    @test state(sim, "counter") === (n = 1,)   # boundary zero is a tick
    @test port(sim, "counter", :n) === 0       # the cell holds what it published
    @test port(sim, "counter", :even) === true
    run!(sim; t_end = 0.5)
    # Six boundaries: zero, then one per step. The store leads the cell by one,
    # the update having run after the output stage at each of them.
    @test state(sim, "counter") === (n = 6,)
    @test port(sim, "counter", :n) === 5
end

@testset "the workspace is scratch, on both tiers (§7.3)" begin
    sim = Simulation(Group((; sm = Smoother(0.5), src = ModedSource(), wg = WorkGain(2.0));
                           wires = ("src/out" => "sm/a",
                                    "src/out" => "sm/b",
                                    "src/out" => "wg/in"));
                     h = 1//10)
    init!(sim)
    @test port(sim, "wg", :out) == 0.0      # 2 × the idle-phase constant
    run!(sim; t_end = 0.3)
    @test state(sim, "sm").v isa SVector{2,Float64}

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
init_s(::WidenedUpdate) = (n = 0,)
output_types(::WidenedUpdate) = (n = Int,)
h_s(::WidenedUpdate, (; s)) = (n = s.n,)
g(::WidenedUpdate, (; s)) = (n = s.n + 0.5,)

@testset "a discrete successor is the store's own type (§7.3)" begin
    @test_throws BuildError build(single(WidenedUpdate()))
end

@testset "the discrete tier is frozen at a non-nominal activation (§7.2)" begin
    # The plant's input is declared `T` and wired to a discrete `Float64` cell:
    # a lawful arrival, embedded as a zero-partial. `ctl`'s stages are outside
    # this activation's executable set, so they do not run at all — boundary
    # zero's wide gate included (D-205 admits entries, and a frozen component
    # has none) — and its cell holds the nominal products §9.4 carried across,
    # pinned for the whole run.
    simd = Simulation(sampled_loop(), D8; h = 1//50)
    @test isempty(walked(simd.exec.bodies.ticks))
    @test length(walked(simd.exec.bodies.sweep_1)) == 1          # plant only; ctl frozen
    @test port(simd, "ctl", :u) isa Float64

    init!(simd, fragment(inputs = (ref = 0.0,)))
    run!(simd; t_end = 0.04)
    @test state(simd, "plant").q isa SVector{2,D8}
    @test port(simd, "ctl", :u) == 0.0              # held, never recomputed
end
