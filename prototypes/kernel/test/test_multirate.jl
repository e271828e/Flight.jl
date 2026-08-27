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
    @test isempty(walked(sim.bodies.sweep_2, :interior))
    @test gated(sim.bodies.sweep_2) == 3
    @test gated(sim.bodies.sweep_1) == 0            # the ramp is continuous
end

@testset "the hyperperiod chart is readable off the cells (§10.5)" begin
    # The ramp makes every sample carry its acquisition time (1 + t), so each
    # cell says when its owner last ticked — the chart's dots, and the
    # deterministic aging of the stagger, in the spec's own numbers.
    sim = Simulation(MultiRate(); h = 1//500)
    init!(sim)                                       # boundary zero: Φ = 0 is due
    @test port(sim, "fcs/inner", :out) == 1.0
    @test port(sim, "gnss", :out) == 1.0

    step!(sim; frames = 2)                           # base tick 2: outer's first tick
    @test port(sim, "fcs/inner", :out) ≈ 1.004       # fresh at every tick
    @test port(sim, "fcs/outer", :out) == 1.0        # gnss's k = 0 sample: two ticks old
    step!(sim; frames = 5)                           # base tick 7
    @test port(sim, "fcs/outer", :out) == 1.0        # the same sample, seven ticks old
    @test port(sim, "gnss", :out) == 1.0             # gnss itself holds until k = 10
    step!(sim; frames = 3)                           # base tick 10
    @test port(sim, "gnss", :out) ≈ 1.02             # its second tick
    run!(sim; t_end = 12 * 0.002)
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
    run!(sim; t_end = 10 * 0.01)
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
    init!(s1); run!(s1; t_end = 0.02)
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
    step!(sim)                                       # base tick 1: the first tick
    @test port(sim, "", :y) ≈ 5.01
    step!(sim)                                       # base tick 2: not due, holds
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
    run!(sim; t_end = N * Δt_ctl)
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
