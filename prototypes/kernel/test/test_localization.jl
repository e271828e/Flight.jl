# --- localization (§10.4, increment 7) ------------------------------------------
# The input epochs these tests rely on are the real ones since increment 9: a
# batch staged between runs is applied by the drain at the next frame top
# (§11.4) — so a slot write is always epoch-aligned here, by construction.

@testset "a localized event fires at t*, within tol of the true crossing (§10.4)" begin
    # Linear trajectory: RK4 and the cubic Hermite are both exact, so the stamp
    # isolates the localization machinery itself — t* within the bracket width.
    m = Group((; src = Sawtooth(1.0), s = Stamper(0.315)),
              ("children/src/q" => "children/s/sig",), (), ())
    sim = Simulation(m; h = 1//10)
    init!(sim)
    run!(sim, 0.5)
    @test modes(sim, "children/s").count == 1
    @test modes(sim, "children/s").t_fired ≈ 0.315 atol = 1e-6
    @test modes(sim, "children/s").t_fired ≥ 0.315          # the holding endpoint
    @test port(sim, "children/s", :armed) === false         # the t* re-sweep published it
    run!(sim, 0.9)                                          # sticky: no further edge
    @test modes(sim, "children/s").count == 1
end

@testset "the timing observable: a discarding reset needs t*, not the boundary" begin
    # Bouncer's handler discards the overshoot, so the trajectory itself moves
    # with the firing instant: localized resets give the exact period
    # level/rate at an h-incommensurate level, where boundary-resolution
    # firings would have accumulated a full overshoot per reset.
    sim = Simulation(single(Bouncer(1.0, 0.315)); h = 1//10)
    init!(sim)
    run!(sim, 2.0)                                          # resets at 0.315·k, k = 1..6
    @test modes(sim, "children/c").count == 6
    @test state(sim, "children/c").q ≈ 2.0 - 6 * 0.315 atol = 1e-5
end

@testset "t* on a nonlinear trajectory is interpolant-accurate, O(h⁴) (§10.4)" begin
    # c(t) = -cos t under the renormalizing project: σ = -cos t + 0.5 crosses
    # at t = π/3, mid-frame. Trials run against the raw interpolated state;
    # the error budget is the Hermite's O(h⁴) plus RK4's own global error.
    m = Group((; src = Rotor(; ω = 1.0, r₀ = SVector(-1.0, 0.0)), s = Stamper(-0.5)),
              ("children/src/c" => "children/s/sig",), (), ())
    sim = Simulation(m; h = 1//10)
    init!(sim)
    run!(sim, 1.5)
    @test modes(sim, "children/s").count == 1
    @test modes(sim, "children/s").t_fired ≈ π / 3 atol = 1e-5
end

@testset "the θ = 0 validation: an epoch-caused edge fires at the frame top, exactly" begin
    # The slot write between runs is the epoch seam: the guard jumps at the
    # frame top without crossing anything, σ₀ re-measured under the frame's own
    # u holds, and not localizing is the action — the event falls through to
    # the frame top's ordinary iteration and stamps the indexed grid time.
    sim = Simulation(fed(Stamper(0.5), "sig"); h = 1//10)
    init!(sim)
    run!(sim, 0.3)
    @test modes(sim, "children/c").count == 0
    stage!(sim, "in" => 1.0)                     # staged, drained at the next frame top (§11.4)
    run!(sim, 0.6)
    @test modes(sim, "children/c").count == 1
    @test modes(sim, "children/c").t_fired == 4 * sim.h     # the grid point itself
end

@testset "t* = tₙ₊₁ exactly degenerates to the grid boundary (§10.4)" begin
    # The ramp reads t, so σ crosses exactly at the frame top: the trigger sees
    # σ₁ = 0, every interior trial is not-holding, and the localization result
    # is discarded — one boundary, one firing, stamping the indexed grid time.
    m = Group((; src = Ramp(0.0), s = Stamper(0.4)),
              ("children/src/out" => "children/s/sig",), (), ())
    sim = Simulation(m; h = 1//10)
    init!(sim)
    run!(sim, 0.6)
    @test modes(sim, "children/s").count == 1
    @test modes(sim, "children/s").t_fired == 4 * sim.h
end

@testset "multiple crossings in one frame: earliest first, re-localized on the remainder" begin
    m2() = Group((; src = Sawtooth(1.0), s1 = Stamper(0.31), s2 = Stamper(0.34)),
                 ("children/src/q" => "children/s1/sig",
                  "children/src/q" => "children/s2/sig"), (), ())
    sim = Simulation(m2(); h = 1//10)
    init!(sim)
    @test_logs run!(sim, 0.5)                               # both localize: no degradation
    @test modes(sim, "children/s1").t_fired ≈ 0.31 atol = 1e-6
    @test modes(sim, "children/s2").t_fired ≈ 0.34 atol = 1e-6
    @test modes(sim, "children/s1").t_fired < modes(sim, "children/s2").t_fired

    # A tie is one localization: both edges stand at the shared t* and fire
    # together inside that boundary's iteration — so a budget of 1 suffices,
    # and no ChatteringBudget degradation is seen.
    mt = Group((; src = Sawtooth(1.0), s1 = Stamper(0.315), s2 = Stamper(0.315)),
               ("children/src/q" => "children/s1/sig",
                "children/src/q" => "children/s2/sig"), (), ())
    simt = Simulation(mt; h = 1//10, localization_budget = 1)
    init!(simt)
    @test_logs run!(simt, 0.5)
    @test modes(simt, "children/s1").t_fired == modes(simt, "children/s2").t_fired
    @test modes(simt, "children/s1").t_fired ≈ 0.315 atol = 1e-6

    # Distinct crossings under budget 1: the earliest localizes, the second
    # spends nothing — it degrades to boundary granularity under the warning,
    # stamping the frame top.
    simb = Simulation(m2(); h = 1//10, localization_budget = 1)
    init!(simb)
    @test_logs (:warn, r"ChatteringBudget: event `cross` at `children/s2`") run!(simb, 0.5)
    @test modes(simb, "children/s1").t_fired ≈ 0.31 atol = 1e-6
    @test modes(simb, "children/s2").t_fired == 4 * simb.h
end

@testset "budget exhaustion degrades to boundary granularity and warns (§10.4)" begin
    # The relaxer re-arms 1 ms below its level, so each remainder re-crosses
    # within the frame: 8 localizations at 50..57 ms spend the default budget,
    # the 9th trigger warns — naming the event and the count — and fires in the
    # frame top's ordinary iteration. The remainder had already completed, so
    # the run proceeds, at boundary granularity for that frame alone.
    sim = Simulation(single(Relaxer(1.0, 0.05, 0.001)); h = 1//10)
    init!(sim)
    @test_logs (:warn, r"ChatteringBudget: event `pop` at `children/c` has localized 8 times") run!(sim, 0.1)
    @test modes(sim, "children/c").count == 9               # 8 at t*, 1 at the frame top
    @test state(sim, "children/c").q ≈ 0.049 atol = 1e-12   # the frame-top firing re-armed it
end

@testset "the gate idiom localizes; a gate flip is an epoch edge (§10.4)" begin
    gated() = Group((; src = Sawtooth(1.0), s = GatedStamper(0.315)),
                    ("children/src/q" => "children/s/sig",),
                    ("gate" => "children/s/gate",), ())
    b = build(gated())
    @test b.policies[index_of(b.flat, "children/s")] === (cross = :localized,)

    # Gate true from the start: the Bool factor is constant over the bracket
    # and the continuous atom localizes as such.
    sim = Simulation(gated(); h = 1//10)
    set_slot!(sim, "gate", true)
    init!(sim)
    run!(sim, 0.5)
    @test modes(sim, "children/s").t_fired ≈ 0.315 atol = 1e-6

    # Gate flipped at a frame top, with σ already past the level: the edge is
    # the u seam's, σ₀ holds under the frame's own u, and the event fires at
    # the frame top exactly — epoch-caused, never root-found.
    sim2 = Simulation(gated(); h = 1//10)
    init!(sim2)
    run!(sim2, 0.5)
    @test modes(sim2, "children/s").count == 0              # gate down: -one(σ) throughout
    stage!(sim2, "gate" => true)                            # the u seam, through the drain (§11.4)
    run!(sim2, 0.8)
    @test modes(sim2, "children/s").count == 1
    @test modes(sim2, "children/s").t_fired == 6 * sim2.h
end

@testset "ticks are never due at t*: the discrete tier holds through it (§10.4)" begin
    # A sampled consumer beside a localized event: the t* boundary runs the
    # full event phase but no g update — a spurious tick there would add the
    # mid-frame sample 0.1·q(t*) to the accumulator.
    m = Group((; src = Sawtooth(1.0), s = Stamper(0.315), ctl = DiscreteIntegrator(1.0)),
              ("children/src/q" => "children/s/sig",
               "children/src/q" => "children/ctl/e"), (), ())
    sim = Simulation(m; h = 1//10)
    init!(sim)
    run!(sim, 0.5)
    @test modes(sim, "children/s").count == 1               # the event did localize
    @test state(sim, "children/ctl").acc ≈ 0.1 * (0.0 + 0.1 + 0.2 + 0.3 + 0.4 + 0.5) rtol = 1e-9
end

@testset "the two deployment keywords are validated with their siblings (§10.4)" begin
    m = single(Bouncer(1.0, 0.315))
    @test occursin("localization_tol",
                   failure(() -> Simulation(m; h = 1//10, localization_tol = 0.0)).msg)
    @test occursin("localization_tol",
                   failure(() -> Simulation(m; h = 1//10, localization_tol = -1e-3)).msg)
    @test occursin("localization_budget",
                   failure(() -> Simulation(m; h = 1//10, localization_budget = 0)).msg)
end

@testset "localization compiles out at a non-nominal activation (§9.4, D-052)" begin
    # Events are outside every non-nominal executable set, so the frame loop's
    # fast-path key is off: the bare step, no arrival machinery, no reset.
    sim = Simulation(single(Bouncer(1.0, 0.315)), D8; h = 1//10)
    @test !sim.has_localized
    init!(sim)
    run!(sim, 0.5)
    @test ForwardDiff.value(state(sim, "children/c").q) ≈ 0.5 rtol = 1e-12
end

@testset "gate 3: localized frames do not allocate (§7.5)" begin
    # A quiet frame pays the arrival sweep and the trigger scan, nothing else.
    mq = Group((; src = Sawtooth(0.1), s = Stamper(100.0)),
               ("children/src/q" => "children/s/sig",), (), ())
    simq = Simulation(mq; h = 1//10)
    init!(simq)
    run!(simq, 0.2)
    @test @ballocated(frame!($simq, 3)) == 0

    # A localizing frame: one crossing, θ = 0 validation, ẋₙ₊₁, the bracketing
    # trials, the t* boundary and the remainder — all against preallocated
    # buffers, re-run from init! each sample.
    siml = Simulation(single(Bouncer(1.0, 0.07)); h = 1//10)
    init!(siml)
    @test @ballocated(frame!($siml, 1), setup = (init!($siml)), evals = 1) == 0
end
