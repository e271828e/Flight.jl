# --- the stepper seam (§10.2, increment 8) --------------------------------------
# Interchangeability is the claim under test: the loop, the localization
# machinery and the deployment surface are method-blind, and the two backends
# differ exactly where their orders say they must.

@testset "the method is a deployment binding, RK4 the default (§10.2)" begin
    sim = Simulation(feedback_model(); h = 1//100)
    @test sim.stepper isa RK4{Float64}
    simh = Simulation(feedback_model(); h = 1//100, method = Heun)
    @test simh.stepper isa Heun{Float64}
    # validated with its siblings: a backend is named by stepper type, and
    # anything else is refused at binding, not deep in a MethodError
    @test occursin("method",
                   failure(() -> Simulation(feedback_model(); h = 1//100, method = 4)).msg)
    @test occursin("method",
                   failure(() -> Simulation(feedback_model(); h = 1//100, method = Int)).msg)
end

@testset "convergence order: each backend is itself — 4 and 2 (§10.2)" begin
    # The sharp form of the interchangeability claim: same model, same loop,
    # same reference, and the fitted order is the method's own. A mislabeled
    # backend cannot sit in both bands, where a shared loose tolerance would
    # let it hide.
    ω, ζ, k, r = 2.0, 0.1, 4.0, 0.7
    A = SMatrix{2,2}(0.0, -ω^2, 1.0, -2ζ * ω)
    B = SVector(0.0, 1.0)
    Acl = A - B * k * SVector(1.0, 0.0)'
    exact(t) = exp(Acl * t) * (Acl \ (B * k * r)) - Acl \ (B * k * r)
    function final_error(method, h)
        sim = Simulation(feedback_model(; k, ω, ζ); h, method)
        init!(sim, fragment(slots = (ref = r,)))
        run!(sim; t_end = 2.0)
        norm(state(sim, "children/plant").q - exact(2.0))
    end
    hs = (1//10, 1//20, 1//40, 1//80)   # errors 5e-8..2e-4: well above float noise
    for (method, p) in ((RK4, 4.0), (Heun, 2.0))
        errs = [final_error(method, h) for h in hs]
        orders = [log2(errs[i] / errs[i+1]) for i in 1:3]
        @test all(o -> abs(o - p) < 0.3, orders)
    end
end

@testset "localization is seam-generic: t* under Heun (§10.2, §10.4)" begin
    # Linear trajectory: Heun and the Hermite are both exact, so the stamp is
    # method-independent down to the bracket width — the machinery, not the
    # method, sets the error.
    m = Group((; src = Sawtooth(1.0), s = Stamper(0.315)),
              ("children/src/q" => "children/s/sig",), (), ())
    sim = Simulation(m; h = 1//10, method = Heun)
    init!(sim)
    run!(sim; t_end = 0.5)
    @test modes(sim, "children/s").count == 1
    @test modes(sim, "children/s").t_fired ≈ 0.315 atol = 1e-6

    # Nonlinear trajectory: the stamp error is now the discrete solution's
    # O(h²), not RK4's O(h⁴) — quartering under h-halving is the method
    # showing through the same machinery. (Boundary-resolution firing would
    # sit at ~h and shrink linearly instead.)
    stamp_error(h) = begin
        mr = Group((; src = Rotor(; ω = 1.0, r₀ = SVector(-1.0, 0.0)), s = Stamper(-0.5)),
                   ("children/src/c" => "children/s/sig",), (), ())
        simr = Simulation(mr; h, method = Heun)
        init!(simr)
        run!(simr; t_end = 1.5)
        @test modes(simr, "children/s").count == 1
        abs(modes(simr, "children/s").t_fired - π / 3)
    end
    e10, e40 = stamp_error(1//10), stamp_error(1//40)
    @test e10 < 5e-3
    @test 10 < e10 / e40 < 24                               # h²: ×16 over two halvings

    # The frame-top claims never depended on the method: an epoch-caused edge
    # falls through to fire at the indexed grid time bitwise under Heun too.
    simf = Simulation(fed(Stamper(0.5), "sig"); h = 1//10, method = Heun)
    init!(simf, fragment(slots = (in = 0.0,)))
    step!(simf; t_plus = 0.3)
    stage!(simf, "in" => 1.0)                   # frame 4's drain, at its frame top
    step!(simf; t_plus = 0.3)
    @test modes(simf, "children/c").t_fired == 4 * simf.h
end

@testset "gate 4: the second backend holds the §7.5 invariant" begin
    sim = Simulation(feedback_model(); h = 1//1000, method = Heun)
    init!(sim, fragment(slots = (ref = 0.0,)))
    step!(sim, 1e-3)
    @test @ballocated(step!($sim, 1e-3)) == 0
    # The localizing frame allocates exactly its t* boundary's publication —
    # the framework-side carve-out (§7.5, §11.2) — as under RK4 (gate 3).
    siml = Simulation(single(Bouncer(1.0, 0.07)); h = 1//10, method = Heun, log = false)
    init!(siml)
    pub = @ballocated publish!($siml)
    @test @ballocated(frame!($siml, 1), setup = (init!($siml)), evals = 1) == pub
end

@testset "the second backend is generic over the scalar (§7.2)" begin
    sim = Simulation(feedback_model(), D8; h = 1//1000, method = Heun)
    init!(sim, fragment(slots = (ref = D8(0.7),)))
    run!(sim; t_end = 0.05)
    @test state(sim, "children/plant").q isa SVector{2,D8}
end
