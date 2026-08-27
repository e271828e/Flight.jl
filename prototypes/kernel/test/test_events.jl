# --- events at boundaries (§2.1, §10.6, increment 6) ----------------------------

# Malformed event and projection declarations, at top level per D-164.
struct HalfEvent <: AbstractComponent end
init_m(::HalfEvent) = (s = :a,)
output_types(::HalfEvent, ::Type{T}) where {T <: Real} = (o = T,)
h_x(::HalfEvent, (; m)) = (o = 1.0,)
halfevent_guard(::HalfEvent, (; m)) = m.s === :a
function halfevent_handler end                   # a name with no method: the missing half
events(::HalfEvent) = (go = Event(halfevent_guard, halfevent_handler),)

struct NotAnEvent <: AbstractComponent end
init_m(::NotAnEvent) = (s = :a,)
output_types(::NotAnEvent, ::Type{T}) where {T <: Real} = (o = T,)
h_x(::NotAnEvent, (; m)) = (o = 1.0,)
events(::NotAnEvent) = (go = 5,)

struct BadGuardForm <: AbstractComponent end
init_m(::BadGuardForm) = (s = :a,)
output_types(::BadGuardForm, ::Type{T}) where {T <: Real} = (o = T,)
h_x(::BadGuardForm, (; m)) = (o = 1.0,)
badguard_guard(::BadGuardForm, (; m)) = "high"
badguard_handler(::BadGuardForm, (; m)) = (m = (; s = :b),)
events(::BadGuardForm) = (go = Event(badguard_guard, badguard_handler),)

struct BadHandlerKey <: AbstractComponent end    # writes `x`, owns only modes
init_m(::BadHandlerKey) = (s = :a,)
output_types(::BadHandlerKey, ::Type{T}) where {T <: Real} = (o = T,)
h_x(::BadHandlerKey, (; m)) = (o = 1.0,)
badkey_guard(::BadHandlerKey, (; m)) = m.s === :a
badkey_handler(::BadHandlerKey, (; m)) = (x = (q = 1.0,),)
events(::BadHandlerKey) = (go = Event(badkey_guard, badkey_handler),)

struct PartialX <: AbstractComponent end         # an incomplete `x` write-back
init_x(::PartialX) = (a = 0.0, b = 0.0)
output_types(::PartialX, ::Type{T}) where {T <: Real} = (a = T,)
h_x(::PartialX, (; x)) = (a = x.a,)
f(::PartialX, (; x)) = (a = 1.0, b = 1.0)
partialx_guard(::PartialX, (; x)) = x.a ≥ 1.0
partialx_handler(::PartialX, (; x)) = (x = (; a = 0.0),)
events(::PartialX) = (reset = Event(partialx_guard, partialx_handler),)

struct EventsOnDiscrete <: AbstractComponent end # `events` is continuous-only
init_s(::EventsOnDiscrete) = (n = 0,)
output_types(::EventsOnDiscrete) = (n = Int,)
h_s(::EventsOnDiscrete, (; s)) = (n = s.n,)
g(::EventsOnDiscrete, (; s)) = (n = s.n + 1,)
eod_guard(::EventsOnDiscrete, (; s)) = s.n > 0
eod_handler(::EventsOnDiscrete, (; s)) = (;)
events(::EventsOnDiscrete) = (go = Event(eod_guard, eod_handler),)

struct ProjectOnDiscrete <: AbstractComponent end
init_s(::ProjectOnDiscrete) = (n = 0,)
output_types(::ProjectOnDiscrete) = (n = Int,)
h_s(::ProjectOnDiscrete, (; s)) = (n = s.n,)
g(::ProjectOnDiscrete, (; s)) = (n = s.n + 1,)
project(::ProjectOnDiscrete, s) = s

struct ProjectNoState <: AbstractComponent end   # nothing to project onto
output_types(::ProjectNoState, ::Type{T}) where {T <: Real} = (o = T,)
h_x(::ProjectNoState, (; t)) = (o = 1.0,)
project(::ProjectNoState, x) = x

struct BadProjectShape <: AbstractComponent end  # wrong fields back
init_x(::BadProjectShape) = (q = 1.0,)
output_types(::BadProjectShape, ::Type{T}) where {T <: Real} = (q = T,)
h_x(::BadProjectShape, (; x)) = (q = x.q,)
f(::BadProjectShape, (; x)) = (q = 0.0,)
project(::BadProjectShape, x) = (v = x.q,)

@testset "the declaration layer and probe reject malformed events (§8.2, §9.3)" begin
    @test occursin("both halves", failure(() -> build(single(HalfEvent()))).msg)
    @test occursin("Event(guard, handler)", failure(() -> build(single(NotAnEvent()))).msg)

    err = failure(() -> build(single(BadGuardForm())))
    @test err isa BuildError
    @test occursin("Bool", err.msg) && occursin("sign value", err.msg)

    err = failure(() -> build(single(BadHandlerKey())))
    @test occursin("stores it writes", err.msg) && occursin("`m`", err.msg)

    @test occursin("complete", failure(() -> build(single(PartialX()))).msg)

    # `events` is continuous-only, beside `init_m` in the tier-agreement check.
    err = failure(() -> classify_tier("c", EventsOnDiscrete()))
    @test err isa BuildError && occursin("events", err.msg)

    @test occursin("continuous-only", failure(() -> build(single(ProjectOnDiscrete()))).msg)
    @test occursin("init_x", failure(() -> build(single(ProjectNoState()))).msg)
    @test occursin("fields", failure(() -> build(single(BadProjectShape()))).msg)
end

@testset "the guard's form is the declared policy (§2.1, §10.4, D-179)" begin
    @test build(fed(Trigger(0.5), "sig")).policies[1] === (fire = :boundary,)
    @test build(single(Sawtooth(0.3))).policies[1] === (wrap = :localized,)
    @test build(single(Rotor())).policies[1] === NamedTuple()
end

@testset "edge semantics: a u-edge fires at its boundary, a sticky predicate once" begin
    sim = Simulation(fed(Trigger(0.5), "sig"); h = 1//10)
    init!(sim, fragment(inputs = (in = 0.0,)))    # authored at 0.0: not holding
    step!(sim; t_plus = 0.3)
    @test modes(sim, "c").count == 0
    stage!(sim, "in" => 1.0)                     # the input epoch seam, staged (§11.4)
    step!(sim; t_plus = 0.3)                     # drained at the frame top: the first boundary sees the edge
    @test modes(sim, "c") === (state = :fired, count = 1)
    @test port(sim, "c", :on) === true  # the settled sweep published the new mode
    step!(sim; t_plus = 0.9)                     # sticky: holding presents no further edge
    @test modes(sim, "c").count == 1

    # Boundary zero establishes every prior as not-holding (§10.6): a predicate
    # already holding in the authored state fires at t₀ — derived, not asserted.
    sim2 = Simulation(fed(Trigger(0.5), "sig"); h = 1//10)
    init!(sim2, fragment(inputs = (in = 1.0,)))
    @test modes(sim2, "c").count == 1
    # A warm restart resets the registers from scratch: it fires again — and
    # from the declared defaults, so `count` restarts at 0 and reaches 1, never 2
    # (D-063).
    init!(sim2, fragment(inputs = (in = 1.0,)))
    @test modes(sim2, "c").count == 1
end

@testset "a cascade settles within one boundary, independently of h (§10.6)" begin
    chain() = Group((; trig = Trigger(0.5), f1 = Follower(), f2 = Follower());
                    wires = ("trig/on" => "f1/go",
                             "f1/on" => "f2/go"),
                    inputs = ("in" => "trig/sig",))
    for h in (1//10, 1//1000)                    # the latency is rounds, never steps
        sim = Simulation(chain(); h)
        init!(sim, fragment(inputs = (in = 1.0,)))    # one boundary: three rounds to quiescence
        @test modes(sim, "trig").state === :fired
        @test modes(sim, "f1").state === :on
        @test modes(sim, "f2").state === :on
        @test port(sim, "f2", :on) === true
    end
end

@testset "priority with re-decision: a blocked edge stands, a falsified premise re-decides (D-191)" begin
    # Both events present the same edge; `first` fires by declaration order and
    # `second` is eligible but blocked — its sample is not overwritten, so the
    # standing edge fires it in the next round.
    sim = Simulation(fed(TwoShot(), "sig"); h = 1//10)
    init!(sim, fragment(inputs = (in = 2.0,)))
    @test modes(sim, "c") === (a = true, b = true)

    # The premise the first transition falsifies: re-decided against the
    # post-transition sweep, `second` never fires on its stale round-1 edge.
    simp = Simulation(fed(Preempted(), "sig"); h = 1//10)
    init!(simp, fragment(inputs = (in = 2.0,)))
    @test modes(simp, "c") === (a = true, b = false)
end

@testset "an x-writing handler's carried overshoot is resolution-invariant" begin
    # The sawtooth's wrap now genuinely localizes, but its carrying handler
    # (`q ← q − 1`) makes the *boundary states* invariant to where the firing
    # lands, so the boundary recursion stays the exact reference — and the
    # second wrap crosses inside an *off-tick* frame under n = 2, which must
    # fire it all the same.
    q_ref(N) = (q = 0.0; for _ in 1:N; q += 0.03; q ≥ 1 && (q -= 1); end; q)
    for n in (1, 2)
        sim = Simulation(single(Sawtooth(0.3)); h = 1//10, n)
        init!(sim)
        step!(sim; t_plus = 6.7)                 # boundary 67: the off-tick wrap
        @test state(sim, "c").q ≈ q_ref(67) rtol = 1e-9
        step!(sim; t_plus = 0.3)
        @test state(sim, "c").q ≈ q_ref(70) rtol = 1e-9
        @test port(sim, "c", :q) ≈ q_ref(70) rtol = 1e-9   # the post-fire re-sweep published it
    end
end

@testset "due updates run after quiescence, from post-transition values (§10.6)" begin
    # At the wrap boundary the handler resets `q`, the re-sweep publishes the
    # post-wrap value, and only then does the integrator's `g` read it: an
    # update-before-quiescence would accumulate 1.02 where the reference has
    # 0.02.
    m = Group((; saw = Sawtooth(0.3), ctl = DiscreteIntegrator(1.0));
              wires = ("saw/q" => "ctl/e",))
    sim = Simulation(m; h = 1//10)
    init!(sim)
    run!(sim; t_end = 4.0)
    q, acc = 0.0, 0.0
    for _ in 1:40
        q += 0.03
        q ≥ 1 && (q -= 1)
        acc += 0.1 * q
    end
    @test state(sim, "ctl").acc ≈ acc rtol = 1e-9
end

@testset "budget exhaustion degrades, reported on the loop's cell (§10.6, §11.8)" begin
    chatty() = Group((; chat = Chatterer(), trig = Trigger(0.5));
                     inputs = ("in" => "trig/sig",))
    sim = Simulation(chatty(); h = 1//10)
    init!(sim, fragment(inputs = (in = 1.0,)))            # exhaustion at boundary zero
    @test modes(sim, "chat").flips == 8         # 2 × the default budget of 4
    @test modes(sim, "trig").count == 1         # the sibling fired normally
    # The report rides the loop's own diagnostic cell (§11.8), folded at the
    # next frame top: boundary zero's snapshot shows nothing yet.
    @test writer_status(latest(sim), "loop").totals.firing == 0
    step!(sim; t_plus = 0.1)                             # one frame: its snapshot carries the delta
    lw = writer_status(latest(sim), "loop")
    fb = only(lw.recent)
    @test fb isa FiringBudget
    @test fb.path == "chat" && fb.event == :up
    @test fb.budget == 4 && fb.count == 4 && fb.t == 0.0
    @test lw.totals.firing == 1
    step!(sim; t_plus = 0.2)                             # the exhausted boundary's samples
    @test modes(sim, "chat").flips == 8         # became honest priors: quiescent now
    lw = writer_status(latest(sim), "loop")              # quiescent frames: the delta has
    @test isempty(lw.recent) && lw.totals.firing == 1    # passed, the session's totals stand

    # Totals count since the run began (§11.8): the warm restart re-exhausts
    # its own boundary zero — the stores restart from the declared defaults
    # (D-063) and the priors reset with them — and the new run's totals carry
    # that one occurrence, never the last trajectory's.
    init!(sim, fragment(inputs = (in = 1.0,)))
    run!(sim; t_end = 0.1)
    @test modes(sim, "chat").flips == 8
    @test writer_status(latest(sim), "loop").totals.firing == 1

    # The budget is a deployment keyword, validated with its siblings.
    sim2 = Simulation(chatty(); h = 1//10, firing_budget = 2)
    init!(sim2, fragment(inputs = (in = 1.0,)))
    @test modes(sim2, "chat").flips == 4
    run!(sim2; t_end = 0.1)
    fb2 = only(writer_status(latest(sim2), "loop").recent)
    @test fb2 isa FiringBudget && fb2.budget == 2 && fb2.count == 2
    @test occursin("firing_budget",
                   failure(() -> Simulation(chatty(); h = 1//10, firing_budget = 0)).msg)
end

@testset "projection runs at every boundary, between write and decode (§5.3)" begin
    sim = Simulation(single(Rotor(; ω = 2.0, r₀ = SVector(2.0, 0.0))); h = 1//100)
    init!(sim)
    # Boundary zero already projects (§14.5), and the sweep decodes the
    # projected state — a wrong order would publish c = 2.0.
    @test state(sim, "c").r ≈ SVector(1.0, 0.0) atol = 1e-15
    @test port(sim, "c", :c) == 1.0
    run!(sim; t_end = 1.0)
    r = state(sim, "c").r
    @test r[1]^2 + r[2]^2 ≈ 1.0 atol = 1e-14             # pinned to the manifold
    @test r ≈ SVector(cos(2.0), sin(2.0)) rtol = 1e-6    # and still the right trajectory
end

@testset "events compile out at a non-nominal activation (§9.4, D-052)" begin
    sim = Simulation(fed(Trigger(0.5), "sig"), D8; h = 1//10)
    @test isempty(sim.events.entries)
    init!(sim, fragment(inputs = (in = D8(1.0),)))
    run!(sim; t_end = 0.3)
    @test modes(sim, "c").count == 0            # the guard never ran

    # Projection is continuous machinery, inside every executable set.
    simr = Simulation(single(Rotor(; r₀ = SVector(2.0, 0.0))), D8; h = 1//100)
    init!(simr)
    @test ForwardDiff.value(state(simr, "c").r[1]) ≈ 1.0 atol = 1e-15
end

@testset "gate 3: a quiet boundary with events does not allocate (§7.5)" begin
    # The first iteration round always runs, so guard evaluation is on the
    # measured path even when nothing ever fires.
    sim = Simulation(fed(Trigger(0.5), "sig"); h = 1//10)
    init!(sim, fragment(inputs = (in = 0.0,)))
    boundary!(sim, 1); offtick_boundary!(sim)
    @test @ballocated(boundary!($sim, 1)) == 0
    @test @ballocated(offtick_boundary!($sim)) == 0
    @test @ballocated(step!($sim, 0.1)) == 0

    simr = Simulation(single(Rotor()); h = 1//100)       # projection on the measured path
    init!(simr)
    boundary!(simr, 1)
    @test @ballocated(boundary!($simr, 1)) == 0
end
