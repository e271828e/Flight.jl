# --- the trim service (§14.7, §14.8; increment 21) --------------------------------
# The problem value, the collected setup diagnostic, the two-half scratch world
# of D-213, the in-house Levenberg–Marquardt behind the backend seam, the
# service's own convergence verdict, and the commit that is literally an
# `init!`. The fixtures live at top level for the README's local-scope reason.

const PEND_G_L = 9.81      # the pendulum's g/l and damping throughout, so the
const PEND_C = 0.5         # analytic solutions below read as what they are

# The smallest model a trim problem can be posed on: one nonlinear continuous
# component whose torque input is the root's (`fed`, test/utils.jl — root input
# `in`, child `c`). The baseline covers that one root input, which is what
# §14.6's totality asks of it.
pend_base() = fragment(inputs = (in = 0.0,))

# The two problem shapes every case below is a variation on: decide the torque
# at a held attitude (linear, one step), or decide the attitude at a held torque
# (nonlinear, iterative). Both read one derivative and balance one equation.
torque_reads() = reads(ω̇ = get_deriv("c", :ω))
torque_only(r, d) = (torque = r.ω̇,)

decide_u(d) = combine(at("c", condition(Pendulum(); θ = 0.5)),
                      fragment(inputs = (in = d.u,)))
decide_θ(d) = combine(at("c", fragment(x = (θ = d.θ, ω = 0.0))),
                      fragment(inputs = (in = 4.0,)))

u_problem() = TrimProblem(guess = (u = 0.0,), lower = (u = -Inf,), upper = (u = Inf,),
                          condition = decide_u, reads = torque_reads(),
                          residuals = torque_only, tolerances = (torque = 1e-9,))

θ_problem() = TrimProblem(guess = (θ = 0.1,), lower = (θ = -π/2,), upper = (θ = π/2,),
                          condition = decide_θ, reads = torque_reads(),
                          residuals = torque_only, tolerances = (torque = 1e-9,))

# The same attitude problem with the torque left to the baseline, for the models
# whose input is not the decision's to write.
decide_θ_alone(d) = at("c", fragment(x = (θ = d.θ, ω = 0.0)))
θ_alone_problem() = TrimProblem(guess = (θ = 0.1,), lower = (θ = -π/2,), upper = (θ = π/2,),
                                condition = decide_θ_alone, reads = torque_reads(),
                                residuals = torque_only, tolerances = (torque = 1e-9,))

# The two-decision shape, for the seam's packing: `torque` balances the moment
# and `hold` pins the attitude, so the residual system is genuinely square and
# both decisions are load-bearing.
decide_both(d) = combine(at("c", fragment(x = (θ = d.θ, ω = 0.0))),
                         fragment(inputs = (in = d.u,)))
both_reads() = reads(ω̇ = get_deriv("c", :ω), θ = get_state("c", :θ))
both_residuals(r, d) = (torque = r.ω̇, hold = r.θ - 0.3)

# D-213's fixture: the pendulum's torque arrives from a *discrete* producer's
# held output, and that producer's `s` is the baseline's. At the seeded
# activation the discrete tier is frozen, so its output cell can only come from
# the nominal half — which is exactly what the ruling is about.
sampled_pend() = Group((; ctl = DiscreteIntegrator(1.0), c = Pendulum());
                       wires = ("ctl/u" => "c/u",), inputs = ("in" => "ctl/e",))
sampled_base(acc = 4.0) = combine(at("ctl", fragment(s = (acc = acc,))),
                                  fragment(inputs = (in = 0.0,)))

# A commit-time mover of the first kind (§14.5): a guard the solved attitude
# already holds, so boundary zero fires it and the report says so.
triggered_pend() = Group((; c = Pendulum(), trig = Trigger(0.3));
                         wires = ("c/θ" => "trig/sig",), inputs = ("in" => "c/u",))

# A commit-time mover of the second kind: a handler that writes `x`, so the
# committed stores land somewhere the solved-point residuals never described.
# A component of its own, because an event is declared per component type.
struct Snapback <: AbstractComponent
    level::Float64
end

init_x(::Snapback) = (θ = 0.0, ω = 0.0)
input_types(::Snapback, ::Type{T}) where {T <: Real} = (u = T,)
output_types(::Snapback, ::Type{T}) where {T <: Real} = (θ = T, ω = T)
h_x(::Snapback, (; x)) = (θ = x.θ, ω = x.ω)
f(::Snapback, (; x, u)) = (θ = x.ω, ω = -PEND_G_L * sin(x.θ) - PEND_C * x.ω + u.u)
snapback_guard(c::Snapback, (; x)) = x.θ > c.level
snapback_handler(::Snapback, (; x)) = (x = (θ = 0.0, ω = x.ω),)
events(::Snapback) = (snap = Event(snapback_guard, snapback_handler),)
condition(::Snapback; θ = 0.0, ω = 0.0) = fragment(x = (θ = θ, ω = ω))

snap_decide_u(d) = combine(at("c", condition(Snapback(0.3); θ = 0.5)),
                           fragment(inputs = (in = d.u,)))

# §14.7's residual return has two observation points, because a lambda may
# answer differently at each: the nominal guess evaluation runs at `Float64`
# and the first seeded one at `Dual`, and this one branches on the scalar it is
# handed — user machinery no shape check can anticipate.
eltype_split(r, d) = r.ω̇ isa Float64 ? (torque = r.ω̇,) : (wrong = r.ω̇,)

@testset "a one-step linear problem solves, commits and reads back (§14.7, §14.8)" begin
    sim = Simulation(fed(Pendulum(), :u); h = 1//10)
    report = trim!(sim, u_problem(); baseline = pend_base(), t₀ = 0.25)

    # The equilibrium torque at θ = 0.5 with ω = 0, which is the whole problem.
    @test report.converged
    @test report.solution.u ≈ PEND_G_L * sin(0.5)
    @test abs(report.residuals.torque) ≤ report.tolerances.torque
    @test report.status === :converged
    @test report.nevals ≤ 10 && report.niters ≤ 10   # exact Jacobians, quadratic (§14.7)
    @test isempty(report.saturated) && isempty(report.fired_events)

    # The commit is an `init!` in every respect: the lifecycle, the anchor, the
    # authored state and the root input the solver decided.
    @test lifecycle(sim) === :initialized
    @test sim.exec.clock.t === 0.25 && sim.exec.clock.t₀ === 0.25
    @test state(sim, "c") === (θ = 0.5, ω = 0.0)
    @test port(sim, "", :in) === report.solution.u

    # And the committed-state residuals are the same numbers here, boundary zero
    # having moved nothing: no projection to run, no guard to fire.
    @test report.committed_residuals !== nothing
    @test abs(report.committed_residuals.torque) ≤ report.tolerances.torque
end

@testset "an iterative nonlinear problem converges, the box picking the branch (§14.7)" begin
    sim = Simulation(fed(Pendulum(), :u); h = 1//10)
    report = trim!(sim, θ_problem(); baseline = pend_base())

    # −(g/l)·sin θ + u = 0 has two branches in [0, π]; the box admits one.
    @test report.converged && report.status === :converged
    @test report.solution.θ ≈ asin(4.0 / PEND_G_L)
    @test !(report.solution.θ ≈ π - asin(4.0 / PEND_G_L))
    @test report.niters ≤ 10
    # The fragment mixes a `Dual` leaf and a held `Float64` one in one payload,
    # which is §14.3's two converter arms meeting inside a single write.
    @test state(sim, "c").ω === 0.0
    @test lifecycle(sim) === :initialized
end

@testset "permuted bounds are a non-event: names pair, order never does (§14.7, §9.5)" begin
    declared = TrimProblem(guess = (θ = 0.1, u = 0.0), lower = (θ = -π/2, u = -Inf),
                           upper = (θ = π/2, u = Inf), condition = decide_both,
                           reads = both_reads(), residuals = both_residuals,
                           tolerances = (torque = 1e-9, hold = 1e-9))
    # Both bound spellings permuted, and the tolerances permuted with them: the
    # canonicalization is by name at both ends of the seam.
    permuted = TrimProblem(guess = (θ = 0.1, u = 0.0), lower = (u = -Inf, θ = -π/2),
                           upper = (u = Inf, θ = π/2), condition = decide_both,
                           reads = both_reads(), residuals = both_residuals,
                           tolerances = (hold = 1e-9, torque = 1e-9))
    a = trim!(Simulation(fed(Pendulum(), :u); h = 1//10), declared; baseline = pend_base())
    b = trim!(Simulation(fed(Pendulum(), :u); h = 1//10), permuted; baseline = pend_base())

    @test a.converged && b.converged
    @test a.solution == b.solution                   # bit for bit, not merely ≈
    @test a.solution.θ ≈ 0.3 && a.solution.u ≈ PEND_G_L * sin(0.3)
    @test keys(a.residuals) === (:torque, :hold) && keys(b.residuals) === (:hold, :torque)
end

@testset "no convergence, no commit: the simulation is untouched (§14.8)" begin
    # −(g/l)·sin θ + 2·(g/l) ≥ g/l > 0 for every θ: the balance is impossible,
    # and hitting the infeasible edge is information, not an exception.
    infeasible = TrimProblem(
        guess = (θ = 0.1,), lower = (θ = -π/2,), upper = (θ = π/2,),
        condition = d -> combine(at("c", fragment(x = (θ = d.θ, ω = 0.0))),
                                 fragment(inputs = (in = 2 * PEND_G_L,))),
        reads = torque_reads(), residuals = torque_only, tolerances = (torque = 1e-9,))

    # On a never-initialized simulation: it stays `built`, and `run!` still says so.
    fresh = Simulation(fed(Pendulum(), :u); h = 1//10)
    report = trim!(fresh, infeasible; baseline = pend_base())
    @test !report.converged
    @test report.status isa Symbol                   # recorded verbatim, decisive of nothing
    @test report.committed_residuals === nothing     # the absence of a commit, not a flag
    @test abs(report.residuals.torque) > report.tolerances.torque
    @test lifecycle(fresh) === :built
    @test failure(() -> run!(fresh; t_end = 0.1)) isa BuildError

    # On an initialized one: every buffer equals its pre-call copy.
    live = Simulation(fed(Pendulum(), :u); h = 1//10)
    init!(live, combine(at("c", condition(Pendulum(); θ = 0.2)),
                        fragment(inputs = (in = 1.0,))))
    before = world(live)
    r2 = trim!(live, infeasible; baseline = pend_base())
    @test !r2.converged && r2.committed_residuals === nothing
    @test world(live) == before
    @test lifecycle(live) === :initialized
end

@testset "a converged solve names the decisions sitting at a bound (§14.8, D-070)" begin
    # The actuator saturates below the equilibrium torque, and the force balance
    # is declared to a tolerance the saturated point still meets: converged, at
    # a bound, with the bound named — the CG-limit diagnostic (§14.8).
    saturating = TrimProblem(guess = (u = 0.0,), lower = (u = -Inf,), upper = (u = 4.0,),
                             condition = decide_u, reads = torque_reads(),
                             residuals = torque_only, tolerances = (torque = 1.0,))
    sim = Simulation(fed(Pendulum(), :u); h = 1//10)
    report = trim!(sim, saturating; baseline = pend_base())

    @test report.converged
    @test report.solution.u === 4.0                  # the projection assigns the bound itself
    @test report.saturated == [(:u, :upper)]
    @test abs(report.residuals.torque) ≤ report.tolerances.torque
    @test port(sim, "", :in) === 4.0                 # and it committed, at the bound

    # An unsaturated solve of the same shape names nothing.
    unbounded = trim!(Simulation(fed(Pendulum(), :u); h = 1//10), u_problem();
                      baseline = pend_base())
    @test unbounded.converged && isempty(unbounded.saturated)
end

@testset "the empty problem is the equilibrium probe, and the solver is bypassed (§14.8)" begin
    probe(u) = TrimProblem(
        guess = (;), lower = (;), upper = (;),
        condition = d -> combine(at("c", condition(Pendulum(); θ = 0.0)),
                                 fragment(inputs = (in = u,))),
        reads = torque_reads(), residuals = torque_only, tolerances = (torque = 1e-9,))

    # θ = 0, ω = 0, u = 0 is an equilibrium: nothing to pack, no seeded
    # activation, no backend call — the establishment round is the evaluation.
    at_rest = Simulation(fed(Pendulum(), :u); h = 1//10)
    yes = trim!(at_rest, probe(0.0); baseline = pend_base())
    @test yes.converged && yes.status === :bypassed
    @test yes.nevals == 1 && yes.niters == 0
    @test yes.solution === (;) && isempty(yes.saturated)
    @test yes.committed_residuals !== nothing && lifecycle(at_rest) === :initialized

    # And the same probe on a baseline that is not an equilibrium answers no,
    # by the ordinary box test, leaving the simulation untouched.
    forced = Simulation(fed(Pendulum(), :u); h = 1//10)
    no = trim!(forced, probe(1.0); baseline = pend_base())
    @test !no.converged && no.status === :bypassed
    @test no.residuals.torque ≈ 1.0 && no.committed_residuals === nothing
    @test lifecycle(forced) === :built
end

@testset "a malformed problem is a collected `TrimProblemInvalid` (§14.7, §13.1)" begin
    sim = Simulation(fed(Pendulum(), :u); h = 1//10)
    before = world(sim)

    # A bounds key-set mismatch, an `Int` guess field and an unresolvable
    # selector, in one throw — the read set's own violation kept verbatim, kind
    # name and all, because the *problem* is what is malformed here (§14.8).
    e = failure(() -> trim!(sim, TrimProblem(
        guess = (u = 0,), lower = (v = -Inf,), upper = (u = Inf,),
        condition = decide_u, reads = reads(ω̇ = get_deriv("c", :ω),
                                            nope = get_state("nope", :q)),
        residuals = torque_only, tolerances = (torque = 1e-9,)); baseline = pend_base()))
    @test e isa BuildError && startswith(e.msg, "TrimProblemInvalid:")
    @test occursin("3 violations", e.msg)
    @test occursin("`lower` names `v` and `guess` names `u`", e.msg)
    @test occursin("`u`::Int64", e.msg)
    @test occursin("ReadResolution: the read labeled `nope`", e.msg)
    @test world(sim) == before

    # The residual key set is the one thing only the setup guess evaluation can
    # observe (§14.7), so it is reported from there — its own collected throw.
    e2 = failure(() -> trim!(sim, TrimProblem(
        guess = (u = 0.0,), lower = (u = -Inf,), upper = (u = Inf,),
        condition = decide_u, reads = torque_reads(),
        residuals = (r, d) -> (wrong = r.ω̇, extra = 1), tolerances = (torque = 1e-9,));
        baseline = pend_base()))
    @test e2 isa BuildError && startswith(e2.msg, "TrimProblemInvalid:")
    @test occursin("`residuals` returns `wrong`, `extra`", e2.msg) &&
          occursin("`tolerances` names `torque`", e2.msg)
    @test world(sim) == before

    # A tolerance that is not a `Float64`, and a read set spelled bare.
    e3 = failure(() -> trim!(sim, TrimProblem(
        guess = (u = 0.0,), lower = (u = -Inf,), upper = (u = Inf,),
        condition = decide_u, reads = (ω̇ = get_deriv("c", :ω),),
        residuals = torque_only, tolerances = (torque = 1,)); baseline = pend_base()))
    @test occursin("`tolerances` field(s) `torque`::Int64", e3.msg)
    @test occursin("`reads` is", e3.msg) && occursin("reads(", e3.msg)
    @test world(sim) == before

    # An inverted box admits no point at all, and the collecting pass says so
    # per decision, with both values in hand — before the projection at the
    # pack site could quietly answer the upper bound instead.
    e4 = failure(() -> trim!(sim, TrimProblem(
        guess = (θ = 0.1, u = 0.0), lower = (θ = -π/2, u = 1.0),
        upper = (θ = π/2, u = -1.0), condition = decide_both, reads = both_reads(),
        residuals = both_residuals, tolerances = (torque = 1e-9, hold = 1e-9));
        baseline = pend_base()))
    @test e4 isa BuildError && startswith(e4.msg, "TrimProblemInvalid:")
    @test occursin("1 violation ", e4.msg)
    @test occursin("`lower` names `u` = 1.0 above `upper`'s -1.0", e4.msg)
    @test occursin("an inverted pair admits no point at all", e4.msg)
    @test world(sim) == before

    # A tolerance is the half-width of a box, so zero and negative name no box
    # at all — and the normalized acceptance test divides by them, which would
    # otherwise reject every trial step and stall the solve at the guess. Both
    # named in one throw.
    e5 = failure(() -> trim!(sim, TrimProblem(
        guess = (θ = 0.1, u = 0.0), lower = (θ = -π/2, u = -Inf),
        upper = (θ = π/2, u = Inf), condition = decide_both, reads = both_reads(),
        residuals = both_residuals, tolerances = (torque = 0.0, hold = -1e-9));
        baseline = pend_base()))
    @test e5 isa BuildError && startswith(e5.msg, "TrimProblemInvalid:")
    @test occursin("2 violations", e5.msg)
    @test occursin("`tolerances` names `torque` = 0.0", e5.msg)
    @test occursin("`tolerances` names `hold` = -1.0e-9", e5.msg)
    @test occursin("finite and strictly positive", e5.msg)
    @test world(sim) == before
end

@testset "the box is honored at every point the backend returns (§14.8, D-070)" begin
    # The guess is projected at the pack site, so the point the backend starts
    # from already lies in the box — including the two returns that never step,
    # the already-within-tolerance one and a stall at the first iteration.
    # Without it an out-of-box guess commits as converged with `saturated`
    # empty, sitting at no bound because it is outside both.
    outside = TrimProblem(guess = (u = 100.0,), lower = (u = -1.0,), upper = (u = 1.0,),
                          condition = decide_u, reads = torque_reads(),
                          residuals = torque_only, tolerances = (torque = 10.0,))
    sim = Simulation(fed(Pendulum(), :u); h = 1//10)
    report = trim!(sim, outside; baseline = pend_base())

    @test report.converged                            # the box point meets the loose tolerance
    @test report.solution.u === 1.0                   # projected, not the 100.0 it was given
    @test report.saturated == [(:u, :upper)]          # and honestly named as saturated
    @test port(sim, "", :in) === 1.0                  # the committed root input is the box's

    # A degenerate box pins the decision outright: `lower == upper`, and the
    # projected guess is that one point. `_saturated` tests `lower` first, so a
    # point that is both bounds is reported as `:lower`.
    degenerate = TrimProblem(guess = (u = 0.0,), lower = (u = 2.0,), upper = (u = 2.0,),
                             condition = decide_u, reads = torque_reads(),
                             residuals = torque_only, tolerances = (torque = 10.0,))
    pinned = Simulation(fed(Pendulum(), :u); h = 1//10)
    r2 = trim!(pinned, degenerate; baseline = pend_base())
    @test r2.converged && r2.solution.u === 2.0
    @test r2.saturated == [(:u, :lower)]
    @test port(pinned, "", :in) === 2.0
end

@testset "the residual return is re-checked at the first seeded point (§14.7, §14.8)" begin
    # The nominal guess evaluation observes the return at `Float64`; a lambda
    # that branches on the scalar answers a different key set at the seeded
    # activation, and that is the problem being malformed rather than the bare
    # `ErrorException` the reorder to `tolerances`' order would raise.
    sim = Simulation(fed(Pendulum(), :u); h = 1//10)
    before = world(sim)
    e = failure(() -> trim!(sim, TrimProblem(
        guess = (u = 0.0,), lower = (u = -Inf,), upper = (u = Inf,),
        condition = decide_u, reads = torque_reads(),
        residuals = eltype_split, tolerances = (torque = 1e-9,)); baseline = pend_base()))
    @test e isa BuildError && startswith(e.msg, "TrimProblemInvalid:")
    @test occursin("`residuals` returns `wrong`", e.msg) &&
          occursin("`tolerances` names `torque`", e.msg)
    @test world(sim) == before && lifecycle(sim) === :built
end

@testset "an incomplete baseline is `UninitializedInputs` at setup (§14.6, §14.8)" begin
    # Setup's own application to the scratch stores is one of §14.6's three
    # sites, and its coverage is a plan-level fact: the check runs before any
    # evaluation, not merely before any write.
    sim = Simulation(fed(Pendulum(), :u); h = 1//10)
    before = world(sim)
    e = failure(() -> trim!(sim, TrimProblem(
        guess = (θ = 0.1,), lower = (θ = -1.0,), upper = (θ = 1.0,),
        condition = d -> at("c", fragment(x = (θ = d.θ, ω = 0.0))),
        reads = torque_reads(), residuals = torque_only, tolerances = (torque = 1e-9,));
        baseline = fragment()))
    @test e isa BuildError && startswith(e.msg, "UninitializedInputs:")
    @test occursin("root input(s) `in`", e.msg) && occursin("`trim!`", e.msg)
    @test world(sim) == before && lifecycle(sim) === :built
end

@testset "the scratch world's frozen cells are established, not probe-seeded (D-213)" begin
    # The torque reaching the pendulum is a discrete producer's held output, and
    # that producer's `s` is the baseline's. At the seeded activation the
    # discrete tier never runs (§9.4), so this cell can only come from the
    # nominal half's establishment round: `asin(4/g_l)` is reachable exactly
    # when it holds the authored `s`'s output, and the probe's synthesized zero
    # would put the solution at 0 instead.
    sim = Simulation(sampled_pend(); h = 1//10)
    report = trim!(sim, θ_alone_problem(); baseline = sampled_base(4.0))
    @test report.converged
    @test report.solution.θ ≈ asin(4.0 / PEND_G_L)
    @test !isapprox(report.solution.θ, 0.0; atol = 1e-3)
    @test state(sim, "ctl") === (acc = 4.0,)          # and the authored `s` committed

    # The same model with a *different* authored `s` moves the solution with it,
    # which is what makes the copy load-bearing rather than incidental.
    other = Simulation(sampled_pend(); h = 1//10)
    @test trim!(other, θ_alone_problem(); baseline = sampled_base(2.0)).solution.θ ≈
          asin(2.0 / PEND_G_L)

    # And a decision variable authored *into* that frozen `s` is refused at
    # resolution, by the seeded half's own `ConditionResolution` — the tier
    # cannot carry partials, and the refusal says so rather than folding into
    # the problem's kind.
    refused = Simulation(sampled_pend(); h = 1//10)
    e = failure(() -> trim!(refused, TrimProblem(
        guess = (θ = 0.1, acc = 4.0), lower = (θ = -π/2, acc = -Inf),
        upper = (θ = π/2, acc = Inf),
        condition = d -> combine(at("c", fragment(x = (θ = d.θ, ω = 0.0))),
                                 at("ctl", fragment(s = (acc = d.acc,)))),
        reads = torque_reads(), residuals = torque_only, tolerances = (torque = 1e-9,));
        baseline = fragment(inputs = (in = 0.0,))))
    @test e isa BuildError && occursin("ConditionResolution: `s.acc` at `ctl`", e.msg)
    @test occursin("a decision variable descends into neither a frozen discrete `s` nor a " *
                   "pinned leaf", e.msg)
    @test lifecycle(refused) === :built               # nothing was written to the sim
end

@testset "the commit is an ordinary boundary zero: fired events are reported (§14.5)" begin
    # A guard the solved attitude already holds fires at `t₀` — derived, not
    # asserted (D-067) — and the report carries the list beside the warning.
    sim = Simulation(triggered_pend(); h = 1//10)
    report = @test_logs (:warn, r"^TrimCommitEvents") trim!(sim, u_problem();
                                                            baseline = pend_base())
    @test report.converged
    @test report.fired_events == [("trig", :fire)]
    @test modes(sim, "trig") === (state = :fired, count = 1)
    # This handler writes only modes, so the committed-state residuals are still
    # the solved ones: a fired event is not by itself a residual mover.
    @test abs(report.committed_residuals.torque) ≤ report.tolerances.torque
end

@testset "a mover at the commit shows in the committed-state residuals (§14.8)" begin
    # The handler resets the attitude the solve chose, so the committed stores
    # sit where the reported solution never described. Both warnings fire, in
    # order, and the numbers that describe where the simulation actually *is*
    # are the committed ones — the verdict is not re-litigated (D-150).
    sim = Simulation(fed(Snapback(0.3), :u); h = 1//10)
    problem = TrimProblem(guess = (u = 0.0,), lower = (u = -Inf,), upper = (u = Inf,),
                          condition = snap_decide_u, reads = torque_reads(),
                          residuals = torque_only, tolerances = (torque = 1e-9,))
    report = @test_logs (:warn, r"^TrimCommitEvents") (:warn, r"^TrimCommitResiduals") begin
        trim!(sim, problem; baseline = pend_base())
    end

    @test report.converged                            # at the solved point, and it stands
    @test abs(report.residuals.torque) ≤ report.tolerances.torque
    @test report.fired_events == [("c", :snap)]
    @test state(sim, "c").θ === 0.0                   # the handler's write, committed
    # ω̇ at θ = 0 is the solved torque itself, not the zero it balanced.
    @test report.committed_residuals.torque ≈ report.solution.u
    @test abs(report.committed_residuals.torque) > report.tolerances.torque
end

@testset "the commit is literally an `init!` over the same composite (§14.8, §14.9)" begin
    sim = Simulation(fed(Pendulum(), :u); h = 1//10)
    twin = Simulation(fed(Pendulum(), :u); h = 1//10)
    report = trim!(sim, u_problem(); baseline = pend_base(), t₀ = 0.75)

    # The hand-written spelling of what `trim!` ran: the same composite over the
    # same baseline at the same anchor, which is why the commit's totality is
    # setup's and the check is structurally unfailable through this path.
    init!(twin, override(pend_base(), decide_u(report.solution)); t₀ = 0.75)
    @test world(twin) == world(sim)
    @test lifecycle(twin) === lifecycle(sim)
end

@testset "a warm restart trims from a capture at its own time (§14.1, §14.8)" begin
    sim = Simulation(fed(Pendulum(), :u); h = 1//10)
    init!(sim, combine(at("c", condition(Pendulum(); θ = 0.2)),
                       fragment(inputs = (in = 1.0,))))
    run!(sim; t_end = 0.4)
    (c, t) = capture(sim)
    @test lifecycle(sim) === :stopped && t === 0.4

    # `trim!(sim, problem; baseline = c, t₀ = t)` is §14.8's resumed spelling:
    # continuity is explicit, and the anchor comes back from the capture.
    report = trim!(sim, u_problem(); baseline = c, t₀ = t)
    @test report.converged && report.committed_residuals !== nothing
    @test lifecycle(sim) === :initialized
    @test sim.exec.clock.t === 0.4 && sim.exec.clock.t₀ === 0.4
    @test state(sim, "c") === (θ = 0.5, ω = 0.0)      # the problem's condition won
    @test port(sim, "", :in) === report.solution.u
end

@testset "the per-iteration write and read are free at the seeded activation (§7.5)" begin
    # The iteration is a raw write → sweep → read cycle (§14.5), and both ends
    # of it are compiled: the shape-compiled plan and the compiled reader, at
    # the service's own seeded scalar. The user's residual lambda is theirs, so
    # the gates are on the two the framework owns.
    sim = Simulation(fed(Pendulum(), :u); h = 1//10)
    b = sim.build
    TD = ForwardDiff.Dual{TrimTag,Float64,1}
    ex = compile(b, activation(b, TD), sim.D, sim.Φ, sim.Δt; chunk_size = sim.chunk_size)
    seeded(v) = (θ = ForwardDiff.Dual{TrimTag}(v, 1.0),)
    plan = compile_plan(override(pend_base(), decide_θ(seeded(0.1))), b, TD)
    reader = compile_reads(torque_reads(), b, TD)
    tree = override(pend_base(), decide_θ(seeded(0.2)))

    apply!(ex, plan, tree)
    evaluate!(ex)
    @test gather(reader, ex).ω̇ isa TD
    @test (@ballocated apply!($ex, $plan, $tree)) == 0
    @test (@ballocated gather($reader, $ex)) == 0

    # And the seeded pass yields value and partials together: one sweep is `r`
    # and the column of `J` beside it (§14.7).
    v = gather(reader, ex).ω̇
    @test ForwardDiff.value(v) ≈ -PEND_G_L * sin(0.2) + 4.0
    @test ForwardDiff.partials(v, 1) ≈ -PEND_G_L * cos(0.2)
end

@testset "`trim!` is a stopped-sim service on a nominal deployment (§14.8, §12.6)" begin
    # A non-nominal `Simulation` is refused outright: the commit runs through
    # boundary zero on the simulation's own stores, and those are the nominal
    # world's — the seeded activation is the service's scratch, never the
    # deployment's.
    dual = Simulation(fed(Pendulum(), :u), D8; h = 1//10)
    e = failure(() -> trim!(dual, u_problem(); baseline = pend_base()))
    @test e isa BuildError && occursin("needs a nominal `Simulation{Float64}`", e.msg)

    # And a value that is not a problem is a directive, not a `MethodError`.
    plain = Simulation(fed(Pendulum(), :u); h = 1//10)
    e2 = failure(() -> trim!(plain, (guess = (u = 0.0,),); baseline = pend_base()))
    @test e2 isa BuildError && occursin("takes a `TrimProblem`", e2.msg)

    # `running` is the §11.3 freeze, as for every other §14 service. Both ends
    # of the run are test-controlled, exactly as in test_readers.
    live = Simulation(armed(); h = 1//100, t_end = 3.0e5, stop_on = ("stop",))
    init!(live, fragment(inputs = (in = 0.0,)))
    task = Threads.@spawn run!(live)
    while lifecycle(live) !== :running
        yield()
    end
    err = failure(() -> trim!(live, u_problem(); baseline = pend_base()))
    stage!(live, "in" => 1.0)
    wait(task)
    @test err isa BuildError && occursin("ServiceLifecycle", err.msg)
    @test occursin("stopped-sim operation and the simulation is running", err.msg)
end
