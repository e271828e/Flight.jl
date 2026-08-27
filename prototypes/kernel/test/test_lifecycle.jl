# --- run lifecycle and termination (§12.6, §13.5, §13.6; increment 16) ----------
# The five-state machine behind `lifecycle(sim)`, `run!`'s keyword policy with
# the constructor defaults, partial advance, the §13.5 termination record and
# §13.6's abnormal entry. The devices below live at top level for the README's
# local-scope reason.

# A monitored ramp: `hit` goes true at the first boundary whose sweep sees the
# ramp at the trigger's level — the boundary-detected stop face.
monitored() = Group((; src = Ramp(0.0), trig = Trigger(0.35)),
                    ("children/src/out" => "children/trig/sig",), (),
                    ("children/trig/on" => "hit",))

# A slot-fed trigger exporting its flag: the boundary-zero stop's model.
armed() = Group((; c = Trigger(0.5)), (), ("in" => "children/c/sig",),
                ("children/c/on" => "stop",))

# The touchdown archetype (§13.5): a sawtooth crossing the overload's level
# mid-frame, so the stop localizes to the crossing's t* boundary.
overloaded() = Group((; src = Sawtooth(1.0), mon = Overload(0.315)),
                     ("children/src/q" => "children/mon/sig",), (),
                     ("children/mon/tripped" => "tripped",))

# A resource-bracket witness for the abnormal tail, and the empty-claim binding
# that rosters it.
mutable struct TailProbe <: AbstractDevice
    log::Vector{Symbol}
end
TailProbe() = TailProbe(Symbol[])
init!(d::TailProbe) = (push!(d.log, :init); nothing)
shutdown!(d::TailProbe) = (push!(d.log, :shutdown); nothing)
function loop(d::TailProbe, h)
    while running(h)
        wait_next_snapshot(h)
    end
    nothing
end
struct NoClaim <: AbstractBinding end
is_input(::NoClaim) = true
claims(::NoClaim) = ()

@testset "the five states, and the gates between them (§12.6)" begin
    sim = Simulation(feedback_model(); h = 1//50, t_end = 1.0)
    @test lifecycle(sim) === :built
    @test termination(sim) === nothing
    e = failure(() -> run!(sim))
    @test e isa BuildError && occursin("`init!` is mandatory", e.msg)
    @test occursin("mandatory", failure(() -> step!(sim)).msg)

    init!(sim, fragment(slots = (ref = 0.0,)))
    @test lifecycle(sim) === :initialized
    init!(sim, fragment(slots = (ref = 0.0,)))  # a warm restart from initialized is legal
    run!(sim)
    @test lifecycle(sim) === :stopped
    e = failure(() -> run!(sim))
    @test e isa BuildError && occursin("stopped → init! → run!", e.msg)
    @test occursin("stopped → init! → run!", failure(() -> step!(sim)).msg)
    init!(sim, fragment(slots = (ref = 0.0,)))  # the supported cycle reopens it
    @test lifecycle(sim) === :initialized
    @test termination(sim) === nothing               # the record cleared with the trajectory
end

@testset "the freeze is the lifecycle's :running — init! and run! refuse it too (§12.6)" begin
    sim = Simulation(fed(Plant(), "u"); h = 1//100000)
    total = fragment(slots = (in = 0.0,))
    init!(sim, total)
    init!(sim, total)                                # warms init!'s path — including this
    #                                                  condition's resolution shape (§14.3):
    #                                                  the mid-run probe below races no JIT
    t = Threads.@spawn run!(sim; t_end = 1.0)        # 100k frames: alive throughout the checks
    while lifecycle(sim) !== :running
        yield()
    end
    err_i = try init!(sim, total) catch e; e end
    err_r = try run!(sim; t_end = 2.0) catch e; e end
    wait(t)
    @test err_i isa BuildError && occursin("ServiceLifecycle", err_i.msg)
    @test err_r isa BuildError && occursin("already running", err_r.msg)
    @test lifecycle(sim) === :stopped
    @test termination(sim).source === EndTimeReached()
end

@testset "t_end: constructor default, run! override, and a bound owed from some site (§13.5)" begin
    sim = Simulation(feedback_model(); h = 1//50, t_end = 1.0)
    init!(sim, fragment(slots = (ref = 0.0,)))
    run!(sim)                                        # the constructor's default
    t = termination(sim)
    @test t isa TerminationRecord{Float64}           # the deployment's own scalar (§7.2, D-203)
    @test t.source === EndTimeReached() && t.t == 1.0
    @test isempty(t.residue)                         # a quiet tail contributes no record
    @test sim.clock.step == 50

    init!(sim, fragment(slots = (ref = 0.0,)))
    run!(sim; t_end = 0.5)                           # this run only
    @test termination(sim).t == 0.5
    init!(sim, fragment(slots = (ref = 0.0,)))
    run!(sim)                                        # nothing was mutated: the default again
    @test termination(sim).t == 1.0

    unbound = Simulation(feedback_model(); h = 1//50)
    init!(unbound, fragment(slots = (ref = 0.0,)))
    e = failure(() -> run!(unbound))
    @test e isa BuildError && occursin("neither at the constructor nor here", e.msg)
    # The override is validated exactly as the constructor validates the default.
    @test failure(() -> Simulation(feedback_model(); h = 1//50, t_end = -1.0)).msg ==
          failure(() -> run!(unbound; t_end = -1.0)).msg
end

@testset "stop_on names root-exported Bool output faces, validated at both sites (§13.5)" begin
    m = feedback_model()                             # "ref" a slot, "y" a Float64 export
    sim = Simulation(m; h = 1//50, t_end = 1.0)
    init!(sim, fragment(slots = (ref = 0.0,)))
    for (bad, frag) in (("nope", "no root face"), ("ref", "root input slot"), ("y", "Bool"))
        ec = failure(() -> Simulation(m; h = 1//50, stop_on = (bad,)))
        er = failure(() -> run!(sim; stop_on = (bad,)))
        @test ec isa BuildError && occursin(frag, ec.msg)
        @test ec.msg == er.msg                       # identical at both binding sites
    end
    @test lifecycle(sim) === :initialized            # a rejected run! bound nothing
end

@testset "a boundary-detected stop face ends the run at its own boundary (§13.5)" begin
    sim = Simulation(monitored(); h = 1//10, t_end = 5.0, stop_on = ("hit",))
    init!(sim)
    run!(sim)
    t = termination(sim)
    @test t.source === ModelRequestedStop(:hit)      # kind + payload, one typed value (D-203)
    @test t.t == 4 * sim.h                           # the sweep at boundary 4 saw 0.4 ≥ 0.35
    @test sim.clock.step == 4                        # the run ended there, not at t_end
    @test latest(sim).t === t.t                      # that snapshot is the final one
end

@testset "an authored condition already terminal ends the run at t₀, integrating nothing (§13.5)" begin
    sim = Simulation(armed(); h = 1//10, t_end = 5.0, stop_on = ("stop",))
    init!(sim, fragment(slots = (in = 1.0,)))        # holds in the authored state:
    #                                                boundary zero derives the firing (§10.6)
    run!(sim)
    t = termination(sim)
    @test t.source === ModelRequestedStop(:stop) && t.t == 0.0
    @test sim.clock.step == 0                        # zero frames: the check precedes the first step
end

@testset "a localized stop ends the run at t*, the crossing state final (§13.5, §10.4)" begin
    sim = Simulation(overloaded(); h = 1//10, t_end = 5.0, stop_on = ("tripped",))
    init!(sim)
    run!(sim)
    t = termination(sim)
    @test t.source === ModelRequestedStop(:tripped)
    @test t.t ≈ 0.315 atol = 1e-6                    # the analytic crossing, not a frame top
    @test t.t == sim.clock.t                         # the frame's remainder was abandoned
    @test latest(sim).t === t.t
    @test logged(sim)[end] === latest(sim)           # the log's terminal endpoint is the t* boundary
end

@testset "stop_on overrides bind per run, like t_end (§13.5)" begin
    sim = Simulation(monitored(); h = 1//10, t_end = 1.0)    # no default stop faces
    init!(sim)
    run!(sim; stop_on = ("hit",))
    @test termination(sim).source isa ModelRequestedStop
    init!(sim)
    run!(sim)                                        # the constructor's (), again
    @test termination(sim).source === EndTimeReached() && termination(sim).t == 1.0
end

@testset "step! advances whole frames and returns the count actually advanced (§12.6)" begin
    sim = Simulation(feedback_model(); h = 1//50, t_end = 1.0)
    init!(sim, fragment(slots = (ref = 0.0,)))
    @test step!(sim) == 1                            # the frames = 1 default
    @test step!(sim; frames = 4) == 4
    @test step!(sim; t_plus = 0.5) == 25             # the duration spelling
    @test lifecycle(sim) === :initialized            # between calls: ready to advance
    @test step!(sim; frames = 100) == 20             # t_end truncates at frame 50
    @test lifecycle(sim) === :stopped
    @test termination(sim).source === EndTimeReached()

    # A stepped trajectory is bit-identical to the same frames under run!.
    ref = Simulation(feedback_model(); h = 1//50, t_end = 1.0)
    init!(ref, fragment(slots = (ref = 0.0,)))
    run!(ref)
    @test port(sim, "children/plant", :y) === port(ref, "children/plant", :y)
    @test state(sim, "children/plant").q === state(ref, "children/plant").q

    sim2 = Simulation(feedback_model(); h = 1//50)
    init!(sim2, fragment(slots = (ref = 0.0,)))
    @test occursin("not both", failure(() -> step!(sim2; frames = 1, t_plus = 0.1)).msg)
    @test occursin("integer ≥ 1", failure(() -> step!(sim2; frames = 0)).msg)
    @test occursin("t_plus", failure(() -> step!(sim2; t_plus = 0.0)).msg)
end

@testset "a stop face inside step! truncates it through the deviceless tail (§12.6, §13.5)" begin
    sim = Simulation(monitored(); h = 1//10, t_end = 5.0, stop_on = ("hit",))
    init!(sim)
    @test step!(sim; frames = 2) == 2                # short of the trigger: an ordinary advance
    @test lifecycle(sim) === :initialized
    @test step!(sim; frames = 10) == 2               # the face holds at boundary 4
    @test lifecycle(sim) === :stopped
    @test termination(sim).source === ModelRequestedStop(:hit)
end

@testset "§13.6: a loop-side throw discards the failed boundary and promotes the last one" begin
    sim = Simulation(fed(Exploder(), "arm"); h = 1//10, t_end = 5.0)
    probe = TailProbe()
    attach!(sim, probe, NoClaim())
    init!(sim, fragment(slots = (in = 0.0,)))
    stage!(sim, "in" => true)                        # armed: frame 1's drain applies it,
    @test_throws Exploded run!(sim)                  # frame 1's integration throws (§13.4's
                                                     # synchronous rethrow, after the tail)
    @test lifecycle(sim) === :errored
    t = termination(sim)
    @test t.source isa LoopError && t.source.exception isa Exploded
    # The failed boundary published nothing: boundary zero is the promoted
    # final snapshot, and the published record ends at it.
    @test t.t == 0.0 && latest(sim).t == 0.0
    @test [s.frame for s in logged(sim)] == [0]      # a post-mortem read, admitted (§13.6)
    @test probe.log == [:init, :shutdown]            # the device took the ordinary tail
    # The stores may hold mid-boundary values — retained for inspection,
    # readable, and worth nothing more than inspection.
    @test state(sim, "children/c").q isa Float64

    # Errored is terminal: never advanced, never re-initialized.
    @test occursin("errored", failure(() -> run!(sim)).msg)
    @test occursin("errored", failure(() -> step!(sim)).msg)
    @test occursin("errored", failure(() -> init!(sim)).msg)
end
