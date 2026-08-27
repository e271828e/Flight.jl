# --- the device contract and its tasks (§11.6, §11.1, §12.4; increment 12) ------
# The handle, the wrapper, the pre-spawn bracket and the tail. Every device
# below lives at top level for the README's local-scope reason, and every
# timing-sensitive assertion is a property check with generous slack, never a
# tight wall-clock equality: single-threaded schedulers must pass.

# A one-shot input device: stages once through its handle, requests the stop,
# and records the contract calls in order on whatever task runs them.
mutable struct OneShot <: AbstractDevice
    v::Float64
    fired::Bool
    log::Vector{Symbol}
end
OneShot(v) = OneShot(v, false, Symbol[])
init!(d::OneShot) = (push!(d.log, :init); nothing)
shutdown!(d::OneShot) = (push!(d.log, :shutdown); nothing)
function loop(d::OneShot, h)
    d.fired && return nothing
    d.fired = true
    push!(d.log, :loop)
    stage!(h, "a" => d.v)
    stop!(h)
    nothing
end

# A boundary-driven consumer: the §12.3 wait as its whole wait structure.
mutable struct Collector <: AbstractDevice
    seen::Vector{Tuple{Float64,Int}}    # (t, boundary ordinal) per observed snapshot
    log::Vector{Symbol}
end
Collector() = Collector(Tuple{Float64,Int}[], Symbol[])
function loop(d::Collector, h)
    while true
        snap = wait_next_snapshot(h)         # returns on publication or on the stop wake
        snap === nothing || push!(d.seen, (Float64(snap.t), snap.boundary))
        running(h) || break                  # record first: the stop wake carries the final world
    end
    push!(d.log, :returned)
    nothing
end

# A crasher: the wrapper's catch path, with shutdown! still guaranteed.
mutable struct Crasher <: AbstractDevice
    log::Vector{Symbol}
end
Crasher() = Crasher(Symbol[])
shutdown!(d::Crasher) = (push!(d.log, :shutdown); nothing)
loop(d::Crasher, h) = (push!(d.log, :loop); error("boom"))

# A failed acquisition: the §12.4 initialization bracket's customer.
mutable struct BadInit <: AbstractDevice
    log::Vector{Symbol}
end
BadInit() = BadInit(Symbol[])
init!(d::BadInit) = (push!(d.log, :init); error("no hardware"))
shutdown!(d::BadInit) = (push!(d.log, :shutdown); nothing)
loop(d::BadInit, h) = (push!(d.log, :loop); nothing)

# The taught obligation violated: no running check, a blocking call with no
# unblock! override — the join-timeout path's customer.
mutable struct Stubborn <: AbstractDevice
    log::Vector{Symbol}
end
Stubborn() = Stubborn(Symbol[])
shutdown!(d::Stubborn) = (push!(d.log, :shutdown); nothing)
loop(d::Stubborn, h) = (sleep(0.8); push!(d.log, :woke); nothing)

# The obligation met: a blocking wait made interruptible by unblock! (§12.4(3)).
mutable struct Blocked <: AbstractDevice
    ch::Channel{Int}
    log::Vector{Symbol}
end
Blocked() = Blocked(Channel{Int}(1), Symbol[])
unblock!(d::Blocked) = (close(d.ch); nothing)
shutdown!(d::Blocked) = (push!(d.log, :shutdown); nothing)
function loop(d::Blocked, h)
    while running(h)
        try
            take!(d.ch)                       # blocks; unblock! closes → raises here
        catch
            break                             # the raise is shutdown, not a crash
        end
    end
    nothing
end

# A calling-task device: records which task ran its body (§11.1's pinning),
# polling between running checks and never blocking across them (§12.4).
mutable struct Inline <: AbstractDevice
    task::Union{Nothing,Task}
    log::Vector{Symbol}
end
Inline() = Inline(nothing, Symbol[])
needs_calling_task(::Inline) = true
function loop(d::Inline, h)
    d.task = current_task()
    while running(h)
        sleep(0.001)
    end
    push!(d.log, :returned)
    nothing
end

# A device with no loop method at all: the error-throwing fallback's customer.
mutable struct Loopless <: AbstractDevice end

@testset "a device stages through its handle from its own task, and departure consults should_abort (§11.6, §12.4)" begin
    sim = Simulation(two_slots(); h = 1//10)
    dev = OneShot(0.7)
    h = attach!(sim, dev, Enumerated("a"); should_abort = true)
    init!(sim)
    run!(sim; t_end = 1000.0)                        # ends by the device's stop, not t_end
    @test sim.clock.step < 10000             # the stop truncated the run
    @test dev.log[1:3] == [:init, :loop, :shutdown]
    @test !running(h)                        # the sticky status, read off the handle
    # The record names the channel and its issuer (§13.5, D-203): the stop
    # rode the departing device's should_abort, so the device is the issuer.
    @test termination(sim).source === ControlRequestedStop("device 1 (OneShot)")
    # The staged batch was applied by a drain the stop did not beat, or still
    # pends in the cell — exactly one of the two, timing's choice — and init!
    # clears whatever pends with the trajectory it predates (§12.6).
    @test (port(sim, "", :a) === 0.7) ⊻ ((@atomic h.writer.cell.pending) !== nothing)
    init!(sim)
    @test (@atomic h.writer.cell.pending) === nothing
end

@testset "wait_next_snapshot observes ordered boundaries and wakes on the stop (§12.3, §12.4)" begin
    sim = Simulation(two_slots(); h = 1//10)
    dev = Collector()
    attach!(sim, dev, Enumerated())          # the may-write-nothing degenerate: a pure reader
    init!(sim)
    run!(sim; t_end = 1.0)
    @test dev.log == [:returned]             # exited through the woken wait, before the join
    @test !isempty(dev.seen)                 # at least one boundary observed in ten frames
    @test issorted([b for (_, b) in dev.seen])       # never a stale wake: newest-wins only
    @test issorted([t for (t, _) in dev.seen])
    # While stopped the wait returns at once instead of parking (§12.3's
    # predicate routes on the sticky status): re-read through the handle.
    snap = wait_next_snapshot(sim.plane.roster[1].handle)
    @test snap === latest(sim)
end

@testset "the boundary ordinal rides in the snapshot (§12.3)" begin
    sim = Simulation(two_slots(); h = 1//10)
    init!(sim)                               # boundary zero: ordinal 0
    run!(sim; t_end = 0.5)
    @test latest(sim).boundary == 5
    @test latest(sim).t ≈ 0.5
    @test logged(sim)[1].boundary == 0
end

@testset "stop! from any task ends the run at a frame top (§12.1, §12.4)" begin
    sim = Simulation(two_slots(); h = 1//10)
    attach!(sim, Pad("p"), Enumerated("a"))  # a rostered device keeps the loop yielding (§12.2)
    init!(sim)
    stopper = Threads.@spawn (sleep(0.05); stop!(sim))
    run!(sim; t_end = 1.0e6)
    wait(stopper)
    @test sim.clock.step < 10^7              # truncated, and stopped is sticky:
    @test !running(sim.plane.roster[1].handle)
    # One channel, and the record names who spoke (§13.5, D-203): stop!(sim)
    # is calling code from any task, issuer :code.
    @test termination(sim).source === ControlRequestedStop(:code)
    # A fresh trajectory owes nothing to this stop: init! clears the word.
    init!(sim)
    @test step!(sim; frames = 3) == 3
    @test sim.clock.step == 3
end

@testset "a crash is caught, shutdown! runs, the run continues, claims persist (§12.4(6))" begin
    sim = Simulation(two_slots(); h = 1//10)
    dev = Crasher()
    attach!(sim, dev, Enumerated("a"))
    init!(sim)
    logs, _ = Test.collect_test_logs() do
        run!(sim; t_end = 0.5)
    end
    @test sim.clock.step == 5                # the run reached t_end regardless
    @test dev.log == [:loop, :shutdown]      # the bracket held on the crash path
    @test crash_accounted(sim, logs, "device 1 (Crasher)")
    # Death is not detach: the claim stands, and the harness cannot take the face.
    stage!(sim, "a" => 9.0)
    cfe = only((@atomic sim.plane.harness_diag.batch).ring)
    @test cfe isa ClaimedFaceEntry && cfe.incumbent == "device 1 (Crasher)"
end

@testset "a crash under should_abort requests the stop (§12.4(6))" begin
    sim = Simulation(two_slots(); h = 1//10)
    attach!(sim, Crasher(), Enumerated("a"); should_abort = true)
    init!(sim)
    logs, _ = Test.collect_test_logs() do
        run!(sim; t_end = 1000.0)
    end
    @test sim.clock.step < 10000             # ended by the crash's stop, not t_end
    @test crash_accounted(sim, logs, "device 1 (Crasher)")
end

@testset "a failed init! is bracketed: shutdown!, no task, claims persist (§12.4)" begin
    sim = Simulation(two_slots(); h = 1//10)
    dev = BadInit()
    attach!(sim, dev, Enumerated("a"))
    init!(sim)
    run!(sim; t_end = 0.5)
    @test dev.log == [:init, :shutdown]      # loop never ran: no task was spawned
    @test sim.clock.step == 5                # flag clear: the run proceeds without it
    # The report was written pre-spawn, addressed by the entry (§12.4), so it
    # deterministically makes the first frame top's fold: the first frame's
    # snapshot carries the delta, the terminal status the totals.
    bw = writer_status(latest(sim), "device 1 (BadInit)")
    @test bw.totals.crash == 1
    dc = only(writer_status(logged(sim)[2], "device 1 (BadInit)").recent)
    @test dc isa DeviceCrash && dc.cause isa ErrorException && dc.abort === false
    # Dead from boundary zero, with no marking machinery (§12.4): the cell was
    # never heartbeated — stale against any clock — and no task ever existed.
    @test stale(bw) && bw.task_state === :none
    stage!(sim, "a" => 9.0)                  # claims persist: death is not detach
    @test only((@atomic sim.plane.harness_diag.batch).ring) isa ClaimedFaceEntry
    # With should_abort set the stop is already pending at the loop's start:
    # the run advances zero frames and ends through the same tail — no frame
    # top ever folds the report, so only the run's-end sweep can present it.
    sim2 = Simulation(two_slots(); h = 1//10)
    attach!(sim2, BadInit(), Enumerated("a"); should_abort = true)
    init!(sim2)
    logs, _ = Test.collect_test_logs() do
        run!(sim2; t_end = 0.5)
    end
    @test sim2.clock.step == 0
    @test any(occursin("DeviceCrash from device 1 (BadInit), past the final", string(l.message))
              for l in logs)
    # The same crash is recorded, not just presented (D-203): the residue
    # carries the device's record, and the abort's stop names it as issuer.
    t = termination(sim2)
    @test t.source === ControlRequestedStop("device 1 (BadInit)")
    rr = only(r for r in t.residue if r.writer == "device 1 (BadInit)")
    @test only(rr.recent) isa DeviceCrash
end

@testset "a body ignoring the predicate is abandoned under join_timeout, by name (§12.4(5))" begin
    sim = Simulation(two_slots(); h = 1//10, join_timeout = 0.2)
    dev = Stubborn()
    attach!(sim, dev, Enumerated())
    init!(sim)
    t0 = time()
    # The abandonment is written to the loop's own cell and presented by the
    # run's-end sweep, the record's renderer (§12.4(5), D-203).
    @test_logs (:warn, r"DeviceJoinTimeout from loop, past the final snapshot's account:.*Stubborn") #=
        =# match_mode=:any run!(sim; t_end = 0.3)
    @test time() - t0 < 0.6                  # abandoned at ~0.2 s, not the sleep's 0.8 s
    @test :woke ∉ dev.log                    # the straggler had not returned when run! did
    # Recorded, not just loud (D-203): the termination record's residue holds
    # the structured kind, by name, with the cap and the final boundary.
    rr = only(r for r in termination(sim).residue if r.writer == "loop")
    jt = only(d for d in rr.recent if d isa DeviceJoinTimeout)
    @test jt.who == "device 1 (Stubborn)" && jt.timeout == 0.2
    @test jt.t == termination(sim).t ≈ 0.3 && jt.boundary == latest(sim).boundary
    # Abandonment is not a kill: let the straggler expire inside this testset —
    # its wrapper still runs shutdown! — rather than leave it parked in the
    # timer wheel across process teardown.
    sleep(1.0)
    @test dev.log == [:woke, :shutdown]
end

@testset "unblock! makes the blocking call return: a clean exit, no timeout (§12.4(3))" begin
    sim = Simulation(two_slots(); h = 1//10, join_timeout = 2.0)
    dev = Blocked()
    attach!(sim, dev, Enumerated())
    init!(sim)
    t0 = time()
    logs, _ = Test.collect_test_logs() do
        run!(sim; t_end = 0.3)
    end
    @test time() - t0 < 1.5                  # joined promptly, well inside the cap
    @test !any(occursin("DeviceJoinTimeout", string(l.message)) for l in logs)
    @test isempty(termination(sim).residue)  # nothing landed past the account (D-203)
    @test dev.log == [:shutdown]
end

@testset "a calling-task device runs inline and the loop moves, trajectory untouched (§11.1)" begin
    sim = Simulation(two_slots(); h = 1//10)
    dev = Inline()
    attach!(sim, dev, Enumerated())
    init!(sim)
    caller = current_task()
    run!(sim; t_end = 0.5)
    @test dev.task === caller                # the pinning: the body ran on run!'s task
    @test dev.log == [:returned]             # and left through the ordinary predicate
    @test sim.clock.step == 5
    ref = Simulation(two_slots(); h = 1//10) # the movable loop moved nothing else
    init!(ref)
    run!(ref; t_end = 0.5)
    @test port(sim, "children/s", :e) === port(ref, "children/s", :e)
end

@testset "a device with no loop method crashes loudly, never idles (§11.6)" begin
    sim = Simulation(two_slots(); h = 1//10)
    attach!(sim, Loopless(), Enumerated())
    init!(sim)
    logs, _ = Test.collect_test_logs() do
        run!(sim; t_end = 0.2)
    end
    @test sim.clock.step == 2                # the crash is the device's alone
    @test crash_accounted(sim, logs, "device 1 (Loopless)")
end

@testset "join_timeout is validated and never trajectory-determining (§12.4, D-198)" begin
    err = failure(() -> Simulation(two_slots(); h = 1//10, join_timeout = 0))
    @test err isa BuildError && occursin("join_timeout", err.msg)
    err = failure(() -> Simulation(two_slots(); h = 1//10, join_timeout = "5"))
    @test err isa BuildError && occursin("join_timeout", err.msg)
    trajectories = map((5.0, 0.01)) do cap
        sim = Simulation(two_slots(); h = 1//10, join_timeout = cap)
        attach!(sim, Pad("p"), Enumerated("a"))
        init!(sim)
        run!(sim; t_end = 0.5)
        port(sim, "children/s", :e)
    end
    @test trajectories[1] === trajectories[2]
end
