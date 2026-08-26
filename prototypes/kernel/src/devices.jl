# The device contract and its task machinery (§11.6, §11.1) with the §12.4
# lifecycle slice: the handle every attached device receives, the authoring
# contract's four functions, the framework wrapper around the author-owned
# loop body, the pre-spawn init bracket and the tail. The control surface here
# is §12.1's stop word and §12.3's counter-plus-condition wait, nothing more:
# pause, pacing, §11.8's diagnostic cells and the operator interrupt are
# absent (README), so every degradation below is a plain @warn carrying the
# spec'd payload.
#
# This file holds the types, the contract and the per-run task mechanics;
# `run!`'s frame anatomy and the Simulation-typed surface (`attach!`,
# `stop!(sim)`) live in sim.jl. The helpers here take `sim` untyped for
# include order only — they run once per run, never inside a frame.

"""
The control surface — §12.1's stop, nothing else. `stop_requested` is the
control-plane stop word: set by `stop!` from any task, consulted by the loop
at frame top, cleared at the top of the next run (which is what lets an init
bracket's `should_abort` failure leave a stop *already pending* at the run's
start, §12.4). `stopped` is §12.4(1)'s sticky status: set only after the
final snapshot is published, and read by `running(handle)`, which is how loop
bodies observe the run's end. The §11.3 freeze (`plane.running`) stays a
separate flag deliberately: it spans the whole of `run!`, tail included,
while `stopped` flips at tail step (1) so device loops exit while the joins
are still ahead.

`cond` and `counter` are §12.3's two artifacts: the counter counts *published
boundaries* — grid, `t*`, boundary zero — mirrored under the lock right
behind each release-store of `latest`, in that normative order, so a waiter
observing `counter > last_seen` can never wake onto a stale snapshot. The
counter is monotonic across runs and never re-armed: its absolute value is
nowhere normative, and monotonicity keeps the predicate sound with no
per-run reset.
"""
mutable struct Control
    @atomic stop_requested::Bool
    @atomic stopped::Bool
    cond::Threads.Condition
    counter::Int
end
Control() = Control(false, true, Threads.Condition(), 0)

"""
The handle (§11.6): the one object every attached device receives, carrying
exactly the primitive capabilities — read (`latest`, `wait_next_snapshot`),
stage (`stage!`) and control access (`running`, `stop!`) — and deliberately
*not* the `Simulation`: what a device may touch is what the handle holds.
`attach!` constructs it and returns it, and the same handle is what the
wrapper passes to `loop(dev, handle)` on the device's task. It also carries
the attachment's binding, read back by `binding(handle)`: the loop's own
`map_input`/`map_output` calls take it from the handle instead of the device
carrying its configuration (§11.6). `last_seen` is the §12.3 waiter's private
register, refreshed at spawn so a run's first wait observes that run's
boundaries. One unguarded edge (README): staging through a handle whose
device was detached lands in an orphaned cell and is silently lost — handles
are run-scoped task equipment, and guarding would put a roster scan back
into `stage!`.
"""
mutable struct DeviceHandle
    const id::Int
    const who::String
    const b::AbstractBinding
    const writer::Writer
    const plane::DataPlane
    const ctl::Control
    const published::Published
    const diag::DiagCell
    const gatherer::Union{Nothing,ReadGather}   # the compiled reads; nothing without an output side
    last_seen::Int
end

# --- the authoring contract (§11.6) --------------------------------------------

"""
The authoring contract: four functions, one optional, one trait
(`needs_calling_task`, roster.jl). `init!` and `shutdown!` are the per-run
resource bracket, with no-op defaults — a stub device holds nothing — and
`shutdown!` must tolerate a partially initialized device: it is guaranteed
on every exit path, the failed-`init!` bracket included. `unblock!` is
§12.4(3)'s optional hook, default no-op: an override makes the device's own
blocking call return. `loop` is the one mandatory function — the
author-owned task body, looping on `running(handle)` with interruptible
blocking — and its fallback throws rather than idling: a device whose loop
was never written would otherwise present as "attached, nothing happens",
the hidden-bug class §11.6 tolerates nowhere. The throw lands in the wrapper
and is reported as the `DeviceCrash` it is.
"""
init!(::AbstractDevice) = nothing
shutdown!(::AbstractDevice) = nothing
unblock!(::AbstractDevice) = nothing
loop(dev::AbstractDevice, handle) = error(
    "$(typeof(dev)) defines no `loop` method — the task body is the authoring " *
    "contract's one mandatory function (§11.6)")

# --- the handle primitives (§11.6) ---------------------------------------------

"""
    running(handle)

§12.4's exit predicate: true from a run's start until tail step (1) sets the
sticky stopped status — after the final snapshot, so a body observing false
may still read a complete final world. The author's loop obligation is to
check it between blocking points (§11.6).

Every handle primitive that a loop pass touches — this one, `latest`,
`stage!`, `wait_next_snapshot`, `gather`, `report!` — stores the liveness
heartbeat on its way through (§11.8, §12.2): the framework observes activity
without owning the loop body, and there is no separate liveness channel to
remember to feed.
"""
running(h::DeviceHandle) = (_beat!(h.diag); !(@atomic h.ctl.stopped))

"""
    stop!(handle)

Request a control-plane stop (§12.1): sets the stop word from any task; the
loop observes it at the next frame top, completes that boundary, publishes,
and enters the tail (§12.4). Idempotent, and inert while already stopped.
"""
stop!(h::DeviceHandle) = (@atomic h.ctl.stop_requested = true; nothing)

"""
    binding(handle)

The attachment's binding (§11.6), for the author-owned loop idiom —
`stage!(handle, map_input(datum, binding(handle))...)`. The framework never
calls `map_input`/`map_output` itself: they are conventions of the loop
body, and the handle carrying the binding is what keeps the device struct
free of its per-deployment configuration (the binding stays an `attach!`
argument, never a device field).
"""
binding(h::DeviceHandle) = h.b

"""
    latest(handle)

The handle's primitive read (§11.6): acquire-load the most recently published
snapshot — exactly `latest(sim)`, through the capability the handle carries.
"""
latest(h::DeviceHandle) = (_beat!(h.diag); @atomic :acquire h.published.latest)

"""
    stage!(handle, "face" => value, ...)

The handle's stage primitive (§11.6, §11.4): the same compiled shim, CAS
merge and per-face newest-wins as every writer's, against this attachment's
claim set — every check at staging, an out-of-claim face discarded under
`OutOfClaimEntry` naming the incumbent, the rest of the batch standing. From
any task, at any wall-clock moment; the batch lands at the top of the next
frame `run!` advances.
"""
function stage!(h::DeviceHandle, pairs::Pair...)
    _beat!(h.diag)
    batch = _normalize(h.writer, pairs, h.plane.claimedby; device = h.who)
    batch === nothing || _stage!(h.writer, batch)
    nothing
end

"""
    gather(handle, snap)

The output side's read (§11.2, §11.6): run the attachment's compiled gather —
`reads(b)`, resolved at attach — over a snapshot, returning the labeled
NamedTuple `map_output` receives. The loop idiom is
`send(dev.socket, map_output(gather(handle, snap), binding(handle)))`, on the
device's own task, against the snapshot §12.3's wait handed it: the compiled
addresses read the frozen store, so no name is resolved per datum and nothing
here touches the running loop. On a handle whose binding declares no output
side the call is a contract misuse, and throws by name.
"""
function gather(h::DeviceHandle, s::Snapshot)
    _beat!(h.diag)
    h.gatherer === nothing && error(
        "$(h.who)'s binding declares no output side — `gather` serves the compiled " *
        "`reads` enumeration (§11.6)")
    _gather(h.gatherer, s)
end

"""
    report!(handle, MalformedDatum(cause))

The bad-datum channel (§11.6, §11.8): the single-writer entry point into the
device's own diagnostic cell, from the author's loop body after catching its
own parser error — catch, stage nothing, report, continue. The tolerance is
bounded by the cell (`DIAG_RING` retained values per frame, the excess a
count), so a peer flooding malformed datagrams costs sixteen retained values
and an integer increment, whatever it does, and no writer can starve another
— the cells are disjoint. This is not a general user-diagnostics channel:
`MalformedDatum` is the one thing it carries, and any exception the loop does
*not* classify as a bad datum propagates to the wrapper as `DeviceCrash`.
The loop drains the cell at frame top and at the run's end, warning each
retained value device-attributed and folding the counts into the per-run
totals `diagnostics(sim)` reads (the framework-status stand-in, README).
"""
report!(h::DeviceHandle, d::MalformedDatum) = (_beat!(h.diag); _report!(h.diag, d))

"""
    wait_next_snapshot(handle)

§12.3's next-snapshot wait: block until a boundary this waiter has not seen
is published, or until the run stops, and return the latest snapshot — at
least the boundary whose publication woke the wait, newest-wins beyond it (a
slow consumer skips boundaries and always receives the current world). The
predicate loop is the canonical idiom: the condition carries no facts, and
the two that matter — the counter and the stopped status — are re-tested at
every wake, which is what makes shutdown work: tail step (2) wakes every
waiter and the predicate routes it out. After a stop return the author's
loop re-checks `running(handle)`, exactly as after any blocking call.
"""
function wait_next_snapshot(h::DeviceHandle)
    _beat!(h.diag)
    ctl = h.ctl
    lock(ctl.cond)
    try
        while ctl.counter <= h.last_seen && !(@atomic ctl.stopped)
            wait(ctl.cond)
        end
        h.last_seen = ctl.counter
    finally
        unlock(ctl.cond)
    end
    latest(h)
end

# --- the wrapper and the run's bracket (§11.6, §12.4) --------------------------

# The guarded release: `shutdown!` is guaranteed on every exit path, and a
# throw out of it must not wreck the bracket or the tail around it (§11.6).
function _shutdown!(e::RosterEntry)
    try
        shutdown!(e.dev)
    catch err
        @warn "shutdown! of $(_who(e)) threw; its resources may leak (§11.6)" #=
            =# exception = (err, catch_backtrace())
    end
    nothing
end

"""
The framework wrapper (§11.6): the author owns the loop body, the framework
owns the bracket. Every exit path — voluntary return, a crash, the
stop-drained predicate — runs `shutdown!` and then consults the attachment's
`should_abort`; a crash is caught here and reported as `DeviceCrash`, the
run continuing with the device's task absent and its claims held to run end
(§12.4(6)). A `needs_calling_task` device runs this identical wrapper inline
on the calling task: the invocation site is its only difference (§11.1).
With §11.8 absent there is no diagnostic cell to mark dead — death is the
task having ended.
"""
function _wrap(e::RosterEntry)
    try
        loop(e.dev, e.handle)
    catch err
        @warn "DeviceCrash: $(_who(e))'s task threw; the run continues with its " *
              "task absent and its claims held to run end (§12.4)" #=
            =# exception = (err, catch_backtrace())
    finally
        _shutdown!(e)
        e.should_abort && stop!(e.handle)
    end
    nothing
end

"""
§12.4's pre-spawn initialization bracket, a step of the shutdown protocol
taken at the top of a run: `init!` once per roster entry, in attachment
order, on the calling task, each call in its own bracket. A throw goes
straight back to `shutdown!` — which is why `shutdown!` owes tolerance of a
partially initialized device — is reported as the ordinary `DeviceCrash`,
and spawns no task: the device is dead from the run's first frame, its
claims persisting to run end and the orphaned slots holding their values.
With `should_abort` set the failure requests a stop, already pending when
the loop would start: the run advances zero frames and ends through the same
tail, every remaining entry still getting its `init!`/`shutdown!` pair
uniformly. Returns the live entries, from which §11.1's topology is derived
— derived *after* initialization, never from the roster alone.
"""
function _init_devices!(sim)
    live = RosterEntry[]
    for e in sim.plane.roster
        ok = try
            init!(e.dev)
            true
        catch err
            _shutdown!(e)
            @warn "DeviceCrash: $(_who(e))'s init! threw; no task is spawned and " *
                  "its claims persist to run end (§12.4)" #=
                =# exception = (err, catch_backtrace())
            e.should_abort && stop!(e.handle)
            false
        end
        ok && push!(live, e)
    end
    live
end

# Spawn the wrapper, one run-scoped task per entry (§11.1): inside `run!`
# only, never at `attach!` — a task exists only once the run it serves does.
# The §12.3 registers are refreshed first, on the calling task.
function _spawn!(entries::Vector{RosterEntry})
    for e in entries
        e.handle.last_seen = e.handle.ctl.counter
    end
    [Threads.@spawn _wrap(e) for e in entries]
end

# Tail steps (1)'s close and (2) (§12.4): the final snapshot is whatever the
# loop last published; the sticky status flips only after it, and the notify
# under the lock wakes every §12.3 waiter, whose predicate routes it out.
# Idempotent, and run even when the loop leaves by a throw — the §13.6 catch
# path is absent (README), but device tasks must never be left parked.
function _finish!(sim)
    ctl = sim.control
    @atomic ctl.stopped = true
    lock(ctl.cond)
    try
        notify(ctl.cond)
    finally
        unlock(ctl.cond)
    end
    nothing
end

"""
Tail steps (3)–(5) (§12.4), on the calling task, after `_finish!` has set
the sticky status and woken the waits: `unblock!` per spawned entry — the
override's own blocking call returns; a throw out of the hook is warned and
the tail proceeds — then the join under one shared `join_timeout` deadline,
D-198's one patience for the whole tail. A task exceeding what remains of
the deadline is reported by name under `DeviceJoinTimeout` and abandoned
rather than left to hang `run!`. The calling-task device sits outside the
join: nothing can abandon the task `run!` stands on.
"""
function _tail!(sim, entries::Vector{RosterEntry}, tasks::Vector{Task})
    for e in entries
        try
            unblock!(e.dev)
        catch err
            @warn "unblock! of $(_who(e)) threw; its task can now exit only " *
                  "through the join timeout (§12.4)" exception = (err, catch_backtrace())
        end
    end
    deadline = time() + sim.join_timeout
    for (e, t) in zip(entries, tasks)
        remaining = deadline - time()
        joined = istaskdone(t) ||
            (remaining > 0 &&
             timedwait(() -> istaskdone(t), remaining; pollint = min(0.01, remaining)) === :ok)
        joined || @warn "DeviceJoinTimeout: $(_who(e)) did not return within " *
                        "$(sim.join_timeout) s of the stop; its task is abandoned (§12.4)"
    end
    nothing
end
