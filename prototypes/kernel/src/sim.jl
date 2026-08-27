# The simulation object with its deployment binding (§9.1, §9.2), the stepper
# seam (§10.2) and the phase-body accessor (§9.7). The framework owns the loop;
# the one delegated operation is "advance the continuous state from `t` by `h`".

"""
The active advance's effective stop policy (§13.5), bound at each `run!` or
`step!` entry from the constructor defaults and the per-run overrides: the
named stop faces in declaration order, their compiled root-cell addresses, and
`hit` — the seam through which a `t*` boundary's stop observation reaches the
loop (localize.jl publishes mid-frame, so the frame reports the holding face
here and abandons its remainder: the `t*` snapshot is the final one).
"""
mutable struct RunPolicy
    faces::Vector{Symbol}
    addrs::Vector{Any}
    hit::Union{Nothing,Symbol}
end

struct Simulation{T,S,B,CL,EV,M}
    store::S
    xbuf::Vector{T}
    ẋbuf::Vector{T}
    clock::CL
    bodies::B
    events::EV                    # the compiled event set with its §10.6 registers
    layout::Layout
    flat::Flat
    h::Float64                    # the continuous step, bound at deployment
    n::Int                        # steps per base tick: Δt_base = n·h (§10.5)
    Δt_base::Float64
    firing_budget::Int            # per-event firings per boundary (§10.6)
    localization_tol::Float64     # relative bracket-width stop (§10.4)
    localization_budget::Int      # t* boundaries permitted per frame (§10.4)
    join_timeout::Float64         # the shutdown tail's join cap, seconds of wall clock (§12.4)
    t_end::Union{Nothing,Float64} # the run's default clock bound (§13.5), run!-overridable
    stop_on::Vector{Symbol}       # the default stop faces (§13.5), run!-overridable
    stop_addrs::Vector{Any}       # their compiled root-cell addresses, validated at binding
    policy::RunPolicy             # the active advance's effective policy, bound per entry
    has_localized::Bool           # any localized event compiled in: the frame loop's fast-path key
    sched::Vector{@NamedTuple{path::String, D::Int, Φ::Int, Δt::Float64}}   # the bound schedule (§9.2)
    sstores::Vector{Any}          # discrete state stores, by component index
    mstores::Vector{Any}          # mode stores, by component index
    stepper::M                    # the seam's backend (§10.2), owning its own scratch
    xnext::Vector{T}              # the retained arrival pair (§10.4): xₙ₊₁ saved before trials clobber
    ẋnext::Vector{T}              # the buffer, ẋₙ₊₁ paid only past a validated trigger
    plane::DataPlane              # the §11.3 roster and harness register, mutable at stopped-sim points
    published::Published          # §11.2's `@atomic latest` holder
    control::Control              # §12.1's stop word, §12.4's sticky status, §12.3's wait (devices.jl)
    log::SnapshotLog              # §11.2's retained snapshots, loop-task bookkeeping
    loop_diag::DiagCell           # the loop's own diagnostic cell (§11.8): the budget degradations
    loop_acct::WriterAccount      # and the account behind it
end

"""
    Simulation(build::Build, T = Float64; h, n = 1, Δt_base = nothing,
               method = RK4, firing_budget = 4, localization_tol = 1e-6,
               localization_budget = 8, join_timeout = 5.0,
               t_end = nothing, stop_on = (), log = true,
               log_every = 1, log_max = 65536, chunk_size = 16)
    Simulation(root, T = Float64; …)

Deployment binding at construction (§9.1, §9.2): `Δt_base` from one of three
cross-validated sources — explicit (a `Rational` or `Period`/`Hz` value), the
`n·h` product (the default path), or GCD derivation over the anchors' constraint
pool, requested as `Δt_base = :derive` and permitted only with every discrete
component anchored. The scalar picks the activation the entries compile over —
the nominal one directly, any other via `activation(b, T)`'s cached Stratum-C
re-run (§9.4). The convenience form is *defined as* `Simulation(build(root), T;
…)`; entry compilation lives behind the binding because `Δt`, `D` and `Φ` are
entry data, and one `Build` backs many `Simulation`s.

`method` selects the integration backend across the stepper seam (§10.2): a
stepper type — `RK4`, the default, or `Heun` — materialized against the flat
buffer at binding like every other deployment product. The method is
trajectory-determining and grid-independent, exactly like the keywords below;
nothing outside the backend's own struct knows which one ran.

`firing_budget` is §10.6's deployment keyword: how many times each declared
event may fire at one boundary, an integer ≥ 1 defaulting to 4 — a legitimate
re-enable is one or two firings deep, a toggling FSM pair chatters without
bound, and 4 separates them without ever binding on a healthy model.

`localization_tol` and `localization_budget` are §10.4's: the relative bracket
width at which root-finding stops (positive, defaulting to 1e-6 — the event
time can only ever be as accurate as the `O(h⁴)` interpolant, so anything
tighter buys nothing while every trial evaluation costs a full sweep), and how
many localizations one frame admits (an integer ≥ 1 defaulting to 8 — a
legitimate multi-event frame needs three or four, chattering needs tens). All
three are trajectory-determining and grid-independent: they stand beside `h`
and `n`, validated here with their siblings, and enter none of the grid
arithmetic.

`join_timeout` is §12.4's: the shutdown tail's join cap in seconds of wall
clock — a positive real defaulting to 5, generous for GUI teardown and socket
closes, short enough that an abandoned join reads as a diagnosed timeout
rather than a hang (D-198). It is the one *operational* keyword: never
trajectory-determining, because the trajectory has ended at the final
snapshot before any join begins, yet not a view policy either — it tunes the
tail's wall-clock patience, nothing more.

`t_end` and `stop_on` are §13.5's termination policy, the `Simulation`'s
defaults with `run!`'s keywords as the per-run overrides: the clock bound — a
finite real ≥ 0, taken to the nearest frame top, with no constructor default
of its own (a run needs one from *some* binding site) — and the model-declared
sibling beside it, a collection naming root-exported `Bool` **output** faces,
OR-combined, sampled at every published boundary. Validation runs here exactly
as at `run!`: an unknown face, a root input slot, or a non-`Bool` face is
refused at whichever site names it. The default `()` is no stop faces — a run
to `t_end`.

`log`, `log_every` and `log_max` are §11.2's retention keywords: the plain
switch (default `true`), the keep-every-kth stride over published boundaries
(an integer ≥ 1, default 1) and the bound on retained snapshot references (an
integer ≥ 1 defaulting to 65536 = 2¹⁶ — about 22 minutes at 50 Hz and full
density before anything is dropped — with `Inf` the explicit opt-out). Unlike
every keyword above they are **view policies, never trajectory-determining**:
two deployments differing only here produce bitwise-identical trajectories,
retention being reference bookkeeping over what publication already built.
"""
function Simulation(b::Build, ::Type{T} = Float64; h = nothing, n = nothing,
                    Δt_base = nothing, method = RK4, firing_budget = 4,
                    localization_tol = 1e-6, localization_budget = 8,
                    join_timeout = 5.0, t_end = nothing, stop_on = (),
                    log = true, log_every = 1,
                    log_max = 65536, chunk_size::Int = 16) where {T}
    method isa Type && method <: AbstractStepper ||
        throw(BuildError("method must be a stepper type — RK4 or Heun — got $method (§10.2)"))
    firing_budget isa Integer && firing_budget ≥ 1 ||
        throw(BuildError("firing_budget must be an integer ≥ 1, got $firing_budget (§10.6)"))
    localization_tol isa Real && localization_tol > 0 ||
        throw(BuildError("localization_tol must be a positive real, got $localization_tol (§10.4)"))
    localization_budget isa Integer && localization_budget ≥ 1 ||
        throw(BuildError("localization_budget must be an integer ≥ 1, got $localization_budget (§10.4)"))
    join_timeout isa Real && join_timeout > 0 ||
        throw(BuildError("join_timeout must be a positive real — the shutdown tail's " *
                         "join cap in seconds of wall clock — got $join_timeout (§12.4)"))
    log isa Bool ||
        throw(BuildError("log must be true or false — the retention switch — got $log (§11.2)"))
    log_every isa Integer && log_every ≥ 1 ||
        throw(BuildError("log_every must be an integer ≥ 1, got $log_every (§11.2)"))
    (log_max isa Integer && log_max ≥ 1) || log_max === Inf ||
        throw(BuildError("log_max must be an integer ≥ 1, or Inf as the explicit " *
                         "opt-out, got $log_max (§11.2)"))
    t_end === nothing || _t_bound(t_end)
    act = activation(b, T)
    (stop_faces, stop_addrs) = _stop_faces(act.layout, stop_on)
    bound = bind_schedule(b, h, n, Δt_base)
    c = compile(b, act, bound.D, bound.Φ, bound.Δt; chunk_size)
    stepper = method(T, length(c.xbuf))
    Simulation{T,typeof(c.store),typeof(c.bodies),typeof(c.clock),typeof(c.events),
               typeof(stepper)}(
        c.store, c.xbuf, c.ẋbuf, c.clock, c.bodies, c.events, act.layout, b.flat,
        bound.h, bound.n, bound.Δt_base, Int(firing_budget), Float64(localization_tol),
        Int(localization_budget), Float64(join_timeout),
        t_end === nothing ? nothing : Float64(t_end), stop_faces, stop_addrs,
        RunPolicy(Symbol[], Any[], nothing), any(c.events.localized), bound.sched,
        c.sstores, c.mstores, stepper, zeros(T, length(c.xbuf)), zeros(T, length(c.xbuf)),
        DataPlane(act.layout, c.store), Published(nothing), Control(),
        SnapshotLog(log, Int(log_every), log_max === Inf ? typemax(Int) : Int(log_max)),
        DiagCell(EMPTY_DIAG), WriterAccount())
end

Simulation(root::AbstractComponent, ::Type{T} = Float64; kw...) where {T} =
    Simulation(build(root), T; kw...)

# §13.5's clock bound, validated identically at both binding sites.
_t_bound(t) = (t isa Real && isfinite(t) && t ≥ 0) ? Float64(t) : throw(BuildError(
    "t_end must be a finite real ≥ 0 — the run's clock bound, taken to the " *
    "nearest frame top — got $t (§13.5)"))

# §13.5's stop-face validation and compilation, run identically at both binding
# sites — the constructor's default and `run!`'s override: each name must be a
# root-exported Bool *output* face. Duplicates collapse; the order kept is the
# declaration's, which is the order the first-holding face is reported in.
function _stop_faces(layout::Layout, stop_on)
    faces, addrs = Symbol[], Any[]
    slot_names = Symbol[f for (f, _) in layout.slots]
    for f in stop_on
        s = Symbol(f)
        haskey(layout.addr, ("", s)) || throw(BuildError(
            "stop_on names `$s`, which is no root face — a stop face is a " *
            "root-exported Bool output face (§13.5)"))
        s in slot_names && throw(BuildError(
            "stop_on names `$s`, a root input slot — a stop face is a root-exported " *
            "*output*: the model detects and exports the condition, and the " *
            "deployment names it (§13.5)"))
        a = layout.addr[("", s)]
        _port_type(a) === Bool || throw(BuildError(
            "stop_on names `$s`, whose declared type is $(_port_type(a)) — stop " *
            "faces are Bool, OR-combined (§13.5)"))
        s in faces || (push!(faces, s); push!(addrs, a))
    end
    (faces, addrs)
end

"""
    lifecycle(sim)

§12.6's five-state lifecycle: `:built` — stores allocated, boundary zero never
run; `:initialized` — boundary-consistent and ready to advance, the state
`init!` establishes and a completed `step!` returns to; `:running` — `run!` or
`step!` holds the loop, and the §11.3 freeze with it; and the two terminal
states, `:stopped` and `:errored` (§13.6). Readable from any task.
"""
lifecycle(sim::Simulation) = @atomic :acquire sim.control.lifecycle

"""
    termination(sim)

§13.5's termination record (devices.jl) — the run's outcome, `nothing` unless
the lifecycle is terminal: which source ended the run (`:t_end`, `:stop_on`
with the face named, `:control_stop`, or `:error` with the cause retained),
and the final snapshot's boundary time.
"""
function termination(sim::Simulation)
    lc = @atomic :acquire sim.control.lifecycle
    lc === :stopped || lc === :errored ? sim.control.termination : nothing
end

# The final snapshot's boundary time, for the record: after `init!` a snapshot
# always exists; the NaN arm covers a failure before boundary zero ever ran.
_t_latest(sim::Simulation) =
    (s = latest(sim); s === nothing ? NaN : s.t)

# §13.5's sampling read, after every publication: the named faces off the
# just-published snapshot, first holding face wins, in declaration order.
function _stop_hit(sim::Simulation, pol::RunPolicy)
    isempty(pol.addrs) && return nothing
    snap = latest(sim)
    for i in eachindex(pol.addrs)
        gather(snap.store, pol.addrs[i]) === true && return pol.faces[i]
    end
    nothing
end

# The shared entry gate of the two advance entries (§12.6): only `:initialized`
# admits an advance, and each refusal names its own way out.
function _assert_advanceable(sim::Simulation, op::String)
    lc = @atomic sim.control.lifecycle
    lc === :initialized && return nothing
    lc === :built && throw(BuildError(
        "`$op` before `init!`: boundary zero has not run — `init!` is mandatory (§12.6)"))
    lc === :running && throw(BuildError(
        "ServiceLifecycle: the simulation is already running (§12.6)"))
    lc === :stopped && throw(BuildError(
        "`$op` on a stopped simulation: re-running is the `stopped → init! → run!` " *
        "cycle — `init!` re-runs boundary zero and opens a fresh trajectory (§12.6)"))
    throw(BuildError(
        "ServiceLifecycle: this simulation ended `errored` — terminally stopped, " *
        "never resumable; reproduction is trace replay, absent here (§13.6)"))
end

"""
    phase_bodies(sim)

The compiled bodies of the nominal activation, bound over this simulation's own
buffers — **these are the bodies the loop runs**, not re-derivations, which is
what makes the §7.5 measurement honest. The roster is fixed and total: a model
with no discrete components still gets `ticks`, empty, compiling to a no-op
whose `@ballocated` assertion passes vacuously, so consumers iterate uniformly
with no per-model branching.
"""
phase_bodies(sim::Simulation) = sim.bodies

# --- evaluation ---------------------------------------------------------------

"""
One RHS evaluation: *evaluating the RHS means running the sweep* (§5.3). The
interior variant of each sweep block, then the `f` block against the complete
fresh table. Leaves `ẋbuf` holding the derivative of whatever `xbuf` holds.
"""
@inline function evaluate!(sim::Simulation)
    sim.bodies.sweep_1()
    sim.bodies.sweep_2()
    sim.bodies.rhs()
    nothing
end

"""
The boundary macro-sequence at a base tick, final form (§5.3, §10.6):

> integrate → project → [sweep → guards → handlers] iterated to quiescence
> (under the firing budget) → all due `g` updates

Integration has just written the state, so projection runs first — between the
write and its decode; the event phase then iterates with the due set fixed for
the whole boundary, and the due updates run last, after quiescence, reading
post-transition values off the settled table. Output stages before updates, so
a discrete component's cells carry `y[k]` computed from `s[k]` while `g`
produces `s[k+1]` — the sampled-data recursion, ordered by construction rather
than by convention.
"""
@inline function boundary!(sim::Simulation, tick::Int)
    _projects!(sim.events)
    event_phase!(sim, tick)
    sim.bodies.ticks(tick)
    nothing
end

"""
The macro-sequence at a step boundary that is *not* a base tick: the tick
counter has not advanced, so the due set is empty — and emptiness is arity
selection, never a sentinel index failing every gate (§10.5, D-185). The
zero-arg interior bodies are exactly the boundary walk with every discrete
entry gated out, so consistency is restored and nothing discrete can move;
projection and the event phase run in full — every step boundary is a boundary
(§10.4).
"""
@inline function offtick_boundary!(sim::Simulation)
    _projects!(sim.events)
    event_phase!(sim, nothing)
    nothing
end

# One iteration round's sweep: the whole gated schedule (§10.6), in the due-set
# arity the boundary fixed — never the update laws, which wait for quiescence.
@inline _round!(sim::Simulation, tick::Int) =
    (sim.bodies.sweep_1(tick); sim.bodies.sweep_2(tick); nothing)
@inline _round!(sim::Simulation, ::Nothing) =
    (sim.bodies.sweep_1(); sim.bodies.sweep_2(); nothing)

"""
The event phase (§10.6): [sweep → guards → handlers] iterated to quiescence,
the fixed point where a round of handlers fires nothing. The first round always
runs — it is the boundary sweep that restores table consistency — so a model
with no events degenerates to exactly that sweep, with no register work at all.

An event fires in a round iff its predicate is observed holding, the sample
observed before it was not-holding, and its firing count for this boundary is
below the budget; at most one event fires per component per round, the first
eligible in declaration order — priority with re-decision. The one register
subtlety is D-191's: an eligible-but-blocked event's last-observed sample is
*not* overwritten, so the edge it presented stands into the next round. Budget
exhaustion degrades rather than throwing: the event's further edges at this
boundary are lost under a `FiringBudget` warning, at most once per event per
boundary, while every other event iterates untouched. At quiescence the prior
is updated unconditionally from the final samples — every prior an honest
observation of a settled boundary.
"""
function event_phase!(sim::Simulation, tick)
    es = sim.events
    _round!(sim, tick)
    n = length(es.prior)
    n == 0 && return nothing
    copyto!(es.last, es.prior)
    fill!(es.count, 0)
    fill!(es.warned, false)
    budget = sim.firing_budget
    while true
        _guards!(es)
        fill!(es.comp_fired, false)
        any_fired = false
        for i in 1:n
            edge = !es.last[i] && es.now[i]
            eligible = edge && es.count[i] < budget
            if edge && !eligible && !es.warned[i]
                es.warned[i] = true       # at most one report per event per boundary
                (path, name) = es.names[i]
                _report!(sim.loop_diag,   # the loop's own cell (§11.8): folded at the next frame top
                         FiringBudget(path, name, Float64(sim.clock.t), budget, es.count[i]))
            end
            firing = eligible && !es.comp_fired[es.owner[i]]
            es.fire[i] = firing
            if firing
                es.comp_fired[es.owner[i]] = true
                es.count[i] += 1
                any_fired = true
            end
            # A blocked edge stays unconsumed (D-191): only the eligible-but-
            # blocked sample stands; every other event takes this round's.
            (eligible && !firing) || (es.last[i] = es.now[i])
        end
        any_fired || break
        _fire!(es)
        _round!(sim, tick)
    end
    copyto!(es.prior, es.last)
    nothing
end

# --- the stepper seam, framework side (§10.2) ----------------------------------

"""
    step!(sim, h)

Advance the continuous state from `t` by `h`, delegated across the stepper
seam to whichever backend the deployment bound (stepper.jl). The seam is never
entered empty: with no continuous state the step degenerates to advancing `t`,
the backend is simply not called, and no backend contract has to say what it
would do at N = 0.
"""
@inline function step!(sim::Simulation, h)
    isempty(sim.xbuf) ? (sim.clock.t += h) : step!(sim.stepper, sim, h)
    nothing
end

"""
Initialize: state at its declared values, table consistent, clock at `t₀`.

Boundary zero is an ordinary boundary with an empty integrate (§10.5, §14.5):
the gate at index 0 admits exactly the components with `Φ = 0`, implemented by
nothing. An offset component holds its probe-populated cells until its first
tick at `Φ·Δt_base`.

Boundary zero also establishes every event prior as not-holding (§10.6), so a
predicate already holding in the authored state fires at `t₀` — derived, not
asserted — and a warm restart (`init!` again) resets all three registers from
scratch: such predicates fire again at the new `t₀`. The log resets with them
(§11.2): a warm restart is a new trajectory, and its boundary zero lands as a
fresh first endpoint.

`init!` is §12.6's one door into `initialized`, and it is mandatory: `run!`
and `step!` refuse a simulation whose boundary zero has not completed. It
opens the fresh trajectory wholesale — the stop word, the §13.5 termination
record *and every staged batch still in a staging cell* clear with the
registers (§12.6: no stale batch survives to clobber the boundary zero it
predates — the pre-run register is `init!` → `stage!` → `run!`, the batch
then waiting for the first frame top as §11.4 says). The diagnostic cells are
deliberately not cleared: a rejection recorded while stopped is a fact about
what happened, not a stale input, and it surfaces in the next run's first
status (§11.8). `init!` is itself a stopped-sim operation: refused while
`running`, and refused on an `errored` simulation, which is terminally
stopped (§13.6) — reproduction is trace replay, not resurrection.
"""
function init!(sim::Simulation{T}; t₀::T = zero(T)) where {T}
    ctl = sim.control
    lc = @atomic ctl.lifecycle
    lc === :running && throw(BuildError(
        "ServiceLifecycle: `init!` is a stopped-sim operation and the simulation " *
        "is running (§11.3, §12.6)"))
    lc === :errored && throw(BuildError(
        "ServiceLifecycle: this simulation ended `errored` — terminally stopped, " *
        "never re-initialized; reproduction is trace replay, absent here (§13.6)"))
    sim.clock.t = t₀
    sim.clock.t₀ = t₀
    sim.clock.step = 0
    fill!(sim.events.prior, false)
    _reset!(sim.log)
    _reset_accounts!(sim)         # a new trajectory opens a fresh account (§11.8)
    for e in sim.plane.roster     # §12.6: no staged batch survives into the
        @atomic e.writer.cell.pending = nothing   # trajectory it predates
    end
    @atomic sim.plane.harness.cell.pending = nothing
    @atomic ctl.stop_requested = false
    ctl.termination = nothing
    boundary!(sim, 0)
    publish!(sim)                 # the boundary-zero snapshot (§11.2, §14.5)
    @atomic :release ctl.lifecycle = :initialized
    nothing
end

"""
    run!(sim; t_end = nothing, stop_on = nothing)

Advance one frame at a time, each frame §11.1's anatomy — drain, integrate,
boundary sequence, publication — under §11.1's task topology and §12.4's
bracket and tail, until a §13.5 termination source ends the run: `t_end`'s
frame reached, a named `stop_on` face observed holding in a published
snapshot, or a control-plane stop observed at frame top. Both keywords bind
for **this run only**, the constructor's values standing where they are not
given (§13.5) — the override validated exactly as the constructor validates
its default — and a run needs a clock bound from one of the two sites. Only
an `initialized` simulation runs (`init!` is mandatory, §12.6), and the run
leaves it terminally `stopped` — the §13.5 record readable through
`termination(sim)` — or `errored` on a loop-side failure (§13.6): the failed
boundary published nothing, so the last published snapshot is already the
promoted final one; the tail runs identically, the cause is retained on the
record, and `run!` rethrows after the tail completes (§13.4's synchronous
rule).

The run's shape, in order: the policy binds and the §11.3 freeze rises (the
lifecycle's `:running`, spanning the tail); the stop word is cleared (a fresh
run owes nothing to the last one's stop); the §12.4 init bracket runs per
roster entry on the calling task; the topology is derived from the *live*
entries — with a `needs_calling_task` holder among them the loop moves to a
spawned task and the calling task runs that device's loop body inline,
otherwise the loop runs here and one task is spawned per live entry — the
loop advances until a termination source fires; and the tail closes the run
(devices.jl): sticky status after the final snapshot, waits woken,
`unblock!`, the join under `join_timeout`. Either way `run!` blocks its
caller until the run ends; what varies is what the calling task spends the
run doing (§11.1).

Inside the loop, `frame!` carries each grid step `[tₖ₋₁, tₖ]` through the
§10.4 localization loop, firing any `t*` boundaries it brackets on the way;
the frame-top boundary (§10.3) is a base tick every `n` frames, where the
gate reads the tick index, and the empty-due-set boundary in between. The
drain runs at the frame top only, never at a `t*` boundary (§10.4), while
publication follows *every* boundary sequence (§11.2) — the frame top's here,
a `t*` boundary's inside the frame loop, before integration resumes — and
every publication is a stop-face sampling point (§13.5), a `t*` hit ending
the run with the `t*` snapshot final. The grid is driven by the step counter,
so `t_end` is taken to the nearest frame top.
"""
function run!(sim::Simulation; t_end = nothing, stop_on = nothing)
    plane, ctl = sim.plane, sim.control
    _assert_advanceable(sim, "run!")
    te = t_end === nothing ? sim.t_end : _t_bound(t_end)
    te === nothing && throw(BuildError(
        "run! has no clock bound: `t_end` was given neither at the constructor " *
        "nor here — the constructor value is the default and the run! keyword " *
        "the per-run override (§13.5)"))
    (faces, addrs) = stop_on === nothing ? (sim.stop_on, sim.stop_addrs) :
                                           _stop_faces(sim.layout, stop_on)
    target = round(Int, te / sim.h)
    pol = sim.policy
    pol.faces, pol.addrs, pol.hit = faces, addrs, nothing
    @atomic :release ctl.lifecycle = :running   # the §11.3 freeze: the roster is fixed for the run
    term = nothing
    try
        @atomic ctl.stop_requested = false
        _reset_accounts!(sim)                 # §11.8: totals count since the run began
        live = _init_devices!(sim)            # §12.4's pre-spawn bracket, attachment order
        @atomic ctl.stopped = false
        ct = findfirst(e -> needs_calling_task(e.dev), live)
        if ct === nothing                     # the unattended register (§11.1)
            tasks = _spawn!(live)
            _register_tasks!(plane, live, tasks)
            try
                term = _advance!(sim, pol, typemax(Int), target)[1]
            finally
                _finish!(sim)                 # tail (1)–(2), even off a loop-side throw
                _tail!(sim, live, tasks)      # tail (3)–(5)
            end
        else                                  # the loop is the movable piece (§11.1)
            e_ct = live[ct]
            others = [live[i] for i in eachindex(live) if i != ct]
            tasks = _spawn!(others)
            _register_tasks!(plane, others, tasks)
            plane.run_tasks[e_ct.id] = current_task()   # the inline body's task (§11.1)
            e_ct.handle.last_seen = ctl.counter
            loop_task = Threads.@spawn try
                _advance!(sim, pol, typemax(Int), target)[1]
            finally
                _finish!(sim)                 # the spawned loop wakes the inline body too
            end
            _wrap(e_ct)                       # the identical wrapper, inline (§11.6)
            try
                term = fetch(loop_task)       # run! blocks until the run ends (§11.1)
            finally
                _tail!(sim, others, tasks)    # the calling-task device sits outside the join
            end
        end
    catch err
        # §13.6's abnormal entry: the failed boundary is discarded by
        # construction — publication is a boundary's last act, so it published
        # nothing and the previous snapshot is already final. The record
        # retains the cause raw (§13.4's wrap is absent, README), unwrapped
        # from the spawned loop's task failure where the topology moved it.
        ctl.termination = Termination(:error, nothing, _t_latest(sim),
                                      err isa TaskFailedException ? err.task.exception : err)
        @atomic :release ctl.lifecycle = :errored
        rethrow()
    finally
        @atomic ctl.stopped = true
        _sweep_tail!(sim)                     # the run's last take (§11.8): what landed past
        empty!(plane.run_tasks)               # the final frame top — presented, never published
        if (@atomic ctl.lifecycle) === :running   # (devices.jl); tasks are run-scoped equipment
            ctl.termination = term
            @atomic :release ctl.lifecycle = :stopped
        end
    end
    nothing
end

# The per-run reset (§11.8): every writer's account starts the run — and, from
# init!, the trajectory — at zero. The cells are deliberately not touched: a
# batch reported or staged while stopped waits for the first frame top's
# drain, exactly as a staged input batch waits (§11.4).
function _reset_accounts!(sim::Simulation)
    for e in sim.plane.roster
        _reset!(e.acct)
    end
    _reset!(sim.plane.harness_acct)
    _reset!(sim.loop_acct)
    nothing
end

_register_tasks!(plane::DataPlane, entries::Vector{RosterEntry}, tasks::Vector{Task}) =
    (for (e, t) in zip(entries, tasks); plane.run_tasks[e.id] = t; end; nothing)

# The frame loop, shared by both advance entries (§12.6: a stepped frame is
# bit-identical to a run frame because this is the same code) — returns
# `(termination, frames advanced)`, the termination `nothing` exactly when the
# frame budget `upto` ran out, which only `step!` binds finitely. The §13.5
# sources in consultation order: the stop word at frame top (§12.1 — the loop
# never stops mid-frame, so a stop observed here leaves the last published
# boundary as the final snapshot, §12.4(1)); `t_end`'s frame completed; and
# the stop faces at every publication — the entry check first (a boundary-zero
# or authored condition already terminal advances nothing, §13.5), then after
# each frame's own publications, where a mid-frame `t*` hit arrives through
# `pol.hit` with the frame's remainder already abandoned. With devices
# rostered every frame yields at least once (§12.2, the unpaced case): the
# explicit yield is the co-resident device tasks' scheduling slot.
function _advance!(sim::Simulation, pol::RunPolicy, upto::Int, t_end_frame::Int)
    plane, ctl = sim.plane, sim.control
    adv = 0
    face = _stop_hit(sim, pol)
    face === nothing ||
        return (Termination(:stop_on, face, _t_latest(sim), nothing), adv)
    while true
        (@atomic ctl.stop_requested) &&
            return (Termination(:control_stop, nothing, _t_latest(sim), nothing), adv)
        sim.clock.step < t_end_frame ||
            return (Termination(:t_end, nothing, _t_latest(sim), nothing), adv)
        sim.clock.step < upto || return (nothing, adv)
        isempty(plane.roster) || yield()
        drain!(sim)
        k = (sim.clock.step += 1)
        adv += 1
        frame!(sim, k)
        if pol.hit === nothing
            k % sim.n == 0 ? boundary!(sim, k ÷ sim.n) : offtick_boundary!(sim)
            publish!(sim)
            face = _stop_hit(sim, pol)
        else
            face = pol.hit        # a t* publication hit (§13.5): that snapshot is final
        end
        face === nothing ||
            return (Termination(:stop_on, face, _t_latest(sim), nothing), adv)
    end
end

"""
    step!(sim; frames = 1)
    step!(sim; t_plus)

§12.6's partial advance: advance whole frames synchronously through the
ordinary frame sequence — drain, integrate, boundaries, publication — and
return the number of frames **actually advanced**. A stepped simulation is
bit-identical to the same frames under `run!`: the two entries share the one
frame loop. `t_plus` is the duration spelling, mutually exclusive with
`frames`: whole frames until the boundary time first covers the duration.

A stepping session is deviceless by construction (§12.6): no task is spawned
and no bracket runs — device tasks are per-`run!` artifacts — while the
frame-top drain still runs, so `stage!(sim, …)` → `step!` → `latest(sim)` is
the advance-assert-advance register, and a batch staged into a rostered
device's cell while stopped is applied exactly as §11.4 says. Between calls
the simulation reports `initialized`, so `attach!` is legal there and `run!`
may follow, continuing from the current boundary.

Termination policy is honored throughout, from the `Simulation`'s own
defaults: `t_end` reached, a `stop_on` face holding, or a pending control
stop ends the run *inside* the call through the deviceless §12.4 tail,
leaves the simulation terminally `stopped` with the §13.5 record set, and
returns the frames advanced before the stop — fewer than requested, which is
how a harness detects the truncation without inspecting the clock. A
loop-side failure ends it `errored` exactly as under `run!` (§13.6).
"""
function step!(sim::Simulation; frames = nothing, t_plus = nothing)
    ctl = sim.control
    _assert_advanceable(sim, "step!")
    frames === nothing || t_plus === nothing || throw(BuildError(
        "step! takes `frames` or `t_plus`, not both — the count and the duration " *
        "are two spellings of one advance (§12.6)"))
    if t_plus === nothing
        nf = frames === nothing ? 1 : frames
        nf isa Integer && nf ≥ 1 || throw(BuildError(
            "frames must be an integer ≥ 1, got $nf (§12.6)"))
        nf = Int(nf)
    else
        t_plus isa Real && isfinite(t_plus) && t_plus > 0 || throw(BuildError(
            "t_plus must be a finite real > 0 — the duration spelling — got " *
            "$t_plus (§12.6)"))
        nf = max(1, ceil(Int, Float64(t_plus) / sim.h - 1e-9))
    end
    pol = sim.policy
    pol.faces, pol.addrs, pol.hit = sim.stop_on, sim.stop_addrs, nothing
    t_end_frame = sim.t_end === nothing ? typemax(Int) : round(Int, sim.t_end / sim.h)
    @atomic :release ctl.lifecycle = :running   # the freeze holds within the call
    term, adv = nothing, 0
    try
        (term, adv) = _advance!(sim, pol, sim.clock.step + nf, t_end_frame)
    catch err
        ctl.termination = Termination(:error, nothing, _t_latest(sim), err)
        @atomic :release ctl.lifecycle = :errored
        rethrow()
    finally
        lc = @atomic ctl.lifecycle
        if lc === :errored                    # §13.6, the stepped register: same tail,
            _finish!(sim)                     # deviceless — waits woken, accounts swept
            _sweep_tail!(sim)
        elseif term === nothing
            @atomic :release ctl.lifecycle = :initialized
        else                                  # a §13.5 source fired inside the call:
            _finish!(sim)                     # the deviceless §12.4 tail, then terminal
            _sweep_tail!(sim)
            ctl.termination = term
            @atomic :release ctl.lifecycle = :stopped
        end
    end
    adv
end

"""
    stop!(sim)

Request a control-plane stop from any task (§12.1) — the calling code's
spelling of the stop a device handle issues with `stop!(handle)`. The loop
observes the word at the next frame top, completes that boundary, publishes,
and enters the tail (§12.4). Idempotent; inert while stopped, the word being
cleared at the top of the next run.
"""
stop!(sim::Simulation) = (@atomic sim.control.stop_requested = true; nothing)

# --- the roster (§11.3): stopped-sim configuration -----------------------------

"""
    attach!(sim, dev::AbstractDevice, b::AbstractBinding;
            should_abort = false) -> DeviceHandle

Roster a device under a binding — a stopped-sim configuration operation
(`ServiceLifecycle` while running: the roster is frozen per run, pause
included). The binding's conformance check runs first (§11.6), then the
three-part admission in spec order (§11.3): identity — this instance already
rostered is `AlreadyAttached`, rebinding being spelled `detach!` then
`attach!` — affinity (`CallerTaskConflict`: at most one `needs_calling_task`
holder), and claims (`ClaimConflict`, which the identity check having run
first always makes two *distinct* devices). On the input side the claim is
staked from its source — the enumeration called once, or the unclaimed
complement computed at this instant and never recomputed, so attaching the
greedy claimant last is the idiom and a second greedy stakes the empty
remainder under an `EmptyGreedyClaim` warning (§11.6) — the entry's writer is
compiled over it, and the harness register's surface is recompiled to the
complement that remains, renormalizing any pending harness batch (§11.4). On
the output side `reads(b)` is called once, resolved against the build and
compiled to the one gather `gather(handle, snap)` runs (§11.2, §14.4) — a
binding that drifted from its model fails here, not with silent garbage on
the wire. An output-only binding stakes no claim: its write surface is empty,
and the harness register keeps every face.

`should_abort` is §11.6's per-attachment failure policy, never a device
property: set, the device's departure — its loop body returning, a crash, or
a failed `init!` — also requests a sim stop (§12.4(6)); clear, the run
continues with the device's task absent and its claims held to run end.

Returns the attachment's handle (devices.jl), carrying the entry's stable
device id and the device's capabilities — read, stage, control access — the
same object the wrapper passes to `loop(dev, handle)` on the device's task.
`attach!` never spawns: it registers, and the task appears at the next
`run!` (§11.1).
"""
function attach!(sim::Simulation, dev::AbstractDevice, b::AbstractBinding;
                 should_abort::Bool = false)
    plane = sim.plane
    assert_stopped(sim.control, "attach!")
    check_binding(b)
    for e in plane.roster                          # identity, before claims (§11.3)
        e.dev === dev && throw(BuildError(
            "AlreadyAttached: this $(typeof(dev)) instance is already rostered as " *
            "$(_who(e)) under $(typeof(e.binding)) — rebinding is spelled `detach!` " *
            "then `attach!` (§11.3)"))
    end
    if needs_calling_task(dev)                     # affinity: a single-slot resource
        i = findfirst(e -> needs_calling_task(e.dev), plane.roster)
        i === nothing || throw(BuildError(
            "CallerTaskConflict: $(typeof(dev)) declares `needs_calling_task`, and " *
            "$(_who(plane.roster[i])) already holds the calling task (§11.1, §11.3)"))
    end
    claim = is_input(b) ? _claim(plane, sim.layout, b) : Symbol[]
    for f in claim                                 # claims: face exclusivity
        haskey(plane.claimedby, f) && throw(BuildError(
            "ClaimConflict: `$f` is claimed by both $(typeof(dev)) and " *
            "$(plane.claimedby[f]) — one writer per slot at any time (§11.3)"))
    end
    rg = is_output(b) ?                            # the output side: reads → one gather,
        _compile_reads(sim.layout, reads(b), typeof(b)) : nothing   # resolved before admission commits
    id = plane.next_id                             # assigned on admission alone: a
    plane.next_id += 1                             # rejected attach consumes no id
    w = Writer(sim.layout, claim)
    diag = DiagCell(EMPTY_DIAG)                    # the device's diagnostic cell (§11.8)
    h = DeviceHandle(id, "device $id ($(typeof(dev)))", b, w, plane, sim.control,
                     sim.published, diag, rg, sim.control.counter)
    push!(plane.roster, RosterEntry(dev, b, id, w, _drain_thunk(sim.store, w),
                                    should_abort, diag, WriterAccount(), h))
    reclaim!(plane, sim.layout)
    is_greedy(b) && isempty(claim) &&
        @warn "EmptyGreedyClaim: device $id ($(typeof(dev))) under $(typeof(b)) staked " *
              "the empty remainder — every root input face is already claimed (§11.6)"
    h
end

"""
    detach!(sim, dev)

Release a rostered device — the same stopped-sim gate as `attach!`. The
claims are released and the harness register's surface regains them; a
pending undrained batch in the entry's cell is discarded with it, detach
being a deliberate reconfiguration, while the slots it fed hold their
last-drained values. The device id retires with the entry, never reused.
"""
function detach!(sim::Simulation, dev::AbstractDevice)
    plane = sim.plane
    assert_stopped(sim.control, "detach!")
    i = findfirst(e -> e.dev === dev, plane.roster)
    i === nothing && throw(BuildError(
        "this $(typeof(dev)) instance is not rostered — `detach!` releases an " *
        "existing attachment (§11.3)"))
    deleteat!(plane.roster, i)
    reclaim!(plane, sim.layout)
    nothing
end

# --- the data plane (§11): staging, the drain, publication ---------------------

"""
    stage!(sim, "face" => value, ...)          # the harness register (§11.3)

Stage a batch of root-slot writes from any task, at any wall-clock moment,
never touching a live slot (§11.1): the entries land in the writer's staging
cell under the one coalescing policy — CAS merge, newest wins per face
(§11.4). Untouched faces survive into the pending batch; re-staged faces take
the newest level. Every check runs here, on the writer's side, its findings
written into the harness register's diagnostic cell (§11.8): a face outside
the harness surface is discarded — `ClaimedFaceEntry` naming the incumbent
when a rostered claim covers the face (a rostered greedy claimant empties the
harness surface outright, D-192) and `OutOfClaimEntry` when nothing does — an
unconvertible value under `EntryTypeMismatch`, and the rest of the batch
stands. The staged batch is applied by the drain at the top
of the next frame `run!` advances. A device stages through the handle its
`attach!` returned — `stage!(handle, pairs...)`, devices.jl — the same shim
against its own claim set.
"""
function stage!(sim::Simulation, pairs::Pair...)
    plane = sim.plane
    w = plane.harness
    batch = _normalize(w, pairs, plane.claimedby, plane.harness_diag)
    batch === nothing || _stage!(w, batch)
    nothing
end

"""
The drain (§11.4): exactly one point in each frame, at its top, where the loop
takes each staged batch with one `atomicswap(cell, nothing)` — an indivisible
take, so there is no lost-write window — and applies it through the entry's
compiled scatter, masked-off positions skipped, every check long since spent
at staging.
Cells drain in attachment order, the harness register's last by convention:
with every surface disjoint the order is unobservable, so the rule exists to
make the record read the same way every time, not to arbitrate (§11.3). Never
at a `t*` boundary (§10.4). Between drains the loop owns its data exclusively,
and the frame's outcome is a pure function of the drained batches.
"""
function drain!(sim::Simulation)
    plane = sim.plane
    for e in plane.roster
        e.drain()
        _fold!(e.acct, e.diag)    # the diagnostic cells drain at the same point (§11.8):
    end                           # retained values into the pending delta, every
    plane.harness_drain()         # occurrence into the totals
    _fold!(plane.harness_acct, plane.harness_diag)
    _fold!(sim.loop_acct, sim.loop_diag)
    nothing
end

"""
Publish this boundary (§11.2): build the snapshot in private memory — one
buffer copy, the framework side of §7.5's scope — then a single release-store
to `latest`. Runs only after the boundary sequence completes, so every
published table is boundary-consistent (§10.3), and the copy is what makes the
binding rule hold: nothing reachable from a published snapshot is ever written
again. The framework status is built here too (§11.8): each writer's record
takes the account's pending delta — so the drain's `recent` rides exactly one
snapshot, the first published after the frame top — beside the totals copy,
the heartbeat acquire-load and the `task_state` read off the run's `Task`
handle (D-193), fresh at every publication. Behind the store the snapshot
enters the log (§11.2) — logging dissolves into publication, retention being
the only thing the log adds — and the §12.3 counter increments under its
lock, *after* the release-store: the normative order, so a waiter observing
the new count finds at least this boundary in `latest`. The counter is read
unlocked to stamp the snapshot's ordinal, the loop task being its only
writer.
"""
function publish!(sim::Simulation)
    ctl = sim.control
    snap = Snapshot(sim.clock.t, sim.clock.step, ctl.counter, capture(sim.store),
                    sim.layout, _status(sim))
    @atomic :release sim.published.latest = snap
    log!(sim.log, snap)
    lock(ctl.cond)
    try
        ctl.counter += 1
        notify(ctl.cond)
    finally
        unlock(ctl.cond)
    end
    nothing
end

# One writer's record (§11.8): the account's pending delta is *taken* — the
# status owns the vector, the account re-arms the shared empty — and the
# totals, isbits, copy by read. `hb` and `ts` are `nothing` for the two
# writers with no task of their own.
function _writer_status(who::String, a::WriterAccount, hb, ts)
    ws = WriterStatus(who, a.recent, a.suppressed, a.totals, hb, ts)
    a.recent = EMPTY_RECENT
    a.suppressed = KindCounts()
    ws
end

# The status assembly (§11.8, §11.2), on the publishing task: per-writer
# records in the drain's order — devices in attachment order, then the
# harness register, then the loop itself.
function _status(sim::Simulation)
    plane = sim.plane
    ws = Vector{WriterStatus}(undef, length(plane.roster) + 2)
    for (i, e) in enumerate(plane.roster)
        t = get(plane.run_tasks, e.id, nothing)
        ws[i] = _writer_status(_who(e), e.acct, _heartbeat(e.diag), _task_state(t))
    end
    ws[end-1] = _writer_status("harness", plane.harness_acct, nothing, nothing)
    ws[end] = _writer_status("loop", sim.loop_acct, nothing, nothing)
    FrameworkStatus(ws)
end

"""
    latest(sim)

Acquire-load the most recently published snapshot — `nothing` before the first
`init!`. A reader on any task observes an immutable, coherent world for as
long as it holds the value, without coordinating with the loop; the calling
task reads the same reference, §12.6's inspection register.
"""
latest(sim::Simulation) = @atomic :acquire sim.published.latest

"""
    logged(sim)

The log's retained snapshots, oldest first (§11.2): the boundary-zero
endpoint, the bounded middle at the effective stride, and the terminal
endpoint — the latest published boundary — deduplicated when retention
already holds it. A stopped-sim read behind the §11.3 gate: the log is
loop-task bookkeeping, so reading it beside a running loop is the same
hazard class as a mid-run `attach!`; a concurrent reader's register is
`latest(sim)`, and it holds snapshots or loses them (§11.2). Empty before
the first `init!`, and empty under `log = false` — the switch gates
retention wholesale.
"""
function logged(sim::Simulation)
    assert_stopped(sim.control, "logged")
    L = sim.log
    out = Snapshot[]
    L.first === nothing && return out
    push!(out, L.first)
    for s in L.snaps
        s === nothing || push!(out, s)
    end
    L.last === out[end] || push!(out, L.last)
    out
end

# --- reading and writing the table outside the loop ---------------------------
# Path-addressed, dictionary-driven, deliberately off the measured path: these
# are boundary-time operations, never called from inside a phase body. Paths are
# §8.6's canonical strings, and an assembly's own path addresses its faces —
# which resolve to the cells they derive from.

port(sim::Simulation, path::String, name::Symbol) = gather(sim.store, sim.layout.addr[(path, name)])

"""
Write a root slot directly, by the root input face's name — a stopped-sim
establishment operation and a stand-in (README): slot initial values are owned
by the §14 init/trim services (§11.3), absent here, and the *running* write
path is staged batches through the drain alone. Kept for authored initial slot
values until the services exist; a write meant to reach a running trajectory
goes through `stage!`.
"""
function set_slot!(sim::Simulation, face::AbstractString, v)
    scatter!(sim.store, sim.layout.addr[("", Symbol(face))], v)
    nothing
end

"""
State at `path`, from whichever home owns it: `x` in the flat buffer on the
continuous tier, `s` in the component's own store on the discrete one (§7.3).
"""
function state(sim::Simulation{T}, path::String) where {T}
    ci = index_of(sim.flat, path)
    sim.sstores[ci] === nothing || return sim.sstores[ci][]
    _tier(i) = classify_tier(sim.flat.paths[i], sim.flat.comps[i])
    _decls(i) = declarations(sim.flat.comps[i], _tier(i), T)
    off = 0
    for i in 1:(ci-1)
        _tier(i) === CONTINUOUS && (off += nleaves(typeof(_decls(i).x)))
    end
    reconstruct(typeof(_decls(ci).x), sim.xbuf, off)
end

"""Modes at `path` (§7.3). Read-only here: modes are written by handlers alone."""
modes(sim::Simulation, path::String) = sim.mstores[index_of(sim.flat, path)][]
