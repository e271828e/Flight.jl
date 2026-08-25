# The simulation object with its deployment binding (§9.1, §9.2), the stepper
# seam (§10.2) and the phase-body accessor (§9.7). The framework owns the loop;
# the one delegated operation is "advance the continuous state from `t` by `h`".

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
    has_localized::Bool           # any localized event compiled in: the frame loop's fast-path key
    sched::Vector{@NamedTuple{path::String, D::Int, Φ::Int, Δt::Float64}}   # the bound schedule (§9.2)
    sstores::Vector{Any}          # discrete state stores, by component index
    mstores::Vector{Any}          # mode stores, by component index
    stepper::M                    # the seam's backend (§10.2), owning its own scratch
    xnext::Vector{T}              # the retained arrival pair (§10.4): xₙ₊₁ saved before trials clobber
    ẋnext::Vector{T}              # the buffer, ẋₙ₊₁ paid only past a validated trigger
    plane::DataPlane              # the §11.3 roster and harness register, mutable at stopped-sim points
    published::Published          # §11.2's `@atomic latest` holder
    log::SnapshotLog              # §11.2's retained snapshots, loop-task bookkeeping
end

"""
    Simulation(build::Build, T = Float64; h, n = 1, Δt_base = nothing,
               method = RK4, firing_budget = 4, localization_tol = 1e-6,
               localization_budget = 8, log = true, log_every = 1,
               log_max = 65536, chunk_size = 16)
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
                    log = true, log_every = 1, log_max = 65536,
                    chunk_size::Int = 16) where {T}
    method isa Type && method <: AbstractStepper ||
        throw(BuildError("method must be a stepper type — RK4 or Heun — got $method (§10.2)"))
    firing_budget isa Integer && firing_budget ≥ 1 ||
        throw(BuildError("firing_budget must be an integer ≥ 1, got $firing_budget (§10.6)"))
    localization_tol isa Real && localization_tol > 0 ||
        throw(BuildError("localization_tol must be a positive real, got $localization_tol (§10.4)"))
    localization_budget isa Integer && localization_budget ≥ 1 ||
        throw(BuildError("localization_budget must be an integer ≥ 1, got $localization_budget (§10.4)"))
    log isa Bool ||
        throw(BuildError("log must be true or false — the retention switch — got $log (§11.2)"))
    log_every isa Integer && log_every ≥ 1 ||
        throw(BuildError("log_every must be an integer ≥ 1, got $log_every (§11.2)"))
    (log_max isa Integer && log_max ≥ 1) || log_max === Inf ||
        throw(BuildError("log_max must be an integer ≥ 1, or Inf as the explicit " *
                         "opt-out, got $log_max (§11.2)"))
    act = activation(b, T)
    bound = bind_schedule(b, h, n, Δt_base)
    c = compile(b, act, bound.D, bound.Φ, bound.Δt; chunk_size)
    stepper = method(T, length(c.xbuf))
    Simulation{T,typeof(c.store),typeof(c.bodies),typeof(c.clock),typeof(c.events),
               typeof(stepper)}(
        c.store, c.xbuf, c.ẋbuf, c.clock, c.bodies, c.events, act.layout, b.flat,
        bound.h, bound.n, bound.Δt_base, Int(firing_budget), Float64(localization_tol),
        Int(localization_budget), any(c.events.localized), bound.sched,
        c.sstores, c.mstores, stepper, zeros(T, length(c.xbuf)), zeros(T, length(c.xbuf)),
        DataPlane(act.layout, c.store), Published(nothing),
        SnapshotLog(log, Int(log_every), log_max === Inf ? typemax(Int) : Int(log_max)))
end

Simulation(root::AbstractComponent, ::Type{T} = Float64; kw...) where {T} =
    Simulation(build(root), T; kw...)

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
                es.warned[i] = true
                (path, name) = es.names[i]
                @warn "FiringBudget: event `$name` at `$path` has fired $budget times at " *
                      "the boundary t = $(sim.clock.t); its further edges here are lost (§10.6)"
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
"""
function init!(sim::Simulation{T}; t₀::T = zero(T)) where {T}
    sim.clock.t = t₀
    sim.clock.t₀ = t₀
    sim.clock.step = 0
    fill!(sim.events.prior, false)
    _reset!(sim.log)
    boundary!(sim, 0)
    publish!(sim)                 # the boundary-zero snapshot (§11.2, §14.5)
    nothing
end

"""
Advance to `t_end` one frame at a time, each frame §11.1's anatomy — drain,
integrate, boundary sequence, publication. `frame!` carries each grid step
`[tₖ₋₁, tₖ]` through the §10.4 localization loop, firing any `t*` boundaries it
brackets on the way; the frame-top boundary (§10.3) is a base tick every `n`
frames, where the gate reads the tick index, and the empty-due-set boundary in
between. The drain runs at the frame top only, never at a `t*` boundary
(§10.4), while publication follows *every* boundary sequence (§11.2) — the
frame top's here, a `t*` boundary's inside the frame loop, before integration
resumes. The grid is driven by the step counter, so `t_end` is taken to the
nearest frame top; partial advance is the periphery's (§12.6) and absent here.
"""
function run!(sim::Simulation{T}, t_end::Real) where {T}
    target = round(Int, Float64(t_end) / sim.h)
    plane = sim.plane
    @atomic plane.running = true      # the §11.3 freeze: the roster is fixed for the run
    try
        while sim.clock.step < target
            drain!(sim)
            k = (sim.clock.step += 1)
            frame!(sim, k)
            k % sim.n == 0 ? boundary!(sim, k ÷ sim.n) : offtick_boundary!(sim)
            publish!(sim)
        end
    finally
        @atomic plane.running = false
    end
    nothing
end

# --- the roster (§11.3): stopped-sim configuration -----------------------------

"""
    attach!(sim, dev::AbstractDevice, b::AbstractBinding) -> device id

Roster a device under a binding — a stopped-sim configuration operation
(`ServiceLifecycle` while running: the roster is frozen per run, pause
included). The binding's conformance check runs first (§11.6), then the
three-part admission in spec order (§11.3): identity — this instance already
rostered is `AlreadyAttached`, rebinding being spelled `detach!` then
`attach!` — affinity (`CallerTaskConflict`: at most one `needs_calling_task`
holder), and claims (`ClaimConflict`, which the identity check having run
first always makes two *distinct* devices). The claim is staked from its
source — the enumeration called once, or the unclaimed complement computed at
this instant and never recomputed, so attaching the greedy claimant last is
the idiom and a second greedy stakes the empty remainder under an
`EmptyGreedyClaim` warning (§11.6) — the entry's writer is compiled over it,
and the harness register's surface is recompiled to the complement that
remains, renormalizing any pending harness batch (§11.4). Returns the entry's
stable device id. With handles absent, staging against the entry is spelled
`stage!(sim, dev, pairs...)` — a stand-in (README).
"""
function attach!(sim::Simulation, dev::AbstractDevice, b::AbstractBinding)
    plane = sim.plane
    assert_stopped(plane, "attach!")
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
    claim = _claim(plane, sim.layout, b)
    for f in claim                                 # claims: face exclusivity
        haskey(plane.claimedby, f) && throw(BuildError(
            "ClaimConflict: `$f` is claimed by both $(typeof(dev)) and " *
            "$(plane.claimedby[f]) — one writer per slot at any time (§11.3)"))
    end
    id = plane.next_id                             # assigned on admission alone: a
    plane.next_id += 1                             # rejected attach consumes no id
    w = Writer(sim.layout, claim)
    push!(plane.roster, RosterEntry(dev, b, id, w, _drain_thunk(sim.store, w)))
    reclaim!(plane, sim.layout)
    is_greedy(b) && isempty(claim) &&
        @warn "EmptyGreedyClaim: device $id ($(typeof(dev))) under $(typeof(b)) staked " *
              "the empty remainder — every root input face is already claimed (§11.6)"
    id
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
    assert_stopped(plane, "detach!")
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
    stage!(sim, dev, "face" => value, ...)     # a rostered device's cell (§11.4)

Stage a batch of root-slot writes from any task, at any wall-clock moment,
never touching a live slot (§11.1): the entries land in the writer's staging
cell under the one coalescing policy — CAS merge, newest wins per face
(§11.4). Untouched faces survive into the pending batch; re-staged faces take
the newest level. Every check runs here, on the writer's side: a face outside
the writer's surface is discarded with a warning — a device's under
`OutOfClaimEntry` against its claim set, the harness register's under
`ClaimedFaceEntry` naming the incumbent when a rostered claim covers the face
(a rostered greedy claimant empties the harness surface outright, D-192) and
`OutOfClaimEntry` when nothing does — an unconvertible value under
`EntryTypeMismatch`, and the rest of the batch stands. The staged batch is
applied by the drain at the top of the next frame `run!` advances. The device
form is the stand-in for the handle's `stage!` (§11.6, README).
"""
function stage!(sim::Simulation, pairs::Pair...)
    plane = sim.plane
    w = plane.harness
    batch = _normalize(w, pairs, plane.claimedby)
    batch === nothing || _stage!(w, batch)
    nothing
end

function stage!(sim::Simulation, dev::AbstractDevice, pairs::Pair...)
    plane = sim.plane
    i = findfirst(e -> e.dev === dev, plane.roster)
    i === nothing && throw(BuildError(
        "this $(typeof(dev)) instance is not rostered — a device stages into its " *
        "own attached cell (§11.4)"))
    e = plane.roster[i]
    batch = _normalize(e.writer, pairs, plane.claimedby; device = _who(e))
    batch === nothing || _stage!(e.writer, batch)
    nothing
end

"""
The drain (§11.4): exactly one point in each frame, at its top, where the loop
takes each staged batch with one `atomicswap(cell, nothing)` — an indivisible
take, so there is no lost-write window — and applies it through the entry's
compiled scatter, `nothing` skipped, every check long since spent at staging.
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
    end
    plane.harness_drain()
    nothing
end

"""
Publish this boundary (§11.2): build the snapshot in private memory — one
buffer copy, the framework side of §7.5's scope — then a single release-store
to `latest`. Runs only after the boundary sequence completes, so every
published table is boundary-consistent (§10.3), and the copy is what makes the
binding rule hold: nothing reachable from a published snapshot is ever written
again. Behind the store the snapshot enters the log (§11.2): logging dissolves
into publication, retention being the only thing the log adds.
"""
function publish!(sim::Simulation)
    snap = Snapshot(sim.clock.t, sim.clock.step, capture(sim.store), sim.layout)
    @atomic :release sim.published.latest = snap
    log!(sim.log, snap)
    nothing
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
    assert_stopped(sim.plane, "logged")
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
