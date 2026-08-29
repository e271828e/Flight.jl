# The input trace (§11.5): the header captured at `init!` — the resolved
# initial state, the root-input values, the writers' schemas and the run's
# deployment block — and behind it one sparse record per drained batch, in the
# one record format D-176 unified on. The trace is *primary* data (D-029,
# D-038): the log is recomputable from it, and what recomputes it is replay
# (§12.7), absent here (README).
#
# This file holds the register and the pure mechanics; the `Simulation`-facing
# surface — the capture's placement in `init!`, the drain's frame stamp,
# `trace(sim)` — lives in sim.jl, beside the loop that runs it. The drain
# thunks are compiled here too (`_install_writers!`), because the writer's
# schema index rides in them and nothing else may fix it.

"""
The trace header (§11.5): the full initial state as **resolved values**, never
the authored overlay (D-038), captured after `apply!` and the root-input
writes and before the boundary-zero sequence — both halves of §14.5's
placement being load-bearing, since replay re-executes boundary zero and a
post-sequence capture would hand it already-transitioned state.

`x` is the flat continuous buffer; `s` and `m` carry each component's store
value (`nothing` where the component owns none). `root_inputs` is what makes
replay possible at all: an unfed `mixture = 0.5` appears in no batch, so the
initial root-input values are recorded beside the stores. `schemas` is the
run's writers as `tag => face-name-by-position` — the positional records are
meaningless without it, and replay does not reconstruct claims (§12.7) — the
tags being §11.8's own writer names, `_who(entry)` and `"harness"`.
`deployment` and `layout` are the two fingerprints: the trajectory depends on
the deployment block exactly as it depends on the stores (the block sits
outside the `Build`, and `t₀` post-dates even deployment), and the structural
half is what a replay target is compared against (§12.7).
"""
struct TraceHeader{T}
    x::Vector{T}                                    # a copy of xbuf after apply!
    s::Vector{Any}                                  # per component: the store value, or nothing
    m::Vector{Any}                                  # likewise for the mode stores
    root_inputs::Vector{Pair{Symbol,Any}}           # face => the resolved value in its cell
    schemas::Vector{Pair{String,Vector{Symbol}}}    # writer tag => face-name-by-position
    deployment::@NamedTuple{t₀::T, Δt_base::Float64, h::Float64, n::Int, method::Symbol,
                            localization_tol::Float64, localization_budget::Int,
                            firing_budget::Int, t_end::Union{Nothing,Float64},
                            stop_on::Vector{Symbol}}
    layout::@NamedTuple{sizes::Vector{Pair{DataType,Int}}, root_faces::Vector{Symbol},
                        paths::Vector{String}, stypes::Vector{Any}, mtypes::Vector{Any}}
end

"""
One drained batch, retained sparse (§11.5, D-176): the masked (touched)
positions as `position ⇒ value` pairs against the writer's schema, which is
`header.schemas[writer]`. `frame` is the frame ordinal the batch was drained
at the top of — the drain runs before the clock's step increments, so a batch
drained at the top of frame `k` carries `frame = k` and replays there.
"""
struct TraceBatch
    frame::Int
    writer::Int                       # index into header.schemas
    entries::Vector{Pair{Int,Any}}
end

"""
What `trace(sim)` hands back (§11.5): header plus batches, the *primary*
record from which the log is derived — a value, detached from the register
that built it. `frames` counts the drains since the header, one per frame, so
a recording whose last frames staged nothing still knows how long it ran.
"""
struct Trace{T}
    header::TraceHeader{T}
    batches::Vector{TraceBatch}       # in drain order: by frame, then by writer index
    frames::Int
end

"""
The `Simulation`'s trace holder: the kill switch fixed at construction
(§11.5's plain switch for memory-constrained marathon sessions, D-029), the
captured header, the records behind it, and the two counters the drain
maintains — `frame`, the ordinal `drain!` stamps the frame's records with,
and `frames`, the recording's length.

`live_writers` names the current writer set's entries in `header.schemas`.
**The schema list only grows**: a roster change is a stopped-sim point that
recompiles the harness writer and may add or remove device writers, while
batches already recorded reference the old indices — so every capture and
every roster change *appends* the current set (the roster in attachment
order, then the harness: the drain's own order) and recompiles the drain
thunks against the new indices. Earlier entries stay, referenced by earlier
batches, and a run's `schemas` may therefore carry superseded ones.
"""
mutable struct TraceRegister
    enabled::Bool
    header::Union{Nothing,TraceHeader}
    batches::Vector{TraceBatch}
    frames::Int
    frame::Int                        # the ordinal this frame's records take
    live_writers::UnitRange{Int}      # the current set's entries in header.schemas
end

# The empty roster leaves the harness register sole writer, so the fresh
# register's provisional set is `1:1`; `_install_writers!` re-fixes it at every
# capture and every roster change.
TraceRegister(enabled::Bool) = TraceRegister(enabled, nothing, TraceBatch[], 0, 0, 1:1)

# §11.5's clearing at `init!`: the header and every record it stood in front of
# go together, the recording's length with them. `live_writers` is left where it
# is — the drain thunks hold those indices, and the capture below re-fixes both
# at once.
function _reset!(reg::TraceRegister)
    reg.header = nothing
    empty!(reg.batches)
    reg.frames = 0
    reg.frame = 0
    nothing
end

# The conversion at the drain (§11.5, D-176), inside the drain thunk so nothing
# on the frame path boxes an argument: an O(surface-width) scan of the mask, the
# touched positions paired with their values. The one allocation per drained
# batch is the `entries` vector and the values boxed into it — the cost D-176
# records rather than argues away. A drained batch always carries at least one
# touched position (`_normalize` returns `nothing` for one that would not), so
# no record here is empty.
function _record!(reg::TraceRegister, widx::Int, batch::Batch)
    mask = batch.mask
    entries = Pair{Int,Any}[]
    for i in 1:length(mask)
        mask[i] && push!(entries, i => batch.vals[i])
    end
    push!(reg.batches, TraceBatch(reg.frame, widx, entries))
    nothing
end

# The run's writers in the drain's own order (§11.3, §11.4): each rostered
# device in attachment order, then the harness register. The tags are §11.8's
# writer names, so a trace and a published status name a writer identically.
function _writer_schemas(plane)
    schemas = Pair{String,Vector{Symbol}}[_who(e) => copy(e.writer.faces) for e in plane.roster]
    push!(schemas, "harness" => copy(plane.harness.faces))
    schemas
end

_reschema(h::TraceHeader{T}, schemas) where {T} =
    TraceHeader{T}(h.x, h.s, h.m, h.root_inputs, schemas, h.deployment, h.layout)

"""
The growth rule (§11.5), and the one site a drain thunk is compiled at: append
the current writer set to the header's schema list, name the appended range
`live_writers`, and recompile every thunk with its writer's new index. The
header is *rebuilt* rather than mutated, so a `Trace` already handed out keeps
the schema list it was given — a valid prefix for its own batches, indices
only ever growing.

Reached from three places, all of them stopped-sim: the header capture below,
which starts from an empty list, and `reclaim!`'s two callers, `attach!` and
`detach!`. With no header — before the first `init!`, or under `trace = false`
— nothing is appended and the indices are provisional; no drain runs before
boundary zero has, so nothing reads them.
"""
function _install_writers!(reg::TraceRegister, plane)
    k = length(plane.roster) + 1
    if reg.header === nothing
        reg.live_writers = 1:k
    else
        h = reg.header
        n = length(h.schemas)
        reg.header = _reschema(h, vcat(h.schemas, _writer_schemas(plane)))
        reg.live_writers = (n + 1):(n + k)
    end
    for (i, e) in enumerate(plane.roster)
        # the entry is immutable and its thunk carries the index, so the
        # recompilation replaces the entry itself (§11.4's stopped-sim compile)
        plane.roster[i] = RosterEntry(
            e.dev, e.binding, e.id, e.writer,
            _drain_thunk(plane.store, e.writer, reg, reg.live_writers[i]),
            e.should_abort, e.diag, e.acct, e.handle)
    end
    plane.harness_drain = _drain_thunk(plane.store, plane.harness, reg, reg.live_writers[k])
    nothing
end

# §11.5's header, read off the simulation at §14.5's placement. `sim` is
# untyped for include order alone — this file precedes sim.jl, the drain
# thunks it compiles being what the data plane is built with.
function _capture_header(sim)
    ex = sim.exec
    layout = ex.act.layout
    T = eltype(ex.xbuf)      # the deployment's scalar, off the buffer that carries it
    s = Any[st === nothing ? nothing : deepcopy(st[]) for st in ex.sstores]
    m = Any[st === nothing ? nothing : deepcopy(st[]) for st in ex.mstores]
    roots = Pair{Symbol,Any}[f => gather(ex.store, layout.addr[("", f)])
                             for (f, _) in layout.root_inputs]
    # the effective termination pair is the one `init!` knows: the constructor's,
    # `run!`'s per-run override post-dating the capture (§13.5)
    deployment = (t₀ = ex.clock.t₀, Δt_base = sim.Δt_base, h = sim.h, n = sim.n,
                  method = nameof(typeof(sim.stepper)),
                  localization_tol = sim.localization_tol,
                  localization_budget = sim.localization_budget,
                  firing_budget = sim.firing_budget, t_end = sim.t_end,
                  stop_on = copy(sim.stop_on))
    structure = (sizes = copy(layout.sizes),
                 root_faces = Symbol[f for (f, _) in layout.root_inputs],
                 paths = copy(sim.build.flat.paths),
                 stypes = Any[st === nothing ? nothing : typeof(st[]) for st in ex.sstores],
                 mtypes = Any[st === nothing ? nothing : typeof(st[]) for st in ex.mstores])
    # parameterized explicitly: the `t_end = nothing` run has a `Nothing` in hand
    # where the field is `Union{Nothing,Float64}`, and only the inner
    # constructor's `convert` bridges a NamedTuple's field types
    TraceHeader{T}(copy(ex.xbuf), s, m, roots, Pair{String,Vector{Symbol}}[],
                   deployment, structure)
end

# The capture, at `init!`'s §14.5 placement: the header first, then the growth
# rule's first application — which is what fills the empty schema list and
# fixes the thunks' indices for the trajectory that is about to open.
function _capture!(sim)
    reg = sim.trace
    reg.enabled || return nothing
    reg.header = _capture_header(sim)
    _install_writers!(reg, sim.plane)
    nothing
end
