# The input trace (§11.5): the header captured at `init!` — the resolved
# initial state, the root-input values, the writers' schemas and the run's
# deployment block — and behind it one sparse record per drained batch, in the
# one record format D-176 unified on. The trace is *primary* data (D-029,
# D-038): the log is recomputable from it, and what recomputes it is replay
# (§12.7), whose entry pass — the validation and normalization a trace is
# admitted through — lives at the tail of this file, beside the capture it
# mirrors.
#
# This file holds the register and the pure mechanics; the `Simulation`-facing
# surface — the capture's placement in `init!`, the drain's frame stamp,
# `trace(sim)`, `_compile_feed`'s scalar gate, `replay!` and the drain
# substitution — lives in sim.jl, beside the loop that runs it. The drain
# thunks are compiled here too
# (`_install_writers!`), because the writer's schema index rides in them and
# nothing else may fix it.

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

# A record is a *value*, and the claim replay makes about a continuation — the
# recording is a bit-identical prefix of what the continued session records
# (§12.7) — is an equality between records built by two different registers.
# The default `==` on a mutable-free struct with a `Vector` field is egal, which
# would make that claim untestable, so the field-wise one is defined here, and
# `hash` with it as Julia's convention requires.
Base.:(==)(a::TraceBatch, b::TraceBatch) =
    a.frame == b.frame && a.writer == b.writer && a.entries == b.entries
Base.hash(b::TraceBatch, h::UInt) = hash(b.entries, hash(b.writer, hash(b.frame, h)))

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
One normalized record (§12.7): the frame ordinal it applies at, the compiled
scatter as a zero-argument thunk, and the `TraceBatch` it came from — carried
beside the thunk so a replayed frame re-records exactly what the recording
held, under the recording's own writer index (§12.7: "the new trace inherits
the old header").
"""
const ReplayRecord = @NamedTuple{frame::Int, thunk::Function, record::TraceBatch}

"""
What `_compile_feed` hands the loop (§12.7): the whole recording normalized to
compiled scatters against the *target* layout, in drain order — by frame, then
by the recording's writer index — with a cursor into it.

The conversion is paid here, once, off the loop: the replay drain applies
compiled scatters exactly as the live drain does, and no face name is resolved
per frame under replay either (D-101). `next` is the first record not yet
applied. The frame budget is not carried here: it is `trc.frames` or
`to_boundary · n`, and `replay!` computes it from the `Trace` in hand.
"""
mutable struct ReplayFeed
    records::Vector{ReplayRecord}
    next::Int
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
    feed::Union{Nothing,ReplayFeed}   # non-`nothing` while a replay drives the loop (§12.7)
end

# The empty roster leaves the harness register sole writer, so the fresh
# register's provisional set is `1:1`; `_install_writers!` re-fixes it at every
# capture and every roster change.
TraceRegister(enabled::Bool) = TraceRegister(enabled, nothing, TraceBatch[], 0, 0, 1:1, nothing)

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
# on the frame path boxes an argument: two O(surface-width) scans of the mask —
# one to count the touched positions, one to fill — so the `entries` vector is
# allocated once at its final length rather than grown. The cost of a drained
# batch is then that vector and the values boxed into it, and nothing
# width-dependent beyond them — the cost D-176 records rather than argues away.
# A drained batch always carries at least one touched position (`_normalize`
# returns `nothing` for one that would not), so no record here is empty.
function _record!(reg::TraceRegister, widx::Int, batch::Batch)
    mask = batch.mask
    k = 0
    for i in 1:length(mask)
        mask[i] && (k += 1)
    end
    entries = Vector{Pair{Int,Any}}(undef, k)
    j = 0
    for i in 1:length(mask)
        mask[i] && (entries[j += 1] = i => batch.vals[i])
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

# The recording's header, detached, for the replay that *inherits* it (§12.7):
# every mutable field copied — the stores deeply, they are user values — so the
# trace the replay goes on to build is a value of its own and nothing the caller
# still holds is reachable from it. `_reschema`'s sharing is the growth rule's,
# within one register; this crosses between two.
_detach(h::TraceHeader{T}) where {T} =
    TraceHeader{T}(copy(h.x), deepcopy(h.s), deepcopy(h.m), copy(h.root_inputs),
                   copy(h.schemas), h.deployment, h.layout)

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

"""
The structural fingerprint (§11.5, §12.7): the layout's cell sizes, the root
input-face list, the flat's component paths and each component store's value
type. One function because it has two sides — the capture writes it into the
header, and replay compares the header's against the target's — and two
spellings of one fingerprint would be a silent way for a replay to pass.
`sim` is untyped for include order alone, as `_capture_header` below is.
"""
function _fingerprint(sim)
    ex = sim.exec
    layout = ex.act.layout
    (sizes = copy(layout.sizes),
     root_faces = Symbol[f for (f, _) in layout.root_inputs],
     paths = copy(sim.build.flat.paths),
     stypes = Any[st === nothing ? nothing : typeof(st[]) for st in ex.sstores],
     mtypes = Any[st === nothing ? nothing : typeof(st[]) for st in ex.mstores])
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
    # parameterized explicitly: the `t_end = nothing` run has a `Nothing` in hand
    # where the field is `Union{Nothing,Float64}`, and only the inner
    # constructor's `convert` bridges a NamedTuple's field types
    TraceHeader{T}(copy(ex.xbuf), s, m, roots, Pair{String,Vector{Symbol}}[],
                   deployment, _fingerprint(sim))
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

# ==============================================================================
# Replay's entry pass (§12.7): validation, then normalization
# ==============================================================================
# "Validation is loud and up front" — the whole trace is checked and converted
# before the first frame, so every refusal precedes every write and the replay
# drain resolves no name. The pass is two-staged: the header's own disagreements
# with the target build fail fast as a collection, and only then are the records
# resolved through the schemas the first stage just validated. `sim` is untyped
# throughout for include order alone; the scalar gate that dispatches into this
# is `_compile_feed` in sim.jl, beside the loop it feeds.

# The header against the target `Build` and its deployment binding (§12.7's
# disposition table): the structural fingerprint compared field for field, then
# the seven trajectory-determining deployment parameters. `t₀` is applied rather
# than compared, and the recorded `t_end`/`stop_on` pair is a fact of the
# recorded session, never a constraint on this one — so neither is looked at.
function _check_header!(diags::Vector{Diagnostic}, sim, h::TraceHeader)
    f = _fingerprint(sim)
    l = h.layout
    l.sizes == f.sizes ||
        push!(diags, ReplayHeaderMismatch(what = :store, name = :sizes,
                                          expected = l.sizes, found = f.sizes))
    l.paths == f.paths ||
        push!(diags, ReplayHeaderMismatch(what = :store, name = :paths,
                                          expected = l.paths, found = f.paths))
    l.root_faces == f.root_faces ||
        push!(diags, ReplayHeaderMismatch(what = :root_input,
                                          expected = l.root_faces, found = f.root_faces))
    # the per-component store types, only where the path lists agree on what a
    # component *index* means — otherwise the comparison would be by position
    # between two different models, and the path mismatch above is the honest fact
    if l.paths == f.paths
        for (i, p) in enumerate(f.paths)
            l.stypes[i] === f.stypes[i] ||
                push!(diags, ReplayHeaderMismatch(what = :store, path = p, name = :s,
                                                  expected = l.stypes[i], found = f.stypes[i]))
            l.mtypes[i] === f.mtypes[i] ||
                push!(diags, ReplayHeaderMismatch(what = :store, path = p, name = :m,
                                                  expected = l.mtypes[i], found = f.mtypes[i]))
        end
    end
    d = h.deployment
    for (name, found) in ((:Δt_base, sim.Δt_base), (:h, sim.h), (:n, sim.n),
                          (:method, nameof(typeof(sim.stepper))),
                          (:localization_tol, sim.localization_tol),
                          (:localization_budget, sim.localization_budget),
                          (:firing_budget, sim.firing_budget))
        getfield(d, name) == found ||
            push!(diags, ReplayHeaderMismatch(what = :deployment, name = name,
                                              expected = getfield(d, name), found = found))
    end
    nothing
end

# Each recorded writer's schema against the target's root-input faces (§12.7):
# a recorded name this model does not export is a replay error, reported per
# writer with the whole schema and the list-in-hand beside it. Superseded schema
# entries are checked with the live ones — a batch may still reference them, and
# the compiled scatters below are built from whatever a batch names.
function _check_schemas!(diags::Vector{Diagnostic}, faces::Vector{Symbol}, h::TraceHeader)
    for (tag, schema) in h.schemas
        unknown = Symbol[s for s in schema if !(s in faces)]
        isempty(unknown) ||
            push!(diags, ReplaySchemaMismatch(writer = tag, schema = schema,
                                              unknown = unknown, faces = faces))
    end
    nothing
end

# The compiled scatter as a zero-argument thunk, the `_drain_thunk` idiom: the
# store, the address tuple and the batch are captured *concretely*, so the
# `@generated` `_apply!` specializes once per writer at this stopped-sim compile
# and the replayed frame calls through a barrier with nothing left to infer.
_apply_thunk(store, addrs::Tuple, batch::Batch) = () -> _apply!(store, addrs, batch)

# One recorded writer's compiled scatter against the *target* layout. The
# schema check above already proved every face addressable, so the guard here
# never fires; it is kept because `Writer` would answer an absent face with a
# `KeyError` rather than with the kind that names it.
function _replay_writer(layout::Layout, diags::Vector{Diagnostic}, tag::String,
                        schema::Vector{Symbol}, frame::Int, faces::Vector{Symbol})
    absent = Symbol[f for f in schema if !haskey(layout.addr, ("", f))]
    isempty(absent) && return Writer(layout, schema)
    for f in absent
        push!(diags, ReplayUnknownFace(face = f, frame = frame, writer = tag, faces = faces))
    end
    nothing
end

"""
The trace-record conversion in reverse (§12.7, D-176), paid once and off the
loop: every sparse record's positions are resolved through the writer's schema,
the values converted to the target's declared types, and the whole batch
rebuilt as the positional `Batch` the compiled scatter takes. The recorded
`TraceBatch` rides beside its thunk, because a replay re-records.

Everything here is *collected* — Appendix C gives `ReplayUnknownFace` the
`collected` policy — so a trace with three bad entries reports three. A batch
that produced a diagnostic contributes no record: a partially applied frame is
not a replay of anything.
"""
function _compile_records!(diags::Vector{Diagnostic}, sim, trc::Trace, faces::Vector{Symbol})
    layout, store, schemas = sim.exec.act.layout, sim.exec.store, trc.header.schemas
    writers = Dict{Int,Writer}()
    recs = ReplayRecord[]
    for b in trc.batches
        if !(1 ≤ b.frame ≤ trc.frames)
            # the batch disagrees with the trace's *own* frame count: the drain
            # visits every frame in `1:frames` exactly once, so an ordinal outside
            # it names a frame that never comes round and the record would silently
            # never apply. Named by the writer's tag where the schema list has one
            push!(diags, ReplayHeaderMismatch(
                what = :frame,
                name = 1 ≤ b.writer ≤ length(schemas) ? Symbol(first(schemas[b.writer])) :
                                                        Symbol("writer #$(b.writer)"),
                expected = 1:trc.frames, found = b.frame))
            continue
        end
        if !(1 ≤ b.writer ≤ length(schemas))
            # no schema to resolve the positions through, and the writer index is
            # what is missing — the tag §11.8 cannot supply is spelled positionally
            for (pos, _) in b.entries
                push!(diags, ReplayUnknownFace(face = pos, frame = b.frame,
                                               writer = "writer #$(b.writer)", faces = faces))
            end
            continue
        end
        (tag, schema) = schemas[b.writer]
        if !haskey(writers, b.writer)
            w = _replay_writer(layout, diags, tag, schema, b.frame, faces)
            w === nothing && continue
            writers[b.writer] = w
        end
        w = writers[b.writer]
        vals, mask, ok = Any[w.blank.vals...], fill(false, length(schema)), true
        for (pos, v) in b.entries
            if !(1 ≤ pos ≤ length(schema))
                push!(diags, ReplayUnknownFace(face = pos, frame = b.frame, writer = tag,
                                               faces = faces))
                ok = false
                continue
            end
            vals[pos] = try
                convert(w.types[pos], v)
            catch
                push!(diags, ReplayHeaderMismatch(what = :root_input, name = schema[pos],
                                                  expected = w.types[pos], found = v))
                ok = false
                continue
            end
            mask[pos] = true
        end
        ok || continue
        batch = Batch(convert(typeof(w.blank.vals), (vals...,)), (mask...,))
        push!(recs, (frame = b.frame, thunk = _apply_thunk(store, w.addrs, batch), record = b))
    end
    # the drain's own order (§11.5): by frame, then by the recording's writer
    # index. Stable, so a trace already in drain order — every trace this
    # register produces — keeps exactly the order it was recorded in.
    sort!(recs; by = r -> (r.frame, r.record.writer), alg = MergeSort)
    recs
end
