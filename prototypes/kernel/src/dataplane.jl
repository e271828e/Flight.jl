# The data plane, minimal slice (§11): the harness register's staging cell
# with its frame-top drain — planes 1 and 2 of §11.1's replacement for the
# shared mutable model — and snapshot publication, plane 3. What stands behind
# them in the spec — the roster, claims, devices, the trace, the log — is
# deliberately absent (README). With no roster built, the harness surface
# (§11.3) is the unclaimed complement in its degenerate form: every root face.
#
# This file holds the types and the pure mechanics; the `Simulation`-facing
# surface — `stage!`, `drain!`, `publish!`, `latest` — lives in sim.jl, beside
# the loop that calls it.

"""
The staging cell (§11.4): one atomic reference holding the pending batch — CAS
merge on the writer's side, `atomicswap` at the drain. The batch rides boxed in
a `Ref` so every handoff is one pointer-width atomic operation (§11.1) whatever
the surface's size, and the CAS compares batch *identity*: a failed replace
means the loaded batch was intercepted — drained, or merged past by another
stager — and the retry re-reads exactly what is pending now. The failure case
is precisely correct: an intercepted batch is already applied (or already
carries the other stager's entries), and must not be re-staged.
"""
mutable struct StagingCell{B}
    @atomic pending::Union{Nothing,Base.RefValue{B}}
end

"""
The harness register (§11.3): the always-present, task-free write path behind
`stage!(sim, "face" => value, …)`. Its surface is *derived* rather than
claimed — the unclaimed complement, which with no roster is the whole root
face set — and its staged representation is compiled against that surface
exactly as a device's would be against its claim (§11.4): a positional tuple
over the face set, `Union{Nothing, Tᵢ}` per face, isbits and pointer-free,
with `nothing` meaning *not touched this time*, never "reset". The face-name →
position schema, the slot types and the per-position cell addresses — the
compiled scatter's data — are all fixed at deployment binding, the stopped-sim
point standing in for attach.
"""
struct Harness{B<:Tuple,A<:Tuple}
    faces::Vector{Symbol}    # position → face name: the schema
    types::Vector{Any}       # position → the slot's declared type Tᵢ
    addrs::A                 # position → the slot's cell address
    cell::StagingCell{B}
end

_port_type(::CellAddr{P,K}) where {P,K} = P

function Harness(layout::Layout)
    faces = Symbol[f for (f, _) in layout.slots]
    addrs = Tuple(layout.addr[("", f)] for f in faces)
    types = Any[_port_type(a) for a in addrs]
    B = Tuple{(Union{Nothing,T} for T in types)...}
    Harness{B,typeof(addrs)}(faces, types, addrs, StagingCell{B}(nothing))
end

"""
§11.4's shim, run at staging on the writer's own side: the author's face ⇒
value pairs become a batch — name → position, convert to the slot's declared
type, fill `nothing`. Every diagnostic site follows the compilation here, to
staging: face validity and value convertibility are static facts of the run,
so a face with no position in the schema is discarded with an `OutOfClaimEntry`
warning and an unconvertible value with `EntryTypeMismatch`, while the rest of
the batch stands. Nothing remains at the drain. Returns `nothing` for a batch
with no surviving entry.
"""
function _normalize(h::Harness{B}, pairs) where {B}
    vals = Any[nothing for _ in h.faces]
    for (face, v) in pairs
        i = findfirst(==(Symbol(face)), h.faces)
        if i === nothing
            @warn "OutOfClaimEntry: `$face` names no face on the staging surface — " *
                  "the root faces are $(join(h.faces, ", ")); the entry is discarded (§11.4)"
            continue
        end
        vals[i] = try
            convert(h.types[i], v)
        catch
            @warn "EntryTypeMismatch: `$face` cannot take $(repr(v))::$(typeof(v)) — " *
                  "its slot declares $(h.types[i]); the entry is discarded (§11.4)"
            continue
        end
    end
    all(isnothing, vals) ? nothing : convert(B, (vals...,))
end

# The one coalescing policy (§11.4): merge, newest wins per face. Untouched
# faces survive; re-staged faces take the newest level — the per-face ZOH.
# Positional, so it compiles straight-line and union-splits. Plain `::Tuple`
# arguments, because a sparse batch's *runtime* type is the narrow tuple of
# what it touched — covariantly inside `B`, but two batches rarely share it.
_merge(pending::Tuple, incoming::Tuple) =
    map((p, i) -> i === nothing ? p : i, pending, incoming)

# The CAS merge loop (§11.4), on the writer's task.
function _stage!(h::Harness{B}, batch::B) where {B}
    cell = h.cell
    while true
        pending = @atomic cell.pending
        merged = pending === nothing ? batch : _merge(pending[], batch)
        (; success) = @atomicreplace cell.pending pending => Base.RefValue{B}(merged)
        success && return nothing
    end
end

# The drain's application (§11.4): the compiled scatter, position → slot cell,
# statically typed, `nothing` skipped — pure application, no checks.
@generated function _apply!(store, addrs::Tuple, batch::Tuple)
    stmts = [:(batch[$i] === nothing || scatter!(store, addrs[$i], batch[$i]))
             for i in 1:fieldcount(batch)]
    quote
        $(stmts...)
        nothing
    end
end

# --- publication (§11.2) -------------------------------------------------------

"""
A snapshot (§11.2): the boundary-consistent signal table — the *whole* table,
root slots riding along as the source cells they are — with `t` and the frame
ordinal. Built in private memory by copying the cell buffers, then handed out
through one release-store; nothing reachable from a published snapshot is ever
written again, which is what makes the lock-free read sound. The state stores
(`x`, `s`, `m`) are deliberately not carried (§11.2) and the framework status
is absent here (README). Read it with `port(snap, path, name)`, addressed
exactly as the live table.
"""
struct Snapshot{T,S<:StoreBundle}
    t::T
    frame::Int
    store::S
    layout::Layout    # build-frozen addressing, shared, never copied
end

port(s::Snapshot, path::String, name::Symbol) = gather(s.store, s.layout.addr[(path, name)])

# One boundary's capture: fresh buffers, one allocation per boundary — the
# framework side of §7.5's scope, which carved publication and logging out.
capture(b::StoreBundle) = StoreBundle(map(cs -> CellStore(copy(cs.buf)), b.stores))

"""
§11.2's `@atomic latest` reference in its own mutable holder, so the
`Simulation` itself stays immutable: release-store at publication,
acquire-load at `latest(sim)`. The exchange is wait-free in both directions —
a wedged reader cannot delay publication, and the loop cannot tear a reader's
view.
"""
mutable struct Published
    @atomic latest::Union{Nothing,Snapshot}
end
