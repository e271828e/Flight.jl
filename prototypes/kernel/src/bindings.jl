# The binding's runtime half (§11.6): the shipped `TableBinding` with its
# generic `map_input` — the one data-driven binding the framework writes, and
# the owner of the shared pure conditioning helper (§11.4) — and the output
# side: the table members of §14.4's read-selector family, the resolution
# `attach!` runs against the build, and the one compiled gather whose labeled
# NamedTuple is what `map_output` receives (§11.2). `map_input`/`map_output`
# are *conventions of the author-owned loop idiom*: the framework never calls
# them, so what is enforced here is only what construction and attach can
# check — the table's own shape, and the enumerated `reads`. Everything the
# mappings touch at runtime is bounded elsewhere: `map_input` by the staging
# checks (§11.4), on the writer's own task; `map_output` by receiving exactly
# the gather's NamedTuple, what it puts on the wire being the peer's business.
#
# The framework-legible half of the binding vocabulary — the roots, the
# declared sides, `claims`/`reads` and the conformance check — lives in
# roster.jl.

"""
    TableBinding(; channel = (face = "...", deadzone = ..., expo = ...), ...)

The shipped data-driven binding (§11.6): the framework writes its `map_input`
once, and a table value — channel name => entry — is constructed per
device × deployment pairing, where configurations are made. Each entry names
the root input face its channel writes (`face`, mandatory) and optionally the
conditioning parameters the generic `map_input` applies (`deadzone`, `expo`);
any other key is a typo caught here, at construction. The entry tuple rides
in the type, so the mapping specializes per table with no dynamic dispatch.

The input side is declared (`is_input`), and the claim is the table's face
set — `claims` derives it, so what the binding may write is exactly what its
entries name. A *code-driven* binding looks identical to the framework: a
JSON telecommand peer whose `claims` returns the vocabulary and whose
`map_input` parses bytes.
"""
struct TableBinding{T<:NamedTuple} <: AbstractBinding
    table::T
end

function TableBinding(; entries...)
    table = NamedTuple(entries)
    for (k, e) in pairs(table)
        e isa NamedTuple || throw(BuildError(
            "TableBinding: entry `$k` must be a NamedTuple — " *
            "(face = ..., deadzone = ..., expo = ...) (§11.6)"))
        haskey(e, :face) || throw(BuildError(
            "TableBinding: entry `$k` names no `face` — the face is what the " *
            "channel writes (§11.6)"))
        e.face isa Union{AbstractString,Symbol} || throw(BuildError(
            "TableBinding: entry `$k`'s face must be a face name, got " *
            "$(repr(e.face)) (§11.6)"))
        for key in keys(e)
            key in (:face, :deadzone, :expo) || throw(BuildError(
                "TableBinding: entry `$k` carries `$key` — the entry vocabulary " *
                "is `face`, `deadzone`, `expo` (§11.6)"))
        end
        dz = get(e, :deadzone, nothing)
        dz === nothing || 0 <= dz < 1 || throw(BuildError(
            "TableBinding: entry `$k`'s deadzone must lie in [0, 1), got $dz (§11.6)"))
        ex = get(e, :expo, nothing)
        ex === nothing || 0 <= ex <= 1 || throw(BuildError(
            "TableBinding: entry `$k`'s expo must lie in [0, 1], got $ex (§11.6)"))
    end
    TableBinding(table)
end

is_input(::TableBinding) = true
claims(b::TableBinding) = String[String(e.face) for e in values(b.table)]

"""
    map_input(datum::NamedTuple, b::TableBinding) -> ("face" => value, ...)

`TableBinding`'s generic mapping — the shared pure conditioning helper, with
an owner (§11.4, §11.6). The datum is whatever the author's loop assembled,
one field per *touched* channel: `map_input` returns face ⇒ value pairs for
exactly those, so a sparse datum stages a sparse batch and merge does the
rest (§11.4). The idiom is `stage!(handle, map_input(datum, binding(handle))...)`,
on the device's own task.

A datum field naming no table channel is configuration drift, not a bad
datum: the mapping and the datum are written by the same author, so the
throw is deliberate — it lands in the wrapper as the `DeviceCrash` it is
(§11.6). Cross-datum state — press counters, edge detection — lives in the
device struct, maintained by the loop, and arrives *inside* the datum:
`map_input` stays pure, and staged values are levels, never deltas (§11.4).
"""
function map_input(datum::NamedTuple, b::TableBinding)
    map(keys(datum)) do k
        haskey(b.table, k) || error(
            "map_input: the datum carries `$k`, which names no channel of this " *
            "TableBinding — its channels are $(_facelist(keys(b.table))) (§11.6)")
        e = b.table[k]
        String(e.face) => _condition(datum[k], e)
    end
end

# --- the output side: selectors, resolution, the compiled gather (§14.4, §11.2) --

# The selectors themselves are §14.4's closed family, declared in readers.jl
# above this file: the three *table* members are what a snapshot-bound reader
# may name, and the source rule is enforced here, at the attach point where the
# source is finally known (§14.4). `get_output` is the inspection register —
# deep paths, zero promises, free access, right for looking at *this* build —
# and `get_face` the integration register: a root-exported output face, named,
# curated, meaning-stable under substitution (§11.2). `get_input` reads a root
# input back, the source cell it is. What this register does not take is depth
# *inside* a cell: a binding read is a whole cell, as every reader of the
# published table is (README).

"""
The compiled gather (§11.2, §14.4): one attachment's `reads`, resolved and
frozen — the labels as a type parameter, the cell addresses as a tuple — so
`gather(handle, snap)` builds its labeled NamedTuple with no name resolved
per read. The exact mirror of the compiled scatter the drain applies
(§11.4), run in the other direction over a published snapshot.
"""
struct ReadGather{L,A<:Tuple}
    addrs::A
end
ReadGather{L}(addrs::A) where {L,A<:Tuple} = ReadGather{L,A}(addrs)

_gather(r::ReadGather{L}, s::Snapshot) where {L} =
    NamedTuple{L}(map(a -> gather(s.store, a), r.addrs))

"""
Resolve one attachment's `reads` against the build and compile the gather —
`attach!`'s output-side work (§11.2), run once at the attach point so a
binding that drifted from its model fails there, not with silent garbage on
the wire. The shape is fixed: `reads` returns a NamedTuple of labeled
selectors, `(; label = get_output(...), ...)`, and the labels are the
NamedTuple `map_output` receives. Every failure names the selector at fault;
the did-you-mean candidate lists are absent (README).
"""
function _compile_reads(layout::Layout, nt, T::Type)
    nt isa NamedTuple || throw(BuildError(
        "BindingContractMismatch: $T's `reads` must return a NamedTuple of labeled " *
        "selectors — (; label = get_output(...), ...) — got $(typeof(nt)) (§11.6, §14.4)"))
    addrs = map(values(nt)) do s
        s isa ReadSelector || throw(BuildError(
            "BindingContractMismatch: $T's `reads` entries must be read selectors — " *
            "get_output, get_input or get_face (§14.4) — got $(repr(s))"))
        _resolve_read(layout, s, T)
    end
    ReadGather{keys(nt)}(addrs)
end

_root_input_names(layout::Layout) = Symbol[f for (f, _) in layout.root_inputs]

# §14.4's source rule, enforced where the source is known: a snapshot carries
# no state stores by construction (§11.2) and `ẋ` is integrator scratch, so a
# snapshot-bound reader naming a store selector is a resolution error at
# attach — in the didactic register, with the remedy named.
_resolve_read(::Layout, s::StoreSelector, T::Type) = throw(BuildError(
    "ReadBindingUnresolved: $T reads $(_spell(s)), a *store* selector — the store " *
    "selectors resolve only against live stores, and a binding reads a published " *
    "snapshot, which deliberately carries none (§14.4, §11.2). The remedy is to " *
    "declare the field public and read the port published from it"))

function _resolve_read(layout::Layout, s::GetOutput, T::Type)
    s.i === nothing || throw(BuildError(
        "ReadBindingUnresolved: $T reads $(_spell(s)) — a binding read is a whole cell, " *
        "and sub-cell index addressing is absent in this register (§14.4, README)"))
    haskey(layout.addr, (s.path, s.name)) || throw(BuildError(
        "ReadBindingUnresolved: $T reads get_output(\"$(s.path)\", :$(s.name)), which " *
        "names no cell — only declared outputs, assembly faces and root inputs are " *
        "addressable (§14.4)"))
    layout.addr[(s.path, s.name)]
end

function _resolve_read(layout::Layout, s::GetInput, T::Type)
    s.face in _root_input_names(layout) || throw(BuildError(
        "ReadBindingUnresolved: $T reads get_input(:$(s.face)), which names no root " *
        "input face — the root inputs are $(_facelist(_root_input_names(layout))) (§14.4)"))
    layout.addr[("", s.face)]
end

function _resolve_read(layout::Layout, s::GetFace, T::Type)
    s.name in _root_input_names(layout) && throw(BuildError(
        "ReadBindingUnresolved: $T reads get_face(:$(s.name)), which names a root " *
        "*input* face — the integration register is the exported output faces, and " *
        "a root input is read back with get_input (§14.4, §11.2)"))
    haskey(layout.addr, ("", s.name)) || throw(BuildError(
        "ReadBindingUnresolved: $T reads get_face(:$(s.name)), which names no " *
        "root-exported output face (§14.4, §11.2)"))
    layout.addr[("", s.name)]
end

"""
    map_output(nt, b) -> wire datum

The output side's convention name (§11.6), declared here for the loop idiom —
`send(dev.socket, map_output(gather(handle, snap), binding(handle)))` — and
never called by the framework: it receives exactly the compiled gather's
labeled NamedTuple, and what it puts on the wire is the peer's business. An
output binding defines its own method; the library's `Readout` returns the
NamedTuple itself.
"""
function map_output end

# The conditioning (§11.4): axis-convention values in [-1, 1], symmetric about
# zero. The deadzone zeroes the band and rescales the remainder so the
# endpoints stay fixed; expo blends linear into cubic — a = (1-e)·a + e·a³ —
# attenuating the midrange with the endpoints again fixed. An entry declaring
# neither passes its value through untouched, which is what carries a
# throttle's [0, 1] level or a press counter (the levels doctrine): faces take
# post-conditioning semantics, and only where conditioning is declared does
# the axis convention bind.
function _condition(v, e)
    dz = get(e, :deadzone, nothing)
    ex = get(e, :expo, nothing)
    dz === nothing && ex === nothing && return v
    x = clamp(float(v), -1, 1)
    a = abs(x)
    dz === nothing || (a = a <= dz ? zero(a) : (a - dz) / (1 - dz))
    ex === nothing || (a = (1 - ex) * a + ex * a^3)
    flipsign(a, x)
end
