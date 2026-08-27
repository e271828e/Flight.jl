# The condition algebra and its resolution (§14.1–§14.3, §14.6): the inert lazy
# tree the fragment functions build, the collecting pass that validates it
# against a `Build`, and the plan §14.4's dynamic-walk register executes.
#
# Two properties carry the section. Composition performs no path arithmetic and
# no validation — `at` stores a prefix and never applies it — so a fragment
# function is stack-only construction and a misaddressed node fails at
# resolution, where the build is in hand. And resolution finishes before
# anything is written, which is what makes `init!` all-or-nothing (§14.6).

# --- the four node kinds (§14.2, §14.6) ---------------------------------------
# Every node is isbits except the prefix strings: a `Fragment`'s payloads are
# the author's own values, and the combinators hold their operands in tuples.

"""
Self-vocabulary payloads at one authoring level (§14.2): `x`, `s` and `m`
fields of the component addressed there, and `slots` naming faces of *that
level's* contract. No paths — addressing children is exclusively `at`'s job.
"""
struct Fragment{X,S,M,L}
    x::X
    s::S
    m::M
    slots::L
end

"""`at(prefix, node)`: scoping. Stores the prefix, never applies it (§14.2)."""
struct Scoped{N}
    prefix::String
    node::N
end

"""
`combine(nodes...)`: symmetric collection. Order is diagnostics only, and a
duplicate leaf is an error at resolution with both provenance chains (§14.2).
"""
struct Combined{T<:Tuple}
    nodes::T
end

"""
`override(base, patches...)`: ordered layering, the fourth node kind (§14.6).
The patch wins on a shared leaf and provenance keeps both sources; collisions
*within* one layer remain errors, and layering is variadic.
"""
struct Override{T<:Tuple}
    layers::T
end

const ConditionNode = Union{Fragment,Scoped,Combined,Override}

"""
    fragment(; x = (;), s = (;), m = (;), slots = (;))

The leaf constructor of the condition tree (§14.2, Appendix B). Payloads are
NamedTuples in the authoring level's own vocabulary; a condition speaks state
(`x`, `s`), modes (`m`) and root slots, never outputs and never workspace
(§14.1).
"""
function fragment(; x = (;), s = (;), m = (;), slots = (;))
    for (name, p) in ((:x, x), (:s, s), (:m, m), (:slots, slots))
        p isa NamedTuple || throw(BuildError(
            "ConditionNodeMisuse: `fragment`'s `$name` payload is $(typeof(p)) — every " *
            "payload is a NamedTuple of the authoring level's own names (§14.2)"))
    end
    Fragment(x, s, m, slots)
end

"""
    at(prefix, node)

Scope `node` under `prefix` (§14.2). Nothing is concatenated here: the prefix
is stored, and flattening at resolution is the one place path strings ever
join (§14.3).
"""
at(prefix::AbstractString, node::ConditionNode) = Scoped(String(prefix), node)
at(::AbstractString, other) = _node_misuse(other, ())

"""
    combine(nodes...)

Symmetric collection of sibling nodes (§14.2, D-204). Collision-intolerant by
design: two entries on one leaf are a resolution error naming both provenance
chains, and the layering spelling is `override`. It is not `Base.merge` and
does not extend it — `Base.merge` is last-writer-wins on NamedTuples, the
exact semantics D-065 rejects here.
"""
combine(nodes::ConditionNode...) = Combined(nodes)

# The blend of a node with a bare NamedTuple, both orders: an error method,
# raised at composition time — before any resolution pass or provenance chain
# exists, which is why it carries its own kind rather than a
# `ConditionResolution` sub-kind (§14.2).
combine(a::ConditionNode, b::NamedTuple) = _node_misuse(b, (nameof(typeof(a)),))
combine(a::NamedTuple, b::ConditionNode) = _node_misuse(a, (nameof(typeof(b)),))
combine(nodes...) = _misuse_in(nodes)

"""
    override(base, patches...)

Ordered layering (§14.6): a leaf present in several layers takes the last
layer's value, with provenance recording both sources. The use case is
"baseline plus tweaks" — the collision with `combine`'s duplicate-leaf error
*is* the intent — and trim commits `override(baseline, solution)`.

A slot's leaf is the *root slot* it resolves to, not the face it was written
with, so a full-coverage baseline authored at the root layers cleanly under a
patch a component's own `condition` method ships against its local face.
"""
override(base::ConditionNode, patches::ConditionNode...) = Override((base, patches...))
override(layers...) = _misuse_in(layers)

_misuse_in(args::Tuple) = _node_misuse(args[findfirst(n -> !(n isa ConditionNode), args)],
                                       Tuple(nameof(typeof(n)) for n in args if n isa ConditionNode))

_node_misuse(v, kinds) = throw(BuildError(
    "ConditionNodeMisuse: $(typeof(v)) is not a condition node" *
    (isempty(kinds) ? "" : ", and the node kind(s) in hand are $(join(kinds, ", "))") *
    " — wrap the NamedTuple in `fragment(…)`, or in `at(prefix, fragment(…))`, and " *
    "combine nodes with nodes (§14.2)"))

# --- flattening (§14.3) --------------------------------------------------------
# The only place path strings are ever concatenated: a recursion with a path
# accumulator, carrying the provenance chain beside it — the `at` prefixes and
# the payload position, which is what a collision diagnostic reports.
#
# A slot entry's *leaf* is the root slot its face resolves to, not the face it
# was written with, so the export chain is walked here, at the same moment the
# paths join, and the resolved face becomes the entry's key. Two spellings of
# one root slot are then one leaf everywhere it matters: `override` layers them
# and `combine` collides on them, which is what §14.6's central use case —
# a root-level baseline under a component-local patch — requires.

struct CEntry
    path::String
    store::Symbol        # :x | :s | :m | :slot
    field::Symbol
    value::Any
    prov::String
    face::Union{Nothing,Symbol}   # slot entries: the root slot the chain lands on
end

_key(e::CEntry) = e.face === nothing ? (e.path, e.store, e.field) : ("", :slot, e.face)
_leaf(e::CEntry) = e.face !== nothing ? "root slot `$(e.face)`" :
                   e.store === :slot ? "slot face `$(e.field)` of $(_at(e.path))" :
                                       "`$(e.store).$(e.field)` at $(_at(e.path))"
_step(prov::String, s::String) = isempty(prov) ? s : prov * " → " * s

function _flat(n::Fragment, path::String, prov::String, flat::Flat, viol::Vector{String})
    out = CEntry[]
    for (store, label, payload) in ((:x, "x", n.x), (:s, "s", n.s),
                                    (:m, "m", n.m), (:slot, "slots", n.slots))
        for (field, v) in pairs(payload)
            e = CEntry(path, store, field, v, _step(prov, "fragment($label).$field"), nothing)
            store === :slot &&
                (e = CEntry(e.path, e.store, e.field, e.value, e.prov,
                            _root_slot(flat, e, viol)))
            push!(out, e)
        end
    end
    out
end

_flat(n::Scoped, path::String, prov::String, flat::Flat, viol::Vector{String}) =
    _flat(n.node, _join(path, n.prefix), _step(prov, "at(\"$(n.prefix)\")"), flat, viol)

_flat(n::Combined, path::String, prov::String, flat::Flat, viol::Vector{String}) =
    reduce(vcat, (_flat(k, path, _step(prov, "combine[$i]"), flat, viol)
                  for (i, k) in enumerate(n.nodes)); init = CEntry[])

# Layering (§14.6): each layer is flattened and checked on its own — a
# within-layer collision is still an error — and then folded onto the
# accumulator, the patch replacing the leaf it overrode and inheriting its
# provenance beside its own.
function _flat(n::Override, path::String, prov::String, flat::Flat, viol::Vector{String})
    acc = CEntry[]
    for (i, layer) in enumerate(n.layers)
        label = i == 1 ? "override[base]" : "override[patch $(i - 1)]"
        es = _flat(layer, path, _step(prov, label), flat, viol)
        _check_duplicates!(es, viol)
        for e in es
            j = findfirst(a -> _key(a) == _key(e), acc)
            j === nothing ? push!(acc, e) :
                (acc[j] = CEntry(e.path, e.store, e.field, e.value,
                                 "$(e.prov) (overrode $(acc[j].prov))", e.face))
        end
    end
    acc
end

# `combine`'s one collision rule, collected rather than thrown (§13.1): both
# provenance chains and the directive naming the layering combinator.
function _check_duplicates!(es::Vector{CEntry}, viol::Vector{String})
    seen = Dict{Tuple{String,Symbol,Symbol},CEntry}()
    for e in es
        k = _key(e)
        if haskey(seen, k)
            push!(viol, "DuplicateConditionLeaf: $(_leaf(e)) is written twice — by " *
                        "$(seen[k].prov), and by $(e.prov). `combine` is " *
                        "collision-intolerant by design — use `override(base, patch)` " *
                        "to layer (§14.2, §14.6)")
        else
            seen[k] = e
        end
    end
    nothing
end

# --- the compiled plan (§14.3) --------------------------------------------------

"""
What a resolved condition compiles to (§14.3): per continuous-state leaf its
`xbuf` offset, per discrete or mode store the whole value the write installs —
`merge(defaults, overlay)`, the genuine last-wins NamedTuple merge baked here
(§14.1's fork, not the condition combinator) — and per root slot its compiled
cell address. Every value has been through its leaf's converter already, so
application is a walk with no decisions left in it.

`faces` is the plan's slot coverage, in entry order: the resolution-time
operand of §14.6's totality check.
"""
struct ConditionPlan
    xs::Vector{Tuple{Int,Any}}             # (xbuf offset, value)
    stores::Vector{Tuple{Symbol,Int,Any}}  # (:s | :m, component index, whole store value)
    slots::Vector{Tuple{Symbol,Any,Any}}   # (root face, cell address, value)
    faces::Vector{Symbol}
end

# --- resolution (§14.3) ---------------------------------------------------------

"""
    resolve(node, b::Build, T = Float64) → ConditionPlan

Flatten the condition tree, validate every entry against `b` in §13.1's
collecting register — full list, violations collected, one `BuildError` — and
compile what survives to a plan.

The checks are §14.3's: the path resolves to a component, the field is
declared in that component's `init_x`/`init_s`/`init_m`, the value converts to
the declared leaf type, a slot face resolves through the export chain to a
root slot, and no leaf is written twice — a slot's leaf being the *root slot*
it resolves to, so two spellings of one slot are one leaf. Schema is the
authority on *may you write this, at what type*; the activation's layout
supplies the destination.
"""
function resolve(node::ConditionNode, b::Build, ::Type{T} = Float64) where {T}
    flat, tiers = b.flat, b.tiers
    act = activation(b, T)
    decls, layout = act.decls, act.layout
    viol = String[]
    entries = _flat(node, "", "", flat, viol)
    _check_duplicates!(entries, viol)

    x_offs = _x_offsets(decls, tiers)
    xs = Tuple{Int,Any}[]
    slots = Tuple{Symbol,Any,Any}[]
    faces = Symbol[]
    overlays = Dict{Tuple{Symbol,Int},Vector{Pair{Symbol,Any}}}()

    for e in entries
        if e.store === :slot
            e.face === nothing && continue     # the chain was reported at flattening
            addr = layout.addr[("", e.face)]
            (ok, v) = _convert(_port_type(addr), e.value)
            ok || (push!(viol, _unconvertible(e, e.value, _port_type(addr))); continue)
            push!(slots, (e.face, addr, v))
            push!(faces, e.face)
            continue
        end
        ci = _component(flat, e, viol)
        ci === nothing && continue
        c, tier, d = flat.comps[ci], tiers[ci], decls[ci]
        declared = e.store === :x ? d.x : e.store === :s ? d.s : init_m(c)
        if isempty(declared)
            push!(viol, _no_store(e, c, tier))
            continue
        end
        haskey(declared, e.field) ||
            (push!(viol, _undeclared(e, c, tier, declared, T)); continue)
        L = typeof(declared[e.field])
        (ok, v) = _convert(L, e.value)
        ok || (push!(viol, _unconvertible(e, e.value, L)); continue)
        if e.store === :x
            push!(xs, (x_offs[ci] + _leaf_offset(d.x, e.field), v))
        else
            push!(get!(() -> Pair{Symbol,Any}[], overlays, (e.store, ci)), e.field => v)
        end
    end
    _report_violations(viol)

    # Overlay partiality baked now (§14.3): the store's write is one whole
    # value, `merge(defaults, overlay)`, so application decides nothing.
    stores = Tuple{Symbol,Int,Any}[]
    for ci in eachindex(flat.comps), store in (:s, :m)
        ov = get(overlays, (store, ci), nothing)
        ov === nothing && continue
        defaults = store === :s ? decls[ci].s : init_m(flat.comps[ci])
        push!(stores, (store, ci, convert(typeof(defaults), merge(defaults, (; ov...)))))
    end
    ConditionPlan(xs, stores, slots, faces)
end

# Anything that is not a node reaching a service entry point is the §14.2
# misuse, not a `MethodError`: the directive is the same one `combine` prints,
# and a bare NamedTuple handed to `init!` is exactly the slip it addresses.
resolve(other, ::Build, ::Type = Float64) = _node_misuse(other, ())

# The `x` destinations, per component then per field: the flat buffer's layout
# is the declaration walk, so the plan bakes the offsets the state accessor
# recomputes (§7.1).
function _x_offsets(decls::Vector{Decls}, tiers::Vector{Tier})
    offs, n = Int[], 0
    for (d, t) in zip(decls, tiers)
        push!(offs, n)
        t === CONTINUOUS && (n += nleaves(typeof(d.x)))
    end
    offs
end

function _leaf_offset(x::NamedTuple, field::Symbol)
    off = 0
    for (k, v) in pairs(x)
        k === field && return off
        off += nleaves(typeof(v))
    end
    off
end

_convert(::Type{P}, v) where {P} =
    try
        (true, convert(P, v))
    catch
        (false, nothing)
    end

# The component a non-slot entry addresses. Assemblies are virtual for
# execution (§10.5) and own no state, so an `at` prefix stopping at one has
# nothing to write — and saying so beats "no such path".
function _component(flat::Flat, e::CEntry, viol::Vector{String})
    i = findfirst(==(e.path), flat.paths)
    i === nothing || return i
    push!(viol, "ConditionResolution: the condition addresses $(_at(e.path)), which is " *
                (any(startswith(p, e.path * "/") for p in flat.paths) ?
                 "an assembly — assemblies own no state, and a condition addresses " *
                 "components and root slots (§14.1, §8.5)" :
                 "no component of this build (§14.3)") * " [$(e.prov)]")
    nothing
end

# A slot payload names a face of the authoring level's contract; resolution
# walks the export chain to the root slot and errors if the face never
# surfaces (§14.2). The chain is what the flat list already holds: an input's
# resolved producer is either a root slot or an internal port, and an
# internally wired input has no slot — writing it would be meaningless,
# because the first sweep overwrites it.
function _root_slot(flat::Flat, e::CEntry, viol::Vector{String})
    if isempty(e.path)
        e.field in flat.slots && return e.field
        push!(viol, "ConditionResolution: `$(e.field)` is no root input face — the root's " *
                    "slots are $(_names(flat.slots)) (§14.2) [$(e.prov)]")
        return nothing
    end
    ci = _component(flat, e, viol)
    ci === nothing && return nothing
    conns = flat.conns[ci]
    k = findfirst(p -> first(p) === e.field, conns)
    if k === nothing
        push!(viol, "ConditionResolution: $(_at(e.path)) declares no input face " *
                    "`$(e.field)` — its input faces are " *
                    "$(_names([first(p) for p in conns])) (§14.2) [$(e.prov)]")
        return nothing
    end
    (path, port) = last(conns[k])
    isempty(path) && return port
    push!(viol, "ConditionResolution: $(_at(e.path))'s input face `$(e.field)` reaches no " *
                "root slot — it is wired internally, to `$path`.$port, and the first sweep " *
                "overwrites it; unexported stays unpokeable (§14.2) [$(e.prov)]")
    nothing
end

_names(ns) = isempty(ns) ? "none" : join(("`$n`" for n in ns), ", ")

# The store a condition names has to exist on the component at all: `x` is the
# continuous tier's state and `s` the discrete one's, disjoint by construction
# (D-195), and `m` is continuous-only (§3.2).
_no_store(e::CEntry, c, tier::Tier) =
    "ConditionResolution: $(_leaf(e)) — $(_at(e.path)) is a $(tier_word(tier)) component " *
    "and declares no `init_$(e.store)`" *
    (e.store === :x && tier === DISCRETE ? "; the discrete tier's state is `s` (D-195)" :
     e.store === :s && tier === CONTINUOUS ? "; the continuous tier's state is `x` (D-195)" :
     e.store === :m ? "; modes are declared by `init_m`, continuous-only (§3.2)" : "") *
    " (§14.3) [$(e.prov)]"

# An undeclared field, discriminated against the component's other name
# families: a condition specifies state, modes and root slots — never outputs,
# which are derived data, and never workspace (§14.1).
function _undeclared(e::CEntry, c, tier::Tier, declared::NamedTuple, ::Type{T}) where {T}
    role = haskey(declared_at(output_types, c, tier), e.field) ? "an output port" :
           haskey(declared_at(input_types, c, tier), e.field) ? "an input face" :
           (_declares_workspace(c, tier) &&
            haskey(_declared_workspace(c, tier, T), e.field)) ? "a workspace entry" : nothing
    "ConditionResolution: $(_leaf(e)) is not declared — `init_$(e.store)` at " *
    "$(_at(e.path)) declares $(_names(keys(declared)))" *
    (role === nothing ? "" :
     "; `$(e.field)` is $role, and a condition specifies state, modes and root " *
     "slots — never outputs, never workspace (§14.1)") * " (§14.3) [$(e.prov)]"
end

_declared_workspace(c, tier::Tier, ::Type{T}) where {T} =
    tier === CONTINUOUS ? workspace(c, T) : workspace(c)

_unconvertible(e::CEntry, v, ::Type{P}) where {P} =
    "ConditionResolution: $(_leaf(e)) takes $(P), and the authored value is " *
    "$(repr(v))::$(typeof(v)), which does not convert (§14.3) [$(e.prov)]"

# §13.1's collecting register: the full list, every violation, one throw.
function _report_violations(viol::Vector{String})
    isempty(viol) && return nothing
    throw(BuildError("the condition does not resolve against this build — " *
                     "$(length(viol)) violation$(length(viol) == 1 ? "" : "s") " *
                     "(§14.3, §13.1):\n  " * join(viol, "\n  ")))
end

# --- slot totality (§14.6) ------------------------------------------------------

"""
§14.6's precondition of starting, checked by the service and *not* a property
of conditions, which are legitimately partial: an application establishing a
complete world over virgin stores covers every root slot. Coverage is a
plan-level fact — both operands are resolution-time data — so the check is one
comparison and runs before any evaluation, not merely before any write. A
shortfall is one collected, declaration-ordered `UninitializedSlots` naming
every uncovered face (D-068).

The services path contains no call to `probe_value`: a slot gets a condition
value or the application errors, and there is no third branch. A fabricated
zero is a fine probe input and a terrible flight condition.
"""
function assert_total(plan::ConditionPlan, flat::Flat, op::String)
    covered = Set(plan.faces)
    uncovered = [f for f in flat.slots if !(f in covered)]
    isempty(uncovered) && return nothing
    throw(BuildError("UninitializedSlots: the condition given to `$op` covers no value for " *
                     "root input slot(s) $(_names(uncovered)) — slots are the one " *
                     "initialized datum with no declared default (§11.3), so an application " *
                     "establishing a complete world authors every one of them; nothing was " *
                     "written (§14.6, D-068)"))
end

# --- the dynamic-walk application register (§14.4) ------------------------------

"""
    apply!(sim, plan)

§14.4's dynamic walk: execute the validated entry list by runtime dispatch per
write — microseconds total, allocation permitted, the stopped-sim path never
having been under §7.5's zero-alloc regime — with no per-shape codegen, so
fifty structurally different scripted conditions cost fifty walks rather than
fifty compiles. Which register a service uses is internal, never user-facing:
the specialized `apply!` the iterating services want is the same plan unrolled
(absent here, README).
"""
function apply!(sim::Simulation, plan::ConditionPlan)
    for (off, v) in plan.xs
        flatten!(sim.xbuf, off, v)
    end
    for (store, ci, v) in plan.stores
        (store === :s ? sim.sstores : sim.mstores)[ci][] = v
    end
    for (_, addr, v) in plan.slots
        scatter!(sim.store, addr, v)
    end
    nothing
end
