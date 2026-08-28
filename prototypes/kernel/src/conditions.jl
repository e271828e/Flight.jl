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
fields of the component addressed there, and `inputs` naming faces of *that
level's* contract. No paths — addressing children is exclusively `at`'s job.
"""
struct Fragment{X,S,M,L}
    x::X
    s::S
    m::M
    inputs::L
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
    fragment(; x = (;), s = (;), m = (;), inputs = (;))

The leaf constructor of the condition tree (§14.2, Appendix B). Payloads are
NamedTuples in the authoring level's own vocabulary; a condition speaks state
(`x`, `s`), modes (`m`) and root inputs, never outputs and never workspace
(§14.1).
"""
function fragment(; x = (;), s = (;), m = (;), inputs = (;))
    for (name, p) in ((:x, x), (:s, s), (:m, m), (:inputs, inputs))
        p isa NamedTuple || throw(BuildError(
            "ConditionNodeMisuse: `fragment`'s `$name` payload is $(typeof(p)) — every " *
            "payload is a NamedTuple of the authoring level's own names (§14.2)"))
    end
    Fragment(x, s, m, inputs)
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

An input entry's leaf is the *root input* it resolves to, not the face it was
written with, so a full-coverage baseline authored at the root layers cleanly
under a patch a component's own `condition` method ships against its local
face.
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
# An input entry's *leaf* is the root input its face resolves to, not the face
# it was written with, so the export chain is walked here, at the same moment
# the paths join, and the resolved face becomes the entry's key. Two spellings
# of one root input are then one leaf everywhere it matters: `override` layers
# them and `combine` collides on them, which is what §14.6's central use case —
# a root-level baseline under a component-local patch — requires.

struct CEntry
    path::String
    store::Symbol        # :x | :s | :m | :input
    field::Symbol
    value::Any
    prov::String
    face::Union{Nothing,Symbol}   # input entries: the root input the chain lands on
end

_key(e::CEntry) = e.face === nothing ? (e.path, e.store, e.field) : ("", :input, e.face)
_leaf(e::CEntry) = e.face !== nothing ? "root input `$(e.face)`" :
                   e.store === :input ? "input face `$(e.field)` of $(_at(e.path))" :
                                       "`$(e.store).$(e.field)` at $(_at(e.path))"
_step(prov::String, s::String) = isempty(prov) ? s : prov * " → " * s

function _flat(n::Fragment, path::String, prov::String, flat::Flat, viol::Vector{String})
    out = CEntry[]
    for (store, label, payload) in ((:x, "x", n.x), (:s, "s", n.s),
                                    (:m, "m", n.m), (:input, "inputs", n.inputs))
        for (field, v) in pairs(payload)
            e = CEntry(path, store, field, v, _step(prov, "fragment($label).$field"), nothing)
            store === :input &&
                (e = CEntry(e.path, e.store, e.field, e.value, e.prov,
                            _root_input(flat, e, viol)))
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
(§14.1's fork, not the condition combinator) — and per root input its compiled
cell address. Every value has been through its leaf's converter already, so
application is a walk with no decisions left in it.

`faces` is the plan's root-input coverage, in entry order: the resolution-time
operand of §14.6's totality check.
"""
struct ConditionPlan
    xs::Vector{Tuple{Int,Any}}             # (xbuf offset, value)
    stores::Vector{Tuple{Symbol,Int,Any}}  # (:s | :m, component index, whole store value)
    inputs::Vector{Tuple{Symbol,Any,Any}}  # (root face, cell address, value)
    faces::Vector{Symbol}
end

# --- resolution (§14.3) ---------------------------------------------------------

"""
    resolve_condition(node, b::Build, T = Float64) → ConditionPlan

Flatten the condition tree, validate every entry against `b` in §13.1's
collecting register — full list, violations collected, one `BuildError` — and
compile what survives to a plan.

The checks are §14.3's: the path resolves to a component, the field is
declared in that component's `init_x`/`init_s`/`init_m`, the value converts to
the declared leaf type, an input face resolves through the export chain to a
root input, and no leaf is written twice — an input entry's leaf being the
*root input* it resolves to, so two spellings of one root input are one leaf.
Schema is the
authority on *may you write this, at what type*; the activation's layout
supplies the destination.
"""
function resolve_condition(node::ConditionNode, b::Build, ::Type{T} = Float64) where {T}
    flat, tiers = b.flat, b.tiers
    act = activation(b, T)
    decls, layout = act.decls, act.layout
    viol = String[]
    entries = _flat(node, "", "", flat, viol)
    _check_duplicates!(entries, viol)

    x_offs = _x_offsets(decls, tiers)
    xs = Tuple{Int,Any}[]
    inputs = Tuple{Symbol,Any,Any}[]
    faces = Symbol[]
    overlays = Dict{Tuple{Symbol,Int},Vector{Pair{Symbol,Any}}}()

    for e in entries
        if e.store === :input
            e.face === nothing && continue     # the chain was reported at flattening
            addr = layout.addr[("", e.face)]
            (ok, v) = _convert(_port_type(addr), e.value)
            ok || (push!(viol, _unconvertible(e, e.value, _port_type(addr))); continue)
            push!(inputs, (e.face, addr, v))
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
    ConditionPlan(xs, stores, inputs, faces)
end

# Anything that is not a node reaching a service entry point is the §14.2
# misuse, not a `MethodError`: the directive is the same one `combine` prints,
# and a bare NamedTuple handed to `init!` is exactly the slip it addresses.
resolve_condition(other, ::Build, ::Type = Float64) = _node_misuse(other, ())

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

# The component a non-input entry addresses. Assemblies are virtual for
# execution (§10.5) and own no state, so an `at` prefix stopping at one has
# nothing to write — and saying so beats "no such path".
function _component(flat::Flat, e::CEntry, viol::Vector{String})
    i = findfirst(==(e.path), flat.paths)
    i === nothing || return i
    push!(viol, "ConditionResolution: the condition addresses $(_at(e.path)), which is " *
                (any(startswith(p, e.path * "/") for p in flat.paths) ?
                 "an assembly — assemblies own no state, and a condition addresses " *
                 "components and root inputs (§14.1, §8.5)" :
                 "no component of this build (§14.3)") * " [$(e.prov)]")
    nothing
end

# An `inputs` payload names a face of the authoring level's contract, whatever
# that level's class; resolution walks the export chain to the root input and
# errors if the face never surfaces (§14.2). The chain is the `Build`'s own
# input-side face graph, total under one-level routing (§9.2, D-207): a face's
# producer is either a root input or an internal port, and a component-fed face
# reaches none — writing it would be meaningless, because the first sweep
# overwrites it.
function _root_input(flat::Flat, e::CEntry, viol::Vector{String})
    if isempty(e.path)
        e.field in flat.root_inputs && return e.field
        push!(viol, "ConditionResolution: `$(e.field)` is no root input face — the root's " *
                    "inputs are $(_names(flat.root_inputs)) (§14.2) [$(e.prov)]")
        return nothing
    end
    k = findfirst(p -> first(p) === (e.path, e.field), flat.in_faces)
    if k === nothing
        here = [f for ((p, f), _) in flat.in_faces if p == e.path]
        push!(viol, isempty(here) && !_addresses_level(flat, e.path) ?
                    "ConditionResolution: the condition addresses $(_at(e.path)), which is " *
                    "no component of this build (§14.3) [$(e.prov)]" :
                    "ConditionResolution: $(_at(e.path)) declares no input face " *
                    "`$(e.field)` — its input faces are $(_names(here)) (§14.2) [$(e.prov)]")
        return nothing
    end
    (path, port) = last(flat.in_faces[k])
    isempty(path) && return port
    push!(viol, "ConditionResolution: $(_at(e.path))'s input face `$(e.field)` reaches no " *
                "root input — it is wired internally, to `$path`.$port, and the first sweep " *
                "overwrites it; unexported stays unpokeable (§14.2) [$(e.prov)]")
    nothing
end

# Does the path name a level of this build at all — a component, or an assembly
# some component sits under? Assemblies leave no row of their own in the flat
# list, so an empty face list alone does not tell a bare typo from an assembly
# that declares no input face.
_addresses_level(flat::Flat, path::String) =
    any(p == path || startswith(p, path * "/") for p in flat.paths)

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
# families: a condition specifies state, modes and root inputs — never outputs,
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
     "inputs — never outputs, never workspace (§14.1)") * " (§14.3) [$(e.prov)]"
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

# --- root-input totality (§14.6) ------------------------------------------------

"""
§14.6's precondition of starting, checked by the service and *not* a property
of conditions, which are legitimately partial: an application establishing a
complete world over virgin stores covers every root input. Coverage is a
plan-level fact — both operands are resolution-time data — so the check is one
comparison and runs before any evaluation, not merely before any write. A
shortfall is one collected, declaration-ordered `UninitializedInputs` naming
every uncovered face (D-068).

The services path contains no call to `probe_value`: a root input gets a
condition value or the application errors, and there is no third branch. A
fabricated zero is a fine probe input and a terrible flight condition.
"""
function assert_total(plan::ConditionPlan, flat::Flat, op::String)
    covered = Set(plan.faces)
    uncovered = [f for f in flat.root_inputs if !(f in covered)]
    isempty(uncovered) && return nothing
    throw(BuildError("UninitializedInputs: the condition given to `$op` covers no value for " *
                     "root input(s) $(_names(uncovered)) — root inputs are the one " *
                     "initialized datum with no declared default (§11.3), so an application " *
                     "establishing a complete world authors every one of them; nothing was " *
                     "written (§14.6, D-068)"))
end

# --- the dynamic-walk application register (§14.4) ------------------------------

"""
    apply!(ex::Executor, plan)
    apply!(sim, plan)

§14.4's dynamic walk: execute the validated entry list by runtime dispatch per
write — microseconds total, allocation permitted, the stopped-sim path never
having been under §7.5's zero-alloc regime — with no per-shape codegen, so
fifty structurally different scripted conditions cost fifty walks rather than
fifty compiles. Which register a service uses is internal, never user-facing:
the specialized `apply!` the iterating services want is the same plan unrolled
(absent here, README).
"""
function apply!(ex::Executor, plan::ConditionPlan)
    for (off, v) in plan.xs
        flatten!(ex.xbuf, off, v)
    end
    for (store, ci, v) in plan.stores
        (store === :s ? ex.sstores : ex.mstores)[ci][] = v
    end
    for (_, addr, v) in plan.inputs
        scatter!(ex.store, addr, v)
    end
    nothing
end

apply!(sim::Simulation, plan::ConditionPlan) = apply!(sim.exec, plan)

# --- capture: the gather twin of the application register (§14.1, §14.10) -------

"""
    capture(sim) → (condition, t)

Read the current committed stores **and root inputs** back as a condition value
(§14.1, glossary): every component's `x` field by field, every discrete `s`,
every mode store and every root input, as one `combine` of per-component
fragments plus the root-input fragment. The pair is capture → tweak → apply:
`init!(sim, c; t₀ = t)` re-establishes exactly the world that was read, and a
`trim!` resumes from it as `trim!(sim, problem; baseline = c, t₀ = t)`
(§14.8).

Root-input coverage is what makes the captured condition **total**, hence
re-applicable under §14.6 — a capture is by construction the one condition that
never needs a baseline under it. The condition is time-free: `t` rides beside
it, because time is not a store of any component (§14.5), and it is passed back
as the `t₀` argument.

Legal in `initialized` and `stopped`, the states whose stores are committed and
boundary-consistent. A `built` simulation has no such stores yet — boundary
zero has never run — and `running` and `errored` are the ordinary service
refusals (§14, §11.3, §13.6); all four are one `ServiceLifecycle`.

Re-applying reproduces the captured world bit for bit, with one caveat that is
boundary zero's rather than capture's: the re-application runs the sequence
(§14.5), so `project` and any guard already holding in the captured state fire
again there — which is exactly what makes a warm restart a *fresh run from
these values* rather than a resumption.
"""
function capture(sim::Simulation{T}) where {T}
    lc = lifecycle(sim)
    lc === :running && throw(BuildError(
        "ServiceLifecycle: `capture` is a stopped-sim operation and the simulation " *
        "is running — the loop owns the stores between drains (§11.3, §14)"))
    lc === :errored && throw(BuildError(
        "ServiceLifecycle: this simulation ended `errored` — terminally stopped; " *
        "post-mortem inspection of its stores stays a diagnostic read, and may not " *
        "become a condition value (§13.6, §14)"))
    lc === :built && throw(BuildError(
        "ServiceLifecycle: `capture` reads committed, boundary-consistent stores, and " *
        "this simulation has never been initialized — `capture` is legal in " *
        "`initialized` and `stopped` (§14, §12.6)"))
    ex, flat, tiers = sim.exec, sim.build.flat, sim.build.tiers
    act = activation(sim.build, T)
    offs = _x_offsets(act.decls, tiers)
    nodes = ConditionNode[]
    for ci in eachindex(flat.comps)
        d = act.decls[ci]
        payload = NamedTuple()
        tiers[ci] === CONTINUOUS && !isempty(d.x) &&
            (payload = merge(payload, (x = _capture_x(d.x, ex.xbuf, offs[ci]),)))
        ex.sstores[ci] === nothing || (payload = merge(payload, (s = ex.sstores[ci][],)))
        ex.mstores[ci] === nothing || (payload = merge(payload, (m = ex.mstores[ci][],)))
        isempty(payload) && continue
        push!(nodes, at(flat.paths[ci], fragment(; payload...)))
    end
    isempty(flat.root_inputs) || push!(nodes, fragment(inputs =
        NamedTuple{Tuple(flat.root_inputs)}(Tuple(gather(ex.store, act.layout.addr[("", f)])
                                                  for f in flat.root_inputs))))
    (combine(nodes...), ex.clock.t)
end

# One component's `x` payload, field by field: the flat buffer is read through
# the same declaration walk the plan writes through, so what comes back is
# exactly what goes in (§7.1).
function _capture_x(x::NamedTuple, buf::Vector, base::Int)
    off, vals = base, Any[]
    for v in values(x)
        push!(vals, reconstruct(typeof(v), buf, off))
        off += nleaves(typeof(v))
    end
    NamedTuple{keys(x)}(Tuple(vals))
end
