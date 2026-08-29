# Hierarchy (§8.5, §8.6, §6.1): class by declaration shape, path resolution, and
# the flatten pass. Assemblies are virtual for execution (§10.5) — what leaves
# this file is a flat list of primitives with their absolute paths and one
# resolved producer per input, and nothing downstream knows the tree existed.

# --- class (§8.5) -------------------------------------------------------------
# Class is not announced either: `child_connections` is the assembly marker, any
# leaf declaration a primitive's, and the rule is total — a component-typed
# struct declaring neither family has no class to read.

@enum Class PRIMITIVE ASSEMBLY

const LEAF_FAMILY = "`init_x`, `init_s`, `init_m`, `workspace`, `input_types`, " *
                    "`output_types`, `events` or a stage (`h_x`, `h_xu`, `h_s`, `h_su`, " *
                    "`f`, `g`, `project`)"

"""The leaf declarations `c` defines, in inventory order (§8.2, §8.5)."""
function leaf_declarations(c)
    found = Symbol[]
    for (name, fn) in ((:init_x, init_x), (:init_s, init_s), (:init_m, init_m))
        _declares(fn, c) && push!(found, name)
    end
    # Either arity is a leaf declaration; which one is lawful is the tier's
    # question, settled by the classifier (§8.2), not this one.
    for (name, fn) in ((:workspace, workspace), (:input_types, input_types),
                       (:output_types, output_types))
        (_declares(fn, c) || _declares(fn, c, Type{Float64})) && push!(found, name)
    end
    _declares(events, c) && push!(found, :events)
    for (name, fn) in ((:h_x, h_x), (:h_xu, h_xu), (:h_s, h_s), (:h_su, h_su),
                       (:f, f), (:g, g), (:project, project))
        has_stage(fn, c) && push!(found, name)
    end
    found
end

"""The class of `c` at `path`, or a `BuildError` naming what makes it unreadable."""
function classify(path::String, c)
    leaves = leaf_declarations(c)
    if _declares(child_connections, c)
        isempty(leaves) ||
            throw(BuildError(ClassMixed(path = path, declarations = leaves)))
        return ASSEMBLY
    end
    isempty(leaves) || return PRIMITIVE
    throw(BuildError(ClassUnreadable(path = path, families = LEAF_FAMILY,
                                     holds_components = _holds_components(c))))
end

_at(path::String) = isempty(path) ? "the root component" : "`$path`"
_join(path::String, seg::String) = isempty(path) ? seg : path * "/" * seg

_holds_components(c) = any(fieldnames(typeof(c))) do name
    v = getfield(c, name)
    v isa AbstractComponent ||
        ((v isa NamedTuple || v isa Tuple) && any(e -> e isa AbstractComponent, v))
end

# --- children (§8.5) ----------------------------------------------------------
# Fields whose type is `<: AbstractComponent` are children, the field name their
# path segment. A `Tuple` or `NamedTuple` field whose elements are *all*
# components is a container: transparent grouping with no contract and no wiring
# of its own, contributing its elements as children *of the parent*, path-named
# `"field/1"…"field/N"` or `"field/key"`. Every other field is inert parameter
# data.
#
# One container field per type may be declared **name-transparent** (D-211), and
# its elements are then contributed under their bare keys — `"key"`, `"1"` —
# everywhere a child name appears. Naming is the only thing the declaration
# changes: the elements are the parent's children exactly as before, in the same
# declaration order, and the container keeps its transparency of contract. What
# pays for it is one declaration check — the declared symbol must name a
# container field — and a collision family in three arms: no two children may
# end up sharing a name, and the two collisions a bare key reaches that no child
# name can, its own field's name (which the rate declaration's sugar already
# spells) and a sibling container field's name, where that field contributes
# children (whose `"field/key"` segment grammar it would shadow).

"""`segment => instance` for every child of `c`, in field order."""
children(path::String, c) = first(_children(path, c))

"""
`(kids, fields)`: the children of `c` as `segment => instance` pairs, and per
child the field name that contributed it — the rate declaration's second key
form, the bare field name applying one entry to every element of a container
(§8.7). The sugar keys on the *field*, so a name-transparent container keeps it
unchanged.
"""
function _children(path::String, c)
    tf = transparent_container(c)
    kids = Pair{String,Any}[]
    fields = Symbol[]
    prov = String[]                            # per child, who contributed it
    # The child-naming pass collects (§13.1): every field is walked, and the
    # whole violation list leaves through one throw at the end. A component with
    # three mixed containers reports three, not the first.
    viol = Diagnostic[]
    # The sibling fields a bare key could shadow: those that currently contribute
    # children, because the grammar a key shadows is the one that currently
    # reaches something. An empty field reserves nothing, and it has to be that
    # way round: an empty `Tuple` or `NamedTuple` is an empty container and empty
    # inert parameter data at once — the value cannot tell them apart, and the
    # walk below already treats them as one case — so reserving every empty
    # field's name would refuse a bare key over inert data, a false positive
    # where no shadow exists. Legality is then per-instantiation, which is the
    # framework's norm: every wire is validated against the instance too.
    # Collected before the walk, because the shadowed field may be declared
    # after the transparent one.
    shadowable = String[String(name) for name in fieldnames(typeof(c))
                        if name !== tf && _is_container(getfield(c, name)) &&
                           !isempty(getfield(c, name))]
    for name in fieldnames(typeof(c))
        v = getfield(c, name)
        if v isa AbstractComponent
            push!(kids, string(name) => v)
            push!(fields, name)
            push!(prov, "field `$name`")
        elseif v isa NamedTuple || v isa Tuple
            n = count(e -> e isa AbstractComponent, v)
            n == 0 && continue                     # inert data; an empty container too
            if n != length(v)
                push!(viol, ContainerMixed(path = path, field = name,
                                           types = unique(Any[typeof(e) for e in v
                                                              if !(e isa AbstractComponent)])))
                continue
            end
            bare = name === tf
            for k in keys(v)
                p = "$(bare ? "name-transparent " : "")container field `$name`, element `$k`"
                # The collision family's other two arms, both reachable only by a
                # bare key and neither of them a duplicate *child* name, so
                # `_check_child_names` below can see neither (§8.5, D-211, D-212).
                hit = false
                if bare && string(k) == string(name)
                    push!(viol, ChildNameCollision(path = path, name = string(k),
                                                   reason = :sample_times_sugar,
                                                   provenance = [p], field = name))
                    hit = true
                end
                if bare && string(k) in shadowable
                    push!(viol, ChildNameCollision(path = path, name = string(k),
                                                   reason = :sibling_field, provenance = [p]))
                    hit = true
                end
                hit && continue                    # a shadowed key names no child
                push!(kids, (bare ? string(k) : string(name, "/", k)) => v[k])
                push!(fields, name)
                push!(prov, p)
            end
        end
    end
    _check_transparent(path, c, tf, viol)
    _check_child_names(path, kids, prov, viol)
    isempty(viol) || throw(BuildError(viol))
    kids, fields
end

# The container form, the empty one included — it contributes zero children, and
# parametric code then needs no special case (§8.5).
_is_container(v) = (v isa NamedTuple || v isa Tuple) && all(e -> e isa AbstractComponent, v)

# The declaration is checked after the walk, so a mixed container reports as one
# rather than as a bad transparency declaration.
function _check_transparent(path::String, c, tf, viol::Vector{Diagnostic})
    tf === nothing && return nothing
    ok = tf in fieldnames(typeof(c)) && _is_container(getfield(c, tf))
    ok || push!(viol, TransparentContainerUnknown(path = path, field = tf,
                                                  component = string(typeof(c))))
    nothing
end

# General, not transparency-specific: whatever mix of fields and containers
# produced them, two children may not share a name. Bare keys make the case
# reachable, but the rule is the older one — a child name is a path segment, and
# a path segment addresses one component.
function _check_child_names(path::String, kids, prov, viol::Vector{Diagnostic})
    for i in eachindex(kids), j in 1:(i - 1)
        first(kids[i]) == first(kids[j]) &&
            push!(viol, ChildNameCollision(path = path, name = first(kids[i]),
                                           reason = :two_children,
                                           provenance = [prov[j], prov[i]]))
    end
    nothing
end

# --- paths (§8.6, §6.1) -------------------------------------------------------
# Slash-separated, relative to the declaring assembly, no leading slash: the one
# canonical form, used verbatim in declarations and in error messages. A terminal
# path's last segment is a port or face name; the prefix names one child. Deep
# structural paths survive on the read side — inspection, provenance, the table
# accessors — and nowhere in the three wiring declarations (D-207).

"""
Resolve terminal `path` against assembly `asm` at `base`, returning the component
it names, that component's absolute path, and the final segment.

§6.1's one-level rule lives here (D-207): a connection endpoint names an
**immediate child and one of its faces** — one child segment, plus the key
segment where the child is a container element (`"units/1/e"` is one level, not
two), and one face name. Anything deeper is a build error whatever the declared
field types along it, which is why the generic-holding question never arises in
this register: an endpoint stops before any field it could traverse past
(§13.3). Faces are the only currency crossing a boundary, so a route through
several levels is declared level by level, each assembly speaking of its own
children alone.
"""
function resolve_terminal(entry::String, base::String, asm, path::AbstractString;
                          owner::String = _at(base))
    segs = String.(split(path, '/'))
    length(segs) > 1 ||
        throw(BuildError(PathResolution(entry = entry, spelling = String(path),
                                        reason = :not_a_terminal, owner = owner)))
    kid, seg = _one_level(entry, base, asm, path, segs, 1; owner)
    kid, _join(base, seg), Symbol(segs[end])
end

# The child half of the rule, shared by both of §13.3's path primitives: `tail`
# is how many segments follow the child — one face name for a terminal path,
# none for `resolve`'s bare child path. Container children match under their
# D-211 naming, which is what the two-segment lookahead serves: an undeclared
# container's element spends two segments on the child, a transparent one's
# spends one, and neither is "deeper".
function _one_level(entry::String, base::String, asm, path::AbstractString,
                    segs::Vector{String}, tail::Int; owner::String = _at(base))
    kids = children(base, asm)
    j = findfirst(kid -> first(kid) == segs[1], kids)
    j === nothing && length(segs) > 1 + tail &&
        (j = findfirst(kid -> first(kid) == segs[1] * "/" * segs[2], kids))
    j === nothing &&
        throw(BuildError(PathResolution(entry = entry, spelling = String(path),
                                        reason = :unknown_child, owner = owner,
                                        segment = segs[1],
                                        candidates = String[first(k) for k in kids])))
    seg, kid = kids[j]
    count(==('/'), seg) + 1 + tail == length(segs) ||
        throw(BuildError(PathResolution(entry = entry, spelling = String(path),
                                        reason = :reaches_past, owner = owner,
                                        segment = seg, level = _join(base, seg), tail = tail)))
    kid, seg
end

# --- §13.3's build primitives -------------------------------------------------
# The four the declaration surface calls: `resolve` and `resolve_terminal` in
# their public, entry-less forms, plus the two face-list accessors. This is the
# *structural* register of §13.3's table — the one-level rule verbatim, the same
# walk wiring resolution runs, entered from a declaration body with no wiring
# entry to attribute the failure to.

"""
    resolve(asm, path) → AbstractComponent

The declared-field walk along `/` segments, container children included under
their D-211 naming (bare keys for a name-transparent container). One level
(§6.1, D-207): `path` names an **immediate child** — one child segment, plus the
key segment where the child is a container element — and anything deeper is a
build error naming the child it reaches past.
"""
function resolve(asm, path::AbstractString)
    who = "`resolve` on `$(nameof(typeof(asm)))`"
    isempty(path) &&
        throw(BuildError(PathResolution(entry = who, spelling = "", reason = :empty_path,
                                        owner = "the component in hand")))
    first(_one_level(who, "", asm, path, String.(split(path, '/')), 0;
                     owner = "the component in hand"))
end

"""
    resolve_terminal(asm, path) → (component, name)

The terminal split: the final segment is the port or face name, the prefix
resolves through `resolve`. The split is unambiguous because face names may
contain dots, never slashes (§8.6).
"""
function resolve_terminal(asm, path::AbstractString)
    comp, _, name = resolve_terminal("`resolve_terminal` on `$(nameof(typeof(asm)))`",
                                     "", asm, path; owner = "the component in hand")
    comp, String(name)
end

# The key set of a contract declaration, whichever arity declares it: the keys are
# a tier-independent fact, and §8.2's classifier is what settles a disagreement.
_contract(fn, c) = _declares(fn, c, Type{Float64}) ? fn(c, Float64) :
                   _declares(fn, c) ? fn(c) : NamedTuple()

"""
    input_faces(c) → Vector{String}

A leaf's `input_types` keys — asked at the nominal `Float64`, the key set being
`T`-independent — or an assembly's `input_connections` face names. Declaration
order is preserved: deterministic printouts, stable diagnostics (§13.3).
"""
input_faces(c) = classify("", c) === PRIMITIVE ?
                 String[String(k) for k in keys(_contract(input_types, c))] :
                 String[String(face) for (face, _) in input_connections(c)]

"""
    output_faces(c) → Vector{String}

`input_faces`' mirror: a leaf's `output_types` keys, or an assembly's
`output_connections` face names, in declaration order (§13.3).
"""
output_faces(c) = classify("", c) === PRIMITIVE ?
                  String[String(k) for k in keys(_contract(output_types, c))] :
                  String[String(face) for (_, face) in output_connections(c)]

# --- §8.8's passthrough helpers -----------------------------------------------
# `input_connections` and `output_connections` are ordinary functions evaluated
# at build against the concrete instance, so they may *compute* entries from
# child contracts. These two are the framework's own sugar over the primitives
# above — the pass-through case, where an assembly re-exports the faces of a
# child it does not itself feed. Two helpers rather than one keyword, because
# after the boundary split a single call cannot emit into two declarations
# (D-209). Nothing else about computed connections is built here: the entries
# they return are ordinary pairs, mixing freely with hand-written ones, and
# every check that meets them is the build's own.

"""
    input_passthrough(asm, child_path; sep = ".",
                      prefix = replace(child_path, "/" => sep),
                      except = (), only = ())

Every input face of `child_path` the assembly does not feed, exposed on its own
boundary under `prefix * sep * face` — splatted into `input_connections` (§8.8).
The default `prefix` folds the path's slash into `sep`, so an undeclared
container element (`"units/1"`) labels its faces `"units.1.…"` — a legal face
name — by default; an explicit `prefix` is used verbatim, and `prefix = ""`
drops the prefixing entirely. `except` and `only` filter face names
within the child's set and are mutually exclusive. `child_path` names an
immediate child (a bare key where the container is name-transparent); a deeper
path meets `resolve`'s one-level rejection like any other wiring endpoint.
"""
function input_passthrough(asm, child_path::AbstractString;
                           sep::AbstractString = ".",
                           prefix::AbstractString = replace(child_path, "/" => sep),
                           except::Tuple = (), only::Tuple = ())
    names = input_faces(resolve(asm, child_path))
    wanted = _passthrough_faces("input_passthrough", child_path, names, except, only)
    Tuple(_labelled(prefix, sep, n) => string(child_path, "/", n) for n in wanted)
end

"""
    output_passthrough(asm, child_path; sep = ".",
                       prefix = replace(child_path, "/" => sep),
                       except = (), only = ())

`input_passthrough`'s sibling on the outward boundary (D-209), splatted into
`output_connections`: the same surface over `output_faces` — the same folded
default `prefix` included — its pairs reading
along the flow — internal source => face name — as every pair in that
declaration does. Its consumer is one-level routing (§6.1): every level
re-exports the outputs it surfaces, so the output side needs the computed
spelling the input side already has.
"""
function output_passthrough(asm, child_path::AbstractString;
                            sep::AbstractString = ".",
                            prefix::AbstractString = replace(child_path, "/" => sep),
                            except::Tuple = (), only::Tuple = ())
    names = output_faces(resolve(asm, child_path))
    wanted = _passthrough_faces("output_passthrough", child_path, names, except, only)
    Tuple(string(child_path, "/", n) => _labelled(prefix, sep, n) for n in wanted)
end

_labelled(prefix, sep, n) = isempty(prefix) ? String(n) : string(prefix, sep, n)

# Exclusivity is enforced, not documented, and a filter naming a face the child
# does not have errors with the list in hand — the same did-you-mean shape every
# declaration-time refusal takes here (§8.8).
function _passthrough_faces(who::String, child_path::AbstractString,
                            names::Vector{String}, except::Tuple, only::Tuple)
    isempty(except) || isempty(only) ||
        throw(BuildError(UnknownFaceSelection(who = who, path = String(child_path),
                                              reason = :both_given)))
    unknown = [String(n) for n in (except..., only...) if !(String(n) in names)]
    isempty(unknown) ||
        throw(BuildError(UnknownFaceSelection(who = who, path = String(child_path),
                                              reason = :unknown_names, names = unknown,
                                              candidates = names)))
    isempty(only) ? setdiff(names, String[String(n) for n in except]) :
                    String[String(n) for n in only]
end

# --- endpoint resolution (§8.6) -----------------------------------------------
# Faces are kind-blind (D-172): an endpoint's final segment resolves to a
# primitive's port or to a sub-assembly's face alike, and a face resolves
# recursively to its own internal endpoint. A face's type and tier are therefore
# derived — they are its ultimate internal endpoint's — never declared.

"""The primitive port a producing endpoint ultimately names, as `(path, port)`."""
function resolve_source(entry::String, base::String, asm, path::AbstractString)
    comp, cpath, name = resolve_terminal(entry, base, asm, path)
    if classify(cpath, comp) === PRIMITIVE
        haskey(_contract(output_types, comp), name) && return (cpath, name)
        _wrong_direction(entry, path, cpath, name, comp, "producer")
    else
        for (src, face) in output_connections(comp)
            String(face) == String(name) && return resolve_source(entry, cpath, comp, src)
        end
        _wrong_direction(entry, path, cpath, name, comp, "producer")
    end
end

"""
The primitive inputs a consuming endpoint ultimately names, as `(path, face)`.
Several, when the endpoint is a sub-assembly's input face fanning out through the
boundary.
"""
function resolve_dest(entry::String, base::String, asm, path::AbstractString)
    comp, cpath, name = resolve_terminal(entry, base, asm, path)
    if classify(cpath, comp) === PRIMITIVE
        haskey(_contract(input_types, comp), name) && return [(cpath, name)]
        _wrong_direction(entry, path, cpath, name, comp, "consumer")
    else
        for (face, inner) in input_connections(comp)
            String(face) == String(name) || continue
            return _fanout(entry, cpath, comp, inner)
        end
        _wrong_direction(entry, path, cpath, name, comp, "consumer")
    end
end

_endpoints(p::AbstractString) = (p,)
_endpoints(ps::Tuple) = ps

_fanout(entry, base, comp, inner) =
    reduce(vcat, (resolve_dest(entry, base, comp, p) for p in _endpoints(inner));
           init = Tuple{String,Symbol}[])

# Direction is declared by the method; the resolved endpoint only cross-checks it.
function _wrong_direction(entry, path, cpath, name, comp, wanted)
    ins, outs = input_faces(comp), output_faces(comp)
    found = String(name) in ins ? "an input" : String(name) in outs ? "an output" : nothing
    found === nothing &&
        throw(BuildError(UnknownPort(entry = entry,
                                     end_ = wanted == "producer" ? :source : :destination,
                                     path = cpath, spelling = String(path), port = name,
                                     candidates = Symbol.(vcat(ins, outs)))))
    throw(BuildError(FaceDirectionConflict(entry = entry, path = cpath,
                                           spelling = String(path),
                                           found = found == "an input" ? :input : :output,
                                           wanted = Symbol(wanted))))
end

# --- the flatten pass ---------------------------------------------------------

"""
What the pipeline consumes: primitives by absolute path, one resolved producer
per declared input, the root's input faces (the [root inputs](§11.3)) and §9.2's
two-sided face table — the assembly faces the periphery may read, aliased onto
the cells they derive from, and beside them every input face at every level with
the producer it routes to. The input side is total: one-level routing gives
every signal crossing a boundary a declared face there (D-207), so a fragment's
`inputs` payload resolves from any authoring level (§14.2).
"""
struct Flat
    paths::Vector{String}
    comps::Vector{Any}
    conns::Vector{Vector{Pair{Symbol,Tuple{String,Symbol}}}}   # face => (producer path, port)
    root_inputs::Vector{Symbol}                                # root input faces, in order
    in_faces::Vector{Pair{Tuple{String,Symbol},Tuple{String,Symbol}}}   # (path, face) => producer
    out_faces::Vector{Pair{Tuple{String,Symbol},Tuple{String,Symbol}}}  # (path, face) => producer
    triples::Vector{NTuple{3,Int}}      # per component: (anchor, m, c), anchor 0 the base grid
    anchors::Vector{NTuple{2,Rational{Int}}}   # anchors 1…K: the exact (T, τ) pairs
    aprov::Vector{String}               # per anchor, the declaring scope and key
end

function index_of(flat::Flat, path::String)
    i = findfirst(==(path), flat.paths)
    i === nothing && throw(InternalInvariant("no component at path `$path`"))
    i
end

struct _Walk
    paths::Vector{String}
    comps::Vector{Any}
    feeds::Dict{Tuple{String,Symbol},Tuple{String,Symbol}}
    claims::Dict{Tuple{String,Symbol},String}                  # who claimed it, for the message
    root_inputs::Vector{Symbol}
    routes::Vector{Tuple{String,Symbol,Vector{Tuple{String,Symbol}}}}   # (path, face, consumers)
    out_faces::Vector{Pair{Tuple{String,Symbol},Tuple{String,Symbol}}}
    triples::Vector{NTuple{3,Int}}
    anchors::Vector{NTuple{2,Rational{Int}}}
    aprov::Vector{String}
    viol::Vector{Diagnostic}            # the walk's collected violations (§13.1)
end

# --- the sample-time fold (§8.7, §9.1, §10.5) -----------------------------------
# Nested rate declarations compile to one `(anchor, m, c)` triple per component,
# folding down the tree beside the wiring walk: the root scope seeds
# `(A₀, 1, 0)`, `Relative(K, φ)` under a scope at `(a, mₛ, cₛ)` steps to
# `(a, K·mₛ, cₛ + φ·mₛ)`, and `Absolute` severs and re-seeds a fresh anchor at
# `(Aₖ, 1, 0)` — a nested anchor simply severs again. Anchor 0 is the base grid,
# symbolic until deployment binds `Δt_base`: final divisors for anchored entries
# do not exist until then (§9.1). The canonical residue `c < m` holds within each
# anchor's subtree by the affine law's one-line induction.

"""
Validate a scope's `sample_times` against §10.5's constraints, with path
attribution (§9.1, §13.1): wrapper-typed values only, `K ≥ 1`, `0 ≤ φ < K`,
`T > 0`, `0 ≤ τ < T`, and every key naming an immediate child.
"""
function _check_sample_times(path::String, st, kids, fields, viol::Vector{Diagnostic})
    _rv(reason; kw...) = push!(viol, RatesViolation(; path = path, reason = reason, kw...))
    if !(st isa NamedTuple)
        _rv(:declaration_shape)                    # nothing further is iterable
        return nothing
    end
    for (k, v) in pairs(st)
        if !(v isa Relative || v isa Absolute)
            _rv(:value_vocabulary; key = k, value = v)
            continue                               # neither residue arm applies
        end
        # The residue bound is stated against the multiplier, so an invalid `K`
        # (or `T`) leaves `φ` (or `τ`) with nothing to be measured against: the
        # dependent check is skipped, not doubled up.
        if v isa Relative
            v.K ≥ 1 ? (0 ≤ v.φ < v.K || _rv(:phase; key = k, value = v.φ)) :
                      _rv(:multiplier; key = k, value = v.K)
        else
            v.T > 0 ? (0 ≤ v.τ < v.T || _rv(:offset; key = k, value = v.τ)) :
                      _rv(:period; key = k, value = v.T)
        end
        any(seg == String(k) for (seg, _) in kids) ||
            any(String(fld) == String(k) for fld in fields) ||
            _rv(:unknown_child; key = k, candidates = String[first(p) for p in kids])
    end
    nothing
end

# The entry scheduling the child at segment `seg`, contributed by field `fld`, or
# `nothing` for the `Relative(1)` default. Exact match first; then the bare field
# name, which applies one declaration to every element of a container (§8.7). The
# sugar keys on the *field*, so a name-transparent container keeps it unchanged:
# `(children = Relative(2),)` is the uniform spelling for a `Group` (D-211).
function _rate_entry(st, seg::String, fld::Symbol)
    st isa NamedTuple || return nothing         # the shape refusal is already collected
    for (k, v) in pairs(st)
        String(k) == seg && return (k, v)
    end
    for (k, v) in pairs(st)
        String(k) == String(fld) && return (k, v)
    end
    nothing
end

"""
The child's scope triple, and whether an explicit key declared it. The unlisted
child continues at the enclosing triple — the `Relative(1)` default is the
affine law at `(K, φ) = (1, 0)`, implemented by nothing.
"""
function _child_scope(w::_Walk, path::String, st, seg::String, fld::Symbol,
                      scope::NTuple{3,Int})
    hit = _rate_entry(st, seg, fld)
    hit === nothing && return scope, false
    (k, v) = hit
    # `_check_sample_times` collects rather than throwing (§13.1), so the fold
    # runs on past a value outside the wrapper vocabulary — which has no affine
    # law to apply. The violation is already recorded; the child continues at the
    # enclosing scope, as an unlisted one would.
    (v isa Relative || v isa Absolute) || return scope, false
    if v isa Relative
        (a, m, c) = scope
        (a, v.K * m, c + v.φ * m), true
    else
        push!(w.anchors, (v.T, v.τ))
        push!(w.aprov, "`sample_times` at $(_at(path)), key `$k`")
        (length(w.anchors), 1, 0), true
    end
end

"""
    flatten(root)

The tree walk of Stratum A (§9.1): components collected by path, classes read,
wiring resolved to absolute leaf terminals, sample times folded to `(anchor, m,
c)` triples, the one-producer-per-input and whole-tree obligation rules
enforced. Any component may be the root (D-208) — a primitive one flattens to
the single leaf at the root path, its `input_types` keys the model's root
inputs.
"""
function flatten(root)
    w = _Walk(String[], Any[], Dict{Tuple{String,Symbol},Tuple{String,Symbol}}(),
              Dict{Tuple{String,Symbol},String}(), Symbol[],
              Tuple{String,Symbol,Vector{Tuple{String,Symbol}}}[],
              Pair{Tuple{String,Symbol},Tuple{String,Symbol}}[],
              NTuple{3,Int}[], NTuple{2,Rational{Int}}[], String[], Diagnostic[])
    _walk!(w, "", root, (0, 1, 0))          # the root scope: anchor 0, the base grid itself

    # The obligation model (§6.1): an input is fed by a wire in some ancestor's
    # `child_connections` or by an `input_connections` chain handing it up level
    # by level, and the chain that never terminates is the error. The one
    # legitimate unfed terminus is the root's own input face.
    conns = Vector{Pair{Symbol,Tuple{String,Symbol}}}[]
    for (path, c) in zip(w.paths, w.comps)
        cs = Pair{Symbol,Tuple{String,Symbol}}[]
        for face in keys(_contract(input_types, c))
            haskey(w.feeds, (path, face)) ||
                (push!(w.viol, UnconnectedInput(path = path, face = face)); continue)
            push!(cs, face => w.feeds[(path, face)])
        end
        push!(conns, cs)
    end

    # Stratum A's barrier (§13.1): the walk's own refusals — face names, sample
    # times, rate keys, empty routes, two-producer claims — and this pass's
    # unfed inputs leave together, in one throw, before anything derived from
    # the wiring is computed. No cascade suppression: a typo'd wire reports its
    # unknown port *and* the input it left unfed.
    isempty(w.viol) || throw(BuildError(w.viol))

    # §9.2's input side, derived once the obligation pass has proved every input
    # fed exactly once: an assembly's face and the leaf entries behind it are
    # claimed together by the one route above them, so the consumers of a face
    # share a producer and `(path, face) => producer` is well defined. A
    # primitive's own entries complete the record, so the graph carries every
    # input face at every level, whatever the level's class.
    in_faces = Pair{Tuple{String,Symbol},Tuple{String,Symbol}}[]
    for (path, face, consumers) in w.routes
        push!(in_faces, (path, face) => w.feeds[first(consumers)])
    end
    for (path, cs) in zip(w.paths, conns), (face, producer) in cs
        push!(in_faces, (path, face) => producer)
    end
    Flat(w.paths, w.comps, conns, w.root_inputs, in_faces, w.out_faces,
         w.triples, w.anchors, w.aprov)
end

function _walk!(w::_Walk, path::String, comp, scope::NTuple{3,Int})
    if classify(path, comp) === PRIMITIVE
        push!(w.paths, path)
        push!(w.comps, comp)
        push!(w.triples, scope)
        # A primitive at the root: its `input_types` keys are the model's root
        # inputs, each face its own consuming entry (§8.6, §11.3, D-208), fed by
        # the same pseudo-producer an assembly root's faces get.
        if isempty(path)
            _check_root_faces(comp, w.viol)
            for face in keys(_contract(input_types, comp))
                push!(w.root_inputs, face)
                _claim!(w, (path, face), ("", face),
                        "the root component's `input_types` entry `$face`")
            end
        end
        return nothing
    end
    _check_face_names(path, comp, w.viol)
    st = sample_times(comp)
    kids, fields = _children(path, comp)
    _check_sample_times(path, st, kids, fields, w.viol)
    for ((seg, kid), fld) in zip(kids, fields)
        kidpath = _join(path, seg)
        kscope, keyed = _child_scope(w, path, st, seg, fld, scope)
        # A key on a continuous child is the Δt-on-continuous error at
        # declaration time (§8.7): keys name discrete or scope children only.
        keyed && classify(kidpath, kid) === PRIMITIVE &&
            classify_tier(kidpath, kid) === CONTINUOUS &&
            push!(w.viol, RatesViolation(path = path, reason = :continuous_child,
                                         key = Symbol(seg)))
        _walk!(w, kidpath, kid, kscope)
    end

    for pair in child_connections(comp)
        entry = _entry("child_connections", path, pair)
        producer = resolve_source(entry, path, comp, first(pair))
        for consumer in resolve_dest(entry, path, comp, last(pair))
            _claim!(w, consumer, producer, entry)
        end
    end

    # Both boundary declarations are resolved wherever they appear, so their
    # entries are checked at every level; only the root's input faces *feed*
    # anything, there being no parent above them to claim the obligation.
    for (face, inner) in input_connections(comp)
        entry = _entry("input_connections", path, face => inner)
        consumers = _fanout(entry, path, comp, inner)
        # Every entry routes to at least one internal endpoint, at every level
        # (D-210): a face feeding nothing declares nothing, and the empty tuple
        # would otherwise reach no consumer, leave no row in §9.2's face graph,
        # and let a condition addressing it misdiagnose as a bare typo.
        if isempty(consumers)
            push!(w.viol, UnknownPort(entry = entry, end_ = :connection, path = path,
                                      port = Symbol(face)))
            continue                       # a route with no consumer registers nothing
        end
        push!(w.routes, (path, Symbol(face), consumers))
        isempty(path) || continue
        push!(w.root_inputs, Symbol(face))
        for consumer in consumers
            _claim!(w, consumer, ("", Symbol(face)), entry)
        end
    end
    for (src, face) in output_connections(comp)
        entry = _entry("output_connections", path, src => face)
        push!(w.out_faces, (path, Symbol(face)) => resolve_source(entry, path, comp, src))
    end
    nothing
end

_entry(method::String, path::String, pair::Pair) =
    "$method at $(_at(path)), entry `$(repr(first(pair))) => $(repr(last(pair)))`"

# Every input takes exactly one connection, and the rule spans levels (§6.1): an
# input fed both by a sibling wire and by an ancestor's route — or handed up while
# also wired — meets its second claim here.
function _claim!(w::_Walk, consumer, producer, entry::String)
    if haskey(w.feeds, consumer)
        push!(w.viol, TwoProducers(path = consumer[1], port = consumer[2],
                                   incumbent = w.claims[consumer], entry = entry))
        return nothing                     # the incumbent keeps the claim
    end
    w.feeds[consumer] = producer
    w.claims[consumer] = entry
    nothing
end

# §8.6's two face-name invariants. Every other naming choice is author convention.
function _check_face_names(path::String, comp, viol::Vector{Diagnostic})
    names = vcat(String[String(face) for (face, _) in input_connections(comp)],
                 String[String(face) for (_, face) in output_connections(comp)])
    for n in names
        occursin('/', n) && push!(viol, FaceNameIllegal(path = path, face = n))
    end
    allunique(names) ||
        push!(viol, FaceNameCollision(path = path, site = :assembly,
                                      faces = unique(n for n in names
                                                     if count(==(n), names) > 1)))
    nothing
end

# The same invariant at a *primitive* root, whose face set is its two contract
# declarations' keys together (§8.6, D-210). The root is where those two first
# share the periphery's address space: a root input places a cell of its own, so
# a shared key would put two cells at one name and the root input's placement
# would silently overwrite the port's. Below the root nothing collides — a
# primitive's input faces alias their producers' cells and place nothing — and
# non-root leaves are left alone.
function _check_root_faces(comp, viol::Vector{Diagnostic})
    outs = String.(keys(_contract(output_types, comp)))
    dup = [n for n in String.(keys(_contract(input_types, comp))) if n in outs]
    isempty(dup) || push!(viol, FaceNameCollision(path = "", faces = dup, site = :root))
    nothing
end
