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
            throw(BuildError("$(_at(path)) declares `child_connections` and the leaf " *
                             "declaration(s) $(join(leaves, ", ")) — an assembly owns no state " *
                             "and no contract of its own (§8.5)"))
        return ASSEMBLY
    end
    isempty(leaves) || return PRIMITIVE
    throw(BuildError("$(_at(path)) declares neither family: `child_connections` would make it " *
                     "an assembly, any of $LEAF_FAMILY a primitive (§8.5)" *
                     (_holds_components(c) ?
                      " — it holds components but declares no `child_connections`" : "")))
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

"""`segment => instance` for every child of `c`, in field order."""
function children(path::String, c)
    kids = Pair{String,Any}[]
    for name in fieldnames(typeof(c))
        v = getfield(c, name)
        if v isa AbstractComponent
            push!(kids, string(name) => v)
        elseif v isa NamedTuple || v isa Tuple
            n = count(e -> e isa AbstractComponent, v)
            n == 0 && continue                     # inert data; an empty container too
            n == length(v) ||
                throw(BuildError("$(_at(path)): container field `$name` mixes components with " *
                                 "$(join(unique(typeof(e) for e in v
                                                if !(e isa AbstractComponent)), ", ")) " *
                                 "— a container holds components only (§8.5)"))
            for k in keys(v)
                push!(kids, string(name, "/", k) => v[k])
            end
        end
    end
    kids
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
function resolve_terminal(entry::String, base::String, asm, path::AbstractString)
    segs = String.(split(path, '/'))
    length(segs) > 1 ||
        throw(BuildError("$entry: `$path` is not an endpoint — a terminal path names a child " *
                         "and one of its ports or faces (§8.6)"))
    kids = children(base, asm)
    j = findfirst(kid -> first(kid) == segs[1], kids)
    j === nothing && length(segs) > 2 &&
        (j = findfirst(kid -> first(kid) == segs[1] * "/" * segs[2], kids))
    j === nothing &&
        throw(BuildError("$entry: `$path` names no child `$(segs[1])` of $(_at(base))" *
                         (isempty(kids) ? " — it has no children" :
                          " — its children are $(join((first(k) for k in kids), ", "))")))
    seg, kid = kids[j]
    count(==('/'), seg) + 2 == length(segs) ||
        throw(BuildError("$entry: `$path` reaches past `$(_join(base, seg))` — a connection " *
                         "endpoint names an immediate child and one of its faces, so an " *
                         "endpoint path is one child segment (plus the key segment where the " *
                         "child is a container element) and one face name; route through " *
                         "`$seg`'s own face instead, declared level by level (§6.1)"))
    kid, _join(base, seg), Symbol(segs[end])
end

# The key set of a contract declaration, whichever arity declares it: the keys are
# a tier-independent fact, and §8.2's classifier is what settles a disagreement.
_contract(fn, c) = _declares(fn, c, Type{Float64}) ? fn(c, Float64) :
                   _declares(fn, c) ? fn(c) : NamedTuple()

input_faces(c) = classify("", c) === PRIMITIVE ? String.(keys(_contract(input_types, c))) :
                 String[String(face) for (face, _) in input_connections(c)]
output_faces(c) = classify("", c) === PRIMITIVE ? String.(keys(_contract(output_types, c))) :
                  String[String(face) for (_, face) in output_connections(c)]

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
        throw(BuildError("$entry: `$path` names no `$name` on $(_at(cpath)) — its faces are " *
                         "$(join(vcat(ins, outs), ", "))"))
    throw(BuildError("$entry: `$path` resolves to $found of $(_at(cpath)), but this entry's " *
                     "endpoint is a $wanted — direction is declared by the method (§8.6)"))
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
    i === nothing && throw(BuildError("no component at path `$path`"))
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
function _check_sample_times(path::String, st, kids)
    st isa NamedTuple ||
        throw(BuildError("$(_at(path)): `sample_times` must return a NamedTuple of child " *
                         "name => `Relative`/`Absolute` entries (§8.7)"))
    for (k, v) in pairs(st)
        entry = "`sample_times` at $(_at(path)), key `$k`"
        v isa Relative || v isa Absolute ||
            throw(BuildError("$entry: $(repr(v)) — the wrappers are the whole value " *
                             "vocabulary; a bare integer or bare quantity is a declaration " *
                             "error (§8.7)"))
        if v isa Relative
            v.K ≥ 1 || throw(BuildError("$entry: K = $(v.K), and K ≥ 1 (§10.5)"))
            0 ≤ v.φ < v.K || throw(BuildError("$entry: φ = $(v.φ), and 0 ≤ φ < K (§10.5)"))
        else
            v.T > 0 || throw(BuildError("$entry: period $(v.T), and T > 0 (§10.5)"))
            0 ≤ v.τ < v.T || throw(BuildError("$entry: offset $(v.τ), and 0 ≤ τ < T (§10.5)"))
        end
        any(seg == String(k) || startswith(seg, String(k) * "/") for (seg, _) in kids) ||
            throw(BuildError("$entry: names no immediate child of $(_at(path)) — keys are " *
                             "immediate child names only; a deep key would edit another " *
                             "type's design from outside (§8.7)" *
                             (isempty(kids) ? "" :
                              "; its children are $(join((first(p) for p in kids), ", "))")))
    end
    nothing
end

# The entry scheduling child segment `seg`, or `nothing` for the `Relative(1)`
# default. Exact match first; a bare container field name applies one
# declaration to every element (§8.7).
function _rate_entry(st, seg::String)
    for (k, v) in pairs(st)
        String(k) == seg && return (k, v)
    end
    for (k, v) in pairs(st)
        startswith(seg, String(k) * "/") && return (k, v)
    end
    nothing
end

"""
The child's scope triple, and whether an explicit key declared it. The unlisted
child continues at the enclosing triple — the `Relative(1)` default is the
affine law at `(K, φ) = (1, 0)`, implemented by nothing.
"""
function _child_scope(w::_Walk, path::String, st, seg::String, scope::NTuple{3,Int})
    hit = _rate_entry(st, seg)
    hit === nothing && return scope, false
    (k, v) = hit
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
              NTuple{3,Int}[], NTuple{2,Rational{Int}}[], String[])
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
                throw(BuildError("`$path`.$face is fed by nothing — every input is fed exactly " *
                                 "once, by a wire or by an `input_connections` chain ending at a " *
                                 "root input face (§6.1)"))
            push!(cs, face => w.feeds[(path, face)])
        end
        push!(conns, cs)
    end

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
            _check_root_faces(comp)
            for face in keys(_contract(input_types, comp))
                push!(w.root_inputs, face)
                _claim!(w, (path, face), ("", face),
                        "the root component's `input_types` entry `$face`")
            end
        end
        return nothing
    end
    _check_face_names(path, comp)
    st = sample_times(comp)
    kids = children(path, comp)
    _check_sample_times(path, st, kids)
    for (seg, kid) in kids
        kidpath = _join(path, seg)
        kscope, keyed = _child_scope(w, path, st, seg, scope)
        # A key on a continuous child is the Δt-on-continuous error at
        # declaration time (§8.7): keys name discrete or scope children only.
        keyed && classify(kidpath, kid) === PRIMITIVE &&
            classify_tier(kidpath, kid) === CONTINUOUS &&
            throw(BuildError("`sample_times` at $(_at(path)) schedules `$seg`, a continuous " *
                             "component — a sample time is a discrete-tier fact; keys name " *
                             "discrete or scope children (§8.7, §10.5)"))
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
        isempty(consumers) &&
            throw(BuildError("$entry: the entry routes to no internal endpoint — every " *
                             "`input_connections` entry routes to at least one, a face " *
                             "feeding nothing declaring nothing (§8.6)"))
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
        throw(BuildError("`$(consumer[1])`.$(consumer[2]) is fed twice: by " *
                         "$(w.claims[consumer]), and by $entry — every input takes exactly one " *
                         "connection, across levels included (§6.1)"))
    end
    w.feeds[consumer] = producer
    w.claims[consumer] = entry
    nothing
end

# §8.6's two face-name invariants. Every other naming choice is author convention.
function _check_face_names(path::String, comp)
    names = vcat(String[String(face) for (face, _) in input_connections(comp)],
                 String[String(face) for (_, face) in output_connections(comp)])
    for n in names
        occursin('/', n) &&
            throw(BuildError("$(_at(path)): face name `$n` contains `/`, which is reserved for " *
                             "structural paths (§8.6)"))
    end
    allunique(names) ||
        throw(BuildError("$(_at(path)): face name(s) " *
                         "$(join(unique(n for n in names if count(==(n), names) > 1), ", ")) " *
                         "appear twice — face names are unique across `input_connections` and " *
                         "`output_connections` together (§8.6)"))
    nothing
end

# The same invariant at a *primitive* root, whose face set is its two contract
# declarations' keys together (§8.6, D-210). The root is where those two first
# share the periphery's address space: a root input places a cell of its own, so
# a shared key would put two cells at one name and the root input's placement
# would silently overwrite the port's. Below the root nothing collides — a
# primitive's input faces alias their producers' cells and place nothing — and
# non-root leaves are left alone.
function _check_root_faces(comp)
    outs = String.(keys(_contract(output_types, comp)))
    dup = [n for n in String.(keys(_contract(input_types, comp))) if n in outs]
    isempty(dup) ||
        throw(BuildError("$(_at("")): face name(s) $(join(dup, ", ")) appear twice — at the " *
                         "root a primitive's faces are its `input_types` and `output_types` " *
                         "keys together, and a key declared in both is the same build error a " *
                         "duplicate assembly face name is (§8.6)"))
    nothing
end
