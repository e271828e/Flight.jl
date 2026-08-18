# Hierarchy (§8.5, §8.6, §6.1): class by declaration shape, path resolution, and
# the flatten pass. Assemblies are virtual for execution (§10.5) — what leaves
# this file is a flat list of primitives with their absolute paths and one
# resolved producer per input, and nothing downstream knows the tree existed.

# --- class (§8.5) -------------------------------------------------------------
# Class is not announced either: `child_connections` is the assembly marker, any
# leaf declaration a primitive's, and the rule is total — a component-typed
# struct declaring neither family has no class to read.

@enum Class PRIMITIVE ASSEMBLY

const LEAF_FAMILY = "`init_x`, `init_m`, `workspace`, `input_types`, " *
                    "`output_types` or a stage (`h_x`, `h_xu`, `f`, `g`)"

"""The leaf declarations `c` defines, in inventory order (§8.2)."""
function leaf_declarations(c)
    found = Symbol[]
    for (name, fn) in ((:init_x, init_x), (:init_m, init_m))
        _declares(fn, c) && push!(found, name)
    end
    # Either arity is a leaf declaration; which one is lawful is the tier's
    # question, settled by the classifier (§8.2), not this one.
    for (name, fn) in ((:workspace, workspace), (:input_types, input_types),
                       (:output_types, output_types))
        (_declares(fn, c) || _declares(fn, c, Type{Float64})) && push!(found, name)
    end
    for (name, fn) in ((:h_x, h_x), (:h_xu, h_xu), (:f, f), (:g, g))
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
    v isa AbstractComponent || (v isa NamedTuple && any(e -> e isa AbstractComponent, v))
end

# --- children (§8.5) ----------------------------------------------------------
# Fields whose type is `<: AbstractComponent` are children, the field name their
# path segment. A `NamedTuple` field whose elements are *all* components is a
# container: transparent grouping with no contract and no wiring of its own,
# contributing its elements as children *of the parent*, path-named `"field/key"`.
# Every other field is inert parameter data.

"""`segment => instance` for every child of `c`, in field order."""
function children(path::String, c)
    kids = Pair{String,Any}[]
    for name in fieldnames(typeof(c))
        v = getfield(c, name)
        if v isa AbstractComponent
            push!(kids, string(name) => v)
        elseif v isa NamedTuple
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

"""
Is the field `name` of `C` **generically held** (§8.8): a `TypeVar` or a
non-concrete type where the *generic struct* declares it? The judgement is about
the declaration's knowledge, not about the instance in hand — a concrete
instantiation resolves deeper than its declaration is allowed to reach.
"""
function generically_held(C::DataType, name::Symbol)
    FT = fieldtype(Base.unwrap_unionall(C.name.wrapper), name)
    !(FT isa Type) || !isconcretetype(FT)
end

# --- paths (§8.6, §6.1) -------------------------------------------------------
# Slash-separated, relative to the declaring assembly, no leading slash: the one
# canonical form, used verbatim in declarations and in error messages. A terminal
# path's last segment is a port or face name; the prefix walks children.

"""
Resolve terminal `path` against assembly `asm` at `base`, returning the component
it names, that component's absolute path, and the final segment.

§6.1's reach rule lives here: a path may traverse any chain of concretely-declared
component fields and **stops at the first generically-held child**, whose faces are
the only things addressable beyond it. Immediate children are never deep, so a
single hop is always legal — which is what lets a container's elements be
addressed by their parent's own declarations while an ancestor sees only faces.
"""
function resolve_terminal(entry::String, base::String, asm, path::AbstractString)
    segs = String.(split(path, '/'))
    length(segs) > 1 ||
        throw(BuildError("$entry: `$path` is not an endpoint — a terminal path names a child " *
                         "and one of its ports or faces (§8.6)"))
    cur, curpath, i = asm, base, 1
    crossed = Tuple{String,DataType,Symbol}[]      # child path, owner type, field crossed
    while i < length(segs)
        kids = children(curpath, cur)
        j = findfirst(kid -> first(kid) == segs[i], kids)
        consumed = 1
        if j === nothing && i + 1 < length(segs)
            j = findfirst(kid -> first(kid) == segs[i] * "/" * segs[i+1], kids)
            consumed = 2
        end
        j === nothing &&
            throw(BuildError("$entry: `$path` names no child `$(segs[i])` of $(_at(curpath))" *
                             (isempty(kids) ? " — it has no children" :
                              " — its children are $(join((first(k) for k in kids), ", "))")))
        seg, kid = kids[j]
        push!(crossed, (_join(curpath, seg), typeof(cur), Symbol(segs[i])))
        cur, curpath, i = kid, _join(curpath, seg), i + consumed
    end
    if length(crossed) > 1
        for (childpath, C, field) in crossed
            generically_held(C, field) &&
                throw(BuildError("$entry: `$path` reaches past `$childpath`, which is held " *
                                 "generically — a path traverses concretely-declared fields only " *
                                 "and stops at the first generically-held child, whose faces are " *
                                 "the only things addressable beyond it (§6.1)"))
        end
    end
    cur, curpath, Symbol(segs[end])
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
per declared input, the root's input faces (the [root slots](§11.3)) and the
assembly faces the periphery may read, aliased onto the cells they derive from.
"""
struct Flat
    paths::Vector{String}
    comps::Vector{Any}
    conns::Vector{Vector{Pair{Symbol,Tuple{String,Symbol}}}}   # face => (producer path, port)
    slots::Vector{Symbol}                                      # root input faces, in order
    faces::Vector{Pair{Tuple{String,Symbol},Tuple{String,Symbol}}}   # (path, face) => producer
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
    slots::Vector{Symbol}
    faces::Vector{Pair{Tuple{String,Symbol},Tuple{String,Symbol}}}
end

"""
    flatten(root)

The tree walk of Stratum A (§9.1): components collected by path, classes read,
wiring resolved to absolute leaf terminals, the one-producer-per-input and
whole-tree obligation rules enforced.
"""
function flatten(root)
    classify("", root) === PRIMITIVE &&
        throw(BuildError("the root component is a primitive — a model's root is an assembly, " *
                         "whose input faces are the root slots (§9.1)"))
    w = _Walk(String[], Any[], Dict{Tuple{String,Symbol},Tuple{String,Symbol}}(),
              Dict{Tuple{String,Symbol},String}(), Symbol[],
              Pair{Tuple{String,Symbol},Tuple{String,Symbol}}[])
    _walk!(w, "", root)

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
    Flat(w.paths, w.comps, conns, w.slots, w.faces)
end

function _walk!(w::_Walk, path::String, comp)
    if classify(path, comp) === PRIMITIVE
        push!(w.paths, path)
        push!(w.comps, comp)
        return nothing
    end
    _check_face_names(path, comp)
    for (seg, kid) in children(path, comp)
        _walk!(w, _join(path, seg), kid)
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
        isempty(path) || continue
        push!(w.slots, Symbol(face))
        for consumer in consumers
            _claim!(w, consumer, ("", Symbol(face)), entry)
        end
    end
    for (src, face) in output_connections(comp)
        entry = _entry("output_connections", path, src => face)
        push!(w.faces, (path, Symbol(face)) => resolve_source(entry, path, comp, src))
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
