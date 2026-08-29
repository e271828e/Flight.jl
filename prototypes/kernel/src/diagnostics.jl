# The diagnostic layer (§13.1, §13.2, Appendix C, D-058, D-214, D-215): the
# closed kind set the build and the services raise, and the one carrier that
# moves it. A diagnostic is a plain value whose *type* is its kind and whose
# fields are Appendix C's payload column — paths and names as strings and
# symbols, never component instances and never model types; the declared and
# observed *port* types are the payload exception, and they are small. Severity
# is a property of the kind, read as `severity(d)` and never stored per
# occurrence (§13.2 as amended, D-214). Where an occurrence surfaces and how it
# is reported are facts of the site, not of the value, so no kind carries them.
#
# Messages are presentation: `message(d)` renders one diagnostic in the didactic
# register (state the fix, show the list-in-hand) and carries no kind name —
# the carrier's `showerror` leads each line with it. Tests match on kind plus
# payload, never on message text (§13.2).
#
# The runtime warning stream's eight kinds live in `src/dataplane.jl` beside the
# cells that carry them; they are parented here so `severity` covers them too.

"""
The root of the closed kind set (§13.2, Appendix C). Every kind is an immutable
struct under it, with three methods: `severity(d)`, `path(d)` and `message(d)`.
"""
abstract type Diagnostic end

"""
The kind's severity (§13.2, D-214): `:error` — an occurrence throws, alone or
within a collection — or `:warning`, which never throws and joins no throw. The
default is `:error`; the warning kinds override it.
"""
severity(::Diagnostic) = :error

"""
The component or assembly path a diagnostic is attributed to, in slash form and
without decoration — the renderer's sort key within a kind group. Kinds that
name no path answer `""`.
"""
path(::Diagnostic) = ""

# `message(d)::String` is the third method; it has no default, so a kind added
# without one fails loudly at its first rendering rather than printing a stub.

# --- the shared renderings ----------------------------------------------------
# The prototype's own spellings, lifted out of the sites so the kinds can render
# without them: `_at`'s path decoration (assembly.jl), `_names`' list-in-hand
# (conditions.jl) and `_facelist`'s face set (dataplane.jl).

_at_path(p::AbstractString) = isempty(p) ? "the root component" : "`$p`"
_namelist(ns) = isempty(ns) ? "none" : join(("`$n`" for n in ns), ", ")
_faceset(ns) = isempty(ns) ? "empty" : "{$(join(ns, ", "))}"
_plainlist(ns) = join(ns, ", ")
_symtuple(ns) = "(" * join((":$n" for n in ns), ", ") * (length(ns) == 1 ? ",)" : ")")

# --- the carrier (§13.1, §13.2, D-058) ----------------------------------------

"""
The one carrier: a fail-fast site throws it holding a single diagnostic, a
stratum barrier holding the whole collection its passes returned (§13.1). Its
`showerror` renders compiler-style — grouped by kind, sorted by path within a
group, the kind name leading each line.
"""
struct BuildError <: Exception
    diagnostics::Vector{Diagnostic}
end

BuildError(d::Diagnostic) = BuildError(Diagnostic[d])
BuildError(ds::AbstractVector{<:Diagnostic}) = BuildError(Vector{Diagnostic}(ds))

"The kinds present in a carrier, in first-appearance order — what a test asks first."
kinds(e::BuildError) = unique(typeof.(e.diagnostics))

# Groups in first-appearance order, each sorted by path; the sort is stable, so
# two diagnostics at one path keep the order the pass produced them in.
function _groups(ds::Vector{Diagnostic})
    order = DataType[]
    for d in ds
        T = typeof(d)
        T in order || push!(order, T)
    end
    [sort(filter(d -> typeof(d) === T, ds); by = path, alg = MergeSort) for T in order]
end

function Base.showerror(io::IO, e::BuildError)
    ds = e.diagnostics
    if length(ds) == 1
        d = only(ds)
        print(io, "BuildError: ", nameof(typeof(d)), ": ", message(d))
        return nothing
    end
    print(io, "BuildError: ", length(ds), " diagnostics")
    for g in _groups(ds), d in g
        print(io, "\n  ", nameof(typeof(d)), ": ", message(d))
    end
    nothing
end

"""
A warning-severity kind's rendering for the log (§13.2, D-214's `logged`
policy): the same kind-name-leading line the carrier prints, minus the carrier
— a logged warning never throws, so it has no `showerror` to lead it. Named
`logline` rather than `logged` because the policy's name is taken here by
`logged(sim)`, the retained-snapshot accessor (§11.4), and one generic function
over both would be two unrelated meanings sharing a name.
"""
logline(d::Diagnostic) = string(nameof(typeof(d)), ": ", message(d))

"""
An internal invariant firing (D-215): not a diagnostic and not a kind, because
it names no failure the user can fix. Its own exception type so the assertions
stay outside the acceptance-test contract.
"""
struct InternalInvariant <: Exception
    msg::String
end

Base.showerror(io::IO, e::InternalInvariant) =
    print(io, "InternalInvariant: internal invariant violated: ", e.msg)

# ==============================================================================
# Stratum A — declaration and wiring (§6.1, §8.2, §8.5–§8.8; collected)
# ==============================================================================

"§6.1, §8.4 w1: a wire end naming no port of the endpoint it resolved to, with that end's port list."
Base.@kwdef struct UnknownPort <: Diagnostic
    entry::String                            # the declaring method and entry: provenance
    end_::Symbol                             # :source | :destination | :connection (D-210)
    path::String = ""                        # the component that end resolved to
    spelling::String = ""                    # the endpoint path as the entry wrote it
    port::Union{Nothing,Symbol} = nothing    # the unknown port or face name
    candidates::Vector{Symbol} = Symbol[]    # that end's port list, the did-you-mean
end
path(d::UnknownPort) = d.path
message(d::UnknownPort) =
    d.end_ === :connection ?
    "$(d.entry): the entry routes to no internal endpoint — every `input_connections` " *
    "entry routes to at least one, a face feeding nothing declaring nothing (§8.6)" :
    "$(d.entry): `$(d.spelling)` names no `$(d.port)` on $(_at_path(d.path)) — its " *
    "faces are $(_plainlist(d.candidates))"

"§6.1, §8.4 w2: an input no wire and no `input_connections` chain feeds."
Base.@kwdef struct UnconnectedInput <: Diagnostic
    path::String
    face::Symbol
end
path(d::UnconnectedInput) = d.path
message(d::UnconnectedInput) =
    "`$(d.path)`.$(d.face) is fed by nothing — every input is fed exactly once, by a " *
    "wire or by an `input_connections` chain ending at a root input face (§6.1)"

"§6.1, §8.8: an input claimed twice, both producers named with their provenance."
Base.@kwdef struct TwoProducers <: Diagnostic
    path::String                             # the destination terminal's component
    port::Symbol                             # the destination terminal's port
    incumbent::String                        # the entry that claimed it first
    entry::String                            # the entry that claimed it second
end
path(d::TwoProducers) = d.path
message(d::TwoProducers) =
    "`$(d.path)`.$(d.port) is fed twice: by $(d.incumbent), and by $(d.entry) — every " *
    "input takes exactly one connection, across levels included (§6.1)"

"§6.1, §8.2, §8.4 w4: a wire whose producer's value does not meet the consumer's declared entry type."
Base.@kwdef struct WireTypeMismatch <: Diagnostic
    path::String                             # the consumer
    face::Symbol
    declared::Any                            # the declared entry type
    producer_path::String                    # "" for a root input
    producer_port::Symbol
    observed::Any                            # the producer face type
    activation::Any = nothing                # the activation scalar, for the pin hint
end
path(d::WireTypeMismatch) = d.path
message(d::WireTypeMismatch) =
    "`$(d.path)`.$(d.face) declared $(d.declared), fed from " *
    (isempty(d.producer_path) ? "root input `$(d.producer_port)`" :
     "`$(d.producer_path)`.$(d.producer_port)") *
    "::$(d.observed)" * _pin(d)

# The embed-accept hint (build.jl's `_pin_hint`) is rendering over the two port
# types and the activation scalar, so it is computed here rather than carried.
_pin(d) = (d.declared isa Type && d.observed isa Type && d.activation isa Type) ?
          _pin_hint(d.declared, d.observed, d.activation) : ""

"""
§8.2: two consumers of one root input face declaring different entry types.
A root input takes its type from the consumer declaration alone, so under
fan-out the concrete declaration has to be unique. The comparison is *at
nominal*, which is what makes a tolerance difference no conflict: `SVector{3,T}`
and `SVector{3,Float64}` both evaluate to `SVector{3,Float64}` there, and their
disagreement about partials is the fan-out meet (D-168), a legitimate model.
"""
Base.@kwdef struct RootInputTypeConflict <: Diagnostic
    face::Symbol                             # the root input face
    paths::Vector{String}                    # the consuming components
    declared::Vector{Any}                    # their entry declarations at nominal
end
message(d::RootInputTypeConflict) =
    "root input `$(d.face)` is declared " *
    join(("$(_at_path(p))::$(P)" for (p, P) in zip(d.paths, d.declared)), ", ") *
    " — a root input is typed by its consumers alone, so its concrete declaration " *
    "has to be unique across the fan-out; agree the entries, or give the face a " *
    "producer (§8.2, D-168)"

"§6.1, §13.3: a path that resolves to nothing, or reaches past the one level a register admits."
Base.@kwdef struct PathResolution <: Diagnostic
    entry::String                            # the register or wiring entry: provenance
    spelling::String                         # the path as written
    reason::Symbol                           # :not_a_terminal|:unknown_child|:reaches_past|:empty_path
    owner::String = ""                       # the component the path was resolved against
    segment::String = ""                     # the offending segment
    level::String = ""                       # the level it stopped at
    candidates::Vector{String} = String[]    # the sibling child names
    tail::Int = 0                            # 1 for a wiring endpoint, 0 for a bare child path
end
path(d::PathResolution) = d.spelling
function message(d::PathResolution)
    d.reason === :not_a_terminal &&
        return "$(d.entry): `$(d.spelling)` is not an endpoint — a terminal path names a " *
               "child and one of its ports or faces (§8.6)"
    d.reason === :empty_path &&
        return "$(d.entry): the empty path names no child — a path in this register names " *
               "an immediate child of the component in hand (§13.3)"
    d.reason === :unknown_child &&
        return "$(d.entry): `$(d.spelling)` names no child `$(d.segment)` of $(d.owner)" *
               (isempty(d.candidates) ? " — it has no children" :
                " — its children are $(_plainlist(d.candidates))")
    "$(d.entry): `$(d.spelling)` reaches past `$(d.level)` — " *
    (d.tail == 1 ?
     "a connection endpoint names an immediate child and one of its faces, so an " *
     "endpoint path is one child segment (plus the key segment where the child is a " *
     "container element) and one face name; route through `$(d.segment)`'s own face " *
     "instead, declared level by level (§6.1)" :
     "a path in this register names an immediate child: one child segment, plus the key " *
     "segment where the child is a container element, and nothing further (§6.1, §13.3)")
end

"§8.2: a declared store with no update — `init_x` without `f`, `init_s` without `g`."
Base.@kwdef struct StoreWithoutUpdate <: Diagnostic
    path::String
    store::Symbol                            # :init_x | :init_s
end
path(d::StoreWithoutUpdate) = d.path
message(d::StoreWithoutUpdate) =
    "`$(d.path)` declares `$(d.store)` but defines neither `f` nor `g` — a store needs " *
    "its update (§8.2)"

"§8.2: an event declared with one half, or an `events` entry that is not an `Event` (D-215)."
Base.@kwdef struct EventHalfMissing <: Diagnostic
    path::String
    event::Symbol
    reason::Symbol                           # :guard | :handler | :not_an_event
    found::Any = nothing                     # the component type, or the entry's type
end
path(d::EventHalfMissing) = d.path
message(d::EventHalfMissing) =
    d.reason === :not_an_event ?
    "`$(d.path)`: events entry `$(d.event)` is $(d.found) — an entry is " *
    "`Event(guard, handler)`, with no detection keyword (§8.2)" :
    "`$(d.path)`: event `$(d.event)`'s $(d.reason) has no method for $(d.found) — an " *
    "event needs both halves (§8.2)"

"§8.5: a component declaring neither family, so its class cannot be read off declaration shape."
Base.@kwdef struct ClassUnreadable <: Diagnostic
    path::String
    families::String                         # the leaf family list, the list-in-hand
    holds_components::Bool = false
end
path(d::ClassUnreadable) = d.path
message(d::ClassUnreadable) =
    "$(_at_path(d.path)) declares neither family: `child_connections` would make it an " *
    "assembly, any of $(d.families) a primitive (§8.5)" *
    (d.holds_components ?
     " — it holds components but declares no `child_connections`" : "")

"§8.5: a component declaring both families — an assembly owns no state and no contract."
Base.@kwdef struct ClassMixed <: Diagnostic
    path::String
    declarations::Vector{Symbol}             # the offending leaf declarations
end
path(d::ClassMixed) = d.path
message(d::ClassMixed) =
    "$(_at_path(d.path)) declares `child_connections` and the leaf declaration(s) " *
    "$(_plainlist(d.declarations)) — an assembly owns no state and no contract of its " *
    "own (§8.5)"

"§8.5: a container field mixing components with plain data."
Base.@kwdef struct ContainerMixed <: Diagnostic
    path::String
    field::Symbol
    types::Vector{Any}                       # the non-component element types
end
path(d::ContainerMixed) = d.path
message(d::ContainerMixed) =
    "$(_at_path(d.path)): container field `$(d.field)` mixes components with " *
    "$(_plainlist(d.types)) — a container holds components only (§8.5)"

"§5.2, §8.2, §8.5: a declaration written in the other tier's form, or `project` off the continuous tier."
Base.@kwdef struct DeclarationOnWrongTier <: Diagnostic
    path::String
    declaration::Symbol                      # the offending declaration
    reason::Symbol                           # :tier_form | :continuous_only | :no_manifold
    found::Union{Nothing,Symbol} = nothing   # the tier the declaration is written in
    announced::Union{Nothing,Symbol} = nothing  # the tier the other declarations announce
end
path(d::DeclarationOnWrongTier) = d.path
message(d::DeclarationOnWrongTier) =
    d.reason === :continuous_only ?
    "`$(d.path)` declares `$(d.declaration)`, which is continuous-only — projection " *
    "normalizes continuous state (§5.2)" :
    d.reason === :no_manifold ?
    "`$(d.path)` declares `$(d.declaration)` but no `init_x` — there is no state " *
    "manifold to project onto (§5.2)" :
    "`$(d.path)`: `$(d.declaration)` is declared in the $(d.found)-tier form, but this " *
    "component's other declarations announce the $(d.announced) tier (§8.2)"

"§8.6: a face name holding `/`, the separator reserved for structural paths."
Base.@kwdef struct FaceNameIllegal <: Diagnostic
    path::String
    face::String
    invariant::Symbol                        # the invariant broken; :contains_slash is §8.6's
end
path(d::FaceNameIllegal) = d.path
message(d::FaceNameIllegal) =
    d.invariant === :contains_slash ?
    "$(_at_path(d.path)): face name `$(d.face)` contains `/`, which is reserved for " *
    "structural paths (§8.6)" :
    "$(_at_path(d.path)): face name `$(d.face)` breaks the `$(d.invariant)` invariant (§8.6)"

"§8.6: one face name declared twice — across an assembly's two methods, or at a primitive root."
Base.@kwdef struct FaceNameCollision <: Diagnostic
    path::String
    faces::Vector{String}
    site::Symbol                             # :assembly | :root
end
path(d::FaceNameCollision) = d.path
message(d::FaceNameCollision) =
    d.site === :root ?
    "$(_at_path("")): face name(s) $(_plainlist(d.faces)) appear twice — at the root a " *
    "primitive's faces are its `input_types` and `output_types` keys together, and a key " *
    "declared in both is the same build error a duplicate assembly face name is (§8.6)" :
    "$(_at_path(d.path)): face name(s) $(_plainlist(d.faces)) appear twice — face names " *
    "are unique across `input_connections` and `output_connections` together (§8.6)"

"§8.6: an entry whose endpoint resolves to a port of the opposite direction."
Base.@kwdef struct FaceDirectionConflict <: Diagnostic
    entry::String                            # the declaring method and entry
    path::String                             # the child the endpoint resolved to
    spelling::String                         # the endpoint path as written
    found::Symbol                            # :input | :output — the resolved direction
    wanted::Symbol                           # :producer | :consumer — the entry's side
end
path(d::FaceDirectionConflict) = d.path
message(d::FaceDirectionConflict) =
    "$(d.entry): `$(d.spelling)` resolves to " *
    (d.found === :input ? "an input" : "an output") * " of $(_at_path(d.path)), but this " *
    "entry's endpoint is a $(d.wanted) — direction is declared by the method (§8.6)"

"§8.8: a passthrough filter naming faces the child does not have, or giving both `except` and `only`."
Base.@kwdef struct UnknownFaceSelection <: Diagnostic
    who::String                              # the calling helper
    path::String                             # the child path
    reason::Symbol                           # :both_given | :unknown_names
    names::Vector{String} = String[]         # the offending names
    candidates::Vector{String} = String[]    # the child's face list
end
path(d::UnknownFaceSelection) = d.path
message(d::UnknownFaceSelection) =
    d.reason === :both_given ?
    "`$(d.who)` at `$(d.path)`: `except` and `only` were both given — one names the faces " *
    "to drop, the other the faces to keep, and they are mutually exclusive (§8.8)" :
    "`$(d.who)` at `$(d.path)`: $(_namelist(d.names)) " *
    "$(length(d.names) == 1 ? "names" : "name") no face of that child — its faces are " *
    "$(_namelist(d.candidates)) (§8.8)"

"§8.7, §10.5: a `sample_times` declaration outside the value vocabulary, the residue bounds or the key rule."
Base.@kwdef struct RatesViolation <: Diagnostic
    path::String
    reason::Symbol   # :declaration_shape|:value_vocabulary|:multiplier|:phase|:period|:offset|:unknown_child|:continuous_child
    key::Union{Nothing,Symbol} = nothing     # the offending key
    value::Any = nothing                     # the offending value
    candidates::Vector{String} = String[]    # the immediate child names
end
path(d::RatesViolation) = d.path
_rates_entry(d) = "`sample_times` at $(_at_path(d.path)), key `$(d.key)`"
function message(d::RatesViolation)
    d.reason === :declaration_shape &&
        return "$(_at_path(d.path)): `sample_times` must return a NamedTuple of child " *
               "name => `Relative`/`Absolute` entries (§8.7)"
    d.reason === :value_vocabulary &&
        return "$(_rates_entry(d)): $(repr(d.value)) — the wrappers are the whole value " *
               "vocabulary; a bare integer or bare quantity is a declaration error (§8.7)"
    d.reason === :multiplier && return "$(_rates_entry(d)): K = $(d.value), and K ≥ 1 (§10.5)"
    d.reason === :phase && return "$(_rates_entry(d)): φ = $(d.value), and 0 ≤ φ < K (§10.5)"
    d.reason === :period && return "$(_rates_entry(d)): period $(d.value), and T > 0 (§10.5)"
    d.reason === :offset && return "$(_rates_entry(d)): offset $(d.value), and 0 ≤ τ < T (§10.5)"
    d.reason === :continuous_child &&
        return "`sample_times` at $(_at_path(d.path)) schedules `$(d.key)`, a continuous " *
               "component — a sample time is a discrete-tier fact; keys name discrete or " *
               "scope children (§8.7, §10.5)"
    "$(_rates_entry(d)): names no immediate child of $(_at_path(d.path)) — keys are " *
    "immediate child names only; a deep key would edit another type's design from " *
    "outside (§8.7)" *
    (isempty(d.candidates) ? "" : "; its children are $(_plainlist(d.candidates))")
end

"§8.5, D-211, D-212: two children under one name — a bare container key against the sugar, a sibling field, or a plain duplicate."
Base.@kwdef struct ChildNameCollision <: Diagnostic
    path::String
    name::String                             # the colliding child name
    reason::Symbol                           # :sample_times_sugar | :sibling_field | :two_children
    provenance::Vector{String} = String[]    # one entry for the bare-key arms, two for the duplicate
    field::Union{Nothing,Symbol} = nothing   # the container field the key came from
end
path(d::ChildNameCollision) = d.path
_prov(d, i) = i ≤ length(d.provenance) ? d.provenance[i] : "an undetermined declaration"
message(d::ChildNameCollision) =
    d.reason === :sample_times_sugar ?
    "$(_at_path(d.path)): the bare key `$(d.name)` — $(_prov(d, 1)) — collides with " *
    "`sample_times`' field-name sugar, which spells one declaration for every element of " *
    "`$(d.field)` under that same name (§8.5, §8.7, D-211)" :
    d.reason === :sibling_field ?
    "$(_at_path(d.path)): the bare key `$(d.name)` — $(_prov(d, 1)) — collides with " *
    "container field `$(d.name)`, whose own children are named `$(d.name)/<key>`: no " *
    "child bears the bare name, but the segment grammar that reaches those children does, " *
    "and the key shadows it — leaving them unreachable behind a diagnostic naming the " *
    "wrong child (§8.5, §6.1, D-212)" :
    "$(_at_path(d.path)): two children are named `$(d.name)` — $(_prov(d, 1)) and " *
    "$(_prov(d, 2)); a child name is a path segment, and a path segment addresses one " *
    "component (§8.5, D-211)"

"§8.5, D-211: `transparent_container` naming no container field of the type."
Base.@kwdef struct TransparentContainerUnknown <: Diagnostic
    path::String
    field::Symbol
    component::String                        # the declaring type, as a string
end
path(d::TransparentContainerUnknown) = d.path
message(d::TransparentContainerUnknown) =
    "$(_at_path(d.path)): `transparent_container` returns `:$(d.field)`, which names no " *
    "container field of `$(d.component)` — a name-transparent declaration names a `Tuple` " *
    "or `NamedTuple` field whose elements are all components, the empty one included " *
    "(§8.5, D-211)"

"§5.2, §8.2, §8.5, D-215: the tier twin of `ClassUnreadable` — no `output_types` and no state."
Base.@kwdef struct TierUnreadable <: Diagnostic
    path::String
    declarations::Vector{Symbol} = Symbol[]  # the tier-announcing declarations found
end
path(d::TierUnreadable) = d.path
message(d::TierUnreadable) =
    "`$(d.path)` declares no `output_types` and owns no state — there is nothing for the " *
    "tier to be read off: declare `output_types`, or give the component an `init_x`/`init_s` " *
    "store with the update law that drives it. Its tier-announcing declarations are " *
    "$(_namelist(d.declarations)) (§8.2)"

"§7.1, §8.2, D-215: the port twin of `IllegalStateLeaf` — a declared type with no numeric leaves."
Base.@kwdef struct IllegalPortType <: Diagnostic
    path::String
    site::Symbol                             # :port | :root_input
    name::Symbol
    declared::Any                            # the offending type
end
path(d::IllegalPortType) = d.path
message(d::IllegalPortType) =
    "$(_at_path(d.path)): $(d.site === :root_input ? "root input" : "port") `$(d.name)` " *
    "declares $(d.declared), which has no leaves"

# ==============================================================================
# Strata B and C — schedule and contract conformance (§5.5, §8.3, §9.3, §9.5)
# ==============================================================================

"§5.5, §5.6: a strongly connected component of the stage-2 port graph."
Base.@kwdef struct AlgebraicCycle <: Diagnostic
    members::Vector{String}                  # the SCC's member terminals, in slash form
end
path(d::AlgebraicCycle) = isempty(d.members) ? "" : first(d.members)
message(d::AlgebraicCycle) =
    "algebraic loop through stage-2 ports: $(join(d.members, " → ")) — break it with a " *
    "stage-1 (`h_x`/`h_s`) port, which carries no input dependence (§5.4/§5.5)"

"§4.3, §8.3: one port written by both stages."
Base.@kwdef struct ProducedByTwoStages <: Diagnostic
    path::String
    ports::Vector{Symbol}
end
path(d::ProducedByTwoStages) = d.path
message(d::ProducedByTwoStages) =
    "`$(d.path)`: $(_plainlist(d.ports)) produced by two stages"

"§8.3: a declared port no stage writes — a cell no one fills, reading as a silent zero."
Base.@kwdef struct DeclaredNotProduced <: Diagnostic
    path::String
    ports::Vector{Symbol}
    products::Vector{Symbol} = Symbol[]      # the stage-product list
end
path(d::DeclaredNotProduced) = d.path
message(d::DeclaredNotProduced) =
    "`$(d.path)`: declared port(s) $(_plainlist(d.ports)) produced by no stage — either a " *
    "stage returns them or `output_types` drops them; the stages return " *
    "$(_namelist(d.products))"

"§8.3, §8.4 w5: a stage returning a field `output_types` does not declare."
Base.@kwdef struct UndeclaredReturnField <: Diagnostic
    path::String
    stage::String
    name::Symbol
    candidates::Vector{Symbol} = Symbol[]    # `output_types`' keys
end
path(d::UndeclaredReturnField) = d.path
message(d::UndeclaredReturnField) =
    "`$(d.path)`: $(d.stage) returns `$(d.name)`, which `output_types` does not declare — " *
    "declare it, or drop it from the return; the declared ports are $(_namelist(d.candidates))"

"""
§9.5: the return laws, probed. `shape` names what the return had to be shaped
like and `reason` which half of the law failed — the return's own type, its
field set, or one field's type.
"""
Base.@kwdef struct ConformanceFailure <: Diagnostic
    path::String
    what::String                             # the function or stage at fault
    reason::Symbol                           # :return_type | :field_set | :field_type
    shape::Symbol                            # :ports|:namedtuple|:state|:init_x|:init_s|:stores|:mode
    observed::Any = nothing
    declared::Any = nothing
    field::Union{Nothing,Symbol} = nothing
    observed_fields::Vector{Symbol} = Symbol[]
    declared_fields::Vector{Symbol} = Symbol[]
    activation::Any = nothing                # the activation scalar, for the pin hint
end
path(d::ConformanceFailure) = d.path

_cf_expect(s::Symbol) =
    s === :ports      ? "must return a NamedTuple of port values" :
    s === :state      ? "must return a NamedTuple shaped like the state" :
    s === :init_x     ? "must return a NamedTuple shaped like `init_x`" :
    s === :init_s     ? "must return a NamedTuple shaped like `init_s`" :
    s === :stores     ? "must return a NamedTuple of the stores it writes" :
    s === :mode       ? "must be a NamedTuple" :
                        "must return a NamedTuple"
_cf_section(s::Symbol) = (s === :stores || s === :mode) ? " (§5.2)" : ""

function message(d::ConformanceFailure)
    d.reason === :return_type &&
        return "`$(d.path)`: $(d.what) $(_cf_expect(d.shape)), got $(d.observed)" *
               _cf_section(d.shape)
    if d.reason === :field_set
        d.shape === :mode &&
            return "`$(d.path)`: $(d.what) writes mode `$(d.field)`, and `init_m` declares " *
                   "$(_symtuple(d.declared_fields)) — `m` is a names-subset write (§5.2)"
        d.shape === :init_s &&
            return "`$(d.path)`: $(d.what) returns $(d.observed), state store is " *
                   "$(d.declared) — a discrete successor is the store's own type exactly (§7.3)"
        d.shape === :init_x &&
            return "`$(d.path)`: $(d.what) returns fields $(_symtuple(d.observed_fields)), " *
                   "state has $(_symtuple(d.declared_fields)) — derivative completeness is " *
                   "structural (§7.1)"
        return "`$(d.path)`: $(d.what) returns fields $(_symtuple(d.observed_fields)), " *
               "state has $(_symtuple(d.declared_fields)) — a state write-back is complete " *
               "against the field set (§9.3, §9.5)"
    end
    d.shape === :ports &&
        return "`$(d.path)`: $(d.what) returns `$(d.field)`::$(d.observed), declared " *
               "$(d.declared)" * _pin(d)
    d.shape === :mode &&
        return "`$(d.path)`: $(d.what) mode `$(d.field)` is $(d.observed), declared " *
               "$(d.declared) (§5.2)"
    d.shape === :init_x &&
        return "`$(d.path)`: derivative field `$(d.field)` is $(d.observed), state field " *
               "is $(d.declared)"
    "`$(d.path)`: $(d.what) field `$(d.field)` is $(d.observed), state field is $(d.declared)"
end

"§9.5: a guard returning neither of the two admissible forms."
Base.@kwdef struct GuardForm <: Diagnostic
    path::String
    event::Symbol
    observed::Any
end
path(d::GuardForm) = d.path
message(d::GuardForm) =
    "`$(d.path)`: event `$(d.event)`'s guard returns $(d.observed) — a guard is " *
    "`Bool`-valued (boundary-detected) or returns the continuous sign value (localized) " *
    "(§2.1, §10.4)"

"§5.2, §9.5: a handler key naming no store the component declares."
Base.@kwdef struct HandlerReturnKey <: Diagnostic
    path::String
    event::Symbol
    key::Symbol
    stores::Vector{Symbol} = Symbol[]        # `{x, m}` narrowed to what exists
end
path(d::HandlerReturnKey) = d.path
message(d::HandlerReturnKey) =
    "`$(d.path)`: event `$(d.event)`'s handler returns `$(d.key)` — a handler's keys name " *
    "the stores it writes, and this component's are " *
    (isempty(d.stores) ? "none" : _namelist(d.stores)) * " (§5.2)"

# ==============================================================================
# Deployment, periphery and services (§9.1, §11, §12, §13.5, §14)
# ==============================================================================

"§12.6: an advance entry called before `init!` has run boundary zero."
Base.@kwdef struct MissingInit <: Diagnostic
    op::Symbol                               # the entry point called
    status::Symbol                           # the simulation's status
end
message(d::MissingInit) =
    "`$(d.op)` before `init!`: boundary zero has not run and this simulation is " *
    "`$(d.status)` — `init!` is mandatory (§12.6)"

"§11.3, §14: a service call against a lifecycle status that does not admit it."
Base.@kwdef struct ServiceLifecycle <: Diagnostic
    op::Symbol
    status::Symbol
    legal::Vector{Symbol} = Symbol[]         # the statuses that admit the operation
end

# `:running` refuses two different operations, and the sentence differs: an
# advance entry is refused *because the loop is already advancing*, while every
# other refusal at this status is a stopped-sim operation meeting a running loop.
_advance_entry(op::Symbol) = op === :run! || op === :step!

message(d::ServiceLifecycle) =
    d.status === :running ?
    (_advance_entry(d.op) ?
     "`$(d.op)`: the simulation is already running — one advance owns the stores between " *
     "drains, and a second one from another thread would race it (§12.6)" :
     "`$(d.op)` is a stopped-sim operation and the simulation is running — the roster is " *
     "frozen per run and the loop owns the stores between drains (§11.3, §12.5, §12.6)") :
    d.status === :errored ?
    "`$(d.op)` on a simulation that ended `errored` — terminally stopped, never resumable " *
    "or re-initialized; reproduction is trace replay, absent here (§13.6)" :
    d.status === :stopped ?
    "`$(d.op)` on a stopped simulation: re-running is the `stopped → init! → run!` cycle " *
    "— `init!` re-runs boundary zero and opens a fresh trajectory (§12.6)" :
    "`$(d.op)` is legal in $(_namelist(d.legal)), and this simulation is `$(d.status)` " *
    "(§12.6, §14)"

"§13.5: a `stop_on` name that is no root-exported `Bool` output face."
Base.@kwdef struct StopFaceInvalid <: Diagnostic
    face::Symbol
    reason::Symbol                           # :unknown | :root_input | :not_bool
    declared::Any = nothing                  # the declared type, for :not_bool
    candidates::Vector{Symbol} = Symbol[]    # the root output-face list
end
message(d::StopFaceInvalid) =
    d.reason === :unknown ?
    "stop_on names `$(d.face)`, which is no root face — a stop face is a root-exported " *
    "Bool output face, and this model exports $(_namelist(d.candidates)) (§13.5)" :
    d.reason === :root_input ?
    "stop_on names `$(d.face)`, a root input — a stop face is a root-exported *output*: " *
    "the model detects and exports the condition, and the deployment names it (§13.5)" :
    "stop_on names `$(d.face)`, whose declared type is $(d.declared) — stop faces are " *
    "Bool, OR-combined (§13.5)"

"§9.1: a deployment parameter outside its constraint, or a grid that does not close."
Base.@kwdef struct DeploymentInvalid <: Diagnostic
    parameter::Symbol
    reason::Symbol   # :range|:inexact|:not_a_quantity|:missing|:unanchored|:no_constraint|:not_harmonic|:disagrees_with_n|:anchor_period|:anchor_offset
    value::Any = nothing
    related::Any = nothing                   # the value it is measured against
    quotient::Union{Nothing,Int} = nothing   # Δt_base / h, where that is the fact
    paths::Vector{String} = String[]         # the unanchored components
    provenance::String = ""                  # the anchor's declaring scope and key
    admissible::Any = nothing                # gcd(pool), the coarsest admissible Δt_base
end

_dep_constraint(p::Symbol) =
    p === :method              ? "must be a stepper type — RK4 or Heun" :
    p === :firing_budget       ? "must be an integer ≥ 1" :
    p === :localization_tol    ? "must be a positive real" :
    p === :localization_budget ? "must be an integer ≥ 1" :
    p === :join_timeout        ? "must be a positive real — the shutdown tail's join cap " *
                                 "in seconds of wall clock" :
    p === :log                 ? "must be true or false — the retention switch" :
    p === :log_every           ? "must be an integer ≥ 1" :
    p === :log_max             ? "must be an integer ≥ 1, or Inf as the explicit opt-out" :
    p === :t_end               ? "must be a finite real ≥ 0 — the run's clock bound, " *
                                 "taken to the nearest frame top" :
    p === :h                   ? "must be positive" :
    p === :n                   ? "must be an integer ≥ 1" :
                                 "is outside its constraint"
_dep_section(p::Symbol) =
    p === :method              ? " (§10.2)" :
    p === :firing_budget       ? " (§10.6)" :
    (p === :localization_tol || p === :localization_budget) ? " (§10.4)" :
    p === :join_timeout        ? " (§12.4)" :
    (p === :log || p === :log_every || p === :log_max) ? " (§11.2)" :
    p === :t_end               ? " (§13.5)" :
    p === :n                   ? " (§9.1)" : ""
_dep_admissible(d) = d.admissible === nothing ? "" :
                     " — an admissible Δt_base divides gcd(pool) = $(d.admissible)"

function message(d::DeploymentInvalid)
    d.reason === :inexact &&
        return "`$(d.parameter)` must be exact — a Rational (`$(d.parameter) = 1//100`) or " *
               "a `Period`/`Hz` value: grid derivation is GCD arithmetic, ill-defined over " *
               "floats (§9.1, §10.5)"
    d.reason === :not_a_quantity &&
        return "`$(d.parameter)` must be a Rational or a `Period`/`Hz` value, got $(d.value)"
    d.reason === :missing &&
        return "deployment needs the continuous step: `Simulation(…; h = 1//100)` — a " *
               "domain rate is not a framework default (§9.1)"
    d.reason === :unanchored &&
        return "Δt_base cannot be derived: `$(join(d.paths, "`, `"))` is/are unanchored, " *
               "with period `m·Δt_base` — an anchor edit anywhere in the tree would " *
               "silently rescale it. Declare the base tick period instead: `Δt_base = …`, " *
               "or `n = …` (§9.1)"
    d.reason === :no_constraint &&
        return "Δt_base cannot be derived: no anchor declares a constraint to derive it " *
               "from (§9.1)"
    d.reason === :not_harmonic &&
        return "harmonic grid: Δt_base = $(d.value) is not an integer multiple of " *
               "h = $(d.related) (Δt_base = n·h, n ≥ 1, §10.5)"
    d.reason === :disagrees_with_n &&
        return "Δt_base = $(d.value) disagrees with n = $(d.related): Δt_base/h = " *
               "$(d.quotient) (§9.1)"
    d.reason === :anchor_period &&
        return "$(d.provenance): period $(d.value) is not an integer multiple of " *
               "Δt_base = $(d.related)$(_dep_admissible(d)) (§9.1)"
    d.reason === :anchor_offset &&
        return "$(d.provenance): offset $(d.value) does not land on the base grid at " *
               "Δt_base = $(d.related)$(_dep_admissible(d)) (§9.1)"
    "`$(d.parameter)` $(_dep_constraint(d.parameter)), got $(d.value)" *
    _dep_section(d.parameter)
end

"§11.3: a claim naming no root input face."
Base.@kwdef struct AttachUnknownFace <: Diagnostic
    binding::String                          # the binding type, as a string
    face::Symbol
    candidates::Vector{Symbol} = Symbol[]    # the root input-face list
end
message(d::AttachUnknownFace) =
    "$(d.binding) claims `$(d.face)`, which names no root input face — the root faces are " *
    "$(_faceset(d.candidates)) (§11.3)"

"§11.3: the same device instance offered to `attach!` twice."
Base.@kwdef struct AlreadyAttached <: Diagnostic
    device::String                           # the candidate's type, as a string
    incumbent::String                        # the existing roster entry's device id
    binding::String                          # the incumbent's binding type
end
message(d::AlreadyAttached) =
    "this $(d.device) instance is already rostered as $(d.incumbent) under $(d.binding) — " *
    "rebinding is spelled `detach!` then `attach!` (§11.3)"

"§11.1, §11.3: two devices claiming the calling task, a single-slot resource."
Base.@kwdef struct CallerTaskConflict <: Diagnostic
    device::String
    incumbent::String
end
message(d::CallerTaskConflict) =
    "$(d.device) declares `needs_calling_task`, and $(d.incumbent) already holds the " *
    "calling task (§11.1, §11.3)"

"§11.3: one root input face claimed by two devices."
Base.@kwdef struct ClaimConflict <: Diagnostic
    face::Symbol
    device::String
    incumbent::String
end
message(d::ClaimConflict) =
    "`$(d.face)` is claimed by both $(d.device) and $(d.incumbent) — one writer per root " *
    "input at any time (§11.3)"

"§11.3, §11.6: a greedy binding whose computed complement was empty — every face already claimed."
Base.@kwdef struct EmptyGreedyClaim <: Diagnostic
    device::String
    binding::String
end
severity(::EmptyGreedyClaim) = :warning
message(d::EmptyGreedyClaim) =
    "$(d.device) under $(d.binding) staked the empty remainder — every root input face is " *
    "already claimed (§11.6)"

"§11.6: a binding whose traits and whose methods disagree, in either direction."
Base.@kwdef struct BindingContractMismatch <: Diagnostic
    binding::String                          # the binding type, as a string
    reason::Symbol   # :claims_missing|:reads_missing|:greedy_without_input|:neither_side|
                     # :greedy_with_claims|:claims_without_input|:reads_without_output|
                     # :reads_not_namedtuple|:reads_not_selectors
    observed::Any = nothing                  # what `reads` returned, where that is the fact
end
function message(d::BindingContractMismatch)
    d.reason === :claims_missing &&
        return "$(d.binding) declares an input side and defines no `claims` enumeration — " *
               "the enumeration is the interface (§11.6)"
    d.reason === :reads_missing &&
        return "$(d.binding) declares an output side and defines no `reads` enumeration — " *
               "the enumeration is the interface (§11.6)"
    d.reason === :greedy_without_input &&
        return "$(d.binding) declares `is_greedy` without `is_input` — greediness is a " *
               "claim source within the input side, and a source without its side is " *
               "meaningless (§11.6)"
    d.reason === :neither_side &&
        return "$(d.binding) declares neither side — a binding that touches nothing is a " *
               "configuration mistake, not a degenerate (§11.6)"
    d.reason === :greedy_with_claims &&
        return "$(d.binding) declares `is_greedy` and defines its own `claims` — the two " *
               "claim sources are alternatives, not layers (§11.6)"
    d.reason === :claims_without_input &&
        return "$(d.binding) defines `claims` while `is_input` reads false — a method " *
               "written and never reached is exactly the drift this check catches (§11.6)"
    d.reason === :reads_without_output &&
        return "$(d.binding) defines `reads` while `is_output` reads false — a method " *
               "written and never reached is exactly the drift this check catches (§11.6)"
    d.reason === :reads_not_namedtuple &&
        return "$(d.binding)'s `reads` must return a NamedTuple of labeled selectors — " *
               "(; label = get_output(...), ...) — got $(d.observed) (§11.6, §14.4)"
    "$(d.binding)'s `reads` entries must be read selectors — get_output, get_input or " *
    "get_face (§14.4) — got $(d.observed)"
end

"§11.6: the device twin of `BindingContractMismatch` — no `loop` method, or `gather` against a binding with no output side."
Base.@kwdef struct DeviceContractMismatch <: Diagnostic
    device::String                           # the device type, or the roster id where that is what the site holds
    reason::Symbol                           # :no_loop | :no_output_side
end
message(d::DeviceContractMismatch) =
    d.reason === :no_loop ?
    "$(d.device) defines no `loop` method — the task body is the authoring " *
    "contract's one mandatory function (§11.6)" :
    "$(d.device)'s binding declares no output side — `gather` serves the compiled " *
    "`reads` enumeration (§11.6)"

"§11.2, §14.4: a binding read that does not resolve against the published snapshot."
Base.@kwdef struct ReadBindingUnresolved <: Diagnostic
    binding::String                          # the binding type, as a string
    selector::String                         # the selector as authored
    reason::Symbol   # :store_selector|:indexed|:unknown_cell|:unknown_root_input|
                     # :root_input_not_output|:unknown_output_face
    path::String = ""
    field::Union{Nothing,Symbol} = nothing
    candidates::Vector{Symbol} = Symbol[]
end
path(d::ReadBindingUnresolved) = d.path
function message(d::ReadBindingUnresolved)
    d.reason === :store_selector &&
        return "$(d.binding) reads $(d.selector), a *store* selector — the store selectors " *
               "resolve only against live stores, and a binding reads a published " *
               "snapshot, which deliberately carries none (§14.4, §11.2). The remedy is to " *
               "declare the field public and read the port published from it"
    d.reason === :indexed &&
        return "$(d.binding) reads $(d.selector) — a binding read is a whole cell, and " *
               "sub-cell index addressing is absent in this register (§14.4, README)"
    d.reason === :unknown_cell &&
        return "$(d.binding) reads $(d.selector), which names no cell — only declared " *
               "outputs, assembly faces and root inputs are addressable (§14.4)"
    d.reason === :unknown_root_input &&
        return "$(d.binding) reads $(d.selector), which names no root input face — the " *
               "root inputs are $(_faceset(d.candidates)) (§14.4)"
    d.reason === :root_input_not_output &&
        return "$(d.binding) reads $(d.selector), which names a root *input* face — the " *
               "integration register is the exported output faces, and a root input is " *
               "read back with get_input (§14.4, §11.2)"
    "$(d.binding) reads $(d.selector): `$(d.field)` is no root-exported output face " *
    "(§14.4, §11.2)"
end

"§14.2, §14.3: a condition leaf that does not resolve against this build."
Base.@kwdef struct ConditionResolution <: Diagnostic
    path::String = ""
    store::Symbol = :input                   # :x | :s | :m | :input
    field::Symbol
    reason::Symbol   # :assembly_path|:unknown_path|:unexported_face|:no_input_face|
                     # :internally_wired|:no_store|:undeclared_field|:unconvertible
    face::Union{Nothing,Symbol} = nothing    # the root input the chain lands on
    provenance::String = ""                  # the tree's chain to this value
    candidates::Vector{Symbol} = Symbol[]
    declared::Any = nothing                  # the declared leaf type
    observed::Any = nothing                  # the authored value's type
    value::Any = nothing
    tier::Union{Nothing,Symbol} = nothing
    role::Union{Nothing,Symbol} = nothing    # :output_port | :input_face | :workspace
    producer::Union{Nothing,Tuple{String,Symbol}} = nothing
    activation::Any = nothing                # set when this is the seeded activation's refusal
end
path(d::ConditionResolution) = d.path

_cleaf(d) = d.face !== nothing ? "root input `$(d.face)`" :
            d.store === :input ? "input face `$(d.field)` of $(_at_path(d.path))" :
                                 "`$(d.store).$(d.field)` at $(_at_path(d.path))"
_ctail(d) = " (§14.3) [$(d.provenance)]"
_role_word(r) = r === :output_port ? "an output port" :
                r === :input_face  ? "an input face" : "a workspace entry"

function message(d::ConditionResolution)
    d.reason === :assembly_path &&
        return "the condition addresses $(_at_path(d.path)), which is an assembly — " *
               "assemblies own no state, and a condition addresses components and root " *
               "inputs (§14.1, §8.5) [$(d.provenance)]"
    d.reason === :unknown_path &&
        return "the condition addresses $(_at_path(d.path)), which is no component of " *
               "this build" * _ctail(d)
    d.reason === :unexported_face &&
        return "`$(d.field)` is no root input face — the root's inputs are " *
               "$(_namelist(d.candidates)) (§14.2) [$(d.provenance)]"
    d.reason === :no_input_face &&
        return "$(_at_path(d.path)) declares no input face `$(d.field)` — its input faces " *
               "are $(_namelist(d.candidates)) (§14.2) [$(d.provenance)]"
    d.reason === :internally_wired &&
        return "$(_at_path(d.path))'s input face `$(d.field)` reaches no root input — it " *
               "is wired internally, to `$(first(d.producer))`.$(last(d.producer)), and " *
               "the first sweep overwrites it; unexported stays unpokeable (§14.2) " *
               "[$(d.provenance)]"
    d.reason === :no_store &&
        return "$(_cleaf(d)) — $(_at_path(d.path)) is a $(d.tier) component and declares " *
               "no `init_$(d.store)`" *
               (d.store === :x && d.tier === :discrete ?
                "; the discrete tier's state is `s` (D-195)" :
                d.store === :s && d.tier === :continuous ?
                "; the continuous tier's state is `x` (D-195)" :
                d.store === :m ?
                "; modes are declared by `init_m`, continuous-only (§3.2)" : "") * _ctail(d)
    d.reason === :undeclared_field &&
        return "$(_cleaf(d)) is not declared — `init_$(d.store)` at $(_at_path(d.path)) " *
               "declares $(_namelist(d.candidates))" *
               (d.role === nothing ? "" :
                "; `$(d.field)` is $(_role_word(d.role)), and a condition specifies state, " *
                "modes and root inputs — never outputs, never workspace (§14.1)") * _ctail(d)
    "$(_cleaf(d)) takes $(d.declared), and the authored value is " *
    "$(repr(d.value))::$(d.observed), which does not convert" *
    (d.activation === nothing ? "" :
     "; this is the seeded activation's own refusal — a value at $(d.activation) is a " *
     "decision variable and this leaf is pinned, and a decision variable descends into " *
     "neither a frozen discrete `s` nor a pinned leaf (§14.3, §9.4)") * _ctail(d)
end

"§14.2: one leaf written by two fragments of a `combine` — collision-intolerant by design."
Base.@kwdef struct DuplicateConditionLeaf <: Diagnostic
    path::String = ""
    store::Symbol = :input
    field::Symbol
    face::Union{Nothing,Symbol} = nothing
    provenance::Vector{String} = String[]    # both chains
end
path(d::DuplicateConditionLeaf) = d.path
message(d::DuplicateConditionLeaf) =
    "$(_cleaf(d)) is written twice — by $(_prov(d, 1)), and by $(_prov(d, 2)). `combine` " *
    "is collision-intolerant by design — use `override(base, patch)` to layer (§14.2, §14.6)"

"§14.2: a value handed to a condition combinator that is not a condition node."
Base.@kwdef struct ConditionNodeMisuse <: Diagnostic
    observed::Any                            # the offending argument's type
    reason::Symbol = :not_a_node             # :not_a_node | :fragment_payload
    payload::Union{Nothing,Symbol} = nothing # which `fragment` payload, for :fragment_payload
    in_hand::Vector{Symbol} = Symbol[]       # the node kinds in hand
end
message(d::ConditionNodeMisuse) =
    d.reason === :fragment_payload ?
    "`fragment`'s `$(d.payload)` payload is $(d.observed) — every payload is a NamedTuple " *
    "of the authoring level's own names (§14.2)" :
    "$(d.observed) is not a condition node" *
    (isempty(d.in_hand) ? "" :
     ", and the node kind(s) in hand are $(_plainlist(d.in_hand))") *
    " — wrap the NamedTuple in `fragment(…)`, or in `at(prefix, fragment(…))`, and " *
    "combine nodes with nodes (§14.2)"

"§14.6: an application whose condition leaves a root input with no value — nothing is written."
Base.@kwdef struct UninitializedInputs <: Diagnostic
    op::Symbol
    faces::Vector{Symbol}                    # every uncovered root face, in declaration order
end
message(d::UninitializedInputs) =
    "the condition given to `$(d.op)` covers no value for root input(s) " *
    "$(_namelist(d.faces)) — root inputs are the one initialized datum with no declared " *
    "default (§11.3), so an application establishing a complete world authors every one of " *
    "them; nothing was written (§14.6, D-068)"

"§14.10: a read selector that does not resolve against this build."
Base.@kwdef struct TapResolution <: Diagnostic
    label::Symbol                            # the read's label in the set
    selector::String                         # the selector as authored
    reason::Symbol   # :assembly_path|:unknown_path|:scalar_index|:undeclared|:discrete_deriv|
                     # :unknown_root_input|:root_input_not_face|:unknown_output_face
    tap::Union{Nothing,Symbol} = nothing     # :x | :u | :y
    path::String = ""
    field::Union{Nothing,Symbol} = nothing
    index::Union{Nothing,Int} = nothing
    declares::Union{Nothing,Symbol} = nothing  # :state_field | :output_port
    declared::Any = nothing
    candidates::Vector{Symbol} = Symbol[]
end
path(d::TapResolution) = d.path

_tap_noun(s) = s === :output_port ? "output port" : "state field"

# The tap set and the index come off the selector's own kind, so every arm has
# them and the shared prefix shows them: which of `x`/`u`/`y` the read addresses
# is §14.10's payload, and the index is the coordinate the author wrote.
_tapviol(d, what) =
    "the read labeled `$(d.label)` is $(d.selector)" *
    (d.tap === nothing ? "" :
     " (tap `$(d.tap)`" * (d.index === nothing ? "" : ", index $(d.index)") * ")") *
    ", and $what (§14.4)"

function message(d::TapResolution)
    d.reason === :assembly_path &&
        return _tapviol(d, "$(_at_path(d.path)) is an assembly — a path selector addresses " *
                           "a component's own declarations, and a root-exported face is " *
                           "read with `get_face`")
    d.reason === :unknown_path &&
        return _tapviol(d, "$(_at_path(d.path)) is no component of this build")
    d.reason === :scalar_index &&
        return _tapviol(d, "the leaf it names is declared $(d.declared) — a scalar has no " *
                           "index, and `i` addresses a component of a vector leaf")
    d.reason === :discrete_deriv &&
        return _tapviol(d, "$(_at_path(d.path)) is a discrete component — a discrete `s` " *
                           "has no derivative, and `ẋ` exists on the continuous tier alone " *
                           "(§7.1, D-195)")
    d.reason === :undeclared &&
        return _tapviol(d, "$(_at_path(d.path)) declares no $(_tap_noun(d.declares)) " *
                           "`$(d.field)` — its $(_tap_noun(d.declares)) names are " *
                           "$(_namelist(d.candidates))")
    d.reason === :unknown_root_input &&
        return _tapviol(d, "`$(d.field)` is no root input face — the root's inputs are " *
                           "$(_namelist(d.candidates))")
    d.reason === :root_input_not_face &&
        return _tapviol(d, "`$(d.field)` is a root *input* face — the integration register " *
                           "is the root-exported output faces, and a root input is read " *
                           "back with `get_input`")
    _tapviol(d, "`$(d.field)` is no root-exported output face — the root exports " *
                "$(_namelist(d.candidates))")
end

"§14.7, §14.8: a `TrimProblem` field that does not meet the problem's closed shape."
Base.@kwdef struct TrimProblemInvalid <: Diagnostic
    field::Symbol
    reason::Symbol   # :not_a_namedtuple|:key_set|:field_types|:inverted_box|
                     # :nonpositive_tolerance|:not_a_read_set
    observed::Any = nothing
    names::Vector{Symbol} = Symbol[]         # the names in hand
    expected::Vector{Symbol} = Symbol[]      # the names it has to match
    key::Union{Nothing,Symbol} = nothing     # the offending decision or residual
    value::Any = nothing
    bound::Any = nothing
    bad::Vector{Pair{Symbol,Any}} = Pair{Symbol,Any}[]   # field => observed type
end

_trim_shape(f::Symbol) =
    f === :tolerances ? "the per-residual convergence test is an all-`Float64` NamedTuple" :
    f === :residuals  ? "the residual system is a NamedTuple of named equations, " *
                        "same-named as `tolerances`" :
                        "the decisions and their two bounds are same-named all-`Float64` " *
                        "NamedTuples"
_trim_floats(f::Symbol) =
    f === :tolerances ? "a tolerance is a `Float64` in its residual's own physical units" :
                        "decisions and bounds are `Float64`"
_trim_verb(f::Symbol) = f === :residuals ? "returned" : "is"
_trim_bad(d) = join(("`$k`::$v" for (k, v) in d.bad), ", ")

function message(d::TrimProblemInvalid)
    d.reason === :not_a_read_set &&
        return "`reads` is $(d.observed) — the declared read set is a `reads(…)` value: " *
               "reads(name = get_deriv(\"path\", :field), …) (§14.4)"
    d.reason === :not_a_namedtuple &&
        return "`$(d.field)` $(_trim_verb(d.field)) $(d.observed) — $(_trim_shape(d.field)) " *
               "(§14.7)"
    d.reason === :key_set &&
        return d.field === :residuals ?
               "`residuals` returns $(_namelist(d.names)) and `tolerances` names " *
               "$(_namelist(d.expected)) — the two share one key set, and order is never a " *
               "mismatch (§14.7)" :
               "`$(d.field)` names $(_namelist(d.names)) and `guess` names " *
               "$(_namelist(d.expected)) — the three share one key set, a permuted " *
               "spelling pairing by name (§14.7)"
    d.reason === :field_types &&
        return d.field === :residuals ?
               "`residuals` field(s) $(_trim_bad(d)) — each residual is a real scalar (§14.7)" :
               "`$(d.field)` field(s) $(_trim_bad(d)) — $(_trim_floats(d.field)) (§14.7)"
    d.reason === :inverted_box &&
        return "`lower` names `$(d.key)` = $(d.value) above `upper`'s $(d.bound) — a " *
               "decision's box is `lower ≤ upper`, and an inverted pair admits no point at " *
               "all (§14.7)"
    "`tolerances` names `$(d.key)` = $(d.value) — a tolerance is finite and strictly " *
    "positive: it is the half-width of the box its residual has to sit in, and the " *
    "normalized acceptance test divides by it (§14.7)"
end

"§14.8: boundary zero fired events at the commit, moving the committed stores off the solved point."
Base.@kwdef struct TrimCommitEvents <: Diagnostic
    events::Vector{Tuple{String,Symbol}}     # component path and event name
end
severity(::TrimCommitEvents) = :warning
message(d::TrimCommitEvents) =
    "boundary zero fired " *
    join(("`$p`.$n" for (p, n) in d.events), ", ") * " at the commit — a handler that " *
    "fires there moves the committed stores off the solved point, and the committed-state " *
    "residuals are where they actually sit (§14.5, §14.8)"

"§14.8: a converged solve whose committed-state residuals leave the box."
Base.@kwdef struct TrimCommitResiduals <: Diagnostic
    residuals::Vector{Tuple{Symbol,Float64,Float64}}   # name, committed value, tolerance
end
severity(::TrimCommitResiduals) = :warning
message(d::TrimCommitResiduals) =
    "this solve converged, and the residuals re-gathered after the commit leave the box: " *
    join(("`$k` = $v against $t" for (k, v, t) in d.residuals), ", ") * " — the mover is " *
    "boundary zero's `project` or a commit-fired handler, and the verdict is not " *
    "re-litigated: it gated the commit, at the solved point (§14.5, §14.8)"

"§14.4: a condition tree whose shape differs from the one its plan was compiled from."
Base.@kwdef struct ConditionShapeDrift <: Diagnostic
    reason::Symbol                           # :tree_type | :prefix
    compiled::Any = nothing                  # the compiled tree type, or the compiled prefix
    observed::Any = nothing
    position::Union{Nothing,Tuple} = nothing # the node position, for :prefix
end

_drift_position(P::Tuple) = "(" * join((s isa Symbol ? ".$s" : "[$s]" for s in P), "") * ")"

message(d::ConditionShapeDrift) =
    d.reason === :prefix ?
    "the `at` prefix at tree position $(_drift_position(d.position)) was " *
    "$(repr(d.compiled)) when this plan was compiled and is $(repr(d.observed)) now — " *
    "prefixes are runtime `String` fields the tree type cannot carry, so the register " *
    "closes the shape with a `===` sweep over them, and a condition function has to return " *
    "one shape for every decision it is evaluated at (§14.4, §9.5, D-066)" :
    "this plan was compiled from a condition tree of type\n    $(d.compiled)\nand the tree " *
    "handed to `apply!` is\n    $(d.observed)\nThe specialized register proves the shape " *
    "by dispatch, so a condition function has to return one shape for every decision it is " *
    "evaluated at — a branch that authors a different field set, a different nesting or a " *
    "different leaf type is a different shape, and needs its own plan (§14.4, §9.5, D-066)"

"§8.7, §11.6, §12.6, §14.7, D-215: an argument outside its constraint — `DeploymentInvalid`'s twin off the deployment surface."
Base.@kwdef struct ArgumentInvalid <: Diagnostic
    call::Symbol                             # :Period|:Hz|:Absolute|:step!|:run!|:trim!|:TableBinding|:selector
    reason::Symbol
    argument::Union{Nothing,Symbol} = nothing
    value::Any = nothing
    entry::Union{Nothing,Symbol} = nothing   # the `TableBinding` entry at fault
    vocabulary::Vector{Symbol} = Symbol[]    # the legal entry keys
end

function message(d::ArgumentInvalid)
    d.reason === :inexact &&
        return d.call === :Period ?
               "a period is an exact Rational — write `Period(1//50)`, not $(repr(d.value)): " *
               "grid derivation is GCD arithmetic (§10.5)" :
               d.call === :Hz ?
               "a frequency is an exact Rational — write `Hz(1//2)` for 0.5 Hz, not " *
               "$(repr(d.value)): grid derivation is GCD arithmetic (§10.5)" :
               "an offset is an exact Rational — write `Absolute(Hz(50), 1//500)`, not " *
               "$(repr(d.value)): grid derivation is GCD arithmetic (§10.5)"
    d.reason === :not_a_quantity &&
        return "`Absolute` takes a quantity value: `Period(1//50)` or `Hz(50)` — got " *
               "$(repr(d.value)) (§10.5)"
    d.reason === :both_given &&
        return "step! takes `frames` or `t_plus`, not both — the count and the duration are " *
               "two spellings of one advance (§12.6)"
    d.reason === :no_clock_bound &&
        return "run! has no clock bound: `t_end` was given neither at the constructor nor " *
               "here — the constructor value is the default and the run! keyword the " *
               "per-run override (§13.5)"
    d.reason === :non_nominal &&
        return "`trim!` needs a nominal `Simulation{Float64}` and this one is $(d.value) — " *
               "trim commits through the nominal world, and the seeded activation it " *
               "iterates on is the service's own scratch, instantiated per invocation " *
               "(§14.8, §9.4)"
    d.reason === :not_a_problem &&
        return "`trim!` takes a `TrimProblem` and was given $(d.value) — the problem is " *
               "one value with a closed field set: TrimProblem(; guess, lower, upper, " *
               "condition, reads, residuals, tolerances) (§14.7)"
    d.reason === :index_not_integer &&
        return "a selector's index must be an integer — the component index of §14.10, " *
               "applied to the read value — got $(repr(d.value)) (§14.4)"
    d.reason === :entry_shape &&
        return "TableBinding: entry `$(d.entry)` must be a NamedTuple — " *
               "(face = ..., deadzone = ..., expo = ...) (§11.6)"
    d.reason === :no_face &&
        return "TableBinding: entry `$(d.entry)` names no `face` — the face is what the " *
               "channel writes (§11.6)"
    d.reason === :face_name &&
        return "TableBinding: entry `$(d.entry)`'s face must be a face name, got " *
               "$(repr(d.value)) (§11.6)"
    d.reason === :vocabulary &&
        return "TableBinding: entry `$(d.entry)` carries `$(d.argument)` — the entry " *
               "vocabulary is $(_namelist(d.vocabulary)) (§11.6)"
    d.reason === :deadzone &&
        return "TableBinding: entry `$(d.entry)`'s deadzone must lie in [0, 1), got " *
               "$(d.value) (§11.6)"
    d.reason === :expo &&
        return "TableBinding: entry `$(d.entry)`'s expo must lie in [0, 1], got " *
               "$(d.value) (§11.6)"
    d.argument === :frames ?
        "frames must be an integer ≥ 1, got $(d.value) (§12.6)" :
        "t_plus must be a finite real > 0 — the duration spelling — got $(d.value) (§12.6)"
end

"§14.4, D-215: a non-selector in a read set, or a non-`Reads` where one is expected — `ConditionNodeMisuse`'s twin."
Base.@kwdef struct ReadSetMisuse <: Diagnostic
    observed::Any                            # the offending argument's type
    reason::Symbol = :not_a_read_set         # :not_a_selector | :not_a_read_set
    label::Union{Nothing,Symbol} = nothing
    in_hand::Vector{Symbol} = Symbol[]       # the selector kinds in hand
end
message(d::ReadSetMisuse) =
    d.reason === :not_a_selector ?
    "the read labeled `$(d.label)` is $(d.observed)" *
    (isempty(d.in_hand) ? "" :
     ", and the selector kind(s) in hand are $(_plainlist(d.in_hand))") *
    " — a read set is labeled *selectors*, " *
    "get_state / get_deriv / get_output / get_input / get_face (§14.4)" :
    "$(d.observed) is not a read set — wrap the labeled selectors in `reads(…)`: " *
    "reads(name = get_deriv(\"path\", :field), …) (§14.4, §14.7)"

"§11.3, D-215: `detach!` offered a device the roster does not hold — `AlreadyAttached`'s mirror."
Base.@kwdef struct NotAttached <: Diagnostic
    device::String
    roster::Vector{String} = String[]        # the roster's device ids
end
message(d::NotAttached) =
    "this $(d.device) instance is not rostered — `detach!` releases an existing " *
    "attachment, and the roster holds " *
    (isempty(d.roster) ? "no device" : join(("`$(r)`" for r in d.roster), ", ")) * " (§11.3)"
