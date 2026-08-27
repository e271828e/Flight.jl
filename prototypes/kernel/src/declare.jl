# The declaration layer (§8.1, §8.2): a plain-Julia trait layer. Declarations
# are ordinary functions of the *instance*, written at concrete `Float64`; the
# framework's activation walk (§7.2) retypes them. No macros, no stage tags.
#
# Both tiers live here, and so does the event declaration surface (§2.1, §8.2).

# --- the one root type (§8.5) -------------------------------------------------
# There is no `AbstractAssembly`: primitives and assemblies share one root, so a
# field declared `E <: AbstractEngine` accepts either, and class is read off
# declaration shape (`src/assembly.jl`).

abstract type AbstractComponent end

"""
The build's own failure: a plain message naming the path or the entry at fault.
§13.2's structured carrier is deliberately absent.
"""
struct BuildError <: Exception
    msg::String
end
Base.showerror(io::IO, e::BuildError) = print(io, "BuildError: ", e.msg)

# --- what an author declares --------------------------------------------------

"""
Continuous state, **by value**, at nominal `Float64`; leaves drawn from §7.1's
closed vocabulary. One argument, and the criterion is D-166's: a by-value
declaration states nominal physics, and its *types* walk by rule — §7.1 admits
no pinned state leaf, so there is no choice for a `T` to record.
"""
init_x(::Any) = NamedTuple()

"""
Discrete state, **by value**, and the discrete tier's own letter (D-195): any
isbits type, pinned wholesale — nothing here walks with the activation, which
is why the declaration takes no `T`. Disjoint from `init_x` by construction, so
a state declaration always carries its tier (§8.2).
"""
init_s(::Any) = NamedTuple()

"""
Modes, **by value**, continuous-only: the event system is continuous-side only,
so declaring `init_m` announces the continuous tier (§8.2). Mode stores are
written by handlers; nothing else may.
"""
init_m(::Any) = NamedTuple()

"""
Mutable scratch, instantiated by the framework and arriving as the bundle's
`ws` field (§7.3). Declaration *is* allocation, and the arity splits the tiers:
`workspace(::C, ::Type{T})` continuous — sizes from the instance, eltypes from
the activation — against plain `workspace(::C)` on the discrete tier. No
fallback: absence is how a component declares no scratch.
"""
function workspace end

"""
Input faces: name => type, **as a function of the activation scalar** (D-167)
on the continuous tier; plain `input_types(::C)` on the discrete tier, where
the declared types are pinned. Entries are read permissively — they state per
leaf what the consumer allows: `T` is tolerant, a literal `Float64` demands a
frozen arrival.
"""
input_types(::Any, ::Type{T}) where {T <: Real} = NamedTuple()

"""
Output ports: name => type, **as a function of the activation scalar** (D-166).
Semantics are literal — a cell's type at an activation *is* this declaration
evaluated at that activation's `T`, with no leaf walk behind it: `T` (alone or
as a parameter) participates, a literal `Float64` is a deliberately pinned
leaf.

One declaration for both stages — there are no stage tags anywhere (§8.2);
which stage produces a port is *discovered* by the build probe (§9.3), and the
declaration is what the probe checks against.
"""
output_types(::Any, ::Type{T}) where {T <: Real} = NamedTuple()

# --- what an assembly declares (§8.6) -----------------------------------------
# All three are ordered collections of `Pair`s of strings, and every arrow reads
# along the signal flow. Direction is declared by the method, never inferred: the
# resolved endpoints only cross-check it.

"""
Wires: `"src/path/port" => "dst/path/port"`, strictly child endpoint => child
endpoint, relative to the declaring assembly.

**Mandatory even when empty**, because defining it *is* the assembly class
marker (§8.5) — which is why it has no fallback to match.
"""
function child_connections end

"""
The boundary, inward: face name => internal endpoint path, or a tuple of paths
for an input face fanning out through the boundary. Absence declares no input
face.
"""
input_connections(::Any) = ()

"""
The boundary, outward: internal source path => face name, so that its pairs, like
every other pair in the three declarations, read along the flow.
"""
output_connections(::Any) = ()

"""
The one optional declaration (§8.5, D-211): a component may declare **at most
one** of its container fields name-transparent by returning that field's
`Symbol`. Its elements are then contributed as children under their bare keys —
`"key"` and `"1"` in place of `"field/key"` and `"field/1"` — everywhere a child
name appears: wiring endpoints, `sample_times` keys, read paths, `at` prefixes,
diagnostics. Naming is the only thing it changes; the container keeps its
transparency of contract. "At most one" is unrepresentable otherwise here: the
return is one `Symbol` or `nothing`, and `nothing` — undeclared containers
keeping their key segment — is the default.
"""
transparent_container(::Any) = nothing

# --- rate scopes (§8.7, §10.5) --------------------------------------------------
# The two registers of the sample-time declaration, and the wrappers are the
# whole value vocabulary — a bare integer or bare quantity in a `sample_times`
# entry is a declaration error. They are plain data carriers (D-185): range
# validation is Stratum A's, with path attribution, in the fold. The one
# constructor-side refusal is floats, because it is a normalization failure
# rather than a range check: grid derivation is GCD arithmetic, ill-defined over
# floats, so every period and offset is an exact `Rational{Int}` and a float
# argument throws the teaching error naming the exact spelling.

"""
A quantity value for `Absolute` and for the `Δt_base` deployment keyword:
`Period(1//50)` or `Hz(50)` — a spelling choice, both normalized to the exact
rational period at construction (§10.5).
"""
struct Period
    T::Rational{Int}
    Period(T::Union{Integer,Rational{<:Integer}}) = new(Rational{Int}(T))
    Period(::AbstractFloat) =
        throw(BuildError("a period is an exact Rational — write `Period(1//50)`, not a " *
                         "float: grid derivation is GCD arithmetic (§10.5)"))
end

Hz(f::Union{Integer,Rational{<:Integer}}) = Period(1 // f)
Hz(::AbstractFloat) =
    throw(BuildError("a frequency is an exact Rational — write `Hz(1//2)` for 0.5 Hz, not " *
                     "a float: grid derivation is GCD arithmetic (§10.5)"))

period(q::Period) = q.T

"""
`Relative(K, φ = 0)`: every `K`-th tick of the enclosing scope, starting from
its `φ`-th — scope ticks as the unit system, composing affinely down the tree
(§10.5). `K = 1` admits no stagger; same-rate siblings stagger one level down.
"""
struct Relative
    K::Int
    φ::Int
end
Relative(K::Integer) = Relative(K, 0)

"""
`Absolute(q, τ = 0)`: tick instants `t = τ + k·period(q)`, in exact rational
seconds — an anchor, severing its child from the enclosing scope's grid; the
final divisor does not exist until deployment binds `Δt_base` (§10.5, §9.1).
"""
struct Absolute
    T::Rational{Int}
    τ::Rational{Int}
    Absolute(q::Period, τ::Union{Integer,Rational{<:Integer}} = 0) = new(q.T, Rational{Int}(τ))
    Absolute(::Period, ::AbstractFloat) =
        throw(BuildError("an offset is an exact Rational — write `Absolute(Hz(50), 1//500)`, " *
                         "not a float: grid derivation is GCD arithmetic (§10.5)"))
    Absolute(::Real, τ...) =
        throw(BuildError("`Absolute` takes a quantity value: `Period(1//50)` or `Hz(50)` (§10.5)"))
end

"""
Rate scope (§8.7): immediate child name => `Relative` or `Absolute`. Optional,
and so is any key — an unlisted discrete child defaults to `Relative(1)`, so
only multiplied, phased or anchored children appear. The declaration belongs to
the assembly type; a deployment-preference multiplier is exposed as a
constructor parameter and read off the instance here (§10.5).
"""
sample_times(::Any) = (;)

# --- events (§2.1, §8.2) ------------------------------------------------------

"""
One guard/handler pair, and nothing else: `Event(guard, handler)` carries no
detection keyword. Detection policy is declared by the guard's *return type* —
a `Bool` guard is boundary-detected, a sign-form guard is localized — read off
the probe it already runs (§10.4, D-179), so the illegal form/policy pairing is
unrepresentable rather than merely diagnosed.
"""
struct Event{G,H}
    guard::G
    handler::H
end

"""
The ordered, named guard/handler collection (§8.2): `name = Event(guard,
handler)` entries. Order is semantics — declaration order is priority with
re-decision at the boundary iteration (§10.6). Continuous-only, like `init_m`:
the event system is continuous-side only (§5.2). Nothing here is inferrable.
"""
events(::Any) = (;)

"""
Manifold projection (§5.2): one store in, the same store out, which is why it
alone stays positional. It runs between a state write and its decode — after
integration and after each handler firing (§5.3) — and its return is written
back to the buffer wholesale, so the probe holds it complete against the state
shape (§9.3).
"""
function project end

"""The §2.1 predicate of a guard's return: the `Bool` form itself, or `σ ≥ 0`."""
_holding(σ::Bool) = σ
_holding(σ) = σ ≥ 0

# --- what an author defines (the §5.2 signatures) -----------------------------
# Every stage takes the component and exactly one NamedTuple bundle of views,
# destructured by name: `h_xu(c::Plant, (; x, u)) = ...`. The framework's call
# is one fixed shape; the bundle law (below) decides what the tuple carries.
#
# The two output stages come in one pair per tier and the pairs are disjoint
# (D-195): `h_x`/`h_xu` continuous, `h_s`/`h_su` discrete. "No `u` in the name"
# is the no-feedthrough marker within either pair, and a stage name on its own
# now carries the tier, which is what makes a wrong-letter declaration
# diagnosable (§8.2).
#
# A component defines the stages it needs. `has_stage` is how the build asks,
# and it is a question about method existence, not a declaration the author
# repeats — the definition site is the single source of truth.

function h_x end
function h_xu end
function h_s end
function h_su end
function f end
function g end

has_stage(fn, c) = hasmethod(fn, Tuple{typeof(c),NamedTuple})

"""
Is `fn` declared *for this component* in the arity `extra` names, as against
matching a framework fallback? `hasmethod` cannot tell the two apart, so the
matched method's own signature is what answers: a fallback carries `Any` in the
component position.
"""
_declares(fn, c, extra...) =
    hasmethod(fn, Tuple{typeof(c),extra...}) &&
    Base.unwrap_unionall(which(fn, Tuple{typeof(c),extra...}).sig).parameters[2] !== Any

# --- tiers (§8.2) -------------------------------------------------------------
# A tier is never announced; it is read off the declaration shape. The enum is
# the build's internal answer, not an authoring surface.

@enum Tier CONTINUOUS DISCRETE

tier_word(t::Tier) = t === CONTINUOUS ? "continuous" : "discrete"

"""The tier form of `fn`'s arity: two-argument continuous, plain discrete."""
declared_at(fn, c, t::Tier) =
    t === CONTINUOUS ? (_declares(fn, c, Type{Float64}) ? fn(c, Float64) : NamedTuple()) :
                       (_declares(fn, c) ? fn(c) : NamedTuple())

# The tier's own name family (D-195). Everything downstream asks for a stage or
# an update law *through* the tier, so no code path ever holds a name that
# serves both.
stage1_of(t::Tier) = t === CONTINUOUS ? h_x : h_s
stage2_of(t::Tier) = t === CONTINUOUS ? h_xu : h_su
update_of(t::Tier) = t === CONTINUOUS ? f : g

# --- the bundle law (§5.2) ----------------------------------------------------
# A name appears in a component's bundle iff the corresponding store or fact
# exists for that component. Undeclared stores are *absent*, never
# `nothing`-filled, so a body destructuring what it does not own fails at the
# destructuring rather than reading a filler.

"""
Bundle field names for `fn` on component `c` at tier `t`, given its discovered
stage-1 ports.

The by-type declarations are asked at nominal `Float64`: presence is a property
of the declaration's key set, which is fixed by the declaration and cannot vary
with the activation scalar, so the bundle shape is the same at every activation.

The per-function, per-tier name sets are closed, and so are the state letters:
`x`/`y_x` with `m` on the continuous tier, `s`/`y_s` with `Δt` on the discrete
(D-195).
"""
function bundle_names(fn, c, t::Tier, stage1_ports::Tuple)
    stage2, update = stage2_of(t), update_of(t)
    names = Symbol[]
    if t === CONTINUOUS
        !isempty(init_x(c)) && push!(names, :x)
        !isempty(init_m(c)) && push!(names, :m)
    else
        !isempty(init_s(c)) && push!(names, :s)
    end
    if fn === stage2 || fn === update
        !isempty(declared_at(input_types, c, t)) && push!(names, :u)
    end
    if fn === stage2
        !isempty(stage1_ports) && push!(names, t === CONTINUOUS ? :y_x : :y_s)
    elseif fn === update
        !isempty(declared_at(output_types, c, t)) && push!(names, :y)
    end
    _declares_workspace(c, t) && push!(names, :ws)
    push!(names, :t)
    t === DISCRETE && push!(names, :Δt)
    tuple(names...)
end

_declares_workspace(c, t::Tier) =
    t === CONTINUOUS ? _declares(workspace, c, Type{Float64}) : _declares(workspace, c)

"""
Bundle field names for a guard or handler (§5.2): the update law's view of the
world — `x, m, y, u, t [, ws]` — one closed set shared by both halves, on the
continuous tier only. The same iff rule as `bundle_names`, without a stage-1
distinction: guards and handlers run against the complete fresh table.
"""
function event_bundle_names(c)
    names = Symbol[]
    !isempty(init_x(c)) && push!(names, :x)
    !isempty(init_m(c)) && push!(names, :m)
    !isempty(declared_at(input_types, c, CONTINUOUS)) && push!(names, :u)
    !isempty(declared_at(output_types, c, CONTINUOUS)) && push!(names, :y)
    _declares_workspace(c, CONTINUOUS) && push!(names, :ws)
    push!(names, :t)
    tuple(names...)
end

# --- probe values (§9.3) -----------------------------------------------------
# Root inputs are the one terminal with no producer; the build synthesizes their
# values here. Overridable per type — that is the seam a constrained type uses
# to state a valid default.

probe_value(::Type{T}) where {T<:Real} = zero(T)
probe_value(::Type{Bool}) = false
probe_value(::Type{P}) where {P<:StaticArray} = zero(P)
probe_value(::Type{P}) where {P} = P()
