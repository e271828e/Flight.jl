# The declaration layer (§8.1, §8.2): a plain-Julia trait layer. Declarations
# are ordinary functions of the *instance*, written at concrete `Float64`; the
# framework's activation walk (§7.2) retypes them. No macros, no stage tags.
#
# Increment 2 covers the continuous tier only. Modes (`init_m`), `workspace`
# and `events` are declaration entries this file deliberately does not have
# yet.

# --- what an author declares --------------------------------------------------

"""
Continuous state, **by value**, at nominal `Float64`; leaves drawn from §7.1's
closed vocabulary. One argument, and the criterion is D-166's: a by-value
declaration states nominal physics, and its *types* walk by rule — §7.1 admits
no pinned state leaf, so there is no choice for a `T` to record.
"""
init_x(::Any) = NamedTuple()

"""
Input faces: name => type, **as a function of the activation scalar** (D-167).
Entries are read permissively — they state per leaf what the consumer allows:
`T` is tolerant, a literal `Float64` demands a frozen arrival.
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

# --- what an author defines (the §5.2 signatures) -----------------------------
# Every stage takes the component and exactly one NamedTuple bundle of views,
# destructured by name: `h_xu(c::Plant, (; x, u)) = ...`. The framework's call
# is one fixed shape; the bundle law (below) decides what the tuple carries.
#
# A component defines the stages it needs. `has_stage` is how the build asks,
# and it is a question about method existence, not a declaration the author
# repeats — the definition site is the single source of truth.

function h_x end
function h_xu end
function f end

has_stage(fn, c) = hasmethod(fn, Tuple{typeof(c),NamedTuple})

# --- the bundle law (§5.2) ----------------------------------------------------
# A name appears in a component's bundle iff the corresponding store or fact
# exists for that component. Undeclared stores are *absent*, never
# `nothing`-filled, so a body destructuring what it does not own fails at the
# destructuring rather than reading a filler.

"""
Bundle field names for `fn` on component `c`, given its discovered stage-1
ports.

The by-type declarations are asked at nominal `Float64`: presence is a property
of the declaration's key set, which is fixed by the declaration and cannot vary
with the activation scalar, so the bundle shape is the same at every activation.
"""
function bundle_names(fn, c, stage1_ports::Tuple)
    names = Symbol[]
    !isempty(init_x(c)) && push!(names, :x)
    if fn === h_xu || fn === f
        !isempty(input_types(c, Float64)) && push!(names, :u)
    end
    if fn === h_xu
        !isempty(stage1_ports) && push!(names, :y_x)
    elseif fn === f
        !isempty(output_types(c, Float64)) && push!(names, :y)
    end
    push!(names, :t)
    tuple(names...)
end

# --- probe values (§9.3) -----------------------------------------------------
# Root slots are the one terminal with no producer; the build synthesizes their
# values here. Overridable per type — that is the seam a constrained type uses
# to state a valid default.

probe_value(::Type{T}) where {T<:Real} = zero(T)
probe_value(::Type{Bool}) = false
probe_value(::Type{P}) where {P<:StaticArray} = zero(P)
probe_value(::Type{P}) where {P} = P()
