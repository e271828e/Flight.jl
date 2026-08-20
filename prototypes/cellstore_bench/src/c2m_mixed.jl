# C2M — C2's representation generalized to mixed-leaf cells (added 2026-08-20,
# after D-162; the C1/C2 rows and their numbers stay frozen). Per-eltype
# homogeneous buffers as before, but a cell may span several of them: its
# address carries one cursor per distinct leaf eltype of the port type, as an
# `NTuple` *field*. `K` is a pure function of `P` — the count of `P`'s distinct
# leaf eltypes — so the address-token type domain is still keyed by `P` alone,
# and two instances of one component type share one `Entry` type and one
# compiled body. The homogeneous cell is the `K = 1` case; C2M uses this
# address for *every* cell, because that is what the kernel would do globally.
#
# The risk this point measures is D-162's, restated for the new shape: the
# `NTuple` offsets must survive as runtime data behind `invoke_sweep`'s
# `@constprop :none` barrier (no per-instance respecialization), and the store
# staying one concrete type per model must keep the compile curve bounded by
# chunk-type count, not model size.

# --- the mixed model ----------------------------------------------------------
# `TaggedPose3` carries an `Int32` leaf beside the `T` leaves, so the mixing is
# real at *both* activations — unlike a pinned `Float64`, which only separates
# from `T` on `Dual8`. `MixedSource`/`MixedWorkhorse` are their `model.jl`
# counterparts with the chained port swapped; bodies keep the same op count so
# the compile columns stay comparable. `Filter` and `Junction` are reused
# verbatim, keeping homogeneous (`K = 1`) cells in the same model.

struct TaggedPose3{T}
    p::SVector{3,T}
    ψ::T
    tag::Int32
end

struct MixedSource
    ω::Float64
end

input_types(::MixedSource, ::Type{T}) where {T} = NamedTuple()
output_types(::MixedSource, ::Type{T}) where {T} =
    (cmd = T, wind_a = SVector{3,T}, wind_b = SVector{3,T}, pose = TaggedPose3{T})
init_x(::MixedSource, ::Type{T}) where {T} = (φ = T(0.1),)

@inline function h_xu(c::MixedSource, x, u)
    s, cs = sincos(c.ω * x.φ)
    cmd = 0.5 * s + 0.25 * cs
    wind_a = SVector(s, cs, s * cs)
    wind_b = SVector(cs - s, s + cs, 0.5 * s)
    pose = TaggedPose3(SVector(x.φ, s, cs), c.ω * cs, Int32(1))
    (; cmd, wind_a, wind_b, pose)
end

struct MixedWorkhorse
    k::Float64
    τ::Float64
end

input_types(::MixedWorkhorse, ::Type{T}) where {T} =
    (cmd = T, wind = SVector{3,T}, pose_in = TaggedPose3{T})
output_types(::MixedWorkhorse, ::Type{T}) where {T} =
    (pose = TaggedPose3{T}, rate = SVector{3,T}, load = T)
init_x(::MixedWorkhorse, ::Type{T}) where {T} =
    (v = SVector{3,T}(1.0, 0.0, -0.5), ω = SVector{3,T}(0.0, 0.1, 0.0), e = T(0.0))

@inline function h_xu(c::MixedWorkhorse, x, u)
    vr = x.v - u.wind
    q = 0.5 * dot(vr, vr)
    ω = x.ω + c.k * cross(vr, u.pose_in.p)
    ψ = u.pose_in.ψ + c.τ * x.e
    p = u.pose_in.p + vr * c.τ + c.k * ω
    load = q * c.k - u.cmd * x.e
    (pose = TaggedPose3(p, ψ, u.pose_in.tag + Int32(1)), rate = ω, load = load)
end

"""`chain_spec(N)` with the chained port mixed: `Int32` tag beside `T` leaves."""
function mixed_chain_spec(N::Int)
    comps = CompSpec[
        CompSpec(:src, MixedSource(0.7), []),
        CompSpec(:fil, Filter(0.25), [:in => (1, :cmd)]),
        CompSpec(:jun, Junction(), [:a => (1, :wind_a), :b => (1, :wind_b)]),
    ]
    for i in 1:N
        pose_src = i == 1 ? (1, :pose) : (3 + i - 1, :pose)
        push!(comps, CompSpec(Symbol(:wrk, i), MixedWorkhorse(0.1 + 0.01i, 0.02),
            [:cmd => (2, :out), :wind => (3, :sum), :pose_in => pose_src]))
    end
    ModelSpec(comps)
end

# --- leaf shape helpers (layout time) ----------------------------------------

"""The element type of each leaf a value of type `P` occupies, in flat order."""
leaf_types(::Type{P}) where {P<:Real} = Type[P]
leaf_types(::Type{P}) where {P<:StaticArray} = repeat(leaf_types(eltype(P)), length(P))
leaf_types(::Type{P}) where {P} =
    reduce(vcat, (leaf_types(FT) for FT in fieldtypes(P)); init = Type[])

"""Distinct leaf eltypes of `P`, in first-appearance order over the flat walk —
the canonical order `MixedAddr`'s cursors follow."""
distinct_leaf_eltypes(::Type{P}) where {P} = unique(leaf_types(P))

# --- the candidate ------------------------------------------------------------

struct C2M <: Candidate end

struct MixedAddr{P,K}
    offs::NTuple{K,Int}
end

struct C2MStore{NT<:NamedTuple}
    bufs::NT       # eltype name => that eltype's contiguous buffer
end

# The expression builders, generalized from one running base to one per eltype:
# a `Real` leaf reads/writes the buffer bound for its own eltype at
# `offs[k] + <static index>`. The `K = 1` expansion is `leaves.jl`'s expression.

function _mreconstruct_expr(::Type{P}, Ls::Vector, bases::Vector{Int}) where {P}
    if P <: Real
        k = findfirst(==(P), Ls)
        bases[k] += 1
        return :(@inbounds $(Symbol(:buf, k))[offs[$k] + $(bases[k])])
    elseif P <: StaticArray
        args = [_mreconstruct_expr(eltype(P), Ls, bases) for _ in 1:length(P)]
        return Expr(:call, P, args...)
    else
        args = [_mreconstruct_expr(FT, Ls, bases) for FT in fieldtypes(P)]
        P <: NamedTuple && return Expr(:call, P, Expr(:tuple, args...))
        return Expr(:call, P, args...)
    end
end

function _mflatten_expr(::Type{P}, v, Ls::Vector, bases::Vector{Int}) where {P}
    stmts = Expr[]
    if P <: Real
        k = findfirst(==(P), Ls)
        bases[k] += 1
        push!(stmts, :(@inbounds $(Symbol(:buf, k))[offs[$k] + $(bases[k])] = $v))
    elseif P <: StaticArray
        for i in 1:length(P)
            push!(stmts, _mflatten_expr(eltype(P), :(@inbounds $v[$i]), Ls, bases))
        end
    else
        for (i, FT) in enumerate(fieldtypes(P))
            push!(stmts, _mflatten_expr(FT, :(getfield($v, $i)), Ls, bases))
        end
    end
    Expr(:block, stmts...)
end

_buf_binds(Ls) = [:($(Symbol(:buf, k)) = getfield(s.bufs, $(QuoteNode(Symbol(L)))))
                  for (k, L) in enumerate(Ls)]

@generated function gather(s::C2MStore, a::MixedAddr{P,K}) where {P,K}
    Ls = distinct_leaf_eltypes(P)
    expr = _mreconstruct_expr(P, Ls, zeros(Int, K))
    quote
        $(Expr(:meta, :inline))
        offs = a.offs
        $(_buf_binds(Ls)...)
        $expr
    end
end

@generated function scatter!(s::C2MStore, a::MixedAddr{P,K}, v) where {P,K}
    Ls = distinct_leaf_eltypes(P)
    block = _mflatten_expr(P, :v, Ls, zeros(Int, K))
    quote
        $(Expr(:meta, :inline))
        offs = a.offs
        $(_buf_binds(Ls)...)
        $block
        nothing
    end
end

function build_store(::Type{C2M}, spec::ModelSpec, ::Type{T}) where {T}
    inv = cell_inventory(spec, T)
    roster = unique(reduce(vcat, (leaf_types(P) for (_, _, P) in inv); init = Type[]))
    sort!(roster; by = string)                 # model-global buffer roster, name-sorted
    cursors = Dict{Type,Int}(L => 0 for L in roster)
    layout = Dict{Tuple{Int,Symbol},Any}()
    for (ci, port, P) in inv
        lts = leaf_types(P)
        Ls = unique(lts)
        layout[(ci, port)] = (Tuple(cursors[L] for L in Ls), P)
        for L in Ls
            cursors[L] += count(==(L), lts)
        end
    end
    bufs = NamedTuple{Tuple(Symbol.(roster))}(Tuple(zeros(L, cursors[L]) for L in roster))
    C2MStore(bufs), layout
end

function cell_addr(::Type{C2M}, layout, ci::Int, port::Symbol)
    offs, P = layout[(ci, port)]
    MixedAddr{P,length(offs)}(offs)
end
