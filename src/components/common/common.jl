module Common

using ComponentArrays, StaticArrays, UnPack, LinearAlgebra
using Flight.Systems

import ControlSystems #avoid bringing ControlSystems.StateSpace in scope


################################################################################
########################### StateSpace Model ###################################

#this yields a non-allocating ComponentVector with an underlying StaticVector
#array, which preserves the original ComponentVector's labels
function StaticArrays.SVector{L}(x::ComponentVector) where {L}
    ComponentVector(SVector{L}(getdata(x)), getaxes(x))
end

const tV = AbstractVector{<:Float64}
const tM = AbstractMatrix{<:Float64}

struct StateSpace{  LX, LU, LY, #state, input and output vector lengths
                tX <: tV, tU <: tV, tY <: tV,
                tA <: tM, tB <: tM, tC <: tM, tD <: tM} <: Component

    ẋ0::tX; x0::tX; u0::tU; y0::tY; #reference values (for linearized systems)
    A::tA; B::tB; C::tC; D::tD; #state-space matrices
    x_cache::tX; y_cache::tY; y_cache_out::tY;
    Δx_cache::tX; Δu_cache::tU

    function StateSpace(ẋ0, x0, u0, y0, A, B, C, D)

        lengths = map(length, (x0, u0, y0))
        types = map(typeof, (x0, u0, y0, A, B, C, D))

        vectors = map(copy, (ẋ0, x0, u0, y0))
        matrices = map(copy, (A, B, C, D))
        caches = map(copy, (x0, y0, y0, x0, u0))

        new{lengths..., types...}(vectors..., matrices..., caches...)

    end

end

StateSpace(; ẋ0, x0, u0, y0, A, B, C, D) = StateSpace(ẋ0, x0, u0, y0, A, B, C, D)

ControlSystems.ss(cmp::StateSpace) = ControlSystems.ss(cmp.A, cmp.B, cmp.C, cmp.D)

Systems.init(::SystemX, cmp::StateSpace) = copy(cmp.x0)
Systems.init(::SystemU, cmp::StateSpace) = copy(cmp.u0)
Systems.init(::SystemY, cmp::StateSpace) = SVector{length(cmp.y0)}(cmp.y0)

function Systems.f_ode!(sys::System{<:StateSpace{LX, LU, LY}}) where {LX, LU, LY}

    @unpack ẋ, x, u, y, params = sys
    @unpack ẋ0, x0, u0, y0, A, B, C, D, x_cache, y_cache, y_cache_out, Δx_cache, Δu_cache = params

    #the following is equivalent to:
    #ẋ = ẋ0 + A * (x - x0) + B * (u - u0)
    #y = y0 + C * (x - x0) + D * (u - u0)
    #... but without allocations

    @. Δx_cache = x - x0
    @. Δu_cache = u - u0

    ẋ .= ẋ0
    mul!(x_cache, A, Δx_cache)
    ẋ .+= x_cache
    mul!(x_cache, B, Δu_cache)
    ẋ .+= x_cache

    y_cache_out .= y0
    mul!(y_cache, C, Δx_cache)
    y_cache_out .+= y_cache
    mul!(y_cache, D, Δu_cache)
    y_cache_out .+= y_cache

    sys.y = SVector{LY}(y_cache_out)

    return nothing

end

function Base.filter(cmp::StateSpace; x = (), u = (), y = ())

    x_ind = (!isempty(x) ? x : keys(cmp.x0))
    u_ind = (!isempty(u) ? u : keys(cmp.u0))
    y_ind = (!isempty(y) ? y : keys(cmp.z0))

    ẋ0 = cmp.ẋ0[x_ind]
    x0 = cmp.x0[x_ind]
    u0 = cmp.u0[u_ind]
    y0 = cmp.y0[y_ind]
    A = cmp.A[x_ind, x_ind]
    B = cmp.B[x_ind, u_ind]
    C = cmp.C[y_ind, x_ind]
    D = cmp.D[y_ind, u_ind]

    return StateSpace(; ẋ0, x0, u0, y0, A, B, C, D)

end

end #module