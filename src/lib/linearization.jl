module Linearization

using ComponentArrays, StaticArrays, UnPack, LinearAlgebra
using FiniteDiff: finite_difference_jacobian! as jacobian!
using ControlSystems: ControlSystemsBase, ControlSystems, ss
using RobustAndOptimalControl
using DataStructures

using Flight.FlightCore

export LinearizedSS, linearize

################################################################################
########################### LinearizedSS ###################################

#represents the linearization
#ẋ = ẋ0 + A(x-x0) + B(u-u0)
#y = y0 + C(x-x0) + D(u-u0)

#of a nonlinear state-space system
#ẋ = f(x,u)
#y = g(x,u)

#around (x0, u0)

const tV = AbstractVector{<:Float64}
const tM = AbstractMatrix{<:Float64}

struct LinearizedSS{ LX, LU, LY, #state, input and output vector lengths
                        tX <: tV, tU <: tV, tY <: tV,
                        tA <: tM, tB <: tM, tC <: tM, tD <: tM} <: ModelDefinition

    ẋ0::tX; x0::tX; u0::tU; y0::tY; #reference values (for linearized systems)
    A::tA; B::tB; C::tC; D::tD; #state-space matrices
    x_cache::tX; y_cache::tY; y_cache_out::tY;
    Δx_cache::tX; Δu_cache::tU

    function LinearizedSS(ẋ0, x0, u0, y0, A, B, C, D)

        lengths = map(length, (x0, u0, y0))
        types = map(typeof, (x0, u0, y0, A, B, C, D))

        vectors = map(copy, (ẋ0, x0, u0, y0))
        matrices = map(copy, (A, B, C, D))
        caches = map(copy, (x0, y0, y0, x0, u0))

        new{lengths..., types...}(vectors..., matrices..., caches...)

    end

end

LinearizedSS(; ẋ0, x0, u0, y0, A, B, C, D) = LinearizedSS(ẋ0, x0, u0, y0, A, B, C, D)


################################## linearize ###################################

function linearize(f::Function, g::Function, x0::AbstractVector{<:Real}, u0::AbstractVector{<:Real})
    ẋ0 = f(x0, u0)
    y0 = g(x0, u0)
    A, B, C, D = ss_matrices(f, g, x0, u0)
    LinearizedSS(; ẋ0, x0, u0, y0, A, B, C, D)
end

# f and g should have the following signatures:
# ẋ = f(x::AbstractVector{<:Real}, u::AbstractVector{<:Real})::FieldVector
# y = g(x::AbstractVector{<:Real}, u::AbstractVector{<:Real})::FieldVector
function linearize(f::Function, g::Function, x0::FieldVector, u0::FieldVector)

    ẋ0 = f(x0, u0)::FieldVector
    y0 = g(x0, u0)::FieldVector
    A, B, C, D = ss_matrices(f, g, x0, u0)

    x_axis = Axis(propertynames(x0))
    u_axis = Axis(propertynames(u0))
    y_axis = Axis(propertynames(y0))

    ẋ0_cv = ComponentVector(Vector(ẋ0), x_axis)
    x0_cv = ComponentVector(Vector(x0), x_axis)
    u0_cv = ComponentVector(Vector(u0), u_axis)
    y0_cv = ComponentVector(Vector(y0), y_axis)

    A_cm = ComponentMatrix(A, x_axis, x_axis)
    B_cm = ComponentMatrix(B, x_axis, u_axis)
    C_cm = ComponentMatrix(C, y_axis, x_axis)
    D_cm = ComponentMatrix(D, y_axis, u_axis)

    LinearizedSS(ẋ0_cv, x0_cv, u0_cv, y0_cv, A_cm, B_cm, C_cm, D_cm)

end

function ss_matrices(f::Function, g::Function, x0::AbstractVector{<:Real}, u0::AbstractVector{<:Real})

    y0 = g(x0, u0)

    f_A!(ẋ, x) = (ẋ .= f(x, u0))
    f_B!(ẋ, u) = (ẋ .= f(x0, u))
    f_C!(y, x) = (y .= g(x, u0))
    f_D!(y, u) = (y .= g(x0, u))

    #preallocate mutable arrays
    A = (x0 * x0') |> Matrix
    B = (x0 * u0') |> Matrix
    C = (y0 * x0') |> Matrix
    D = (y0 * u0') |> Matrix

    jacobian!(A, f_A!, Vector(x0))
    jacobian!(B, f_B!, Vector(u0))
    jacobian!(C, f_C!, Vector(x0))
    jacobian!(D, f_D!, Vector(u0))

    return (A, B, C, D)

end

function subsystem(cmp::LinearizedSS; x = keys(cmp.x0), u = keys(cmp.u0), y = keys(cmp.y0))

    #to do: generalize for scalars

    x_ind = x
    u_ind = u
    y_ind = y

    ẋ0 = cmp.ẋ0[x_ind]
    x0 = cmp.x0[x_ind]
    u0 = cmp.u0[u_ind]
    y0 = cmp.y0[y_ind]

    A = cmp.A[x_ind, x_ind]
    B = cmp.B[x_ind, u_ind]
    C = cmp.C[y_ind, x_ind]
    D = cmp.D[y_ind, u_ind]

    return LinearizedSS(; ẋ0, x0, u0, y0, A, B, C, D)

end


############################# Update Functions #################################

Modeling.X(cmp::LinearizedSS) = copy(cmp.x0)
Modeling.U(cmp::LinearizedSS) = copy(cmp.u0)
Modeling.Y(cmp::LinearizedSS) = SVector{length(cmp.y0)}(cmp.y0)

@no_periodic LinearizedSS
@no_step LinearizedSS

function Modeling.f_ode!(mdl::Model{<:LinearizedSS{LX, LU, LY}}) where {LX, LU, LY}

    @unpack ẋ, x, u, y, parameters = mdl
    @unpack ẋ0, x0, u0, y0, A, B, C, D, x_cache, y_cache, y_cache_out, Δx_cache, Δu_cache = parameters

    #non-allocating equivalent of:
    #ẋ = ẋ0 + A * (x - x0) + B * (u - u0)
    #y = y0 + C * (x - x0) + D * (u - u0)

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

    mdl.y = SVector{LY}(y_cache_out)

    return nothing

end


############################### ControlSystems #################################

ControlSystems.ss(cmp::LinearizedSS) = ControlSystems.ss(cmp.A, cmp.B, cmp.C, cmp.D)

function RobustAndOptimalControl.named_ss(lss::LinearizedSS)
    x_labels, u_labels, y_labels = map(collect ∘ propertynames, (lss.x0, lss.u0, lss.y0))
    named_ss(ss(lss), x = x_labels, u = u_labels, y = y_labels)
end


function subsystem(nss::RobustAndOptimalControl.NamedStateSpace;
                  x = nss.x, u = nss.u, y = nss.y)

    #to do: generalize for scalars

    x_axis = Axis(nss.x)
    u_axis = Axis(nss.u)
    y_axis = Axis(nss.y)

    A_nss = ComponentMatrix(nss.A, x_axis, x_axis)
    B_nss = ComponentMatrix(nss.B, x_axis, u_axis)
    C_nss = ComponentMatrix(nss.C, y_axis, x_axis)
    D_nss = ComponentMatrix(nss.D, y_axis, u_axis)

    A_sub = A_nss[x, x]
    B_sub = B_nss[x, u]
    C_sub = C_nss[y, x]
    D_sub = D_nss[y, u]

    ss_sub = ss(A_sub, B_sub, C_sub, D_sub)

    named_ss(ss_sub; x, u, y)

end

function LinearizedSS(mdl::ControlSystemsBase.StateSpace{ControlSystemsBase.Continuous, <:AbstractFloat})
    @unpack A, B, C, D, nx, nu, ny = mdl
    ẋ0 = zeros(nx); x0 = zeros(nx); u0 = zeros(nu); y0 = zeros(ny)
    LinearizedSS(; ẋ0, x0, u0, y0, A, B, C, D)
end


end #module