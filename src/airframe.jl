module Airframe

using StaticArrays: SVector, SMatrix
using LinearAlgebra
using UnPack
using ComponentArrays

using Flight.WGS84
using Flight.Attitude
using Flight.Kinematics
using Flight.Airdata
using Flight.System
import Flight.System: X, U, D, typeof_y, f_output!

export Wrench, MassData, AbstractComponent, ComponentFrame, ComponentGroup
export v2skew, inertia_wrench, gravity_wrench, f_vel!

Base.@kwdef struct MassData
    m::Float64 = 1.0
    J_Ob_b::SMatrix{3, 3, Float64, 9} = SMatrix{3,3,Float64}(I)
    r_ObG_b::SVector{3, Float64} = zeros(SVector{3})
end


const WrenchAxes = getaxes(ComponentVector(F = zeros(3), M = zeros(3)))

struct Wrench{D <: AbstractVector{Float64}} <: AbstractVector{Float64}
    data::ComponentVector{Float64, D, typeof(WrenchAxes)}
end

#avoid ComponentVector(F = F, M = M), which forwards to ComponentArray{T}
#(::NamedTuple), which needs to create the axes from scratch and is very slow.
#much faster to provide the axes directly

function Wrench(input::AbstractVector{<:Real})
    if length(input) != 6
        throw(ArgumentError("Got input length $(length(input)), expected 6"))
    end
    data = ComponentVector{Float64}(undef, WrenchAxes)
    data .= input
    Wrench(data)
end

function Wrench(; F = SVector(0,0,0), M = SVector(0,0,0))
    data = ComponentVector{Float64}(undef, WrenchAxes)
    data.F = F; data.M = M
    Wrench(data)
end
ComponentArrays.ComponentVector(wr::Wrench) = copy(wr.data)

Base.size(::Wrench) = (6,)
Base.length(::Wrench) = 6

Base.getproperty(wr::Wrench, i::Symbol) = getproperty(wr, Val(i))
Base.getproperty(wr::Wrench, ::Val{:data}) = getfield(wr, :data)
Base.getproperty(wr::Wrench, i::Val{S} where {S}) = getproperty(wr.data, i)

Base.setproperty!(wr::Wrench, i::Symbol, v) = setproperty!(wr, Val(i), v)
Base.setproperty!(wr::Wrench, ::Val{:F}, v) = (wr.F .= v) #gets the F block, then broadcasts
Base.setproperty!(wr::Wrench, ::Val{:M}, v) = (wr.M .= v) #gets the M block, then broadcasts

Base.getindex(wr::Wrench, i) = getindex(wr.data, i)
Base.setindex!(wr::Wrench, v, i) = setindex!(wr.data, v, i)

Base.eltype(::Wrench) = Float64 #helps with allocation efficiency
Base.similar(::Wrench) = Wrench(ComponentVector{Float64}(undef, WrenchAxes))
Base.similar(::Type{<:Wrench}) = Wrench(ComponentVector{Float64}(undef, WrenchAxes))

Base.show(io::IO, wr::Wrench) = print(io, "Wrench(F = $(wr.F), M = $(wr.M))")
Base.show(io::IO, ::MIME"text/plain", wr::Wrench) = print(io, "Wrench(F = $(wr.F), M = $(wr.M))")

#since Wrench <: AbstractVector, broadcasting works out of the box (and with
#it the non-broadcast operators +, -, etc), but it falls back to the default
#broadcast implementation in Base, which returns a generic Vector. we want
#broadcasted operations to return a Wrench. so we do this:

struct WrenchStyle{D} <: Broadcast.AbstractArrayStyle{1} end

WrenchStyle{D}(::Val{1}) where {D} = WrenchStyle{D}()
Base.BroadcastStyle(::Type{Wrench{D}}) where {D} = WrenchStyle{D}()
function Base.similar(::Broadcast.Broadcasted{WrenchStyle{D}}, ::Type{ElType}) where {D,ElType}
    similar(Wrench{D})
end
function Base.BroadcastStyle(::WrenchStyle{D1}, ::WrenchStyle{D2}) where {D1,D2}
    WrenchStyle{promote_type(D1, D2)}()
end


abstract type AbstractComponent <: AbstractSystem end

"""
#Specifies a local ComponentFrame fc(Oc, Ɛc) relative to the airframe reference
frame fb(Ob, Ɛb) by:
#a) the position vector of the local frame origin Oc relative to the reference
#frame origin Ob, projected in the reference frame axes
# b) the attitude of the local frame axes relative to the reference
#frame axes, given by rotation quaternion q_bc
"""
Base.@kwdef struct ComponentFrame
    r_ObOc_b::SVector{3,Float64} = zeros(SVector{3})
    q_bc::RQuat = RQuat()
end

"""
Translate a Wrench specified on a local ComponentFrame fc(Oc, εc) to the
airframe reference frame fb(Ob, εb) given the relative ComponentFrame
specification f_bc
"""
function Base.:*(f_bc::ComponentFrame, wr_Oc_c::Wrench)

    F_Oc_c = wr_Oc_c.F
    M_Oc_c = wr_Oc_c.M

    #project on the reference axes
    F_Oc_b = f_bc.q_bc * F_Oc_c
    M_Oc_b = f_bc.q_bc * M_Oc_c

    #translate them to airframe origin
    F_Ob_b = F_Oc_b
    M_Ob_b = M_Oc_b + f_bc.r_ObOc_b × F_Oc_b
    Wrench(F = F_Ob_b, M = M_Ob_b) #wr_Ob_b

end

struct ComponentGroup{N,T,L,C} <: AbstractComponent
    function ComponentGroup(nt::NamedTuple{L, NTuple{N, T}}) where {L, N, T <: AbstractComponent} #Dicts are not ordered, so they won't do
        new{N,T,L,values(nt)}()
    end
end

X(::ComponentGroup{N,T,L,C}) where {N,T,L,C} = ComponentVector(NamedTuple{L}(X.(C)))
U(::ComponentGroup{N,T,L,C}) where {N,T,L,C} = ComponentVector(NamedTuple{L}(U.(C)))
D(::ComponentGroup{N,T,L,C}) where {N,T,L,C} = NamedTuple{L}(D.(C))
typeof_y(::ComponentGroup{N,T,L,C}) where {N,T,L,C} = NamedTuple{L, NTuple{N, typeof_y(T)}}

@inline function f_output!(ẋ::Any, x::Any, u::Any, t::Real, data::Any, ::ComponentGroup{N,T,L,C}) where {N,T,L,C}

    v = Vector{typeof_y(T)}(undef, N)
    for (i, k) in enumerate(valkeys(x)) #valkeys is the only way to avoid allocations
        v[i] = f_output!(getproperty(ẋ, k), getproperty(x, k), getproperty(u, k), t, data[i], C[i])
    end
    tup = Tuple(v)::NTuple{N, typeof_y(T)}
    return NamedTuple{L}(tup)::NamedTuple{L, NTuple{N, typeof_y(T)}}
end


function inertia_wrench(mass::MassData, vel::VelY, h_rot_b::AbstractVector{<:Real})

    @unpack m, J_Ob_b, r_ObG_b = mass
    @unpack ω_ie_b, ω_eb_b, ω_ib_b, v_eOb_b = vel #these are already SVectors

    h_rot_b = SVector{3,Float64}(h_rot_b)

    #angular momentum of the overall airframe as a rigid body
    h_rbd_b = J_Ob_b * ω_ib_b

    #total angular momentum
    h_all_b = h_rbd_b + h_rot_b

    #exact
    a_1_b = (ω_eb_b + 2 * ω_ie_b) × v_eOb_b
    F_in_Ob_b = -m * (a_1_b + ω_ib_b × (ω_ib_b × r_ObG_b) + r_ObG_b × (ω_eb_b × ω_ie_b ))
    M_in_Ob_b = - ( J_Ob_b * (ω_ie_b × ω_eb_b) + ω_ib_b × h_all_b + m * r_ObG_b × a_1_b)

    Wrench(F = F_in_Ob_b, M = M_in_Ob_b)

end

function gravity_wrench(mass::MassData, pos::PosY)

    #strictly, the gravity vector should be evaluated at G, with its direction
    #given by the z-axis of LTF(G). however, since g(G) ≈ g(Ob) and LTF(G) ≈
    #LTF(Ob), we can instead evaluate g at Ob, assuming its direction given by
    #LTF(Ob), and then apply it at G.
    g_G_l = gravity(pos.Ob)

    #the resultant consists of the force of gravity acting on G along the local
    #vertical and a null torque
    F_G_l = mass.m * g_G_l
    M_G_l = zeros(SVector{3})
    wr_G_l = Wrench(F = F_G_l, M = M_G_l)

    #with the previous assumption, the transformation from body frame to local
    #gravity frame is given by the translation r_ObG_b and the (passive)
    #rotation from b to LTF(Ob) (instead of LTF(G)), which is given by pos.l_b'
    wr_Oc_c = wr_G_l
    f_bc = ComponentFrame(r_ObOc_b = mass.r_ObG_b, q_bc = pos.q_lb')
    return f_bc * wr_Oc_c #wr_Ob_b

end

function f_vel!(ẋ_vel, wr_ext_Ob_b::Wrench, h_rot_b::AbstractVector{<:Real},
                mass::MassData, y_pv::PosVelY)

    #wr_ext_Ob_b: External Wrench on the airframe due to aircraft components

    #h_rot_b: Additional angular momentum due to rotating aircraft components
    #(computed using their angular velocity wrt the airframe, not the inertial
    #frame)

    #wr_ext_Ob_b and h_rot_b, as well as mass data, are produced by aircraft
    #components, so they must be computed by the aircraft's x_dot method. and,
    #since y_pv is needed by those components, it must be called from the
    #aircraft's kinematic state vector

    #clearly, r_ObG_b cannot be arbitrarily large, because J_Ob_b is larger than
    #J_G_b (Steiner). therefore, at some point J_G_b would become zero (or at
    #least singular)!

    @unpack m, J_Ob_b, r_ObG_b = mass
    @unpack ω_eb_b, ω_ie_b, v_eOb_b = y_pv.vel
    @unpack Ob, q_eb = y_pv.pos

    #preallocating is faster than directly concatenating the blocks
    A = Array{Float64}(undef, (6,6))

    r_ObG_b_sk = v2skew(r_ObG_b)
    A[1:3, 1:3] .= J_Ob_b
    A[1:3, 4:6] .= m * r_ObG_b_sk
    A[4:6, 1:3] .= -m * r_ObG_b_sk
    A[4:6, 4:6] .= m * SMatrix{3,3,Float64}(I)

    A = SMatrix{6,6}(A)

    wr_g_Ob_b = gravity_wrench(mass, y_pv.pos)
    wr_in_Ob_b = inertia_wrench(mass, y_pv.vel, h_rot_b)
    wr_Ob_b = wr_ext_Ob_b + wr_g_Ob_b + wr_in_Ob_b
    b = SVector{6}([wr_Ob_b.M ; wr_Ob_b.F])

    ẋ_vel .= A\b

    α_eb_b = SVector{3}(ẋ_vel.ω_eb_b) #α_eb_b == ω_eb_b_dot
    α_ib_b = α_eb_b - ω_eb_b × ω_ie_b

    v̇_eOb_b = SVector{3}(ẋ_vel.v_eOb_b)
    r_eO_e = rECEF(Ob)
    r_eO_b = q_eb' * r_eO_e
    a_eOb_b = v̇_eOb_b + ω_eb_b × v_eOb_b
    a_iOb_b = v̇_eOb_b + (ω_eb_b + 2ω_ie_b) × v_eOb_b + ω_ie_b × (ω_ie_b × r_eO_b)

    AccY(α_eb_b, α_ib_b, a_eOb_b, a_iOb_b)

end


"""
Computes the skew-symmetric matrix corresponding to 3-element vector v.
"""
function v2skew(v::AbstractVector{T}) where {T<:Real}
    #much slower, each indexing operation yields an allocation
    # [0. -v[3] v[2]; v[3] 0. -v[1]; -v[2] v[1] 0.]
    M = zeros(T, 3, 3)
                    M[1,2] = -v[3];  M[1,3] = v[2]
    M[2,1] = v[3];                   M[2,3] = -v[1]
    M[3,1] = -v[2]; M[3,2] = v[1]

    SMatrix{3,3}(M)
end


end #module