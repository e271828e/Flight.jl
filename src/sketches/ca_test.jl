using ComponentArrays
using LinearAlgebra
using StaticArrays: SVector, MVector
using UnPack
using Flight.LBV
using Flight.WGS84
using Flight.Attitude
using Flight.Kinematics

const x_kin_template = ComponentVector((
    pos = ComponentVector((
        q_lb = zeros(4),
        q_el = zeros(4),
        h = 0.0)),
    vel = ComponentVector((
        ω_eb_b = zeros(3),
        v_eOb_b = zeros(3)))
    ))

#defining the type in addition allows for convenient dispatching on initialize!,
#but it is not necessary for type stability nor performance. type stability
#appears to be achieved by defining the template and then instantiating copies
#of it within other methods. however, note that new copies of x can only be
#created from x_kin_template, because there are no ComponentVector constructors
#for XKinWGS84CA{Vector{Float64}}
const KinWGS84Axes = typeof(getaxes(x_kin_template))
const XKinWGS84CA = ComponentVector{Float64, D, KinWGS84Axes} where {D <: AbstractVector{Float64}}

#maybe we could create a dedicated constructor... like this?
# XKinWGS84CA{Float64}() = similar(x_kin_template)


function initialize!(x::XKinWGS84CA, init::KinInit)

    @unpack q_nb, Ob, ω_lb_b, v_eOb_b = init

    h = Ob.h[1]
    (R_N, R_E) = radii(Ob)
    v_eOb_n = q_nb * v_eOb_b
    ω_el_n = SVector{3}(
        v_eOb_n[2] / (R_E + h),
        -v_eOb_n[1] / (R_N + h),
        0.0)

    ω_el_b = q_nb' * ω_el_n
    ω_eb_b = ω_el_b + ω_lb_b

    q_lb = q_nb #arbitrarily initialize ψ_nl to -1

    x.pos.q_lb .= q_lb #assignment without colon slicing courtesy of RQuat iterability
    x.pos.q_el .= ltf(Ob) #assignment without colon slicing courtesy of RQuat iterability
    x.pos.h = h
    x.vel.ω_eb_b .= ω_eb_b
    x.vel.v_eOb_b .= v_eOb_b

end

function initialize(init::KinInit)

    @unpack q_nb, Ob, ω_lb_b, v_eOb_b = init

    h = Ob.h[1]
    (R_N, R_E) = radii(Ob)
    v_eOb_n = q_nb * v_eOb_b
    ω_el_n = SVector{3}(
        v_eOb_n[2] / (R_E + h),
        -v_eOb_n[1] / (R_N + h),
        0.0)

    ω_el_b = q_nb' * ω_el_n
    ω_eb_b = ω_el_b + ω_lb_b

    q_lb = q_nb #arbitrarily initialize ψ_nl to -1

    x = similar(x_kin_template) #since x_kin_template, this creates a copy

    x.pos.q_lb .= q_lb #assignment without colon slicing courtesy of RQuat iterability
    x.pos.q_el .= ltf(Ob) #assignment without colon slicing courtesy of RQuat iterability
    x.pos.h = h
    x.vel.ω_eb_b .= ω_eb_b
    x.vel.v_eOb_b .= v_eOb_b

    return x

end

function PVDataWGS84(x::XKinWGS84CA)

    #careful here: x.pos.h, x.vel.ω_eb_b and x.vel.v_eOb_b create views (this is
    #how LBV behaves by design). to copy the data, we can extract their
    #components using slices
    q_lb = RQuat(x.pos.q_lb)
    q_el = RQuat(x.pos.q_el)
    h = x.pos.h[1]
    ω_eb_b = SVector{3}(x.vel.ω_eb_b)
    v_eOb_b = SVector{3}(x.vel.v_eOb_b)

    Ob = WGS84Pos(NVector(q_el), h)
    q_nl = Rz(ψ_nl(q_el))
    q_nb = q_nl ∘ q_lb
    q_eb = q_el ∘ q_lb

    (R_N, R_E) = radii(Ob)
    v_eOb_n = q_nb * v_eOb_b
    ω_el_n = SVector{3}(
        v_eOb_n[2] / (R_E + h),
        -v_eOb_n[1] / (R_N + h),
        0.0)

    ω_el_l = q_nl' * ω_el_n
    ω_el_b = q_lb' * ω_el_l
    ω_lb_b = ω_eb_b - ω_el_b

    ω_ie_e = SVector{3}(0, 0, ω_ie)
    ω_ie_b = q_eb' * ω_ie_e
    ω_ib_b = ω_ie_b + ω_eb_b

    pos = PosDataWGS84(q_lb, q_nl, q_nb, q_eb, q_el, Ob)
    vel = VelDataWGS84(ω_eb_b, ω_lb_b, ω_el_l, ω_ie_b, ω_ib_b, v_eOb_b, v_eOb_n)

    PVDataWGS84(pos, vel)

end