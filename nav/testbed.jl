module TestBed

using LinearAlgebra, StaticArrays, ComponentArrays, UnPack, HDF5, Interpolations
using Logging
using Reexport
using NLopt

using Flight.FlightCore
using Flight.FlightLib


################################################################################
############################# DynamicInverter ##################################

    #to interact with this vehicle we need to dig into its components and assign
    #the DynamicInverter inputs directly. is there a still simple but more
    #elegant way? maybe not. moreover, later on we may want to add another
    #component in parallel with the dynamic inverter that implements a
    #proportional controller to track v and ω instead of their derivatives.

    #a further question is whether we could command ω_wb_b instead of ω_eb_b,
    #since that one seems more useful. in theory, we could do this by modifying
    #the dynamic inversion layer to command ω̇_wb_b directly. but this requires
    #expressing ω̇_wb_b in terms of ω̇_eb_b. but this is a nightmare. we could
    #instead try to drive the error in ω_wb_b to zero by commanding ω̇_eb_b and
    #cancelling the discrepancy between them (which is the transport rate) by
    #adding an integral term, as we would for an external disturbance.

    #but wait, wait. what we can plug into our ω_eb_b control loop as a
    #reference input is ω_eb_b_ref. however, we would like to command ω_wb_b_ref
    #(or maybe ω_nb_b_ref). but the kinematic relation between them is perfectly
    #known. we just need to compute the required ω_eb_b!!!

    #now, we might want to add an outer loop to explicitly control the attitude
    #with respect to the NED frame. we might do so by designing a set of
    #independent scalar controllers that track ψ, θ and ϕ. we model each one of
    #them as a simple integrator. that is, for example, u = θ̇, x = θ. with this
    #open loop system, a simple proportional controller on each channel would be
    #enough (the plant is already type 1). once the required ψ_dot_ref,
    #θ_dot_ref and ϕ_dot_ref are determined, we would find the corresponding
    #ω_nb_b_ref = (p_ref, q_ref, r_ref), which is given by the kinematic
    #equations:

    #p_ref = ϕ_dot_ref - q * sin(ϕ)*tan(θ) - r*cos(ϕ)*tan(θ)
    #q_ref = 1/cos(ϕ) * θ_dot_ref + r * tan(ϕ)
    #r_ref = cos(θ)/cos(ϕ) * ψ_dot_ref - q * tan(ϕ)

    #another question is, since we probably want to be able to
    #control both ω_wb_b directly and the attitude, we are probably dealing with
    #a multi-mode control system. this means it might be better to implement it
    #as a discrete controller within Avionics.

@kwdef mutable struct DynamicInverterU
    ω̇_eb_b::MVector{3,Float64} = zeros(3)
    v̇_eb_b::MVector{3,Float64} = zeros(3)
end

@kwdef struct DynamicInverterY
    wr_b::SVector{3,Float64} = zeros(SVector{3})
end

struct DynamicInverter <: AbstractComponentSet end

Modeling.U(::DynamicInverter) = DynamicInverterU()
Modeling.Y(::DynamicInverter) = DynamicInverterY()

Dynamics.get_mp_b(::System{DynamicInverter}) = MassProperties(RigidBodyDistribution())
Dynamics.get_hr_b(::System{DynamicInverter}) = zeros(SVector{3})
Dynamics.get_wr_b(sys::System{DynamicInverter}) = sys.y.wr_b


############################# Update Methods ###################################

#the following computations are valid for any mp_Σ_b and ho_Σ_b returned by the
#methods above. since for this vehicle these values are arbitrary, in principle
#we could ensure that m_Σ = 1, J_Σ_b = I, r_bc = zeros(3) and ho_Σ = zeros(3),
#and exploit this to simplify the computation. however, sticking to the general
#case keeps the commonality with VehicleDynamics

function Modeling.f_ode!(sys::System{DynamicInverter},
                        kin_data::KinData,
                        ::AirflowData,
                        ::System{<:AbstractTerrain})

    @unpack ω_eb_b, v_eb_b = kin_data

    ω̇_eb_b = SVector(sys.u.ω̇_eb_b)
    v̇_eb_b = SVector(sys.u.v̇_eb_b)

    #translate derivatives to c
    ω̇_ec_c = ω̇_eb_b
    v̇_ec_c = v̇_eb_b + ω̇_eb_b × r_bc_b

    ω_ie_e = SVector{3, Float64}(0, 0, ω_ie) #use WGS84 constant
    ω_ie_b = q_eb'(ω_ie_e)

    mp_Σ_b = get_mp_b(components) #mass properties at b projected in b
    ho_Σ_b = get_hr_b(components) #internal angular momentum projected in b

    #frame transform from c (CoM) to b (body)
    r_bc_b = mp_Σ_b.r_OG
    t_cb = FrameTransform(r = -r_bc_b) #pure translation
    t_bc = t_cb'

    #translate data to frame c
    mp_Σ_c = t_cb(mp_Σ_b)
    ho_Σ_c = ho_Σ_b

    m_Σ = mp_Σ_c.m; J_Σ_c = mp_Σ_c.J

    ω_ec_c = ω_eb_b
    v_ec_c = v_eb_b + ω_ec_c × r_bc_b

    ω_ie_c = ω_ie_b
    ω_ic_c = ω_ie_c + ω_ec_c

    #compute geographic position of Oc
    r_bc_e = q_eb(r_bc_b)
    r_ec_e = r_eb_e + r_bc_e
    Oc = Cartesian(r_ec_e)

    #define auxiliary local-level frame l with Ol = Oc
    q_el = ltf(Oc)
    q_be = q_eb'
    q_ce = q_be
    q_cl = q_ce ∘ q_el

    #compute gravity at c
    g_c_l = SVector{3,Float64}(0, 0, gravity(Oc)) #gravity at c, l coordinates\
    g_c_c = q_cl(g_c_l) #gravity at c, c coordinates

    #solve for F_Σ_c and τ_Σ_c
    hc_Σ_c = J_Σ_c * ω_ic_c + ho_Σ_c
    F_Σ_c = m_Σ * (v̇_ec_c - g_c_c + (ω_ec_c + 2ω_ie_c) × v_ec_c)
    τ_Σ_c = J_Σ_c * (ω̇_ec_c + ω_ie_c × ω_ec_c) + ω_ic_c × hc_Σ_c

    wr_Σ_c = Wrench(F_Σ_c, τ_Σ_c)
    wr_Σ_b = t_bc(wr_Σ_c)

    sys.y = DynamicInverterY(; wr_Σ_b)

end

@no_step DynamicInverter
@no_disc DynamicInverter

################################################################################
################################# Templates ####################################

const Vehicle{F, K, T} = AircraftBase.Vehicle{F, K, T} where {F <: Components, K <: AbstractKinematicDescriptor, T <: AbstractTerrain}


end