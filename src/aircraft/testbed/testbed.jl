module TestBed

using LinearAlgebra, StaticArrays, ComponentArrays, UnPack, HDF5, Interpolations
using Logging
using Reexport
using NLopt

using Flight.FlightCore
using Flight.FlightLib


################################################################################
################################ Airframe #####################################

struct Airframe <: SystemDefinition end

# This component represents the platform's structure, together with any
# components rigidly attached to it, such as powerplant or landing gear, but not
# payload or fuel contents. Its mass corresponds roughly to the aircraft's
# Standard Empty Weight

#Airframe mass properties computed in the vehicle reference frame b
const mp_b_afm = let
    #define the airframe as a RigidBodyDistribution
    afm_c = RigidBodyDistribution(767.0, SA[820.0 0 0; 0 1164.0 0; 0 0 1702.0])
    #a RigidBodyDistribution is always specified in a reference frame c with
    #origin at its center of mass. now, define the transform from vehicle
    #reference frame b to reference frame c (pure translation)
    t_bc = FrameTransform(r = SVector{3}(0.056, 0, 0.582))
    #translate the airframe's mass properties to frame b
    MassProperties(afm_c, t_bc)
end

#the airframe itself receives no external actions. these are considered to act
###upon the vehicle's aerodynamics, power plant and landing gear. the same goes
#for rotational angular momentum.


################################################################################
################################ Components ######################################

struct DynamicInverter <: AbstractComponents end

Dynamics.get_mp_b(::System{DynamicInverter}) = MassProperties(RigidBodyDistribution(1.0, SA[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]))
Dynamics.get_wr_b(::System{DynamicInverter}) = Wrench()
Dynamics.get_hr_b(::System{DynamicInverter}) = zeros(SVector{3})


############################# Update Methods ###################################

function Systems.f_ode!(inverter::System{DynamicInverter},
                        kin::KinData,
                        air::AirData,
                        trn::AbstractTerrain)

    Systems.update_y!(components)

end

function Systems.f_step!(components::System{<:Components})
    @unpack aero, ldg, pwp, fuel = components

    f_step!(aero)
    f_step!(ldg)
    f_step!(pwp, is_fuel_available(fuel))

end


#################################### GUI #######################################

function GUI.draw!( components::System{<:Components}, ::System{A},
                    p_open::Ref{Bool} = Ref(true),
                    label::String = "Cessna 172 Components") where {A<:AbstractAvionics}

    @unpack act, pwp, ldg, aero, fuel, pld = components

    CImGui.Begin(label, p_open)

        @cstatic(c_act=false, c_aero=false, c_ldg=false, c_pwp=false, c_fuel=false, c_pld=false,
        begin
            @c CImGui.Checkbox("Actuation", &c_act)
            @c CImGui.Checkbox("Aerodynamics", &c_aero)
            @c CImGui.Checkbox("Landing Gear", &c_ldg)
            @c CImGui.Checkbox("Power Plant", &c_pwp)
            @c CImGui.Checkbox("Fuel", &c_fuel)
            @c CImGui.Checkbox("Payload", &c_pld)
            if c_act
                if A === NoAvionics
                    @c GUI.draw!(act, &c_act)
                else
                    @c GUI.draw(act, &c_act)
                end
            end
            c_aero && @c GUI.draw(aero, &c_aero)
            c_ldg && @c GUI.draw(ldg, &c_ldg)
            c_pwp && @c GUI.draw(pwp, &c_pwp)
            c_fuel && @c GUI.draw(fuel, &c_fuel)
            c_pld && @c GUI.draw!(pld, &c_pld)
        end)

    CImGui.End()

end


################################################################################
################################# Templates ####################################

const Vehicle{F, K, T} = AircraftBase.Vehicle{F, K, T} where {F <: Components, K <: AbstractKinematicDescriptor, T <: AbstractTerrain}

############################### Trimming #######################################
################################################################################

#first 2 are aircraft-agnostic
@kwdef struct TrimState <: AbstractTrimState{7}
    α_a::Float64 = 0.1 #angle of attack, aerodynamic axes
    φ_nb::Float64 = 0.0 #bank angle
    n_eng::Float64 = 0.75 #normalized engine speed (ω/ω_rated)
    throttle::Float64 = 0.47
    aileron::Float64 = 0.014
    elevator::Float64 = -0.0015
    rudder::Float64 = 0.02 #rudder↑ -> aero.u.r↓ -> right yaw
end

@kwdef struct TrimParameters <: AbstractTrimParameters
    Ob::Geographic{NVector, Ellipsoidal} = Geographic(NVector(), HEllip(1050)) #3D location of vehicle frame origin
    ψ_nb::Float64 = 0.0 #geographic heading
    EAS::Float64 = 50.0 #equivalent airspeed
    γ_wb_n::Float64 = 0.0 #wind-relative flight path angle
    ψ_wb_dot::Float64 = 0.0 #WA-relative turn rate
    θ_wb_dot::Float64 = 0.0 #WA-relative pitch rate
    β_a::Float64 = 0.0 #sideslip angle measured in the aerodynamic reference frame
    x_fuel::Ranged{Float64, 0., 1.} = 0.5 #normalized fuel load
    mixture::Ranged{Float64, 0., 1.} = 0.5 #engine mixture control
    flaps::Ranged{Float64, 0., 1.} = 0.0 #flap setting
    payload::C172.PayloadU = C172.PayloadU(m_pilot = 75, m_copilot = 75, m_baggage = 50)
end


function Kinematics.Initializer(trim_state::TrimState,
                                trim_params::TrimParameters,
                                atm_data::AtmData)

    @unpack EAS, β_a, γ_wb_n, ψ_nb, ψ_wb_dot, θ_wb_dot, Ob = trim_params
    @unpack α_a, φ_nb = trim_state

    TAS = Air.EAS2TAS(EAS; ρ = ISAData(Ob, atm_data).ρ)
    v_wb_a = Air.get_velocity_vector(TAS, α_a, β_a)
    v_wb_b = C172.f_ba.q(v_wb_a) #wind-relative aircraft velocity, body frame

    θ_nb = AircraftBase.θ_constraint(; v_wb_b, γ_wb_n, φ_nb)
    e_nb = REuler(ψ_nb, θ_nb, φ_nb)
    q_nb = RQuat(e_nb)

    e_wb = e_nb #initialize WA arbitrarily to NED
    ė_wb = SVector(ψ_wb_dot, θ_wb_dot, 0.0)
    ω_wb_b = Attitude.ω(e_wb, ė_wb)

    loc = NVector(Ob)
    h = HEllip(Ob)

    v_wb_n = q_nb(v_wb_b) #wind-relative aircraft velocity, NED frame
    v_ew_n = atm_data.v_ew_n
    v_eb_n = v_ew_n + v_wb_n

    Kinematics.Initializer(; q_nb, loc, h, ω_wb_b, v_eb_n, Δx = 0.0, Δy = 0.0)

end


function cost(vehicle::System{<:C172.Vehicle})

    @unpack ẋ, y = vehicle

    v_nd_dot = SVector{3}(ẋ.dynamics.v_eb_b) / norm(y.kinematics.data.v_eb_b)
    ω_dot = SVector{3}(ẋ.dynamics.ω_eb_b) #ω should already of order 1
    n_eng_dot = ẋ.components.pwp.engine.ω / vehicle.components.pwp.engine.constants.ω_rated

    sum(v_nd_dot.^2) + sum(ω_dot.^2) + n_eng_dot^2

end

function get_f_target(vehicle::System{<:C172.Vehicle},
                      trim_params::TrimParameters)

    let vehicle = vehicle, trim_params = trim_params
        function (x::TrimState)
            AircraftBase.assign!(vehicle, trim_params, x)
            return cost(vehicle)
        end
    end

end

function Systems.init!(vehicle::System{<:C172.Vehicle}, trim_params::TrimParameters)

    trim_state = TrimState() #could provide initial condition as an optional input

    f_target = get_f_target(vehicle, trim_params)

    #wrapper with the interface expected by NLopt
    f_opt(x::Vector{Float64}, ::Vector{Float64}) = f_target(TrimState(x))

    n = length(trim_state)
    x0 = zeros(n); lower_bounds = similar(x0); upper_bounds = similar(x0); initial_step = similar(x0)

    x0[:] .= trim_state

    lower_bounds[:] .= TrimState(
        α_a = -π/12,
        φ_nb = -π/3,
        n_eng = 0.4,
        throttle = 0,
        aileron = -1,
        elevator = -1,
        rudder = -1)

    upper_bounds[:] .= TrimState(
        α_a = vehicle.components.aero.constants.α_stall[2], #critical AoA is 0.28 < 0.36
        φ_nb = π/3,
        n_eng = 1.1,
        throttle = 1,
        aileron = 1,
        elevator = 1,
        rudder = 1)

    initial_step[:] .= 0.05 #safe value for all optimization variables

    #any of these three algorithms works
    # opt = Opt(:LN_NELDERMEAD, length(x0))
    opt = Opt(:LN_BOBYQA, length(x0))
    # opt = Opt(:GN_CRS2_LM, length(x0))
    opt.min_objective = f_opt
    opt.maxeval = 100000
    opt.stopval = 1e-16
    opt.lower_bounds = lower_bounds
    opt.upper_bounds = upper_bounds
    opt.initial_step = initial_step

    # @btime optimize($opt, $x0)

    (minf, minx, exit_flag) = optimize(opt, x0)

    success = (exit_flag === :STOPVAL_REACHED)
    if !success
        @warn("Trimming optimization failed with exit_flag $exit_flag")
    end
    trim_state_opt = TrimState(minx)
    AircraftBase.assign!(vehicle, trim_params, trim_state_opt)
    return (success = success, trim_state = trim_state_opt)

end

function AircraftBase.trim!(ac::System{<:AircraftBase.Aircraft{<:C172.Vehicle}},
                            trim_params::TrimParameters = TrimParameters())
    Systems.init!(ac, trim_params)
end

function AircraftBase.linearize!(ac::System{<:AircraftBase.Aircraft{<:C172.Vehicle}},
                                trim_params::TrimParameters = TrimParameters())
    linearize!(ac.vehicle, trim_params)
end

################################################################################
############################### C172 Variants ##################################

include(normpath("c172r/c172r.jl")); @reexport using .C172R
include(normpath("c172rpa/c172rpa.jl")); @reexport using .C172RPA

end