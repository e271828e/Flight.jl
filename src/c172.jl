module C172

using LinearAlgebra
using StaticArrays
using ComponentArrays
using UnPack
using Unitful

using Flight.Modeling
using Flight.Plotting
using Flight.Misc
using Flight.Attitude
using Flight.Terrain
using Flight.Airdata
using Flight.Kinematics
using Flight.Dynamics
using Flight.Components
using Flight.Propulsion
using Flight.Aerodynamics
using Flight.LandingGear
using Flight.Aircraft: AircraftBase, AbstractAircraftID
using Flight.Input

import Flight.Modeling: init_x0, init_y0, init_u0, init_d0, f_cont!, f_disc!
import Flight.Plotting: plots
import Flight.Components: MassTrait, WrenchTrait, AngularMomentumTrait, get_wr_b, get_mass_properties
import Flight.Aircraft: assign_joystick_inputs!

export C172Aircraft

struct ID <: AbstractAircraftID end


############################## Powerplant ################################

struct Pwp <: SystemGroupDescriptor
    left::EThruster
    right::EThruster
end

WrenchTrait(::System{<:Pwp}) = HasWrench()
AngularMomentumTrait(::System{<:Pwp}) = HasAngularMomentum()

function Pwp()

    prop = SimpleProp(kF = 4e-3, J = 0.25)

    left = EThruster(propeller = prop, motor = ElectricMotor(α = CW))
    right = EThruster(propeller = prop, motor = ElectricMotor(α = CCW))

    Pwp(left, right)

end

############################ LandingGear ############################

struct Ldg{L <: LandingGearUnit, R <: LandingGearUnit,
    C <: LandingGearUnit} <: SystemGroupDescriptor
    left::L
    right::R
    center::C
end

WrenchTrait(::System{<:Ldg}) = HasWrench()
AngularMomentumTrait(::System{<:Ldg}) = HasNoAngularMomentum()

function Ldg()

    mlg_damper = SimpleDamper(k_s = 25000, k_d_ext = 1000, k_d_cmp = 1000)
    nlg_damper = SimpleDamper(k_s = 25000, k_d_ext = 1000, k_d_cmp = 1000)

    left = LandingGearUnit(
        strut = Strut(
            t_bs = FrameTransform(r = [-1, -1.25, 1], q = RQuat() ),
            l_0 = 0.0,
            damper = mlg_damper),
        braking = DirectBraking())

    right = LandingGearUnit(
        strut = Strut(
            t_bs = FrameTransform(r = [-1, 1.25, 1], q = RQuat() ),
            l_0 = 0.0,
            damper = mlg_damper),
        braking = DirectBraking())

    center = LandingGearUnit(
        strut = Strut(
            t_bs = FrameTransform(r = [2, 0, 1] , q = RQuat()),
            l_0 = 0.0,
            damper = nlg_damper),
        steering = DirectSteering())

    Ldg(left, right, center)

end

###############################  Aerodynamics ###########################

#we need a generic aerodynamics system, which computes all the input values to
#the aerodynamic coefficient tables, including filtered alphadot and betadot,
#and total thrust (which can be retrieved by calling get_wr_b on the pwp model)
#its interface should be f_cont!(aero, air, pwp, ldg, kin, trn). This System's
#descriptor can be parameterized with a type AerodynamicsDataset that allows
#dispatching on a method get_aerodynamic_coefficients(::AeroDataset, alpha,
#alpha_dot, beta, beta_dot, etc), which returns C_X, C_Y, C_Z, C_l, C_m, etc. to
#the main Aerodynamics System. the AerodynamicsDataset and its associated
#get_coefficients method has all the specifics required by each aircraft. the
#AerodynamicsDataset type itself can be empty (a dummy type) or have actual
#fields required for the computation

#should probably define a struct AerodynamicCoefficients with fields Cx, Cy, Cz,
#Cl, Cm, Cn to clearly define the interface for an AeroDataset. since the
#scaling factors for transforming that into component axes are always the same
#(dynamic pressure, etc) maybe we should keep that in the generic root
#Aerodynamics System. If it turns out the non-dimensionalization depends on the
#specific dataset, we should return the dimensional quantities instead.

#if the specific AeroDataset requires an intermediate expression in stability or
#wind axes, these can always be constructed alpha and beta by the standard
#methods

#nope... realistically, all computations should be probably grouped under a
#single Aerodynamics System, because the input arguments to the AeroDataset may
#change greatly from one aircraft to another

Base.@kwdef struct Aero <: AbstractAerodynamics
    S::Float64 = 23.23 #wing area
    b::Float64 = 14.63 #wingspan
    c::Float64 = 1.5875 #mean aerodynamic chord
    δe_max::Float64 = 20 |> deg2rad #maximum elevator deflection (rad)
    δa_max::Float64 = 20 |> deg2rad #maximum aileron deflection (rad)
    δr_max::Float64 = 20 |> deg2rad #maximum rudder deflection (rad)
    δf_max::Float64 = 30 |> deg2rad #maximum flap deflection (rad)
    τ::Float64 = 0.1 #time constant for airflow angle filtering
end

Base.@kwdef mutable struct AeroU
    δe::Bounded{Float64, -1, 1} = 0.0 #(+ pitch down)
    δa::Bounded{Float64, -1, 1} = 0.0 #(+ roll left)
    δr::Bounded{Float64, -1, 1} = 0.0 #(+ yaw left)
    δf::Bounded{Float64, 0, 1} = 0.0 #(+ flap down)
end

Base.@kwdef struct AeroY
    α::Float64 = 0.0
    α_filt::Float64 = 0.0
    α_filt_dot::Float64 = 0.0
    β::Float64 = 0.0
    β_filt::Float64 = 0.0
    β_filt_dot::Float64 = 0.0
    δe::Float64 = 0.0
    δa::Float64 = 0.0
    δr::Float64 = 0.0
    δf::Float64 = 0.0
    wr_b::Wrench = Wrench()
end

init_x0(::Aero) = ComponentVector(α_filt = 0.0, β_filt = 0.0) #filtered airflow angles
init_y0(::Aero) = AeroY()
init_u0(::Aero) = AeroU()


############################## Controls #################################

struct Controls <: SystemDescriptor end

Base.@kwdef mutable struct ControlsU
    throttle::Bounded{Float64, 0, 1} = 0.0
    yoke_x::Bounded{Float64, -1, 1} = 0.0 #ailerons (+ bank right)
    yoke_x_trim::Bounded{Float64, -1, 1} = 0.0 #ailerons (+ bank right)
    yoke_y::Bounded{Float64, -1, 1} = 0.0 #elevator (+ pitch up)
    yoke_y_trim::Bounded{Float64, -1, 1} = 0.0 #elevator (+ pitch up)
    pedals::Bounded{Float64, -1, 1} = 0.0 #rudder and nose wheel (+ yaw right)
    brake_left::Bounded{Float64, 0, 1} = 0.0 #[0, 1]
    brake_right::Bounded{Float64, 0, 1} = 0.0 #[0, 1]
    flaps::Bounded{Float64, 0, 1} = 0.0 #[0, 1]
end

#const is essential when declaring type aliases!
Base.@kwdef struct ControlsY
    throttle::Float64
    yoke_x::Float64
    yoke_x_trim::Float64
    yoke_y::Float64
    yoke_y_trim::Float64
    pedals::Float64
    brake_left::Float64
    brake_right::Float64
    flaps::Float64
end

init_u0(::Controls) = ControlsU()
init_y0(::Controls) = ControlsY(zeros(SVector{9})...)


############################## Airframe #################################

Base.@kwdef struct Airframe{ Aero <: SystemDescriptor, Pwp <: SystemDescriptor,
                        Ldg <: SystemDescriptor} <: SystemGroupDescriptor
    aero::Aero = Aero()
    pwp::Pwp = Pwp()
    ldg::Ldg = Ldg()
end

MassTrait(::System{<:Airframe}) = HasMass()
WrenchTrait(::System{<:Airframe}) = HasWrench()
AngularMomentumTrait(::System{<:Airframe}) = HasAngularMomentum()


################## Aero Update Functions #######################

function f_cont!(sys::System{Aero}, pwp::System{Pwp},
    air::AirData, kinematics::KinData, terrain::AbstractTerrain)

    #USING BEAVER AERODYNAMICS

    #in this aircraft, the aerodynamics' frame is the airframe itself (b), so
    #we can just use the airflow angles computed by the air data module for the
    #airframe axes. no need to recompute them on a local component frame

    #for near-zero TAS, the airflow angles are likely to chatter between 0, -π
    #and π. this introduces noise in airflow angle derivatives and general
    #unpleasantness. to fix it we fade them in from zero to a safe threshold. as
    #for V, it's only used for non-dimensionalization. we just need to avoid
    #dividing by zero. for TAS < TAS_min, dynamic pressure will be close to zero
    #and therefore forces and moments will vanish anyway.

    @unpack ẋ, x, u, params = sys
    @unpack α_filt, β_filt = x
    @unpack S, b, c, δe_max, δa_max, δr_max, δf_max, τ = params
    @unpack TAS, q, α_b, β_b = air
    ω_lb_b = kinematics.vel.ω_lb_b

    TAS_min = 1.0
    χ_TAS = min(TAS/TAS_min, 1.0) #linear fade-in for TAS < TAS_thr

    α = α_b * χ_TAS
    β = β_b * χ_TAS
    V = max(TAS, TAS_min)

    α_filt_dot = 1/τ * (α - α_filt)
    β_filt_dot = 1/τ * (β - β_filt)

    p_nd = ω_lb_b[1] * b / (2V) #non-dimensional roll rate
    q_nd = ω_lb_b[2] * c / V #non-dimensional pitch rate
    r_nd = ω_lb_b[3] * b / (2V) #non-dimensional yaw rate

    α_dot_nd = α_filt_dot * c / (2V) #not used
    β_dot_nd = β_filt_dot * b / (2V)

    # T = get_wr_b(pwp).F[1]
    # C_T = T / (dyn_p * S) #thrust coefficient, not used
    δe = Float64(u.δe) * δe_max
    δa = Float64(u.δa) * δa_max
    δr = Float64(u.δr) * δr_max
    δf = Float64(u.δf) * δf_max

    C_X, C_Y, C_Z, C_l, C_m, C_n = get_coefficients(; α, β, p_nd, q_nd, r_nd,
        δa, δr, δe, δf, α_dot_nd, β_dot_nd)

    F_x = C_X * q * S
    F_y = C_Y * q * S
    F_z = C_Z * q * S
    F_aero_b = SVector{3,Float64}(F_x, F_y, F_z)

    M_x = C_l * q * S * b
    M_y = C_m * q * S * c
    M_z = C_n * q * S * b
    M_aero_b = SVector{3,Float64}(M_x, M_y, M_z)

    wr_b = Wrench(F_aero_b, M_aero_b)

    # @show wr_b

    ẋ.α_filt = α_filt_dot
    ẋ.β_filt = β_filt_dot

    sys.y = AeroY(; α, α_filt, α_filt_dot, β, β_filt, β_filt_dot,
        δe, δa, δr, δf, wr_b)

end

function get_coefficients(; α, β, p_nd, q_nd, r_nd, δa, δr, δe, δf, α_dot_nd, β_dot_nd)

    #stall is modeled, so +/-1 seem like sensible limits for airflow angles
    α = clamp(α, -1.0, 1.0)
    β = clamp(β, -1.0, 1.0)

    α² = α^2; α³ = α^3; β² = β^2; β³ = β^3

    C_X = -0.03554 + 0.00292α + 5.459α² - 5.162α³  - 0.6748q_nd + 0.03412δr - 0.09447δf + 1.106(δf*α)
    C_Y = -0.002226 - 0.7678β - 0.1240p_nd + 0.3666r_nd -0.02956δa + 0.1158δr + 0.5238(δr*α) - 0.16β_dot_nd
    C_Z = -0.05504 - 5.578α + 3.442α³ - 2.988q_nd - 0.3980δe - 1.377δf - 1.261(δf*α) - 15.93(δe*β²)
    C_l = 5.91e-4 - 0.0618β - 0.5045p_nd + 0.1695r_nd - 0.09917δa + 6.934e-3δr - 0.08269(δa*α)
    C_m = 0.09448 - 0.6028α - 2.14α² -15.56q_nd - 1.921δe + 0.6921β² - 0.3118r_nd + 0.4072δf
    C_n = -3.117e-3 + 6.719e-3β + 0.1373β³ - 0.1585p_nd + 0.1595q_nd - 0.1112r_nd - 3.872e-3δa - 0.08265δr

    return (C_X, C_Y, C_Z, C_m, C_l, C_n)

end

function f_disc!(::System{<:Aero})
    return false
end

# get_wr_b(sys::System{Aero}) = sys.y.wr_b
get_wr_b(sys::System{Aero}) = Wrench()


######################## Controls Update Functions ###########################

function f_cont!(ctl::System{<:Controls}, ::System{<:Airframe},
                ::KinData, ::AirData, ::AbstractTerrain)

    #here, controls do nothing but update their output state. for a more complex
    #aircraft a continuous state-space autopilot implementation could go here
    @unpack throttle, yoke_x, yoke_x_trim, yoke_y, yoke_y_trim,
            pedals, brake_left, brake_right, flaps = ctl.u

    return ControlsY(; throttle, yoke_x, yoke_x_trim, yoke_y, yoke_y_trim,
                            pedals, brake_left, brake_right, flaps)

end

function f_disc!(::System{<:Controls}, ::System{<:Airframe})
    #this currently does nothing, but it may be used to implement open loop or
    #closed loop control laws, predefined maneuvers, etc. by calling an
    #externally overloadable function
    return false
end


####################### Airframe Update Functions ###########################


function f_cont!(afr::System{<:Airframe}, ctl::System{<:Controls},
                kin::KinData, air::AirData, trn::AbstractTerrain)

    @unpack aero, pwp, ldg = afr.subsystems

    #could this go in the main aircraft f_cont!?
    assign_component_inputs!(afr, ctl)
    # f_cont!(srf, air) #update surface actuator continuous state & outputs
    f_cont!(ldg, kin, trn) #update landing gear continuous state & outputs
    f_cont!(pwp, kin, air) #update powerplant continuous state & outputs
    f_cont!(aero, pwp, air, kin, trn)
    # f_cont!(aero, air, kin, srf, trn) #requires previous srf update

    afr.y = (aero = aero.y, pwp = pwp.y, ldg = ldg.y)

end

function f_disc!(afr::System{<:Airframe}, ::System{<:Controls})
    #fall back to the default SystemGroup implementation, the f_disc! for the
    # components don't have to deal with the Controls
    return f_disc!(afr)
end

function get_mass_properties(::System{<:Airframe})

    #for an aircraft implementing a fuel system the current mass properties are
    #computed here by querying the fuel system for the contributions of the
    #different fuel tanks

    MassProperties(
        #upreferred(2650u"lb") |> ustrip,
        m = 1202.0197805,
        #upreferred.([948, 1346, 1967]u"m") |> ustrip |> diagm |> SMatrix{3,3,Float64},
        J_Ob_b = SA[948.0 0 0; 0 1346.0 0; 0 0 1967.0],
        r_ObG_b = zeros(SVector{3}))
end

#get_wr_b and get_hr_b use the fallback for SystemGroups, which will in turn
#call get_wr_b and get_hr_b on aero, pwp and ldg

function assign_component_inputs!(afr::System{<:Airframe},
    ctl::System{<:Controls})

    @unpack throttle, yoke_x, yoke_x_trim, yoke_y, yoke_y_trim,
            pedals, brake_left, brake_right, flaps = ctl.u
    @unpack aero, pwp, ldg = afr.subsystems

    #here, we interpret yoke_y as the offset with respect to the force-free yoke
    #y position. the absolute yoke y position is the sum of yoke_y and
    #yoke_y_trim, from which it is measured

    pwp.u.left.throttle = throttle
    pwp.u.right.throttle = throttle
    ldg.u.center.steering[] = steering_curve(pedals)
    ldg.u.left.braking[] = brake_curve(brake_left)
    ldg.u.right.braking[] = brake_curve(brake_right)
    aero.u.δe = -elevator_curve(yoke_y_trim + yoke_y) #+yoke_y and +yoke_y_trim are back and +δe is pitch down, need to invert it
    aero.u.δa = -aileron_curve(yoke_x_trim + yoke_x) #+yoke_x and +yoke_x_trim are right and +δa is roll left, need to invert it
    aero.u.δr = -rudder_curve(pedals) # +pedals is right and +δr is yaw left
    aero.u.δf = flaps # +flaps is flaps down and +δf is flaps down

    return nothing
end

elevator_curve(x) = exp_axis_curve(x, strength = 0.5, deadzone = 0.05)
aileron_curve(x) = exp_axis_curve(x, strength = 0.5, deadzone = 0.05)
rudder_curve(x) = exp_axis_curve(x, strength = 1.5, deadzone = 0.1)
steering_curve(x) = exp_axis_curve(x, strength = 1.5, deadzone = 0.1)
brake_curve(x) = exp_axis_curve(x, strength = 0, deadzone = 0.05)

function exp_axis_curve(x::Bounded{T}, args...; kwargs...) where {T}
    exp_axis_curve(T(x), args...; kwargs...)
end

function exp_axis_curve(x::Real; strength::Real = 0.0, deadzone::Real = 0.0)

    a = strength
    x0 = deadzone

    abs(x) <= 1 || throw(ArgumentError("Input to exponential curve must be within [-1, 1]"))
    (x0 >= 0 && x0 <= 1) || throw(ArgumentError("Exponential curve deadzone must be within [0, 1]"))

    if x > 0
        y = max(0, (x - x0)/(1 - x0)) * exp( a * (abs(x) -1) )
    else
        y = min(0, (x + x0)/(1 - x0)) * exp( a * (abs(x) -1) )
    end
end

######################## XBoxController Input Interface ########################

function assign_joystick_inputs!(ac::System{<:AircraftBase{ID}}, joystick::XBoxController)

    ac.u.yoke_x = get_axis_value(joystick, :right_analog_x)
    ac.u.yoke_x_trim += 0.01 * is_released(joystick, :dpad_right)
    ac.u.yoke_x_trim -= 0.01 * is_released(joystick, :dpad_left)

    ac.u.yoke_y = get_axis_value(joystick, :right_analog_y)
    ac.u.yoke_y_trim -= 0.01 * is_released(joystick, :dpad_up)
    ac.u.yoke_y_trim += 0.01 * is_released(joystick, :dpad_down)

    ac.u.pedals = get_axis_value(joystick, :left_analog_x)

    ac.u.brake_left = get_axis_value(joystick, :left_trigger)
    ac.u.brake_right = get_axis_value(joystick, :right_trigger)

    ac.u.throttle += 0.1 * is_released(joystick, :button_Y) #manifold pressure
    ac.u.throttle -= 0.1 * is_released(joystick, :button_A)

    # ac.u.propeller_speed += 0.1 * is_released(joystick, :button_X) #rpms
    # ac.u.propeller_speed -= 0.1 * is_released(joystick, :button_B)

    ac.u.flaps += 0.5 * is_released(joystick, :right_bumper)
    ac.u.flaps -= 0.5 * is_released(joystick, :left_bumper)

    # Y si quisiera landing gear up y down, podria usar option como
    #modifier

end


#Aircraft constructor override keyword inputs to customize

function C172Aircraft(; id = ID(), kin = KinLTF(), afm = Airframe(), ctl = Controls())
    AircraftBase( id, kin, afm, ctl)
end



end #module