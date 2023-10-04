module C172FBWCAS

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays

using Flight.FlightCore.Systems
using Flight.FlightCore.GUI
using Flight.FlightCore.IODevices
using Flight.FlightCore.Joysticks
using Flight.FlightCore.Utils: Ranged, saturation, wrap_to_π

using Flight.FlightPhysics.Attitude
using Flight.FlightPhysics.Kinematics
using Flight.FlightPhysics.RigidBody
using Flight.FlightPhysics.Environment

using Flight.FlightComponents.Control
using Flight.FlightComponents.Piston
using Flight.FlightComponents.Aircraft
using Flight.FlightComponents.World
using Flight.FlightComponents.Control: PIDDiscreteY

using ...C172
using ..C172FBW

export Cessna172FBWCAS


################################################################################
########################## AbstractControlChannel ##############################

abstract type AbstractControlChannel <: SystemDefinition end

################################################################################
################################## ThrottleControl #############################

@enum ThrottleMode begin
    direct_throttle_mode = 0
    airspeed_throttle_mode = 1
end

@kwdef struct ThrottleControl <: AbstractControlChannel
    TAS_comp::PIDDiscrete{1} = PIDDiscrete{1}(k_p = 1.8, k_i = 0.6, k_d = 0.2, τ_d = 0.04)
end

@kwdef mutable struct ThrottleControlU
    mode::ThrottleMode = direct_throttle_mode #throttle control mode
    thr_dmd::Ranged{Float64, 0., 1.} = 0.0 #throttle actuation demand
    TAS_dmd::Float64 = 0.0 #airspeed demand
end

@kwdef struct ThrottleControlY
    mode::ThrottleMode = direct_throttle_mode
    thr_dmd::Ranged{Float64, 0., 1.} = 0.0
    TAS_dmd::Float64 = 0.0
    thr_cmd::Ranged{Float64, 0., 1.} = 0.0 #throttle actuation command
    thr_sat::Int64 = 0 #throttle saturation state
    TAS_comp::PIDDiscreteY{1} = PIDDiscreteY{1}()
end

Systems.init(::SystemU, ::ThrottleControl) = ThrottleControlU()
Systems.init(::SystemY, ::ThrottleControl) = ThrottleControlY()

function Control.reset!(sys::System{<:ThrottleControl})
    #set default inputs and states
    sys.u.mode = direct_throttle_mode
    sys.u.thr_dmd = 0
    sys.u.TAS_dmd = 0
    #reset compensators
    Control.reset!.(values(sys.subsystems))
    #set default outputs
    sys.y = ThrottleControlY()
end

function Systems.f_disc!(sys::System{ThrottleControl}, air::AirData, Δt::Real)

    @unpack mode, thr_dmd, TAS_dmd = sys.u
    @unpack TAS_comp = sys.subsystems

    TAS_comp.u.feedback .= air.TAS

    if mode === direct_throttle_mode
        thr_cmd = thr_dmd
    else #airspeed mode
        TAS_comp.u.setpoint .= TAS_dmd
        f_disc!(TAS_comp, Δt)
        thr_cmd = Ranged(TAS_comp.y.out[1], 0., 1.) #compensator sign inversion
    end

    thr_sat = saturation(thr_cmd)
    TAS_comp.u.sat_ext .= thr_sat #saturation must also be inverted!

    sys.y = ThrottleControlY(; mode, thr_dmd, TAS_dmd,
                            thr_cmd, thr_sat, TAS_comp = TAS_comp.y)

end

function GUI.draw(sys::System{<:ThrottleControl})

    @unpack TAS_comp = sys.subsystems
    @unpack mode, thr_dmd, TAS_dmd, thr_cmd, thr_sat = sys.y

    CImGui.Begin("Throttle Control")

    CImGui.Text("Mode: $mode")
    CImGui.Text(@sprintf("Throttle demand: %.3f", Float64(thr_dmd)))
    CImGui.Text(@sprintf("Airspeed demand: %.3f m/s", rad2deg(TAS_dmd)))

    CImGui.Text(@sprintf("Throttle command: %.3f", Float64(thr_cmd)))
    CImGui.Text("Throttle saturation: $r_sat")

    if CImGui.TreeNode("Airspeed Compensator")
        GUI.draw(TAS_comp)
        CImGui.TreePop()
    end

    CImGui.End()

end

################################################################################
############################### PitchControl ################################

@enum PitchMode begin
    direct_elevator_mode = 0
    pitch_rate_mode = 1
    pitch_angle_mode = 2
    climb_rate_mode = 3
end

############################### PitchRateCmp ##################################

@kwdef struct PitchRateCmp <: SystemDefinition
    c1::PIDDiscrete{1} = PIDDiscrete{1}(k_p = 0, k_i = 1, k_d = 0) #pure integrator
    c2::PIDDiscrete{1} = PIDDiscrete{1}(k_p = 4.8, k_i = 20, k_d = 0.2, τ_d = 0.04) #see notebook
end

#overrides the default NamedTuple built from subsystem u's
@kwdef mutable struct PitchRateCmpU
    setpoint::MVector{1,Float64} = zeros(MVector{1})
    feedback::MVector{1,Float64} = zeros(MVector{1})
    sat_ext::MVector{1,Int64} = zeros(MVector{1, Int64})
end

@kwdef struct PitchRateCmpY
    setpoint::SVector{1,Float64} = zeros(SVector{1})
    feedback::SVector{1,Float64} = zeros(SVector{1})
    sat_ext::SVector{1,Int64} = zeros(SVector{1, Int64})
    out::SVector{1,Float64} = zeros(SVector{1})
    c1::PIDDiscreteY{1} = PIDDiscreteY{1}()
    c2::PIDDiscreteY{1} = PIDDiscreteY{1}()
end

Systems.init(::SystemU, ::PitchRateCmp) = PitchRateCmpU()
Systems.init(::SystemY, ::PitchRateCmp) = PitchRateCmpY()

function Control.reset!(sys::System{<:PitchRateCmp})
    sys.u.setpoint .= 0
    sys.u.feedback .= 0
    sys.u.sat_ext .= 0
    Control.reset!.(values(sys.subsystems))
end

#we leave the compensators' outputs unbounded (the default at initialization),
#the integrators will halt only when required by sat_ext
function Systems.f_disc!(sys::System{PitchRateCmp}, Δt::Real)
    @unpack setpoint, feedback, sat_ext = sys.u
    @unpack c1, c2 = sys.subsystems

    c1.u.setpoint .= setpoint
    c1.u.feedback .= feedback
    c1.u.sat_ext .= sat_ext
    f_disc!(c1, Δt)

    c2.u.setpoint .= c1.y.out #connected to c1's output
    c2.u.feedback .= 0.0 #no feedback, just feedforward path
    c2.u.sat_ext .= sat_ext
    f_disc!(c2, Δt)

    out = c2.y.out

    sys.y = PitchRateCmpY(; setpoint, feedback, sat_ext, out, c1 = c1.y, c2 = c2.y)

end

function GUI.draw(pitch_rate_comp::System{<:PitchRateCmp})
    if CImGui.TreeNode("Integrator")
        GUI.draw(pitch_rate_comp.c1)
        CImGui.TreePop()
    end
    if CImGui.TreeNode("PID")
        GUI.draw(pitch_rate_comp.c2)
        CImGui.TreePop()
    end
end

################################################################################

@kwdef struct PitchControl <: AbstractControlChannel
    q_comp::PitchRateCmp = PitchRateCmp()
    θ_comp::PIDDiscrete{1} = PIDDiscrete{1}(k_p = 1.8, k_i = 0.6, k_d = 0.12, τ_d = 0.04) #replace design with pure pitch rate feedback
    c_comp::PIDDiscrete{1} = PIDDiscrete{1}() ##########################TO DO
end

#overrides the default NamedTuple built from subsystem u's
@kwdef mutable struct PitchControlU
    mode::PitchMode = direct_elevator_mode #pitch control mode
    e_dmd::Ranged{Float64, -1., 1.} = 0.0 #elevator actuation demand
    q_dmd::Float64 = 0.0 #pitch rate demand
    θ_dmd::Float64 = 0.0 #pitch angle demand
    c_dmd::Float64 = 0.0 #climb rate demand
end

@kwdef struct PitchControlY
    mode::PitchMode = direct_elevator_mode
    e_dmd::Ranged{Float64, -1., 1.} = 0.0
    q_dmd::Float64 = 0.0
    θ_dmd::Float64 = 0.0
    c_dmd::Float64 = 0.0
    e_cmd::Ranged{Float64, -1., 1.} = 0.0 #elevator actuation command
    e_sat::Int64 = 0 #elevator saturation state
    q_comp::PitchRateCmpY = PitchRateCmpY()
    θ_comp::PIDDiscreteY{1} = PIDDiscreteY{1}()
    c_comp::PIDDiscreteY{1} = PIDDiscreteY{1}()
end

Systems.init(::SystemU, ::PitchControl) = PitchControlU()
Systems.init(::SystemY, ::PitchControl) = PitchControlY()

function Control.reset!(sys::System{<:PitchControl})
    #set default inputs and states
    sys.u.mode = direct_elevator_mode
    sys.u.e_dmd = 0
    sys.u.q_dmd = 0
    sys.u.θ_dmd = 0
    sys.u.c_dmd = 0
    #reset compensators
    Control.reset!.(values(sys.subsystems))
    #set default outputs
    sys.y = PitchControlY()
end

function Systems.f_disc!(sys::System{PitchControl}, kin::KinematicData, Δt::Real)

    @unpack mode, e_dmd, q_dmd, θ_dmd, c_dmd = sys.u
    @unpack q_comp, θ_comp, c_comp = sys.subsystems

    c_comp.u.feedback .= -kin.v_eOb_n[3]
    θ_comp.u.feedback .= kin.e_nb.θ
    q_comp.u.feedback .= kin.ω_lb_b[2]

    if mode === direct_elevator_mode
        e_cmd = e_dmd
    else
        if mode === pitch_rate_mode
            q_comp.u.setpoint .= q_dmd
        else #pitch_angle, climb_rate
            if mode === pitch_angle_mode
                θ_comp.u.setpoint .= θ_dmd
            else #climb rate
                c_comp.u.setpoint .= c_dmd
                f_disc!(c_comp, Δt)
                θ_comp.u.setpoint .= c_comp.y.out[1]
            end
            f_disc!(θ_comp, Δt)
            q_comp.u.setpoint .= θ_comp.y.out[1]
        end
        f_disc!(q_comp, Δt)
        e_cmd = Ranged(q_comp.y.out[1], -1., 1.)
    end

    e_sat = saturation(e_cmd)

    #will take effect on the next call
    q_comp.u.sat_ext .= e_sat
    θ_comp.u.sat_ext .= e_sat
    c_comp.u.sat_ext .= e_sat

    sys.y = PitchControlY(; mode, e_dmd, q_dmd, θ_dmd, c_dmd,
                            e_cmd, e_sat, q_comp = q_comp.y,
                            θ_comp = θ_comp.y, c_comp = c_comp.y)

end


function GUI.draw(sys::System{<:PitchControl})

    @unpack q_comp, θ_comp, c_comp = sys.subsystems
    @unpack mode, e_dmd, q_dmd, θ_dmd, c_dmd, e_cmd, e_sat = sys.y

    CImGui.Begin("Pitch Control")

    CImGui.Text("Mode: $mode")
    CImGui.Text(@sprintf("Elevator demand: %.3f", Float64(e_dmd)))
    CImGui.Text(@sprintf("Pitch rate demand: %.3f deg/s", rad2deg(q_dmd)))
    CImGui.Text(@sprintf("Pitch angle demand: %.3f deg", rad2deg(θ_dmd)))
    CImGui.Text(@sprintf("Climb rate demand: %.3f m/s", c_dmd))

    CImGui.Text(@sprintf("Elevator command: %.3f", Float64(e_cmd)))
    CImGui.Text("Elevator saturation: $e_sat")

    if CImGui.TreeNode("Pitch Rate Compensator")
        GUI.draw(q_comp)
        CImGui.TreePop()
    end
    if CImGui.TreeNode("Pitch Angle Compensator")
        GUI.draw(θ_comp)
        CImGui.TreePop()
    end
    if CImGui.TreeNode("Climb Rate Compensator")
        GUI.draw(c_comp)
        CImGui.TreePop()
    end

    CImGui.End()

end


################################################################################
################################## RollControl #################################

@enum RollMode begin
    direct_aileron_mode = 0
    roll_rate_mode = 1
    bank_angle_mode = 2
    course_angle_mode = 3
end

#TO DO: redesign φ compensator with actual p feedback
@kwdef struct RollControl <: AbstractControlChannel
    p_comp::PIDDiscrete{1} = PIDDiscrete{1}(k_p = 0.6, k_i = 8.0, k_d = 0.04, τ_d = 0.04)
    φ_comp::PIDDiscrete{1} = PIDDiscrete{1}(k_p = 6.0, k_i = 1.5, k_d = 0.6, τ_d = 0.04)
    χ_comp::PIDDiscrete{1} = PIDDiscrete{1}() ################ TO DO
end

@kwdef mutable struct RollControlU
    mode::RollMode = direct_aileron_mode
    a_dmd::Ranged{Float64, -1., 1.} = 0.0 #aileron actuation demand
    p_dmd::Float64 = 0.0 #roll rate demand
    φ_dmd::Float64 = 0.0 #bank angle demand
    χ_dmd::Float64 = 0.0 #course angle demand
end

@kwdef struct RollControlY
    mode::RollMode = direct_aileron_mode
    a_dmd::Ranged{Float64, -1., 1.} = Ranged(0.0, -1.0, 1.0) #aileron actuation demand
    p_dmd::Float64 = 0.0
    φ_dmd::Float64 = 0.0
    χ_dmd::Float64 = 0.0
    a_cmd::Ranged{Float64, -1., 1.} = Ranged(0.0, -1.0, 1.0) #aileron actuation command
    a_sat::Int64 = 0 #aileron saturation state
    p_comp::PIDDiscreteY{1} = PIDDiscreteY{1}()
    φ_comp::PIDDiscreteY{1} = PIDDiscreteY{1}()
    χ_comp::PIDDiscreteY{1} = PIDDiscreteY{1}()
end

Systems.init(::SystemU, ::RollControl) = RollControlU()
Systems.init(::SystemY, ::RollControl) = RollControlY()

function Systems.init!(sys::System{RollControl})
    #initialize course angle compensator bank angle demand output limits
    sys.χ_comp.u.bound_lo .= -π/4
    sys.χ_comp.u.bound_hi .= π/4
end

function Control.reset!(sys::System{<:RollControl})
    #set default inputs and states
    sys.u.mode = direct_aileron_mode
    sys.u.a_dmd = 0
    sys.u.p_dmd = 0
    sys.u.φ_dmd = 0
    sys.u.χ_dmd = 0
    #reset compensators
    Control.reset!.(values(sys.subsystems))
    #set default outputs
    sys.y = RollControlY()
end

function Systems.f_disc!(sys::System{RollControl}, kin::KinematicData, Δt::Real)

    @unpack mode, a_dmd, p_dmd, φ_dmd, χ_dmd = sys.u
    @unpack p_comp, φ_comp, χ_comp = sys.subsystems

    #χ_err = χ_dmd - χ must be wrapped between -π and π. so instead of letting
    #the PID subtract them internally, we set the feedback to zero and pass
    #χ_err directly as the setpoint
    p_comp.u.feedback .= kin.ω_lb_b[1]
    φ_comp.u.feedback .= kin.e_nb.φ
    χ_comp.u.feedback .= 0

    if mode === direct_aileron_mode
        a_cmd = a_dmd
    else
        if mode === roll_rate_mode
            p_comp.u.setpoint .= p_dmd
        else
            if mode === bank_angle_mode
                φ_comp.u.setpoint .= φ_dmd
            else #course angle
                χ = Attitude.azimuth(kin.v_eOb_n)
                χ_err = wrap_to_π(χ_dmd - χ)
                χ_comp.u.setpoint .= χ_err #use error as setpoint
                f_disc!(χ_comp, Δt)
                φ_comp.u.setpoint .= χ_comp.y.out[1]
            end
            f_disc!(φ_comp, Δt)
            p_comp.u.setpoint .= φ_comp.y.out[1]
        end
        f_disc!(p_comp, Δt)
        a_cmd = Ranged(p_comp.y.out[1], -1., 1.)
    end

    a_sat = saturation(a_cmd)

    #will take effect on the next call
    p_comp.u.sat_ext .= a_sat
    φ_comp.u.sat_ext .= a_sat
    χ_comp.u.sat_ext .= a_sat

    sys.y = RollControlY(; mode, a_dmd, p_dmd, φ_dmd, χ_dmd,
                        a_cmd, a_sat, p_comp = p_comp.y,
                        φ_comp = φ_comp.y, χ_comp = χ_comp.y)

end


function GUI.draw(sys::System{<:RollControl})

    @unpack p_comp, φ_comp, χ_comp = sys.subsystems
    @unpack mode, a_dmd, p_dmd, φ_dmd, χ_dmd, a_cmd, a_sat = sys.y

    CImGui.Begin("Roll Control")

    CImGui.Text("Mode: $mode")
    CImGui.Text(@sprintf("Aileron demand: %.3f", Float64(a_dmd)))
    CImGui.Text(@sprintf("Roll rate demand: %.3f deg/s", rad2deg(p_dmd)))
    CImGui.Text(@sprintf("Bank angle demand: %.3f deg", rad2deg(φ_dmd)))
    CImGui.Text(@sprintf("Track angle demand: %.3f deg", rad2deg(χ_dmd)))

    CImGui.Text(@sprintf("Aileron command: %.3f", Float64(a_cmd)))
    CImGui.Text("Aileron saturation: $a_sat")

    if CImGui.TreeNode("Roll Rate Compensator")
        GUI.draw(p_comp)
        CImGui.TreePop()
    end
    if CImGui.TreeNode("Bank Angle Compensator")
        GUI.draw(φ_comp)
        CImGui.TreePop()
    end
    if CImGui.TreeNode("Track Angle Compensator")
        GUI.draw(χ_comp)
        CImGui.TreePop()
    end

    CImGui.End()

end


################################################################################
#################################### YawControl ################################

@enum YawMode begin
    direct_rudder_mode = 0
    sideslip_mode = 1
end

@kwdef struct YawControl <: AbstractControlChannel
    β_comp::PIDDiscrete{1} = PIDDiscrete{1}(k_p = 12, k_i = 30, k_d = 5, τ_d = 0.04)
end

@kwdef mutable struct YawControlU
    mode::YawMode = direct_rudder_mode #yaw control mode
    r_dmd::Ranged{Float64, -1., 1.} = 0.0 #aileron actuation demand
    β_dmd::Float64 = 0.0 #sideslip angle demand
end

@kwdef struct YawControlY
    mode::YawMode = direct_rudder_mode
    r_dmd::Ranged{Float64, -1., 1.} = 0.0
    β_dmd::Float64 = 0.0
    r_cmd::Ranged{Float64, -1., 1.} = 0.0 #rudder actuation command
    r_sat::Int64 = 0 #rudder saturation state
    β_comp::PIDDiscreteY{1} = PIDDiscreteY{1}()
end

Systems.init(::SystemU, ::YawControl) = YawControlU()
Systems.init(::SystemY, ::YawControl) = YawControlY()

function Control.reset!(sys::System{<:YawControl})
    #set default inputs and states
    sys.u.mode = direct_rudder_mode
    sys.u.r_dmd = 0
    sys.u.β_dmd = 0
    #reset compensators
    Control.reset!.(values(sys.subsystems))
    #set default outputs
    sys.y = YawControlY()
end

function Systems.f_disc!(sys::System{YawControl}, air::AirData, Δt::Real)

    @unpack mode, r_dmd, β_dmd = sys.u
    @unpack β_comp = sys.subsystems

    β_comp.u.feedback .= air.β_b

    if mode === direct_rudder_mode
        r_cmd = r_dmd
    else #sideslip mode
        β_comp.u.setpoint .= β_dmd
        f_disc!(β_comp, Δt)
        r_cmd = Ranged(-β_comp.y.out[1], -1., 1.) #compensator sign inversion
    end

    r_sat = saturation(r_cmd)
    β_comp.u.sat_ext .= -r_sat #saturation must also be inverted!

    sys.y = YawControlY(; mode, r_dmd, β_dmd, r_cmd, r_sat, β_comp = β_comp.y)

end


function GUI.draw(sys::System{<:YawControl})

    @unpack β_comp = sys.subsystems
    @unpack mode, r_dmd, β_dmd, r_cmd, r_sat = sys.y

    CImGui.Begin("Yaw Control")

    CImGui.Text("Mode: $mode")
    CImGui.Text(@sprintf("Rudder demand: %.3f", Float64(r_dmd)))
    CImGui.Text(@sprintf("Sideslip angle demand: %.3f deg", rad2deg(β_dmd)))

    CImGui.Text(@sprintf("Rudder command: %.3f", Float64(r_cmd)))
    CImGui.Text("Rudder saturation: $r_sat")

    if CImGui.TreeNode("Sideslip Angle Compensator")
        GUI.draw(β_comp)
        CImGui.TreePop()
    end

    CImGui.End()

end

################################################################################
########################## Longitudinal Control ################################

#with θ-based airspeed control there is no risk of stalling. θ will be reduced
#as required to hold the requested airspeed. when the aircraft cannot climb any
#further, it will simply remain in acquire mode at constant altitude.

#but what if we are descending at some given speed and then when we reach the
#target altitude and switch to altitude hold we cannot hold that altitude at
#maximum throttle above the stall speed? this is theoretically possible, but it
#should not happen in practice, because how would we have climbed above it in
#the first place?

# we need e_sat to halt integration in case elevator saturates after enabling
#acquire mode, but we don't need throttle saturation, since in acquire mode we
#are setting it manually, and in hold mode we have no control over v_dmd, which
#is what ultimately affects thr_cmd via autothrottle

#do we need to reset the compensators on each state change? airspeed2theta: if
#we are at altitude hold and h_dmd changes, we will switch to airspeed2theta. if
#airspeed2theta retains its previous states, what theta will it command by
#default? that which achieved the previous airspeed at the previous altitude
#altitude2climbrate: what climbrate will it command? in principle zero, because
#the last time we disabled it, we were in altitude_hold, and therefore altitude
#was already at the commanded value. what we should probably do is design these
#compensators to have a relatively soft response

@enum LonControlAutoState begin
    altitude_acquire = 0
    altitude_hold = 1
end

@enum AltitudeDatum begin
    ellipsoidal = 0
    orthometric = 1
end

@kwdef struct LonControlAuto <: AbstractControlChannel
    h_comp::PIDDiscrete{1} = PIDDiscrete{1}() #altitude to climb rate compensator ########### TO DO
    TAS_comp::PIDDiscrete{1} = PIDDiscrete{1}() #airspeed to θ compensator ############TO DO
end

@kwdef mutable struct LonControlAutoU
    h_dmd::Tuple{Float64, AltitudeDatum} = (0.0, ellipsoidal) #altitude demand
    TAS_dmd::Float64 = 0.0 #airspeed demand
end

@kwdef struct LonControlAutoY
    state::LonControlAutoState = altitude_acquire
    throttle_mode::ThrottleMode = direct_throttle_mode
    pitch_mode::PitchMode = pitch_angle_mode
    thr_dmd::Float64 = 0.0
    θ_dmd::Float64 = 0.0
    TAS_dmd::Float64 = 0.0
    c_dmd::Float64 = 0.0
    h_comp::PIDDiscreteY{1} = PIDDiscreteY{1}()
    TAS_comp::PIDDiscreteY{1} = PIDDiscreteY{1}()
end

Systems.init(::SystemU, ::LonControlAuto) = LonControlAutoU()
Systems.init(::SystemY, ::LonControlAuto) = LonControlAutoY()


function Systems.f_disc!(sys::System{LonControlAuto}, kin::KinematicData, air::AirData, Δt::Real)

    @unpack h_dmd, TAS_dmd = sys.u
    @unpack h_comp, TAS_comp = sys.subsystems

    h_threshold = 20 #within 20 m we switch to altitude_hold
    h = (h_dmd[2] === ellipsoidal) ? Float64(kin.h_e) : Float64(kin.h_o)
    state = abs(h_dmd[1] - h) > h_threshold ? altitude_acquire : altitude_hold

    if state === altitude_acquire

        throttle_mode = direct_throttle_mode
        pitch_mode = pitch_angle_mode

        TAS_comp.u.setpoint .= TAS_dmd
        TAS_comp.u.feedback .= air.TAS
        f_disc!(TAS_comp, Δt)
        print("Remember sign inversion in θ_dmd from TAS, DC gain will be negative")

        thr_dmd = h_dmd[1] > h ? 1.0 : 0.0 #full throttle to climb, idle to descend
        θ_dmd = TAS_comp.y.out[1]
        c_dmd = 0.0 #no effect

    else #altitude_hold

        throttle_mode = airspeed_throttle_mode
        pitch_mode = climb_rate_mode

        h_comp.u.setpoint .= h_dmd[1]
        h_comp.u.feedback .= h
        f_disc!(h_comp, Δt)

        thr_dmd = 0.0 #no effect
        θ_dmd = 0.0 #no effect
        c_dmd = h_comp.y.out[1]

    end

    #note: in altitude_hold mode, TAS_dmd just passes through to ThrottleControl
    sys.y = LonControlAutoY(; state, throttle_mode, pitch_mode,
                            thr_dmd, θ_dmd, TAS_dmd, c_dmd,
                            h_comp = h_comp.y, TAS_comp = TAS_comp.y)

end

##################################################################################
################################## Avionics ######################################

@enum FlightPhase begin
    phase_gnd = 0
    phase_air = 1
end

@enum LonMode begin
    lon_mode_semi = 0
    lon_mode_auto = 1
end

@enum LatMode begin
    lat_mode_semi = 0
end

@kwdef struct Avionics <: AbstractAvionics
    throttle_ctl::ThrottleControl = ThrottleControl()
    roll_ctl::RollControl = RollControl()
    pitch_ctl::PitchControl = PitchControl()
    yaw_ctl::YawControl = YawControl()
    lon_ctl::LonControlAuto = LonControlAuto()
end

@kwdef mutable struct PhysicalInputs
    eng_start::Bool = false
    eng_stop::Bool = false
    mixture::Ranged{Float64, 0., 1.} = 0.5
    throttle::Ranged{Float64, 0., 1.} = 0.0 #used in direct_throttle_mode
    roll_input::Ranged{Float64, -1., 1.} = 0.0 #used in aileron_mode and roll_rate_mode
    pitch_input::Ranged{Float64, -1., 1.} = 0.0 #used in direct_elevator_mode and pitch_rate_mode
    yaw_input::Ranged{Float64, -1., 1.} = 0.0 #used in rudder_mode and sideslip_mode
    aileron_cmd_offset::Ranged{Float64, -1., 1.} = 0.0
    elevator_cmd_offset::Ranged{Float64, -1., 1.} = 0.0
    rudder_cmd_offset::Ranged{Float64, -1., 1.} = 0.0
    flaps::Ranged{Float64, 0., 1.} = 0.0
    brake_left::Ranged{Float64, 0., 1.} = 0.0
    brake_right::Ranged{Float64, 0., 1.} = 0.0
end

#the β control loop tracks the β_dmd input. a positive β_dmd increment initially
#produces a negative yaw rate. the sign inversion β_dmd_sf keeps consistency in
#the perceived behaviour between direct rudder and β control modes.

@kwdef mutable struct DigitalInputs
    throttle_mode_sel::ThrottleMode = direct_throttle_mode #selected throttle channel mode
    roll_mode_sel::RollMode = direct_aileron_mode #selected roll channel mode
    pitch_mode_sel::PitchMode = direct_elevator_mode #selected pitch channel mode
    yaw_mode_sel::YawMode = direct_rudder_mode #selected yaw channel mode
    lon_mode_sel::LonMode = lon_mode_semi #selected longitudinal control mode
    lat_mode_sel::LatMode = lat_mode_semi #selected lateral control mode
    TAS_dmd::Float64 = 40.0 #airspeed demand
    θ_dmd::Float64 = 0.0 #pitch angle demand
    c_dmd::Float64 = 0.0 #climb rate demand
    φ_dmd::Float64 = 0.0 #bank angle demand
    χ_dmd::Float64 = 0.0 #course angle demand
    h_dmd::Tuple{Float64, AltitudeDatum} = (0.0, ellipsoidal) #altitude demand
    p_dmd_sf::Float64 = 0.2 #roll_input to p_dmd scale factor
    q_dmd_sf::Float64 = 0.2 #pitch_input to q_dmd scale factor
    β_dmd_sf::Float64 = -deg2rad(10) #yaw_input β_dmd scale factor, sign inverted
end

@kwdef struct AvionicsU
    physical::PhysicalInputs = PhysicalInputs()
    digital::DigitalInputs = DigitalInputs()
end

@kwdef struct AvionicsModing
    flight_phase::FlightPhase = phase_gnd
    throttle_mode::ThrottleMode = direct_throttle_mode
    roll_mode::RollMode = direct_aileron_mode
    pitch_mode::PitchMode = direct_elevator_mode
    yaw_mode::YawMode = direct_rudder_mode
    lon_mode::LonMode = lon_mode_semi
    lat_mode::LatMode = lat_mode_semi
end

@kwdef struct ActuationCommands
    eng_start::Bool = false
    eng_stop::Bool = false
    mixture::Ranged{Float64, 0., 1.} = 0.5
    throttle_cmd::Ranged{Float64, 0., 1.} = 0.0
    aileron_cmd::Ranged{Float64, -1., 1.} = 0.0
    elevator_cmd::Ranged{Float64, -1., 1.} = 0.0
    rudder_cmd::Ranged{Float64, -1., 1.} = 0.0
    aileron_cmd_offset::Ranged{Float64, -1., 1.} = 0.0
    elevator_cmd_offset::Ranged{Float64, -1., 1.} = 0.0
    rudder_cmd_offset::Ranged{Float64, -1., 1.} = 0.0
    flaps::Ranged{Float64, 0., 1.} = 0.0
    brake_left::Ranged{Float64, 0., 1.} = 0.0
    brake_right::Ranged{Float64, 0., 1.} = 0.0
end

@kwdef struct AvionicsY
    moding::AvionicsModing = AvionicsModing()
    actuation::ActuationCommands = ActuationCommands()
    throttle_ctl::ThrottleControlY = ThrottleControlY()
    roll_ctl::RollControlY = RollControlY()
    pitch_ctl::PitchControlY = PitchControlY()
    yaw_ctl::YawControlY = YawControlY()
    lon_ctl::LonControlAutoY = LonControlAutoY()
end

Systems.init(::SystemU, ::Avionics) = AvionicsU()
Systems.init(::SystemY, ::Avionics) = AvionicsY()
Systems.init(::SystemS, ::Avionics) = nothing #keep subsystems local


# ########################### Update Methods #####################################

function Systems.f_disc!(avionics::System{<:Avionics}, Δt::Real,
                        physics::System{<:AircraftPhysics},
                        ::System{<:AbstractEnvironment})

    @unpack eng_start, eng_stop, mixture, throttle,
            roll_input, pitch_input, yaw_input,
            aileron_cmd_offset, elevator_cmd_offset, rudder_cmd_offset,
            flaps, brake_left, brake_right = avionics.u.physical

    @unpack throttle_mode_sel, roll_mode_sel, pitch_mode_sel, yaw_mode_sel,
            lon_mode_sel, lat_mode_sel, TAS_dmd, θ_dmd, c_dmd, φ_dmd, χ_dmd, h_dmd,
            p_dmd_sf, q_dmd_sf, β_dmd_sf = avionics.u.digital

    @unpack throttle_ctl, roll_ctl, pitch_ctl, yaw_ctl, lon_ctl = avionics.subsystems

    @unpack airframe, kinematics, air = physics.y

    #direct surface and low level demands always come from physical inputs
    roll_ctl.u.a_dmd = roll_input
    pitch_ctl.u.e_dmd = pitch_input
    yaw_ctl.u.r_dmd = yaw_input

    roll_ctl.u.p_dmd = p_dmd_sf * Float64(roll_input)
    pitch_ctl.u.q_dmd = q_dmd_sf * Float64(pitch_input)
    yaw_ctl.u.β_dmd = β_dmd_sf * Float64(yaw_input)

    any_wow = any(SVector{3}(leg.strut.wow for leg in airframe.ldg))
    flight_phase = any_wow ? phase_gnd : phase_air

    if flight_phase === phase_gnd

        #these are irrelevant on ground, but must be defined in all paths
        lon_mode = lon_mode_semi
        lat_mode = lat_mode_semi

        throttle_ctl.u.mode = direct_throttle_mode
        roll_ctl.u.mode = direct_aileron_mode
        pitch_ctl.u.mode = direct_elevator_mode
        yaw_ctl.u.mode = direct_rudder_mode

    else #air

        lon_mode = lon_mode_sel
        lat_mode = lat_mode_sel

        if lon_mode === lon_mode_semi

            throttle_ctl.u.mode = throttle_mode_sel
            throttle_ctl.u.thr_dmd = throttle
            throttle_ctl.u.TAS_dmd = TAS_dmd

            pitch_ctl.u.mode = pitch_mode_sel
            pitch_ctl.u.θ_dmd = θ_dmd
            pitch_ctl.u.c_dmd = c_dmd

        else #lon_mode === lon_auto

            lon_ctl.u.h_dmd = h_dmd
            lon_ctl.u.TAS_dmd = TAS_dmd
            f_disc!(lon_ctl, kinematics, air, Δt)

            throttle_ctl.u.mode = lon_ctl.y.throttle_mode
            throttle_ctl.u.thr_dmd = lon_ctl.y.thr_dmd
            throttle_ctl.u.TAS_dmd = lon_ctl.y.TAS_dmd

            pitch_ctl.u.mode = lon_ctl.y.pitch_mode
            pitch_ctl.u.θ_dmd = lon_ctl.y.θ_dmd
            pitch_ctl.u.c_dmd = lon_ctl.y.c_dmd

        end

        if lat_mode === lat_mode_semi

            roll_ctl.u.mode = roll_mode_sel
            roll_ctl.u.φ_dmd = φ_dmd

            yaw_ctl.u.mode = yaw_mode_sel

        end

    end

    f_disc!(throttle_ctl, air, Δt)
    f_disc!(roll_ctl, kinematics, Δt)
    f_disc!(pitch_ctl, kinematics, Δt)
    f_disc!(yaw_ctl, air, Δt)

    throttle_cmd = throttle_ctl.y.thr_cmd
    aileron_cmd = roll_ctl.y.a_cmd
    elevator_cmd = pitch_ctl.y.e_cmd
    rudder_cmd = yaw_ctl.y.r_cmd

    moding = AvionicsModing(;
        flight_phase,
        throttle_mode = throttle_ctl.y.mode,
        roll_mode = roll_ctl.y.mode,
        pitch_mode = pitch_ctl.y.mode,
        yaw_mode = yaw_ctl.y.mode,
        lon_mode,
        lat_mode
      )

    #all signal  except for throttle, roll_input, pitch_input and yaw_input pass through
    actuation = ActuationCommands(; eng_start, eng_stop, mixture,
                throttle_cmd, aileron_cmd, elevator_cmd, rudder_cmd,
                aileron_cmd_offset, elevator_cmd_offset, rudder_cmd_offset,
                flaps, brake_left, brake_right)

    avionics.y = AvionicsY(; moding, actuation,
                            throttle_ctl = throttle_ctl.y,
                            roll_ctl = roll_ctl.y,
                            pitch_ctl = pitch_ctl.y,
                            yaw_ctl = yaw_ctl.y,
                            lon_ctl = lon_ctl.y
                            )

    return false

end

function Aircraft.map_controls!(airframe::System{<:C172.Airframe},
                                avionics::System{Avionics})

    @unpack eng_start, eng_stop, mixture, throttle_cmd, aileron_cmd,
            elevator_cmd, rudder_cmd, aileron_cmd_offset, elevator_cmd_offset,
            rudder_cmd_offset, flaps, brake_left, brake_right = avionics.y.actuation

    @pack! airframe.act.u = eng_start, eng_stop, mixture, throttle_cmd, aileron_cmd,
           elevator_cmd, rudder_cmd, aileron_cmd_offset, elevator_cmd_offset,
           rudder_cmd_offset, flaps, brake_left, brake_right

end


##################################### Tools ####################################

#makes Avionics inputs consistent with the trim solution obtained for the
#aircraft physics. this enables the trim condition to be preserved during
#simulation
function assign!(sys::System{<:Avionics},
                trim_params::C172FBW.TrimParameters,
                trim_state::C172FBW.TrimState)

    @unpack mixture, flaps = trim_params
    @unpack throttle, aileron, elevator, rudder = trim_state

    sys.u.physical.throttle = throttle
    sys.u.physical.roll_input = 0
    sys.u.physical.pitch_input = 0
    sys.u.physical.yaw_input = 0
    sys.u.physical.aileron_cmd_offset = aileron
    sys.u.physical.elevator_cmd_offset = elevator
    sys.u.physical.rudder_cmd_offset = rudder
    sys.u.physical.mixture = mixture
    sys.u.physical.flaps = flaps

    sys.u.digital.lon_mode_sel = C172FBWCAS.lon_mode_semi
    sys.u.digital.lat_mode_sel = C172FBWCAS.lat_mode_semi
    sys.u.digital.throttle_mode_sel = C172FBWCAS.direct_throttle_mode
    sys.u.digital.roll_mode_sel = C172FBWCAS.direct_aileron_mode
    sys.u.digital.pitch_mode_sel = C172FBWCAS.direct_elevator_mode
    sys.u.digital.yaw_mode_sel = C172FBWCAS.direct_rudder_mode

end


# # ################################## GUI #########################################

# # function mode_button_color(button_mode, selected_mode, active_mode)
# #     if active_mode === mode
# #         return HSV_green
# #     elseif selected_mode === mode
# #         return HSV_amber
# #     else
# #         return HSV_gray
# #     end
# # end

# function GUI.draw!(avionics::System{<:Avionics}, airframe::System{<:C172.Airframe},
#                     label::String = "Cessna 172 FBW CAS Avionics")

#     u = avionics.u

#     CImGui.Begin(label)

#     CImGui.PushItemWidth(-60)

#     if airframe.y.pwp.engine.state === Piston.eng_off
#         eng_start_HSV = HSV_gray
#     elseif airframe.y.pwp.engine.state === Piston.eng_starting
#         eng_start_HSV = HSV_amber
#     else
#         eng_start_HSV = HSV_green
#     end
#     dynamic_button("Engine Start", eng_start_HSV, 0.1, 0.2)
#     u.eng_start = CImGui.IsItemActive()
#     CImGui.SameLine()
#     dynamic_button("Engine Stop", HSV_gray, (HSV_gray[1], HSV_gray[2], HSV_gray[3] + 0.1), (0.0, 0.8, 0.8))
#     u.eng_stop = CImGui.IsItemActive()
#     CImGui.SameLine()
#     CImGui.Text(@sprintf("%.3f RPM", Piston.radpersec2RPM(airframe.y.pwp.engine.ω)))
#     CImGui.Separator()

#     if avionics.y.interface.CAS_state === CAS_disabled
#         CAS_HSV = HSV_gray
#     elseif avionics.y.interface.CAS_state === CAS_standby
#         CAS_HSV = HSV_amber
#     else
#         CAS_HSV = HSV_green
#     end
#     dynamic_button("CAS", CAS_HSV, 0.1, 0.1)
#     CImGui.IsItemActivated() ? u.CAS_enable = !u.CAS_enable : nothing
#     CImGui.SameLine()
#     CImGui.Text("Flight Phase: $(avionics.y.logic.flight_phase)")

#     @unpack roll_mode, pitch_mode, yaw_mode = avionics.y.interface

#     CImGui.Text("Roll Control Mode: "); CImGui.SameLine()
#     dynamic_button("Aileron", control_mode_HSV(aileron_mode, u.roll_mode_select, roll_mode), 0.1, 0.1); CImGui.SameLine()
#     CImGui.IsItemActive() ? u.roll_mode_select = aileron_mode : nothing
#     dynamic_button("Roll Rate", control_mode_HSV(roll_rate_mode, u.roll_mode_select, roll_mode), 0.1, 0.1); CImGui.SameLine()
#     CImGui.IsItemActive() ? u.roll_mode_select = roll_rate_mode : nothing
#     dynamic_button("Roll Angle", control_mode_HSV(roll_angle_mode, u.roll_mode_select, roll_mode), 0.1, 0.1); CImGui.SameLine()
#     CImGui.IsItemActive() ? u.roll_mode_select = roll_angle_mode : nothing

#     CImGui.Separator()
#     CImGui.Text("Pitch Control Mode: "); CImGui.SameLine()
#     dynamic_button("Elevator", control_mode_HSV(elevator_mode, u.pitch_mode_select, pitch_mode), 0.1, 0.1); CImGui.SameLine()
#     CImGui.IsItemActive() ? u.pitch_mode_select = elevator_mode : nothing
#     dynamic_button("Pitch Rate", control_mode_HSV(pitch_rate_mode, u.pitch_mode_select, pitch_mode), 0.1, 0.1); CImGui.SameLine()
#     CImGui.IsItemActive() ? u.pitch_mode_select = pitch_rate_mode : nothing
#     dynamic_button("Pitch Angle", control_mode_HSV(pitch_angle_mode, u.pitch_mode_select, pitch_mode), 0.1, 0.1); CImGui.SameLine()
#     CImGui.IsItemActive() ? u.pitch_mode_select = pitch_angle_mode : nothing

#     CImGui.Separator()
#     CImGui.Text("Yaw Control Mode: "); CImGui.SameLine()
#     dynamic_button("Rudder", control_mode_HSV(rudder_mode, u.yaw_mode_select, yaw_mode), 0.1, 0.1); CImGui.SameLine()
#     CImGui.IsItemActive() ? u.yaw_mode_select = rudder_mode : nothing
#     dynamic_button("Sideslip", control_mode_HSV(sideslip_mode, u.yaw_mode_select, yaw_mode), 0.1, 0.1); CImGui.SameLine()
#     CImGui.IsItemActive() ? u.yaw_mode_select = sideslip_mode : nothing

#     CImGui.Separator()
#     u.throttle = safe_slider("Throttle", u.throttle, "%.6f")
#     u.roll_input = safe_slider("Roll Input", u.roll_input, "%.6f")
#     u.pitch_input = safe_slider("Pitch Input", u.pitch_input, "%.6f")
#     u.yaw_input = safe_slider("Yaw Input", u.yaw_input, "%.6f")
#     u.aileron_offset = safe_input("Aileron Offset", u.aileron_offset, 0.001, 0.1, "%.6f")
#     u.elevator_offset = safe_input("Elevator Offset", u.elevator_offset, 0.001, 0.1, "%.6f")
#     u.rudder_offset = safe_input("Rudder Offset", u.rudder_offset, 0.001, 0.1, "%.6f")
#     u.flaps = safe_slider("Flaps", u.flaps, "%.6f")
#     u.mixture = safe_slider("Mixture", u.mixture, "%.6f")
#     u.brake_left = safe_slider("Left Brake", u.brake_left, "%.6f")
#     u.brake_right = safe_slider("Right Brake", u.brake_right, "%.6f")

#     #Internals
#     CImGui.Separator()
#     @unpack roll_control, pitch_control, yaw_control = avionics.subsystems

#     if CImGui.TreeNode("Internals")
#         show_roll_control = @cstatic check=false @c CImGui.Checkbox("Roll Control", &check); CImGui.SameLine()
#         show_roll_control && GUI.draw(roll_control)
#         show_pitch_control = @cstatic check=false @c CImGui.Checkbox("Pitch Control", &check); CImGui.SameLine()
#         show_pitch_control && GUI.draw(pitch_control)
#         show_yaw_control = @cstatic check=false @c CImGui.Checkbox("Yaw Control", &check); CImGui.SameLine()
#         show_yaw_control && GUI.draw(yaw_control)
#         CImGui.TreePop()
#     end


#     CImGui.PopItemWidth()

#     CImGui.End()

# end

################################################################################
############################# Cessna172RCAS #####################################

#Cessna172R with control augmenting Avionics
const Cessna172FBWCAS{K} = C172FBW.Template{K, Avionics} where {K}
Cessna172FBWCAS(kinematics = LTF()) = C172FBW.Template(kinematics, Avionics())


# # ############################ Joystick Mappings #################################

# function IODevices.assign!(sys::System{<:Cessna172RCAS}, joystick::Joystick,
#                            mapping::InputMapping)
#     IODevices.assign!(sys.avionics, joystick, mapping)
# end

# elevator_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
# aileron_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
# rudder_curve(x) = exp_axis_curve(x, strength = 1.5, deadzone = 0.05)
# brake_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)

# function IODevices.assign!(sys::System{Avionics},
#                            joystick::XBoxController,
#                            ::DefaultMapping)

#     u = sys.u

#     u.roll_input = get_axis_value(joystick, :right_analog_x) |> aileron_curve
#     u.pitch_input = get_axis_value(joystick, :right_analog_y) |> elevator_curve
#     u.yaw_input = get_axis_value(joystick, :left_analog_x) |> rudder_curve
#     u.brake_left = get_axis_value(joystick, :left_trigger) |> brake_curve
#     u.brake_right = get_axis_value(joystick, :right_trigger) |> brake_curve

#     u.aileron_offset -= 0.01 * was_released(joystick, :dpad_left)
#     u.aileron_offset += 0.01 * was_released(joystick, :dpad_right)
#     u.elevator_offset += 0.01 * was_released(joystick, :dpad_down)
#     u.elevator_offset -= 0.01 * was_released(joystick, :dpad_up)

#     u.throttle += 0.1 * was_released(joystick, :button_Y)
#     u.throttle -= 0.1 * was_released(joystick, :button_A)

#     u.flaps += 0.3333 * was_released(joystick, :right_bumper)
#     u.flaps -= 0.3333 * was_released(joystick, :left_bumper)

# end

# function IODevices.assign!(sys::System{Avionics},
#                            joystick::T16000M,
#                            ::DefaultMapping)

#     u = sys.u

#     u.throttle = get_axis_value(joystick, :throttle)
#     u.roll_input = get_axis_value(joystick, :stick_x) |> aileron_curve
#     u.pitch_input = get_axis_value(joystick, :stick_y) |> elevator_curve
#     u.yaw_input = get_axis_value(joystick, :stick_z) |> rudder_curve

#     u.brake_left = is_pressed(joystick, :button_1)
#     u.brake_right = is_pressed(joystick, :button_1)

#     u.aileron_offset -= 2e-4 * is_pressed(joystick, :hat_left)
#     u.aileron_offset += 2e-4 * is_pressed(joystick, :hat_right)
#     u.elevator_offset += 2e-4 * is_pressed(joystick, :hat_down)
#     u.elevator_offset -= 2e-4 * is_pressed(joystick, :hat_up)

#     u.flaps += 0.3333 * was_released(joystick, :button_3)
#     u.flaps -= 0.3333 * was_released(joystick, :button_2)

# end

# function IODevices.assign!(sys::System{Avionics},
#                            joystick::GladiatorNXTEvo,
#                            ::DefaultMapping)

#     u = sys.u

#     u.throttle = get_axis_value(joystick, :throttle)
#     u.roll_input = get_axis_value(joystick, :stick_x) |> aileron_curve
#     u.pitch_input = get_axis_value(joystick, :stick_y) |> elevator_curve
#     u.yaw_input = get_axis_value(joystick, :stick_z) |> rudder_curve

#     u.brake_left = is_pressed(joystick, :red_trigger_half)
#     u.brake_right = is_pressed(joystick, :red_trigger_half)

#     u.aileron_offset -= 2e-4 * is_pressed(joystick, :A3_left)
#     u.aileron_offset += 2e-4 * is_pressed(joystick, :A3_right)
#     u.elevator_offset += 2e-4 * is_pressed(joystick, :A3_down)
#     u.elevator_offset -= 2e-4 * is_pressed(joystick, :A3_up)

#     if is_pressed(joystick, :A3_press)
#         u.aileron_offset = 0
#         u.elevator_offset = 0
#     end

#     u.flaps += 0.3333 * was_released(joystick, :switch_down)
#     u.flaps -= 0.3333 * was_released(joystick, :switch_up)

# end



end #module