module C172Xv1

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays
using StructTypes

using Flight.FlightCore
using Flight.FlightLib

using ...C172
using ..C172X
using ..C172X.C172XControl: Controller, ControllerY

export Cessna172Xv1


################################################################################
############################### Avionics #######################################

@kwdef struct Avionics{C} <: AbstractAvionics
    ctl::C = Subsampled(Controller(), 2)
end

############################### Update methods #################################

@no_ode C172Xv1.Avionics
@no_step C172Xv1.Avionics
@ss_disc C172Xv1.Avionics

function AircraftBase.assign!(components::System{<:C172X.Components},
                          avionics::System{<:C172Xv1.Avionics})
    AircraftBase.assign!(components, avionics.ctl)
end

################################# Trimming #####################################

function Systems.init!(avionics::System{<:C172Xv1.Avionics},
                            vehicle::System{<:C172X.Vehicle})

    Systems.reset!(avionics)
    Systems.init!(avionics.ctl, vehicle)
    update_output!(avionics)

end

################################## GUI #########################################

function GUI.draw!(avionics::System{<:C172Xv1.Avionics},
                    vehicle::System{<:C172X.Vehicle},
                    p_open::Ref{Bool} = Ref(true),
                    label::String = "Cessna172Xv1 Avionics")

    CImGui.Begin(label, p_open)

    @cstatic c_ctl=false begin
        @c CImGui.Checkbox("Controller", &c_ctl)
        c_ctl && @c GUI.draw!(avionics.ctl, vehicle, &c_ctl)
    end

    CImGui.End()

end


################################################################################
############################# Cessna172Xv1 ###################################

const Cessna172Xv1{K, A} = Cessna172X{K, A} where {
    K <: AbstractKinematicDescriptor, A <: C172Xv1.Avionics}

function Cessna172Xv1(kinematics = WA())
    AircraftBase.Aircraft(C172X.Vehicle(kinematics), C172Xv1.Avionics())
end


############################ Joystick Mappings #################################

pitch_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
roll_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
yaw_curve(x) = exp_axis_curve(x, strength = 1.5, deadzone = 0.05)
brake_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)


function Systems.assign_input!(sys::System{<:Cessna172Xv1},
                           data::XBoxControllerData, ::IOMapping)

    u = sys.avionics.ctl.u

    p_sf = 0.5 #roll rate sensitivity
    q_sf = 0.5 #pitch rate sensitivity

    roll_input = get_axis_value(data, :right_stick_x) |> roll_curve
    pitch_input = get_axis_value(data, :right_stick_y) |> pitch_curve
    yaw_input = get_axis_value(data, :left_stick_x) |> yaw_curve

    u.throttle_sp_input += 0.1 * was_released(data, :button_Y)
    u.throttle_sp_input -= 0.1 * was_released(data, :button_A)

    u.aileron_sp_input = roll_input
    u.elevator_sp_input = pitch_input
    u.rudder_sp_input = yaw_input

    u.p_sp = p_sf * roll_input
    u.q_sp = q_sf * pitch_input

    u.steering = get_axis_value(data, :left_stick_x) |> yaw_curve
    u.brake_left = get_axis_value(data, :left_trigger) |> brake_curve
    u.brake_right = get_axis_value(data, :right_trigger) |> brake_curve

    u.aileron_sp_offset -= 1e-3 * was_released(data, :dpad_left)
    u.aileron_sp_offset += 1e-3 * was_released(data, :dpad_right)
    u.elevator_sp_offset += 1e-3 * was_released(data, :dpad_down)
    u.elevator_sp_offset -= 1e-3 * was_released(data, :dpad_up)

    u.flaps += 0.3333 * was_released(data, :right_bumper)
    u.flaps -= 0.3333 * was_released(data, :left_bumper)

end

function Systems.assign_input!(sys::System{<:Cessna172Xv1},
                           data::T16000MData, ::IOMapping)

    u = sys.avionics.ctl.u

    p_sf = 0.5 #roll rate sensitivity
    q_sf = 0.5 #pitch rate sensitivity

    throttle_input = get_axis_value(data, :throttle)
    roll_input = get_axis_value(data, :stick_x) |> roll_curve
    pitch_input = get_axis_value(data, :stick_y) |> pitch_curve
    yaw_input = get_axis_value(data, :stick_z) |> yaw_curve

    u.throttle_sp_input = throttle_input
    u.aileron_sp_input = roll_input
    u.elevator_sp_input = pitch_input
    u.rudder_sp_input = yaw_input

    u.p_sp = p_sf * roll_input
    u.q_sp = q_sf * pitch_input

    u.steering = get_axis_value(data, :stick_z) |> yaw_curve
    u.brake_left = is_pressed(data, :button_1)
    u.brake_right = is_pressed(data, :button_1)

    u.aileron_sp_offset -= 1e-3 * is_pressed(data, :hat_left)
    u.aileron_sp_offset += 1e-3 * is_pressed(data, :hat_right)
    u.elevator_sp_offset += 1e-3 * is_pressed(data, :hat_down)
    u.elevator_sp_offset -= 1e-3 * is_pressed(data, :hat_up)

    u.flaps += 0.3333 * was_released(data, :button_3)
    u.flaps -= 0.3333 * was_released(data, :button_2)

end

function Systems.assign_input!(sys::System{<:Cessna172Xv1},
                           data::GladiatorNXTEvoData, ::IOMapping)

    u = sys.avionics.ctl.u

    p_sf = 0.5 #roll rate sensitivity
    q_sf = 0.5 #pitch rate sensitivity

    throttle_input = get_axis_value(data, :throttle)
    roll_input = get_axis_value(data, :stick_x) |> roll_curve
    pitch_input = get_axis_value(data, :stick_y) |> pitch_curve
    yaw_input = get_axis_value(data, :stick_z) |> yaw_curve

    u.throttle_sp_input = throttle_input
    u.aileron_sp_input = roll_input
    u.elevator_sp_input = pitch_input
    u.rudder_sp_input = yaw_input

    u.p_sp = p_sf * roll_input
    u.q_sp = q_sf * pitch_input

    u.steering = get_axis_value(data, :stick_z) |> yaw_curve
    u.brake_left = is_pressed(data, :red_trigger_half)
    u.brake_right = is_pressed(data, :red_trigger_half)

    u.aileron_sp_offset -= 1e-3 * is_pressed(data, :A3_left)
    u.aileron_sp_offset += 1e-3 * is_pressed(data, :A3_right)
    u.elevator_sp_offset += 1e-3 * is_pressed(data, :A3_down)
    u.elevator_sp_offset -= 1e-3 * is_pressed(data, :A3_up)

    if is_pressed(data, :A3_press)
        u.aileron_sp_offset = 0
        u.elevator_sp_offset = 0
    end

    u.flaps += 0.3333 * was_released(data, :switch_down)
    u.flaps -= 0.3333 * was_released(data, :switch_up)

end


end #module