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


function IODevices.assign_input!(sys::System{<:Cessna172Xv1},
                           ::GenericInputMapping, data::T16000MData)

    u = sys.avionics.ctl.u

    p_sf = 0.5 #roll rate sensitivity
    q_sf = 0.5 #pitch rate sensitivity

    @unpack axes, buttons, hat = data

    throttle_input = axes.throttle
    roll_input = roll_curve(axes.stick_x)
    pitch_input =  pitch_curve(axes.stick_y)
    yaw_input = yaw_curve(axes.stick_z)

    u.throttle_sp_input = throttle_input
    u.aileron_sp_input = roll_input
    u.elevator_sp_input = pitch_input
    u.rudder_sp_input = yaw_input

    u.p_sp = p_sf * roll_input
    u.q_sp = q_sf * pitch_input

    u.steering = yaw_input
    u.brake_left = is_pressed(buttons.button_1)
    u.brake_right = is_pressed(buttons.button_1)

    u.aileron_sp_offset -= 5e-3 * was_released(hat.left)
    u.aileron_sp_offset += 5e-3 * was_released(hat.right)
    u.elevator_sp_offset += 5e-3 * was_released(hat.down)
    u.elevator_sp_offset -= 5e-3 * was_released(hat.up)

    u.flaps += 0.3333 * was_released(buttons.button_3)
    u.flaps -= 0.3333 * was_released(buttons.button_2)

end

end #module