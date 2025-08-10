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
@sm_periodic C172Xv1.Avionics

function AircraftBase.assign!(systems::Model{<:C172X.Systems},
                          avionics::Model{<:C172Xv1.Avionics})
    AircraftBase.assign!(systems, avionics.ctl)
end

################################# Trimming #####################################

function Modeling.init!(avionics::Model{<:C172Xv1.Avionics},
                            vehicle::Model{<:C172X.Vehicle})

    Control.reset!(avionics.ctl)
    Modeling.init!(avionics.ctl, vehicle)
    update_output!(avionics)

end

################################## GUI #########################################

function GUI.draw!(avionics::Model{<:C172Xv1.Avionics},
                    vehicle::Model{<:C172X.Vehicle},
                    p_open::Ref{Bool} = Ref(true),
                    label::String = "Cessna172Xv1 Avionics")

    CImGui.Begin(label, p_open)

    @cstatic c_ctl=false begin
        @c CImGui.Checkbox("Flight Control", &c_ctl)
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


function IODevices.assign_input!(mdl::Model{<:Cessna172Xv1},
                           ::GenericInputMapping, data::T16000MData)

    u = mdl.avionics.ctl.u

    p_sf = 0.5 #roll rate sensitivity
    q_sf = 0.5 #pitch rate sensitivity

    @unpack axes, buttons, hat = data

    throttle_axis = axes.throttle
    roll_axis = roll_curve(axes.stick_x)
    pitch_axis =  pitch_curve(axes.stick_y)
    yaw_axis = yaw_curve(axes.stick_z)

    u.throttle_axis = throttle_axis
    u.aileron_axis = roll_axis
    u.elevator_axis = pitch_axis
    u.rudder_axis = yaw_axis

    u.p_ref = p_sf * roll_axis
    u.q_ref = q_sf * pitch_axis

    u.brake_left = is_pressed(buttons.button_1)
    u.brake_right = is_pressed(buttons.button_1)

    u.aileron_offset -= 5e-3 * was_released(hat.left)
    u.aileron_offset += 5e-3 * was_released(hat.right)
    u.elevator_offset += 5e-3 * was_released(hat.down)
    u.elevator_offset -= 5e-3 * was_released(hat.up)

    u.flaps += 0.3333 * was_released(buttons.button_3)
    u.flaps -= 0.3333 * was_released(buttons.button_2)

end

end #module