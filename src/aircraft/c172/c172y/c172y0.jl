module C172Yv0

using UnPack

using Flight.FlightCore
using Flight.FlightLib

using ...C172
using ..C172Y

export Cessna172Yv0

################################################################################
################################# Cessna172Yv0 ################################

const Cessna172Yv0{K} = Cessna172Y{K, NoAvionics} where { K <: AbstractKinematicDescriptor}

function Cessna172Yv0(kinematics = WA())
    AircraftBase.Aircraft(C172Y.Vehicle(kinematics), NoAvionics())
end

############################ Joystick Mappings #################################

pitch_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
roll_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
yaw_curve(x) = exp_axis_curve(x, strength = 1.5, deadzone = 0.05)
brake_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)

function IODevices.assign_input!(mdl::Model{<:Cessna172Yv0},
                           ::GenericInputMapping, data::T16000MData)

    #mixture not assigned
    @unpack throttle, mixture, aileron, elevator, rudder, flaps,
            brake_left, brake_right = mdl.vehicle.systems.act.submodels

    @unpack axes, buttons, hat = data

    throttle.u[] = axes.throttle
    aileron.u[] = axes.stick_x |> roll_curve
    elevator.u[] = axes.stick_y |> pitch_curve
    rudder.u[] = axes.stick_z |> yaw_curve
    brake_left.u[] = is_pressed(buttons.button_1)
    brake_right.u[] = is_pressed(buttons.button_1)

    flaps.u[] += 0.3333 * was_released(buttons.button_3)
    flaps.u[] -= 0.3333 * was_released(buttons.button_2)

end


end #module