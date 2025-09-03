module C172Zv0

using UnPack

using Flight.FlightCore
using Flight.FlightLib

using ...C172
using ..C172Z

export Cessna172Zv0

################################################################################
################################# Cessna172Zv0 ################################

const Cessna172Zv0{K} = Cessna172Z{K, NoAvionics} where { K <: AbstractKinematicDescriptor}

function Cessna172Zv0(kinematics = WA())
    AircraftBase.Aircraft(C172Z.Vehicle(kinematics), NoAvionics())
end

############################ Joystick Mappings #################################

pitch_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
roll_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
yaw_curve(x) = exp_axis_curve(x, strength = 1.5, deadzone = 0.05)
brake_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)

function IODevices.assign_input!(mdl::Model{<:Cessna172Zv0},
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