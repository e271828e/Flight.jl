module C172Sv0

using Flight.FlightCore
using Flight.FlightLib

using ...C172
using ..C172S

export Cessna172Sv0

################################################################################
############################## Cessna172Sv0 ####################################

const Cessna172Sv0{K} = Cessna172S{K, NoAvionics} where { K <: AbstractKinematicDescriptor}

function Cessna172Sv0(kinematics = WA())
    AircraftBase.Aircraft(C172S.Vehicle(kinematics), NoAvionics())
end

############################ Joystick Mappings #################################

elevator_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
aileron_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
rudder_curve(x) = exp_axis_curve(x, strength = 1.5, deadzone = 0.05)
brake_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)

function IODevices.assign_input!(mdl::Model{<:Cessna172Sv0},
                                ::GenericInputMapping, data::T16000MData)

    u = mdl.vehicle.systems.act.u

    (; axes, buttons, hat) = data

    u.throttle = axes.throttle
    u.aileron = axes.stick_x |> aileron_curve
    u.elevator = axes.stick_y |> elevator_curve
    u.rudder = axes.stick_z |> rudder_curve

    u.brake_left = is_pressed(buttons.button_1)
    u.brake_right = is_pressed(buttons.button_1)

    u.aileron_offset -= 2e-4 * was_released(hat.left)
    u.aileron_offset += 2e-4 * was_released(hat.right)
    u.elevator_offset += 2e-4 * was_released(hat.down)
    u.elevator_offset -= 2e-4 * was_released(hat.up)

    u.flaps += 0.3333 * was_released(buttons.button_3)
    u.flaps -= 0.3333 * was_released(buttons.button_2)

end


end #module
