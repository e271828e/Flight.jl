module C172Yv2

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays
using StructTypes

using CImGui: Begin, End, PushItemWidth, PopItemWidth, AlignTextToFramePadding,
    Dummy, SameLine, NewLine, IsItemActive, IsItemActivated, Separator, Text,
    Checkbox, RadioButton, TableNextColumn, TableNextRow, BeginTable, EndTable

using Flight.FlightCore
using Flight.FlightLib

using ...C172
using ..C172Y: Cessna172Y, Vehicle, Systems
using ..C172Y.C172YControl: ControlLaws
using ..C172Y.C172YGuidance: GuidanceLaws

export Cessna172Yv2


################################################################################
############################# Cessna172Yv2 ###################################

@kwdef struct Avionics{C <: ControlLaws} <: AbstractAvionics
    ctl::C = ControlLaws()
    gdc::GuidanceLaws = GuidanceLaws()
end

@no_ode Avionics
@no_step Avionics

########################### Update Methods #####################################

function Modeling.f_periodic!(::NoScheduling, avionics::Model{<:Avionics},
                                vehicle::Model{<:Vehicle})

    @unpack ctl, gdc = avionics

    f_periodic!(gdc, ctl, vehicle)
    f_periodic!(ctl, vehicle)
    f_output!(avionics)

end

function AircraftBase.assign!(systems::Model{<:Systems}, avionics::Model{<:Avionics})
    AircraftBase.assign!(systems, avionics.ctl)
end


############################# Initialization ###################################

function Modeling.init!(avionics::Model{<:Avionics}, vehicle::Model{<:Vehicle})
    @unpack ctl, gdc = avionics
    Modeling.init!(ctl, vehicle)
    f_output!(avionics)
end


############################# ParametricTypes ##################################

const Cessna172Yv2{K, A} = Cessna172Y{K, A} where {
    K <: AbstractKinematicDescriptor, A <: Avionics}

function Cessna172Yv2(kinematics = WA())
    AircraftBase.Aircraft(Vehicle(kinematics), Avionics())
end


################################################################################

function GUI.draw!(avionics::Model{<:Avionics}, vehicle::Model{<:Vehicle},
                    p_open::Ref{Bool} = Ref(true))

    Begin("Cessna172Yv2 Avionics", p_open)

    @cstatic c_gdc=false c_ctl=false begin
        @c CImGui.Checkbox("Guidance", &c_gdc)
        c_gdc && @c GUI.draw!(avionics.gdc, vehicle, &c_gdc)
        @c CImGui.Checkbox("Control", &c_ctl)
        c_ctl && @c GUI.draw!(avionics.ctl, vehicle, &c_ctl)
    end

    End()

end

############################ Joystick Mappings #################################

pitch_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
roll_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
yaw_curve(x) = exp_axis_curve(x, strength = 1.5, deadzone = 0.05)
brake_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)

function IODevices.assign_input!(mdl::Model{<:Cessna172Yv2},
                           ::GenericInputMapping, data::T16000MData)

    @unpack act = mdl.vehicle.systems
    u_ctl = mdl.avionics.u

    q_sf = 0.5 #pitch rate sensitivity
    p_sf = 0.5 #roll rate sensitivity

    @unpack axes, buttons, hat = data

    throttle_axis = axes.throttle
    roll_axis = roll_curve(axes.stick_x)
    pitch_axis =  pitch_curve(axes.stick_y)
    yaw_axis = yaw_curve(axes.stick_z)

    u_ctl.lon.throttle_axis = throttle_axis
    u_ctl.lon.elevator_axis = pitch_axis
    u_ctl.lon.elevator_offset += 5e-3 * was_released(hat.down)
    u_ctl.lon.elevator_offset -= 5e-3 * was_released(hat.up)
    u_ctl.lon.q_ref = q_sf * pitch_axis

    u_ctl.lat.aileron_axis = roll_axis
    u_ctl.lat.rudder_axis = yaw_axis
    u_ctl.lat.aileron_offset -= 5e-3 * was_released(hat.left)
    u_ctl.lat.aileron_offset += 5e-3 * was_released(hat.right)
    u_ctl.lat.p_ref = p_sf * roll_axis

    act.brake_left.u = is_pressed(buttons.button_1)
    act.brake_right.u = is_pressed(buttons.button_1)

    act.flaps.u += 0.3333 * was_released(buttons.button_3)
    act.flaps.u -= 0.3333 * was_released(buttons.button_2)

end


function IODevices.assign_input!(mdl::Model{<:Cessna172Yv2},
                           ::GenericInputMapping, data::GladiatorNXTEvoData)

    @unpack act = mdl.vehicle.systems
    u_ctl = mdl.avionics.u

    q_sf = 0.5 #pitch rate sensitivity
    p_sf = 0.5 #roll rate sensitivity

    @unpack axes, buttons, hat = data

    throttle_axis = axes.throttle
    roll_axis = roll_curve(axes.stick_x)
    pitch_axis =  pitch_curve(axes.stick_y)
    yaw_axis = yaw_curve(axes.stick_z)

    u_ctl.lon.throttle_axis = throttle_axis
    u_ctl.lon.elevator_axis = pitch_axis
    u_ctl.lon.elevator_offset += 5e-3 * was_released(buttons.A4_down)
    u_ctl.lon.elevator_offset -= 5e-3 * was_released(buttons.A4_up)
    u_ctl.lon.q_ref = q_sf * pitch_axis

    u_ctl.lat.aileron_axis = roll_axis
    u_ctl.lat.rudder_axis = yaw_axis
    u_ctl.lat.aileron_offset -= 5e-3 * was_released(buttons.A4_left)
    u_ctl.lat.aileron_offset += 5e-3 * was_released(buttons.A4_right)
    u_ctl.lat.p_ref = p_sf * roll_axis

    act.brake_left.u[] = is_pressed(buttons.F2)
    act.brake_right.u[] = is_pressed(buttons.F3)

    act.flaps.u[] += 0.3333 * was_released(buttons.switch_down)
    act.flaps.u[] -= 0.3333 * was_released(buttons.switch_up)

end


end #module