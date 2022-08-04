module Input

using StaticArrays
using UnPack

export AbstractInputInterface

using GLFW: GLFW, Joystick as JoystickSlot, DeviceConfigEvent, JoystickPresent,
        GetJoystickAxes, GetJoystickButtons, GetJoystickName
export connected_joysticks, init_joysticks, joystick_callback
export get_axis_value, get_button_state, get_button_change, is_pressed, is_released
export XBoxController

abstract type AbstractInputInterface end

update!(input::AbstractInputInterface) = throw(MethodError(update!, (input,)))
assign!(target::Any, input::AbstractInputInterface) = throw(MethodError(assign!, (target, input)))
#to be extended by users for concrete target and interface

#############################################################
#################### AbstractJoystick #######################

abstract type AbstractJoystick <: AbstractInputInterface end

const connected_joysticks = Dict{JoystickSlot, AbstractJoystick}()

get_slot(joystick::AbstractJoystick) = throw(MethodError(get_slot, (joystick,)))

function update!(joystick::AbstractJoystick)
    slot = get_slot(joystick)
    JoystickPresent(slot) ? _update!(joystick) : error(
        "$(typeof(joystick)) not found at slot $slot"
    )
end

#to be overridden by concrete joystick subtypes
_update!(j::AbstractJoystick) = throw(MethodError(update!, (j,)))



##################### Joystick configuration ################

#adds any joysticks already connected when GLFW is first imported. from that
#moment on, any connections and disconnections can only be detected via
#callback, and then retrieved by calling GLFW.PollEvents()
function init_joysticks()

    for slot in instances(JoystickSlot)
        delete!(connected_joysticks, slot)
        if JoystickPresent(slot)
            add_joystick(slot)
        end
    end

    return connected_joysticks

end

function add_joystick(slot::JoystickSlot)

    if !JoystickPresent(slot)
        println("Error while adding device at slot $slot: not found")
        return
    end

    joystick_model = GetJoystickName(slot)

    if joystick_model === "Xbox Controller"
        joystick = XBoxController(slot)
        println("$joystick_model active at slot $slot")
    else
        error("Interface for joystick $joystick_model not implemented")
    end

    connected_joysticks[slot] = joystick

end

function remove_joystick(slot::JoystickSlot)

    if JoystickPresent(slot)
        joystick_model = GetJoystickName(slot)
        println("Can't remove $joystick_model from slot $slot: device still connected")
        return
    end

    delete!(connected_joysticks, slot)
    println("$slot is no longer active")
end

function joystick_callback(slot::JoystickSlot, event::DeviceConfigEvent)

    print("Joystick callback called")
    if event === GLFW.CONNECTED
        add_joystick(slot)
    elseif event === GLFW.DISCONNECTED
        remove_joystick(slot)
    end

end

#these only work when run from the module that uses the joystick(s)
# init_joysticks()
# GLFW.SetJoystickCallback(joystick_callback)
# GLFW.PollEvents() #check for newly connected joysticks

###################

@enum ButtonChange begin
    unchanged = 0
    pressed = 1
    released = 2
end
Base.convert(::Type{ButtonChange}, n::Integer) = ButtonChange(n)
Base.zero(::Type{ButtonChange}) = unchanged


##################################################################
######################### XBox Controller ########################

########################## Axes

const XBoxAxisLabels = (
    :left_analog_x, :left_analog_y, :right_analog_x, :right_analog_y,
    :left_trigger, :right_trigger
)
const default_axis_mapping = NamedTuple{XBoxAxisLabels}(1:6)

Base.@kwdef struct XBoxAxes
    mapping::NamedTuple{XBoxAxisLabels, NTuple{6,Int}} = default_axis_mapping
    data::MVector{6,Float32} = zeros(MVector{6, Float32})
end

Base.getindex(axes::XBoxAxes, s::Symbol) = axes.data[axes.mapping[s]]
Base.setindex!(axes::XBoxAxes, v, s::Symbol) = (axes.data[axes.mapping[s]] = v)

######################## Buttons

const XBoxButtonLabels = (
    :button_A, :button_B, :button_X, :button_Y, :left_bumper, :right_bumper,
    :view, :menu, :left_analog, :right_analog,
    :dpad_up, :dpad_right, :dpad_down, :dpad_left
)

const default_button_mapping = NamedTuple{XBoxButtonLabels}(1:14)

Base.@kwdef struct XBoxButtons
    mapping::NamedTuple{XBoxButtonLabels, NTuple{14,Int}} = default_button_mapping
    state::MVector{14,Bool} = zeros(MVector{14,Bool})
    change::MVector{14,ButtonChange} = zeros(MVector{14,ButtonChange})
end

####################### Main

Base.@kwdef struct XBoxController <: AbstractJoystick
    slot::JoystickSlot = GLFW.JOYSTICK_1
    axes::XBoxAxes = XBoxAxes()
    buttons::XBoxButtons = XBoxButtons()
end

XBoxController(slot::JoystickSlot) = XBoxController(; slot)
get_slot(joystick::XBoxController) = joystick.slot

function _update!(joystick::XBoxController)

    @unpack slot, axes, buttons = joystick

    axes.data .= GetJoystickAxes(slot)
    axes[:left_trigger] = 0.5*(1 + axes[:left_trigger])
    axes[:right_trigger] = 0.5*(1 + axes[:right_trigger])

    buttons_state_new = SVector{14,Bool}(GetJoystickButtons(slot))

    for i in 1:length(buttons.state)
        if !buttons.state[i] && buttons_state_new[i]
            buttons.change[i] = pressed
        elseif buttons.state[i] && !buttons_state_new[i]
            buttons.change[i] = released
        else
            buttons.change[i] = unchanged
        end
    end

    buttons.state .= buttons_state_new

end

get_axis_value(joy::XBoxController) = NamedTuple{XBoxAxisLabels}(Tuple(joy.axes.data))
get_button_state(joy::XBoxController) = NamedTuple{XBoxButtonLabels}(Tuple(joy.buttons.state))
get_button_change(joy::XBoxController) = NamedTuple{XBoxButtonLabels}(Tuple(joy.buttons.change))

function get_axis_value(joystick::XBoxController, s::Symbol)
    @unpack data, mapping = joystick.axes
    data[mapping[s]]
end

function get_button_state(joystick::XBoxController, s::Symbol)
    @unpack state, mapping = joystick.buttons
    state[mapping[s]]
end

function get_button_change(joystick::XBoxController, s::Symbol)
    @unpack change, mapping = joystick.buttons
    change[mapping[s]]
end

function is_pressed(joystick::XBoxController, s::Symbol)
    @unpack change, mapping = joystick.buttons
    change[mapping[s]] === pressed
end

function is_released(joystick::XBoxController, s::Symbol)
    @unpack change, mapping = joystick.buttons
    change[mapping[s]] === released
end


end #module