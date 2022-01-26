module Input

using StaticArrays
using UnPack

using GLFW: GLFW, Joystick as JoystickSlot, DeviceConfigEvent, JoystickPresent,
        GetJoystickAxes, GetJoystickButtons, GetJoystickName

export connected_joysticks, init_joysticks, update_joystick, joystick_callback
export get_axis_value, get_button_state, get_button_change, is_pressed, is_released
export XBoxController

abstract type AbstractJoystick end

assign_joystick_inputs!(args...) = throw(MethodError(assign_joystick_inputs!, args))

const connected_joysticks = Dict{JoystickSlot, AbstractJoystick}()

###################

@enum ButtonChange begin
    unchanged = 0
    pressed = 1
    released = 2
end
Base.convert(::Type{ButtonChange}, n::Integer) = ButtonChange(n)
Base.zero(::Type{ButtonChange}) = unchanged


######################### XBox Controller ########################

const XBoxAxisLabels = (
    :left_analog_x, :left_analog_y, :right_analog_x, :right_analog_y,
    :left_trigger, :right_trigger
)
const XBoxButtonLabels = (
    :button_A, :button_B, :button_X, :button_Y, :left_bumper, :right_bumper,
    :view, :menu, :left_analog, :right_analog,
    :dpad_up, :dpad_right, :dpad_down, :dpad_left
)

const default_axis_mapping = NamedTuple{XBoxAxisLabels}(1:6)
const default_button_mapping = NamedTuple{XBoxButtonLabels}(1:14)

Base.@kwdef struct XBoxAxes
    mapping::NamedTuple{XBoxAxisLabels, NTuple{6,Int}} = default_axis_mapping
    data::MVector{6,Float32} = zeros(MVector{6, Float32})
end

function Base.getindex(axes::XBoxAxes, s::Symbol)
    data = getproperty(axes, :data); mapping = getproperty(axes, :mapping)
    data[mapping[s]]
end

function Base.setindex!(axes::XBoxAxes, v, s::Symbol)
    data = getproperty(axes, :data); mapping = getproperty(axes, :mapping)
    data[mapping[s]] = v
end

Base.@kwdef struct XBoxButtons
    mapping::NamedTuple{XBoxButtonLabels, NTuple{14,Int}} = default_button_mapping
    state::MVector{14,Bool} = zeros(MVector{14,Bool})
    change::MVector{14,ButtonChange} = zeros(MVector{14,ButtonChange})
end

Base.@kwdef struct XBoxController <: AbstractJoystick
    slot::JoystickSlot = GLFW.JOYSTICK_1
    axes::XBoxAxes = XBoxAxes()
    buttons::XBoxButtons = XBoxButtons()
end

XBoxController(slot::JoystickSlot) = XBoxController(; slot)

function update_joystick(joystick::XBoxController)

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
##################### Joystick configuration ################

#adds any joysticks already connected when GLFW is first imported. from this
#moment, connections and disconnections will need to be caught by the callback
#by calling GLFW.PollEvents()
function init_joysticks()

    for slot in instances(JoystickSlot)
        delete!(connected_joysticks, slot)
        if JoystickPresent(slot)
            add_joystick(slot)
        end
    end

end

function joystick_callback(slot::JoystickSlot, event::DeviceConfigEvent)

    if event === GLFW.CONNECTED
        add_joystick(slot)
    elseif event === GLFW.DISCONNECTED
        remove_joystick(slot)
    end

end

function add_joystick(slot::JoystickSlot)

    if !JoystickPresent(slot)
        println("Error while adding device at slot $slot: not found")
        return
    end

    joystick_model = GetJoystickName(slot)

    if joystick_model === "Xbox Controller"
        joystick = XBoxController(slot)
        println("$joystick_model now active at slot $slot")
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


end #module