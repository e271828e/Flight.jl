module Input

using StaticArrays
using UnPack


using GLFW: GLFW, Joystick as JoystickSlot, DeviceConfigEvent, JoystickPresent,
        GetJoystickAxes, GetJoystickButtons, GetJoystickName, SetJoystickCallback

export AbstractInputInterface, AbstractInputMapping, DefaultInputMapping
export AbstractJoystick, XBoxController
export connected_joysticks, init_joysticks
export get_axis_value, get_button_state, get_button_change, is_pressed, is_released


################################################################################
########################### AbstractInputInterface #############################

abstract type AbstractInputInterface end

#for every input interface it supports, a target should extend the three
#argument assign! method below, with a ::DefaultInputMapping argument. for
#alternative mappings, it can define a subtype of AbstractInputMapping, and then
#extend assign! using that type as third argument
abstract type AbstractInputMapping end
struct DefaultInputMapping <: AbstractInputMapping end

#to be extended for concrete targets and interfaces
update!(input::AbstractInputInterface) = throw(MethodError(update!, (input,)))
assign!(target::Any, input::AbstractInputInterface, ::DefaultInputMapping) = throw(MethodError(assign!, (target, input)))


assign!(target::Any, input::AbstractInputInterface) = assign!(target, input, DefaultInputMapping())


################################################################################
############################## AbstractJoystick ################################


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


############################## ButtonChange ####################################

@enum ButtonChange begin
    unchanged = 0
    pressed = 1
    released = 2
end
Base.convert(::Type{ButtonChange}, n::Integer) = ButtonChange(n)
Base.zero(::Type{ButtonChange}) = unchanged



############################# Joystick initialization ##########################

function init_joysticks()

    #simply calling PollEvents refreshes all the joystick slots, so when
    #JoystickPresent is then called, it returns their updated states. there is
    #no need to explicitly handle the addition or removal of joysticks directly
    #from a joystick_callback
    GLFW.PollEvents()
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
        println("Could not add joystick at slot $slot: not found")
        return
    end

    joystick_model = GetJoystickName(slot)

    if joystick_model === "Xbox Controller"
        joystick = XBoxController(slot)
        println("XBoxController active at slot $slot")
    else
        println("$joystick_model not supported")
    end

    connected_joysticks[slot] = joystick

end

# #let's see if this
# function is_connected(joy::AbstractJoystick)

# struct GenericJoystick end


default_axes_labels(N::Integer) = Symbol.("axis_".*string.(Tuple(1:N)))
default_button_labels(N::Integer) = Symbol.("button_".*string.(Tuple(1:N)))

struct AxisSet{N, L}
    mapping::NamedTuple{L, NTuple{N,Int}}
    data::MVector{6,Float32}
end

function AxisSet(labels::NTuple{N,Symbol} = default_axes_labels(N)) where {N}
    mapping = NamedTuple{labels}(Tuple(1:N))
    data = zeros(MVector{6, Float32})
    AxisSet{N, labels}(mapping, data)
end

AxisSet{N}() where {N} = AxisSet(default_axes_labels(N))
AxisSet(N::Integer) = AxisSet{N}()


struct ButtonSet5{N, L}
    mapping::NamedTuple{L, NTuple{N,Int}}
    state::MVector{N,Bool}
    change::MVector{N,ButtonChange}
end

function ButtonSet5(labels::NTuple{N,Symbol} = default_button_labels(N)) where {N}
    mapping = NamedTuple{labels}(Tuple(1:N))
    state = zeros(MVector{N,Bool})
    change = zeros(MVector{N,ButtonChange})
    ButtonSet5{N, labels}(mapping, state, change)
end

ButtonSet5{N}() where {N} = ButtonSet5(default_button_labels(N))
ButtonSet5(N::Integer) = ButtonSet5{N}()


################################################################################
############################# XBox Controller ##################################

################################ Axes ##########################################


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


############################### Buttons ########################################

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

############################### Controller #####################################

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

#to show values in the REPL
function Base.show(::IO, ::MIME"text/plain", joystick::XBoxController)
    println("Hi")
end

#to print
# function Base.show(::IO, joy::XBoxController)
#     println("Hi from print")
# end

end #module