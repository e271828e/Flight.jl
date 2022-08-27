module Joysticks

using StaticArrays
using UnPack
using GLFW: GLFW, Joystick as JoystickSlot, DeviceConfigEvent, JoystickPresent,
        GetJoystickAxes, GetJoystickButtons, GetJoystickName, SetJoystickCallback

using Flight.IODevices
using Flight.Utils: Ranged

export JoystickSlot, Joystick
export get_connected_joysticks
export get_axis_value, exp_axis_curve
export get_button_state, get_button_change, was_pressed, was_released

export XBoxControllerID, XBoxController


################################################################################
################################ AxisSet #######################################

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

Base.getindex(axes::AxisSet, s::Symbol) = axes.data[axes.mapping[s]]
Base.setindex!(axes::AxisSet, v, s::Symbol) = (axes.data[axes.mapping[s]] = v)

default_axes_labels(N::Integer) = Symbol.("axis_".*string.(Tuple(1:N)))

update!(axes::AxisSet, slot::JoystickSlot) = (axes.data .= GetJoystickAxes(slot))

function exp_axis_curve(x::Ranged{T}, args...; kwargs...) where {T}
    exp_axis_curve(T(x), args...; kwargs...)
end

function exp_axis_curve(x::Real; strength::Real = 0.0, deadzone::Real = 0.0)

    a = strength
    x0 = deadzone

    abs(x) <= 1 || throw(ArgumentError("Input to exponential curve must be within [-1, 1]"))
    (x0 >= 0 && x0 <= 1) || throw(ArgumentError("Exponential curve deadzone must be within [0, 1]"))

    if x > 0
        y = max(0, (x - x0)/(1 - x0)) * exp( a * (abs(x) -1) )
    else
        y = min(0, (x + x0)/(1 - x0)) * exp( a * (abs(x) -1) )
    end
end


################################################################################
################################ ButtonSet #####################################

@enum ButtonChange begin
    unchanged = 0
    pressed = 1
    released = 2
end
Base.convert(::Type{ButtonChange}, n::Integer) = ButtonChange(n)
Base.zero(::Type{ButtonChange}) = unchanged


struct ButtonSet{N, L}
    mapping::NamedTuple{L, NTuple{N,Int}}
    state::MVector{N,Bool}
    change::MVector{N,ButtonChange}
end

function ButtonSet(labels::NTuple{N,Symbol} = default_button_labels(N)) where {N}
    mapping = NamedTuple{labels}(Tuple(1:N))
    state = zeros(MVector{N,Bool})
    change = zeros(MVector{N,ButtonChange})
    ButtonSet{N, labels}(mapping, state, change)
end

ButtonSet{N}() where {N} = ButtonSet(default_button_labels(N))
ButtonSet(N::Integer) = ButtonSet{N}()

default_button_labels(N::Integer) = Symbol.("axis_".*string.(Tuple(1:N)))

function update!(buttons::ButtonSet{N, L}, slot::JoystickSlot) where {N, L}

    buttons_state_new = SVector{N,Bool}(GetJoystickButtons(slot))

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


################################################################################
################################# Joystick #####################################

abstract type AbstractJoystickID end

#update_interval: number of display updates per input device update:
#T_update = T_display * update_interval (where typically T_display =
#16.67ms). update_interval=0 uncaps the update rate (not recommended!)

mutable struct Joystick{T <: AbstractJoystickID, A <: AxisSet, B <: ButtonSet} <: IODevice
    id::T
    axes::A
    buttons::B
    slot::JoystickSlot
    update_interval::Int64
    window::GLFW.Window
    function Joystick(  id::AbstractJoystickID,
                        slot::JoystickSlot = GLFW.JOYSTICK_1,
                        update_interval::Integer = 1)
        axes = AxisSet(id)
        buttons = ButtonSet(id)
        new{typeof(id), typeof(axes), typeof(buttons)}(
            id, axes, buttons, slot, update_interval) #window uninitialized
    end
end

function is_connected(joystick::Joystick)
    #does the reference in our supposed slot still point to us?
    slot = joystick.slot
    return (slot in keys(active_slots) && active_slots[slot][] === joystick)
end

#override as required by each joystick ID
rescale!(::AxisSet, ::AbstractJoystickID) = nothing

function get_axis_value(joystick::Joystick, s::Symbol)
    @unpack data, mapping = joystick.axes
    data[mapping[s]]
end

function get_button_state(joystick::Joystick, s::Symbol)
    @unpack state, mapping = joystick.buttons
    state[mapping[s]]
end

function get_button_change(joystick::Joystick, s::Symbol)
    @unpack change, mapping = joystick.buttons
    change[mapping[s]]
end

function was_pressed(joystick::Joystick, s::Symbol)
    @unpack change, mapping = joystick.buttons
    change[mapping[s]] === pressed
end

function was_released(joystick::Joystick, s::Symbol)
    @unpack change, mapping = joystick.buttons
    change[mapping[s]] === released
end

function Base.show(::IO, joystick::Joystick)

    @unpack id, slot, axes, buttons = joystick
    println()
    println("$(typeof(id)) at slot $(slot), connected = $(is_connected(joystick))")
    println("Axes:")
    for label in keys(axes.mapping)
        println("\t", label, ": ", get_axis_value(joystick, label))
    end
    println("Buttons:")
    for label in keys(buttons.mapping)
        println("\t", label, ": ", get_button_state(joystick, label), ", ", get_button_change(joystick,label))
    end

end

## REPL-specific:
# function Base.show(::IO, ::MIME"text/plain", joystick::Joystick)
# end

########################### IODevices extensions ###############################

function IODevices.init!(joystick::Joystick)
    joystick.window = GLFW.CreateWindow(640, 480, "$(string(typeof(joystick.id)))")
    @unpack window, update_interval = joystick
    # GLFW.HideWindow(window)
    GLFW.MakeContextCurrent(window)
    GLFW.SwapInterval(update_interval)
end

function IODevices.update!(joystick::Joystick)

    @unpack id, slot, axes, buttons, window = joystick

    if !is_connected(joystick)
        println("Can't update $(joystick.id) at slot $slot, no longer connected")
        return
    end

    GLFW.SwapBuffers(window) #honor the requested update_interval
    update!(axes, slot)
    update!(buttons, slot)
    rescale!(axes, id)
    GLFW.PollEvents() #see if we got a shutdown request

end

IODevices.should_close(joystick::Joystick) = GLFW.WindowShouldClose(joystick.window)
IODevices.shutdown!(joystick::Joystick) = GLFW.DestroyWindow(joystick.window)


############################# Initialization ###################################

const active_slots = Dict{JoystickSlot, Ref{<:Joystick}}()

function refresh_joystick_slots()

    #simply calling PollEvents refreshes all the joystick slots, so when
    #JoystickPresent is then called, it returns their updated states. there is
    #no need to explicitly handle the addition or removal of joysticks directly
    #from a joystick_callback
    GLFW.PollEvents()
    for slot in instances(JoystickSlot)
        delete!(active_slots, slot)
        if JoystickPresent(slot)
            add_joystick(slot)
        end
    end

    return active_slots
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

    active_slots[slot] = Ref(joystick)

end

function get_connected_joysticks()

    active_slots = refresh_joystick_slots()
    return Tuple(jref[] for jref in values(active_slots))

end


################################################################################
############################# XBox Controller ##################################

############################# Axes & Buttons ###################################

const XBoxAxisLabels = (
    :left_analog_x, :left_analog_y, :right_analog_x, :right_analog_y,
    :left_trigger, :right_trigger
)

const XBoxButtonLabels = (
    :button_A, :button_B, :button_X, :button_Y, :left_bumper, :right_bumper,
    :view, :menu, :left_analog, :right_analog,
    :dpad_up, :dpad_right, :dpad_down, :dpad_left
)

############################### Controller #####################################

struct XBoxControllerID <: AbstractJoystickID end

AxisSet(::XBoxControllerID) = AxisSet(XBoxAxisLabels)
ButtonSet(::XBoxControllerID) = ButtonSet(XBoxButtonLabels)

XBoxController(slot::JoystickSlot = GLFW.JOYSTICK_1) = Joystick(XBoxControllerID(), slot)

function rescale!(axes::AxisSet, ::XBoxControllerID)
    axes[:left_trigger] = 0.5*(1 + axes[:left_trigger])
    axes[:right_trigger] = 0.5*(1 + axes[:right_trigger])
end


end #module