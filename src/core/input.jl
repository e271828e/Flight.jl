module Input

using StaticArrays
using UnPack

using GLFW: GLFW, Joystick as JoystickSlot, DeviceConfigEvent, JoystickPresent,
        GetJoystickAxes, GetJoystickButtons, GetJoystickName, SetJoystickCallback

export AbstractInputInterface, AbstractInputMapping, DefaultInputMapping
export JoystickSlot, Joystick, XBoxController
export get_connected_joysticks
export get_axis_value, get_button_state, get_button_change, was_pressed, was_released


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
################################# Joystick #####################################


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


################################ Joystick ######################################

abstract type AbstractJoystickID end

struct Joystick{T <: AbstractJoystickID, A <: AxisSet, B <: ButtonSet} <: AbstractInputInterface
    id::T
    axes::A
    buttons::B
    slot::JoystickSlot
end

function is_connected(joystick::Joystick)
    #does the reference in our supposed slot still point to us?
    return active_slots[joystick.slot][] === joystick
end

function update!(joystick::Joystick)

    @unpack id, slot, axes, buttons = joystick

    if !is_connected(joystick)
        println("Can't update $(joystick.id) at slot $slot, it's no longer connected")
        return
    end

    update!(axes, slot)
    update!(buttons, slot)
    rescale!(axes, id)

end

#to override as required by each joystick ID
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


############################# Joystick initialization ##########################

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
        joystick = xboxcontroller(slot)
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

################################ Axes ##########################################

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

struct XBoxController <: AbstractJoystickID end

function xboxcontroller(slot::JoystickSlot = GLFW.JOYSTICK_1)
    axes = AxisSet(XBoxAxisLabels)
    buttons = ButtonSet(XBoxButtonLabels)
    Joystick(XBoxController(), axes, buttons, slot)
end

function rescale!(axes::AxisSet, ::XBoxController)
    axes[:left_trigger] = 0.5*(1 + axes[:left_trigger])
    axes[:right_trigger] = 0.5*(1 + axes[:right_trigger])
end

#to print
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

# to show in the REPL
# function Base.show(::IO, ::MIME"text/plain", joystick::Joystick)
# end



end #module