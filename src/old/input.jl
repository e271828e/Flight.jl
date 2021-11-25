# module Input

using StaticArrays
using ComponentArrays
using UnPack

using GLFW: GLFW, Joystick, DeviceConfigEvent, JoystickPresent,
        GetJoystickAxes, GetJoystickButtons, GetJoystickName

export XBoxInterface
export get_axis_data, get_button_state, get_button_change

abstract type JoystickInterface end

const connected_joysticks = Dict{Joystick, JoystickInterface}()

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
const default_axis_ordering = XBoxAxisLabels
const default_button_ordering = XBoxButtonLabels

# Base.@kwdef struct XBoxAxes
#     mapping::NamedTuple{XBoxAxisLabels, NTuple{6,Int}} = default_axis_mapping
#     data::MVector{6,Float32} = zeros(MVector{6, Float32})
# end

# function Base.getindex(axes::XBoxAxes, s::Symbol)
#     data = getproperty(axes, :data); mapping = getproperty(axes, :mapping)
#     data[mapping[s]]
# end

# function Base.setindex!(axes::XBoxAxes, v, s::Symbol)
#     data = getproperty(axes, :data); mapping = getproperty(axes, :mapping)
#     data[mapping[s]] = v
# end

struct XBoxVector{N,L,T}
    data::MVector{N,T}
    mapping::NamedTuple{L, NTuple{N,Int}}
    function XBoxVector{T}(mapping::NamedTuple{L, NTuple{N,Int}}) where {N,L,T}
        new{N,L,T}(zeros(MVector{N, T}), mapping)
    end
end

function Base.getindex(buttons::XBoxVector, s::Symbol)
    data = getproperty(buttons, :data); mapping = getproperty(buttons, :mapping)
    data[mapping[s]]
end

function Base.setindex!(buttons::XBoxVector, v, s::Symbol)
    state = getproperty(buttons, :data); mapping = getproperty(buttons, :mapping)
    state[mapping[s]] = v
end

Base.@kwdef struct XBoxInterface{A,B,C} <: JoystickInterface
    id::Joystick = GLFW.JOYSTICK_1
    axes::A
    button_states::B
    button_changes::C
    function XBoxInterface(axis_labels::NTuple{NAx,Symbol}, button_labels::NTuple{NBt,Symbols}) where {NAx, NBt}
        axes = ComponentVector(NamedTuple{axis_labels}(zeros(Float32, NAx)))
        button_states = ComponentVector(NamedTuple{button_labels}(zeros(Bool, NBt)))
        button_changes = ComponentVector(NamedTuple{button_labels}(zeros(Integer, NBt)))
        new{typeof()}
end

# Base.@kwdef struct XBoxInterface <: JoystickInterface
#     id::Joystick = GLFW.JOYSTICK_1
#     axes::XBoxAxes = XBoxAxes()
#     button_state::XBoxButtons{Bool} = XBoxButtons{Bool}()
#     button_change::XBoxButtons{ButtonChange} = XBoxButtons{ButtonChange}()
# end


# XBoxInterface(id::Joystick) = XBoxInterface(; id)

# function update(interface::XBoxInterface)

#     @unpack id, axes, button_state, button_change = interface

#     axes.data .= GetJoystickAxes(id)
#     left_trigger = axes[:left_trigger]
#     left_trigger = axes[:left_trigger]
#     axes[:left_trigger] += 1
#     axes[:right_trigger] += 1

#     buttons_state_new = SVector{14,Bool}(GetJoystickButtons(id))

#     for i in 1:length(buttons.state)
#         if !buttons.state[i] && buttons_state_new[i]
#             buttons.change[i] = pressed
#         elseif buttons.state[i] && !buttons_state_new[i]
#             buttons.change[i] = released
#         else
#             buttons.change[i] = unchanged
#         end
#     end

#     buttons.state .= buttons_state_new

# end

# function get_axis_data(interface::XBoxInterface, s::Symbol)
#     @unpack state, mapping = interface.axes
#     state[mapping[s]]
# end

# function get_button_state(interface::XBoxInterface, s::Symbol)
#     @unpack state, mapping = interface.buttons
#     state[mapping[s]]
# end

# function get_button_change(interface::XBoxInterface, s::Symbol)
#     @unpack change, mapping = interface.buttons
#     change[mapping[s]]
# end


# ##################### Joystick configuration ################

# #adds any joysticks already connected when GLFW is first imported. from this
# #moment, connections and disconnections will need to be caught by the callback
# #by calling GLFW.PollEvents()
# function init_joysticks() #used to initialize any joysticks connected upon GLFW
#     for joy_id in instances(Joystick)
#         delete!(connected_joysticks, joy_id)
#         if JoystickPresent(joy_id)
#             add_joystick(joy_id)
#         end
#     end
# end

# function joystick_callback(joy_id::Joystick, event::DeviceConfigEvent)
#     if event === GLFW.CONNECTED
#         add_joystick(joy_id)
#     elseif event === GLFW.DISCONNECTED
#         remove_joystick(joy_id)
#     end
# end

# function set_joystick_callback()
#     GLFW.SetJoystickCallback(joystick_callback)
# end

# function add_joystick(joy_id::Joystick)
#     if !JoystickPresent(joy_id)
#         println("Failed to add $joy_id")
#         return
#     end
#     joystick_model = GetJoystickName(joy_id)
#     if joystick_model === "Xbox Controller"
#         interface = XBoxInterface(joy_id)
#         println("Successfully added $joystick_model as $joy_id")
#     else
#         error("Interface for joystick $joystick_model not implemented")
#     end
#     connected_joysticks[joy_id] = interface
# end

# function remove_joystick(joy_id::Joystick)
#     if JoystickPresent(joy_id)
#         joystick_model = GetJoystickName(joy_id)
#         println("Failed to remove $joystick_model as $joy_id, still connected")
#         return
#     end
#     delete!(connected_joysticks, joy_id)
# end











# end #module