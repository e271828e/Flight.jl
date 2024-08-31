module Joysticks

using StaticArrays
using UnPack
using GLFW: GLFW, Joystick as JoystickSlot, DeviceConfigEvent, JoystickPresent,
        GetJoystickAxes, GetJoystickButtons, GetJoystickName, SetJoystickCallback

using ..IODevices
using ..Utils

export JoystickSlot, Joystick
export get_connected_joysticks
export get_axis_value, exp_axis_curve
export get_button_state, is_pressed, is_released
export get_button_change, was_pressed, was_released

export AbstractJoystickID, XBoxController, T16000M, GladiatorNXTEvo


################################################################################
################################ AxisSet #######################################

@enum ButtonChange begin
    unchanged = 0
    pressed = 1
    released = 2
end
Base.convert(::Type{ButtonChange}, n::Integer) = ButtonChange(n)
Base.zero(::Type{ButtonChange}) = unchanged

abstract type AbstractAxisSet{N} <: FieldVector{N, Float32} end

abstract type AbstractButtonSet{N, T <: Union{Bool, ButtonChange}} <: FieldVector{N, T} end


########################################

@kwdef struct XBox360ControllerAxes <: AbstractAxisSet{6}
    left_stick_x::Float32 = 0.0
    left_stick_y::Float32 = 0.0
    right_stick_x::Float32 = 0.0
    right_stick_y::Float32 = 0.0
    left_trigger::Float32 = 0.0
    right_trigger::Float32 = 0.0
end


@kwdef struct XBox360ControllerButtons{T} <: AbstractButtonSet{14, T}
    button_A::T = false
    button_B::T = false
    button_X::T = false
    button_Y::T = false
    left_bumper::T = false
    right_bumper::T = false
    view::T = false
    menu::T = false
    left_stick_press::T = false
    right_stick_press::T = false
    dpad_up::T = false
    dpad_right::T = false
    dpad_down::T = false
    dpad_left::T = false
end



########################################

@kwdef struct T16000Axes <: AbstractAxisSet{4}
    stick_x::Float32 = 0.0
    stick_y::Float32 = 0.0
    stick_z::Float32 = 0.0
    throttle::Float32 = 0.0
end


@kwdef struct T16000Buttons{T} <: AbstractButtonSet{20, T}
    button_0::T = false
    button_1::T = false
    button_2::T = false
    button_3::T = false
    button_4::T = false
    button_5::T = false
    button_6::T = false
    button_7::T = false
    button_8::T = false
    button_9::T = false
    button_10::T = false
    button_11::T = false
    button_12::T = false
    button_13::T = false
    button_14::T = false
    button_15::T = false
    hat_up::T = false
    hat_right::T = false
    hat_down::T = false
    hat_left::T = false
end

#ideally we would like JoystickData to be immutable. to compute the values of
#button_change, we need to construct a button_change based on the previous
#values of

#we could store an instance of JoystickData within Joystick, along with the
#slot.
#when get_data is called, we compare

struct JoystickData{NA, NB, A <: AbstractAxisSet{NA}, BS <: AbstractButtonSet{NB, Bool}, BC <: AbstractButtonSet{NB, ButtonChange}}
    axes::A
    button_state::BS
    button_change::BC
end

T16000Data() = JoystickData(T16000Axes(), T16000Buttons{Bool}(), T16000Buttons{ButtonChange}())

# mutable struct Joystick{T <: JoystickData, BS <: AbstractButtonSet}
#     slot::JoystickSlot
#     button_state_last::BS
#     function Joystick(data::T, slot::JoystickSlot = GLFW.JOYSTICK_1) where {T<:JoystickData{NA, NB, A, BS, BC}} where {NA, NB, A, BS, BC}
#         new{T, BS}(slot, data.button_state)
#     end
# end


@kwdef mutable struct Joystick{T <: JoystickData}
    data::T
    slot::JoystickSlot = GLFW.JOYSTICK_1
    function Joystick(data::T, slot::JoystickSlot = GLFW.JOYSTICK_1) where {T<:JoystickData}
        new{T}(data, slot)
    end
end

# function IODevices.get_data!(joystick::Joystick{T}) where {T <: JoystickData{NA, NB, A, BS, BC}} where {NA, NB, A, BS, BC}
#     button_state_current = BS()
#     button_state_last = joystick.data.button_state
#     #using SAcollect, construct a tuple of all the button changes
#     button_change = map(button_state_current, button_state_last) do current, last
#         current
#     end |> BS
#     joystick

# end


# struct Joystick{T, A <: AbstractAxisSet, BS <: AbstractButtonSet{Bool}, BC <: AbstractButtonSet{Float64}}
#     slot::JoystickSlot
#     axes::A
#     button_state_last::BS
#     button_change::BC
#     function Joystick{A, B}() where {A <: AbstractAxisSet, B <: AbstractButtonSet}
#         data_type =

#     end
# end



# XBoxController = Joystick{JoystickData{XBox360ControllerAxes, XBox360ControllerButtons}}

# function JoystickData(joystick::Joystick{A, B}) where {A <: AbstractAxisSet, B <: AbstractButtonSet}
#     button_state = A{Bool}()

# end

###################

#a

# struct AxisSet{N, L}
#     mapping::NamedTuple{L, NTuple{N,Int}}
#     data::MVector{N,Float32}
# end

# function AxisSet(labels::NTuple{N,Symbol} = default_axes_labels(N)) where {N}
#     mapping = NamedTuple{labels}(Tuple(1:N))
#     data = zeros(MVector{N, Float32})
#     AxisSet{N, labels}(mapping, data)
# end

# AxisSet{N}() where {N} = AxisSet(default_axes_labels(N))
# AxisSet(N::Integer) = AxisSet{N}()

# Base.getindex(axes::AxisSet, s::Symbol) = axes.data[axes.mapping[s]]
# Base.setindex!(axes::AxisSet, v, s::Symbol) = (axes.data[axes.mapping[s]] = v)

# default_axes_labels(N::Integer) = Symbol.("axis_".*string.(Tuple(1:N)))

# update!(axes::AxisSet, slot::JoystickSlot) = (axes.data .= GetJoystickAxes(slot))

# function exp_axis_curve(x::Ranged{T}, args...; kwargs...) where {T}
#     exp_axis_curve(T(x), args...; kwargs...)
# end

# function exp_axis_curve(x::Real; strength::Real = 0.0, deadzone::Real = 0.0)

#     a = strength
#     x0 = deadzone

#     abs(x) <= 1 || throw(ArgumentError("Input to exponential curve must be within [-1, 1]"))
#     (x0 >= 0 && x0 <= 1) || throw(ArgumentError("Exponential curve deadzone must be within [0, 1]"))

#     if x > 0
#         y = max(0, (x - x0)/(1 - x0)) * exp( a * (abs(x) -1) )
#     else
#         y = min(0, (x + x0)/(1 - x0)) * exp( a * (abs(x) -1) )
#     end
# end


################################################################################
################################ ButtonSet #####################################



# struct ButtonSet{N, L}
#     mapping::NamedTuple{L, NTuple{N,Int}}
#     state::MVector{N,Bool}
#     change::MVector{N,ButtonChange}
# end

# function ButtonSet(labels::NTuple{N,Symbol} = default_button_labels(N)) where {N}
#     mapping = NamedTuple{labels}(Tuple(1:N))
#     state = zeros(MVector{N,Bool})
#     change = zeros(MVector{N,ButtonChange})
#     ButtonSet{N, labels}(mapping, state, change)
# end

# ButtonSet{N}() where {N} = ButtonSet(default_button_labels(N))
# ButtonSet(N::Integer) = ButtonSet{N}()

# default_button_labels(N::Integer) = Symbol.("button_".*string.(Tuple(1:N)))

# function update!(buttons::ButtonSet{N, L}, slot::JoystickSlot) where {N, L}

#     buttons_state_new = SVector{N,Bool}(GetJoystickButtons(slot))

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


# ################################################################################
# ################################# Joystick #####################################

# abstract type AbstractJoystickID end

# #update_interval: number of display updates per input device update:
# #T_update = T_display * update_interval (where typically T_display =
# #16.67ms). update_interval=0 uncaps the update rate (not recommended!)

# mutable struct Joystick{T <: AbstractJoystickID, A <: AxisSet, B <: ButtonSet} <: InputDevice
#     id::T
#     axes::A
#     buttons::B
#     slot::JoystickSlot
#     update_interval::Int64
#     window::GLFW.Window
#     function Joystick(  id::AbstractJoystickID,
#                         slot::JoystickSlot = GLFW.JOYSTICK_1,
#                         update_interval::Integer = 1)
#         axes = AxisSet(id)
#         buttons = ButtonSet(id)
#         new{typeof(id), typeof(axes), typeof(buttons)}(
#             id, axes, buttons, slot, update_interval) #window uninitialized
#     end
# end

# function Joystick{ID}(args...) where {ID <: AbstractJoystickID}
#     Joystick(ID(), args...)
# end

# function is_connected(joystick::Joystick)
#     #does the reference in our supposed slot still point to us?
#     slot = joystick.slot
#     return (slot in keys(active_slots) && active_slots[slot][] === joystick)
# end

# #override as required by each joystick ID
# rescale!(::AxisSet, ::AbstractJoystickID) = nothing

# function get_axis_value(joystick::Joystick, s::Symbol)
#     @unpack data, mapping = joystick.axes
#     data[mapping[s]]
# end

# function get_button_state(joystick::Joystick, s::Symbol)
#     @unpack state, mapping = joystick.buttons
#     state[mapping[s]]
# end

# function is_pressed(joystick::Joystick, s::Symbol)
#     @unpack state, mapping = joystick.buttons
#     state[mapping[s]] === true
# end

# function is_released(joystick::Joystick, s::Symbol)
#     @unpack state, mapping = joystick.buttons
#     state[mapping[s]] === false
# end

# function get_button_change(joystick::Joystick, s::Symbol)
#     @unpack change, mapping = joystick.buttons
#     change[mapping[s]]
# end

# function was_pressed(joystick::Joystick, s::Symbol)
#     @unpack change, mapping = joystick.buttons
#     change[mapping[s]] === pressed
# end

# function was_released(joystick::Joystick, s::Symbol)
#     @unpack change, mapping = joystick.buttons
#     change[mapping[s]] === released
# end

# function Base.show(::IO, joystick::Joystick)

#     @unpack id, slot, axes, buttons = joystick
#     println()
#     println("$(typeof(id)) at slot $(slot), connected = $(is_connected(joystick))")
#     println("Axes:")
#     for label in keys(axes.mapping)
#         println("\t", label, ": ", get_axis_value(joystick, label))
#     end
#     println("Buttons:")
#     for label in keys(buttons.mapping)
#         println("\t", label, ": ", get_button_state(joystick, label), ", ", get_button_change(joystick,label))
#     end

# end

# function update!(joystick::Joystick)

#     @unpack id, slot, axes, buttons = joystick

#     if !is_connected(joystick)
#         println("Can't update $(joystick.id) at slot $slot, no longer connected")
#         return
#     end
#     update!(axes, slot)
#     update!(buttons, slot)
#     rescale!(axes, id)
#     return nothing

# end
# ## REPL-specific:
# # function Base.show(::IO, ::MIME"text/plain", joystick::Joystick)
# # end

# ########################### IODevices extensions ###############################

# function IODevices.init!(joystick::Joystick)
#     joystick.window = GLFW.CreateWindow(640, 480, "$(string(typeof(joystick.id)))")
#     @unpack window, update_interval = joystick
#     GLFW.HideWindow(window)
#     GLFW.MakeContextCurrent(window)
#     GLFW.SwapInterval(update_interval)
# end

# function IODevices.update_input!(joystick::Joystick, args...)
#     GLFW.SwapBuffers(joystick.window) #blocks to achieve the requested update_interval
#     update!(joystick)
#     GLFW.PollEvents() #see if we got a shutdown request
# end

# IODevices.should_close(joystick::Joystick) = GLFW.WindowShouldClose(joystick.window)
# IODevices.shutdown!(joystick::Joystick) = GLFW.DestroyWindow(joystick.window)


# ############################# Initialization ###################################

# const active_slots = Dict{JoystickSlot, Ref{<:Joystick}}()

# function refresh_joystick_slots()

#     #simply calling PollEvents refreshes all the joystick slots, so when
#     #JoystickPresent is then called, it returns their updated states. there is
#     #no need to explicitly handle the addition or removal of joysticks directly
#     #from a joystick_callback
#     GLFW.PollEvents()
#     for slot in instances(JoystickSlot)
#         delete!(active_slots, slot)
#         if JoystickPresent(slot)
#             add_joystick(slot)
#         end
#     end

#     return active_slots
# end

# function add_joystick(slot::JoystickSlot)

#     if !JoystickPresent(slot)
#         @warn("Could not add joystick at slot $slot: not found")
#         return
#     end

#     joystick_model = GetJoystickName(slot) |> strip

#     if joystick_model == "Xbox Controller"
#         joystick = XBoxController(slot)
#         @info("Adding XBoxController at slot $slot")
#     elseif joystick_model == "T.16000M"
#         joystick = T16000M(slot)
#         @info("Adding Thrustmaster T16000M at slot $slot")
#     elseif joystick_model == "VKBsim Gladiator EVO  R"
#         joystick = GladiatorNXTEvo(slot)
#         @info("Adding VKBSim Gladiator NXT Evo at slot $slot")
#     else
#         @warn("Failed to add $joystick_model at slot $slot, device not supported")
#         return
#     end

#     active_slots[slot] = Ref(joystick)

# end

# function get_connected_joysticks()

#     active_slots = refresh_joystick_slots()
#     return Tuple(jref[] for jref in values(active_slots))

# end


# ################################################################################
# ############################# XBox Controller ##################################

# ############################# Axes & Buttons ###################################

# const XBoxAxisLabels = (
#     :left_stick_x, :left_stick_y, :right_analog_x, :right_analog_y,
#     :left_trigger, :right_trigger
# )

# const XBoxButtonLabels = (
#     :button_A, :button_B, :button_X, :button_Y, :left_bumper, :right_bumper,
#     :view, :menu, :left_stick, :right_analog,
#     :dpad_up, :dpad_right, :dpad_down, :dpad_left
# )

# ############################### Controller #####################################

# struct XBoxControllerID <: AbstractJoystickID end
# const XBoxController = Joystick{XBoxControllerID}

# AxisSet(::XBoxControllerID) = AxisSet(XBoxAxisLabels)
# ButtonSet(::XBoxControllerID) = ButtonSet(XBoxButtonLabels)

# function rescale!(axes::AxisSet, ::XBoxControllerID)
#     axes[:left_trigger] = 0.5*(1 + axes[:left_trigger])
#     axes[:right_trigger] = 0.5*(1 + axes[:right_trigger])
# end

# ################################################################################
# ########################### Thrustmaster T.16000M ##############################

# ############################# Axes & Buttons ###################################

# const T16000MAxisLabels = (
#     :stick_x, :stick_y, :stick_z, :throttle
# )

# const T16000MButtonLabels = (
#     :button_0, :button_1, :button_2, :button_3, :button_4, :button_5,
#     :button_6, :button_7, :button_8, :button_9, :button_10, :button_11,
#     :button_12, :button_13, :button_14, :button_15,
#     :hat_up, :hat_right, :hat_down, :hat_left
# )

# ############################### Controller #####################################

# struct T16000M_ID <: AbstractJoystickID end
# const T16000M = Joystick{T16000M_ID}

# AxisSet(::T16000M_ID) = AxisSet(T16000MAxisLabels)
# ButtonSet(::T16000M_ID) = ButtonSet(T16000MButtonLabels)

# function rescale!(axes::AxisSet, ::T16000M_ID)
#     axes[:throttle] = 0.5*(1 - axes[:throttle])
# end

# ################################################################################
# ######################## VKBSim Gladiator NXT Evo ##############################

# ############################# Axes & Buttons ###################################

# #factory defaults
# const GladiatorNXTEvoAxisLabels = (
#     :stick_x, :stick_y, :throttle, :analog_hat_x, :analog_hat_y, :stick_z,
# )

# function gen_GladiatorNXTEvoButtonLabels()
#     labels = Symbol.("unmapped_".*string.(1:132))
#     labels[1:29] .= (:red_trigger_half, :red_trigger_full, :A2, :B1, :D1,
#     :A3_up, :A3_right, :A3_down, :A3_left, :A3_press,
#     :A4_up, :A4_right, :A4_down, :A4_left, :A4_press,
#     :C1_up, :C1_right, :C1_down, :C1_left, :C1_press,
#     :black_trigger_up, :black_trigger_down, :encoder_up, :encoder_down,
#     :switch_up, :switch_down, :F1, :F2, :F3)
#     return Tuple(labels)
# end

# #factory defaults
# const GladiatorNXTEvoButtonLabels = gen_GladiatorNXTEvoButtonLabels()

# ############################## Controller #####################################

# struct GladiatorNXTEvo_ID <: AbstractJoystickID end
# const GladiatorNXTEvo = Joystick{GladiatorNXTEvo_ID}

# AxisSet(::GladiatorNXTEvo_ID) = AxisSet(GladiatorNXTEvoAxisLabels)
# ButtonSet(::GladiatorNXTEvo_ID) = ButtonSet(GladiatorNXTEvoButtonLabels)

# function rescale!(axes::AxisSet, ::GladiatorNXTEvo_ID)
#     axes[:throttle] = 0.5*(1 - axes[:throttle])
# end

end #module