module Joysticks

using StaticArrays
using UnPack
using GLFW: GLFW, Joystick as JoystickSlot, JoystickPresent,
        GetJoystickAxes, GetJoystickButtons, GetJoystickName

using ..IODevices
using ..Types

export JoystickData, Joystick
export get_connected_joysticks
export get_axis_value, exp_axis_curve
export get_button_state, is_pressed
export get_button_change, was_pressed, was_released

export XBoxController, XBoxControllerData
export T16000M, T16000MData
export GladiatorNXTEvo, GladiatorNXTEvoData


################################################################################
############################# ButtonChange #####################################

@enum ButtonChange begin
    ButtonUnchanged = 0
    ButtonPressed = 1
    ButtonReleased = 2
end

Base.convert(::Type{ButtonChange}, n::Integer) = ButtonChange(n)


################################################################################
############################# AbstractAxisSet ##################################

abstract type AbstractAxisSet{N} <: FieldVector{N, Float32} end

rescale(axes::AbstractAxisSet) = axes


################################################################################
############################# AbstractButtonSet ################################

abstract type AbstractButtonSet{N, T <: Union{Bool, ButtonChange}} <: FieldVector{N, T} end


################################################################################
############################### JoystickData ###################################

struct JoystickData{A <: AbstractAxisSet, BS <: AbstractButtonSet, BC <: AbstractButtonSet}
    axes::A
    button_state::BS
    button_change::BC
end

get_axis_value(data::JoystickData, s::Symbol) = getproperty(data.axes, s)
get_button_state(data::JoystickData, s::Symbol) = getproperty(data.button_state, s)
is_pressed(data::JoystickData, s::Symbol) = get_button_state(data, s) === true

get_button_change(data::JoystickData, s::Symbol) = getproperty(data.button_change, s)
was_pressed(data::JoystickData, s::Symbol) = get_button_change(data, s) === ButtonPressed
was_released(data::JoystickData, s::Symbol) = get_button_change(data, s) === ButtonReleased


################################################################################
################################# Joystick #####################################

mutable struct Joystick{T <: JoystickData} <: InputDevice
    slot::JoystickSlot
    cache::T
end

function IODevices.get_data!(joystick::Joystick{T}) where {T <: JoystickData{A, BS, BC}} where {A, BS, BC}

    #from the axis values returned by GetJoystickAxes we only keep as many
    #as defined by our AbstractButtonSet
    axes = view(GetJoystickAxes(joystick.slot), 1:length(A)) |> A |> rescale

    #from the button values returned by GetJoystickButtons we only keep as many
    #as defined by our AbstractButtonSet
    button_state = view(GetJoystickButtons(joystick.slot), 1:length(BS)) |> BS

    button_state_last = joystick.cache.button_state

    button_change = map(button_state, button_state_last) do current, last
        (current && !last) && return ButtonPressed
        (!current && last) && return ButtonReleased
        return ButtonUnchanged
    end |> BC

    data = JoystickData(axes, button_state, button_change)

    joystick.cache = data

    return data
end

################################################################################
########################### Thrustmaster T.16000M ##############################

@kwdef struct T16000MAxes <: AbstractAxisSet{4}
    stick_x::Float32 = 0.0
    stick_y::Float32 = 0.0
    stick_z::Float32 = 0.0
    throttle::Float32 = 0.0
end

@kwdef struct T16000MButtons{T} <: AbstractButtonSet{20, T}
    button_0::T = 0
    button_1::T = 0
    button_2::T = 0
    button_3::T = 0
    button_4::T = 0
    button_5::T = 0
    button_6::T = 0
    button_7::T = 0
    button_8::T = 0
    button_9::T = 0
    button_10::T = 0
    button_11::T = 0
    button_12::T = 0
    button_13::T = 0
    button_14::T = 0
    button_15::T = 0
    hat_up::T = 0
    hat_right::T = 0
    hat_down::T = 0
    hat_left::T = 0
end

const T16000MData = JoystickData{T16000MAxes,
                                T16000MButtons{Bool},
                                T16000MButtons{ButtonChange}}

const T16000M = Joystick{T16000MData}

get_T16000MData() = JoystickData(T16000MAxes(),
                                T16000MButtons{Bool}(),
                                T16000MButtons{ButtonChange}())

get_T16000M(slot::JoystickSlot) = Joystick(slot, get_T16000MData())

function rescale(axes::T16000MAxes)
    @unpack stick_x, stick_y, stick_z, throttle = axes
    T16000MAxes(; stick_x, stick_y, stick_z, throttle = 0.5*(1 - throttle))
end

################################################################################
############################ XBox Controller ###################################

@kwdef struct XBoxControllerAxes <: AbstractAxisSet{6}
    left_stick_x::Float32 = 0.0
    left_stick_y::Float32 = 0.0
    right_stick_x::Float32 = 0.0
    right_stick_y::Float32 = 0.0
    left_trigger::Float32 = 0.0
    right_trigger::Float32 = 0.0
end

@kwdef struct XBoxControllerButtons{T} <: AbstractButtonSet{14, T}
    button_A::T = 0
    button_B::T = 0
    button_X::T = 0
    button_Y::T = 0
    left_bumper::T = 0
    right_bumper::T = 0
    view::T = 0
    menu::T = 0
    left_stick_press::T = 0
    right_stick_press::T = 0
    dpad_up::T = 0
    dpad_right::T = 0
    dpad_down::T = 0
    dpad_left::T = 0
end

const XBoxControllerData = JoystickData{XBoxControllerAxes,
                                        XBoxControllerButtons{Bool},
                                        XBoxControllerButtons{ButtonChange}}

const XBoxController = Joystick{XBoxControllerData}

get_XBoxControllerData() = JoystickData(XBoxControllerAxes(),
                                        XBoxControllerButtons{Bool}(),
                                        XBoxControllerButtons{ButtonChange}())

get_XBoxController(slot::JoystickSlot) = Joystick(slot, get_XBoxControllerData())

function rescale(axes::XBoxControllerAxes)
    @unpack left_stick_x, left_stick_y, right_stick_x, right_stick_y,
            left_trigger, right_trigger = axes
    XBoxControllerAxes(; left_stick_x, left_stick_y, right_stick_x, right_stick_y,
            left_trigger = 0.5*(1 + left_trigger),
            right_trigger = 0.5*(1 + right_trigger))
end


################################################################################
######################## VKBSim Gladiator NXT Evo ##############################

@kwdef struct GladiatorNXTEvoAxes <: AbstractAxisSet{6}
    stick_x::Float32 = 0.0
    stick_y::Float32 = 0.0
    throttle::Float32 = 0.0
    analog_hat_x::Float32 = 0.0
    analog_hat_y::Float32 = 0.0
    stick_z::Float32 = 0.0
end

#when queried with GetJoystickButtons, the Gladiator NXT Evo returns an array of
#132 values, from which the first 29 actually correspond to physical buttons. we
#only keep those here
@kwdef struct GladiatorNXTEvoButtons{T} <: AbstractButtonSet{29, T}
    red_trigger_half::T = 0
    red_trigger_full::T = 0
    A2::T = 0
    B1::T = 0
    D1::T = 0
    A3_up::T = 0
    A3_right::T = 0
    A3_down::T = 0
    A3_left::T = 0
    A3_press::T = 0
    A4_up::T = 0
    A4_right::T = 0
    A4_down::T = 0
    A4_left::T = 0
    A4_press::T = 0
    C1_up::T = 0
    C1_right::T = 0
    C1_down::T = 0
    C1_left::T = 0
    C1_press::T = 0
    black_trigger_up::T = 0
    black_trigger_down::T = 0
    encoder_up::T = 0
    encoder_down::T = 0
    switch_up::T = 0
    switch_down::T = 0
    F1::T = 0
    F2::T = 0
    F3::T = 0
end

const GladiatorNXTEvoData = JoystickData{GladiatorNXTEvoAxes,
                                        GladiatorNXTEvoButtons{Bool},
                                        GladiatorNXTEvoButtons{ButtonChange}}

const GladiatorNXTEvo = Joystick{GladiatorNXTEvoData}

get_GladiatorNXTEvoData() = JoystickData(GladiatorNXTEvoAxes(),
                                        GladiatorNXTEvoButtons{Bool}(),
                                        GladiatorNXTEvoButtons{ButtonChange}())

get_GladiatorNXTEvo(slot::JoystickSlot) = Joystick(slot, get_GladiatorNXTEvoData())

function rescale(axes::GladiatorNXTEvoAxes)
    @unpack stick_x, stick_y, throttle, analog_hat_x, analog_hat_y, stick_z = axes
    GladiatorNXTEvoAxes(; stick_x, stick_y, throttle = 0.5*(1 - throttle),
                        analog_hat_x, analog_hat_y, stick_z)
end


################################################################################
############################# Initialization ###################################

const connected_slots = Dict{JoystickSlot, Joystick}()

#specify a Joystick constructor for each joystick model
const supported_joysticks = Dict{String, Function}(
    "Xbox Controller" => get_XBoxController, #XBox 360 Controller
    "T.16000M" => get_T16000M, #Thrustmaster T16000M
    "VKBsim Gladiator EVO  R" => get_GladiatorNXTEvo, #VKBSim Gladiator NXT Evo
)

function refresh_joysticks()

    #calling PollEvents refreshes all the joystick slots, so JoystickPresent
    #then returns their updated states. no need to explicitly handle the
    #addition or removal of joysticks directly from a joystick_callback (which
    #would also require a PollEvents call to work)
    GLFW.PollEvents()
    empty!(connected_slots)
    for slot ∈ instances(JoystickSlot)
        JoystickPresent(slot) && add_joystick(slot)
    end

    return connected_slots
end

get_connected_joysticks() = Tuple(values(refresh_joysticks()))

function add_joystick(slot::JoystickSlot)

    if !JoystickPresent(slot)
        @warn("Failed to add joystick at slot $slot: not found")
        return
    end

    joystick_model = GetJoystickName(slot) |> strip

    if !(joystick_model ∈ keys(supported_joysticks))
        @warn("Failed to add $joystick_model at slot $slot: not supported")
        return
    end

    connected_slots[slot] = supported_joysticks[joystick_model](slot)

    return connected_slots[slot]

end

function is_connected(joystick::Joystick)
    slot = joystick.slot
    #this joystick's slot must be listed in connected_slots, and the
    #corresponding key must actually point to this joystick instance
    return (slot ∈ keys(connected_slots) && connected_slots[slot] === joystick)
end


################################################################################
############################ Misc Stuff ########################################

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



end #module