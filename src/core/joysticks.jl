module Joysticks

using StaticArrays
using StaticArrays: sacollect
using SDL2_jll

using ..IODevices
using ..Types

export AbstractJoystickData, Joystick
export T16000MData, TWCSData, GladiatorNXTEvoData

export connected_joysticks, update_connected_joysticks
export is_pressed, was_pressed, was_released
export exp_axis_curve


################################################################################
############################ SDL Initialization ################################

#observations:
#0) if a joystick is disconnected, update_joysticks() must be called for
#is_connected to return false. if we close it explicitly, it also returns false
#1) we unplug a joystick, call update_joysticks, plug it again and call
#update_joysticks. even if the second time its device_index is the same, its
#instance id will have increased. so an instance id is never reused even if it
#becomes available again, unless we restart the REPL

const SDL_INIT_JOYSTICK = Cuint(0x00000200)

const SDL_HAT_UP = 0x01
const SDL_HAT_RIGHT = 0x02
const SDL_HAT_DOWN = 0x04
const SDL_HAT_LEFT = 0x08
const SDL_HAT_4POS = (SDL_HAT_UP, SDL_HAT_RIGHT, SDL_HAT_DOWN, SDL_HAT_LEFT)

const SDL_JoystickInstanceID = Int32

struct SDL_JoystickGUID
    data::NTuple{16, UInt8}
end

mutable struct SDL_Joystick end

#initialize SDL2's joystick functionality
init_joysticks() = ccall((:SDL_Init, libsdl2), Cint, (UInt32,), SDL_INIT_JOYSTICK)

num_joysticks() = ccall((:SDL_NumJoysticks, libsdl2), Cint, ())

#update joystick status and data (thread-safe)
function update_joysticks()
    ccall((:SDL_LockJoysticks, libsdl2), Cvoid, ())
    ccall((:SDL_JoystickUpdate, libsdl2), Cvoid, ())
    ccall((:SDL_UnlockJoysticks, libsdl2), Cvoid, ())
end

get_instance_id(idx::Integer) = ccall(
    (:SDL_JoystickGetDeviceInstanceID, libsdl2), SDL_JoystickInstanceID, (Cint,), idx)

get_guid(idx::Integer) = ccall(
    (:SDL_JoystickGetDeviceGUID, libsdl2), SDL_JoystickGUID, (Cint,), idx)

get_name(idx::Integer) = unsafe_string(ccall(
    (:SDL_JoystickNameForIndex, libsdl2), Ptr{Cchar}, (Cint,), idx))

open(idx::Integer)::Ptr{SDL_Joystick} = ccall(
    (:SDL_JoystickOpen, libsdl2), Ptr{SDL_Joystick}, (Cint,), idx)


################################ post-open #####################################

get_instance_id(ptr::Ptr{SDL_Joystick}) = ccall(
    (:SDL_JoystickInstanceID, libsdl2), SDL_JoystickInstanceID, (Ptr{SDL_Joystick},), ptr)

get_guid(ptr::Ptr{SDL_Joystick}) = ccall(
    (:SDL_JoystickGetGUID, libsdl2), SDL_JoystickGUID, (Ptr{SDL_Joystick},), ptr)

get_name(ptr::Ptr{SDL_Joystick}) = ccall(
    (:SDL_JoystickName, libsdl2), Ptr{Cchar}, (Ptr{SDL_Joystick},), ptr) |> unsafe_string

close(ptr::Ptr{SDL_Joystick}) = ccall(
    (:SDL_JoystickClose, libsdl2), Cvoid, (Ptr{SDL_Joystick},), ptr)

is_connected(ptr::Ptr{SDL_Joystick}) = ccall(
    (:SDL_JoystickGetAttached, libsdl2), Bool, (Ptr{SDL_Joystick},), ptr)

num_axes(ptr::Ptr{SDL_Joystick}) = ccall(
    (:SDL_JoystickNumAxes, libsdl2), Cint, (Ptr{SDL_Joystick},), ptr)

num_balls(ptr::Ptr{SDL_Joystick}) = ccall(
    (:SDL_JoystickNumBalls, libsdl2), Cint, (Ptr{SDL_Joystick},), ptr)

num_hats(ptr::Ptr{SDL_Joystick}) = ccall(
    (:SDL_JoystickNumHats, libsdl2), Cint, (Ptr{SDL_Joystick},), ptr)

num_buttons(ptr::Ptr{SDL_Joystick}) = ccall(
    (:SDL_JoystickNumButtons, libsdl2), Cint, (Ptr{SDL_Joystick},), ptr)

get_axis(ptr::Ptr{SDL_Joystick}, axis::Integer) = ccall(
    (:SDL_JoystickGetAxis, libsdl2), Int16, (Ptr{SDL_Joystick}, Cint), ptr, axis)

get_hat(ptr::Ptr{SDL_Joystick}, hat::Integer) = ccall(
    (:SDL_JoystickGetHat, libsdl2), UInt8, (Ptr{SDL_Joystick}, Cint), ptr, hat)

get_button(ptr::Ptr{SDL_Joystick}, button::Integer) = ccall(
    (:SDL_JoystickGetButton, libsdl2), UInt8, (Ptr{SDL_Joystick}, Cint), ptr, button)

function warn_not_connected(ptr::Ptr{SDL_Joystick})
    if !is_connected(ptr)
        id = get_instance_id(ptr)
        name = get_name(ptr)
        @warn "Can't retrieve data from device $name with instance ID $id, device is not connected"
    end
end


################################################################################
############################# ButtonChange #####################################

@enum ButtonChange begin
    ButtonUnchanged = 0
    ButtonPressed = 1
    ButtonReleased = 2
end

Base.convert(::Type{ButtonChange}, n::Integer) = ButtonChange(n)

function get_change(current::Bool, last::Bool)
    (current && !last) && return ButtonPressed
    (!current && last) && return ButtonReleased
    return ButtonUnchanged
end

################################################################################
############################# AbstractAxisSet ##################################

abstract type AbstractAxisSet{N, T <: Real} <: FieldVector{N, T} end

################################################################################
############################# AbstractButtonSet ################################

@kwdef struct ButtonData
    state::Bool = false
    change::ButtonChange = ButtonUnchanged
end

Base.zero(::Type{ButtonData}) = ButtonData()

is_pressed(data::ButtonData) = data.state === true
was_pressed(data::ButtonData) = data.change === ButtonPressed
was_released(data::ButtonData) = data.change === ButtonReleased

abstract type AbstractButtonSet{N, T <: Union{Bool, ButtonChange, ButtonData}} <: FieldVector{N, T} end

@kwdef struct HatButtons{T} <: AbstractButtonSet{4, T}
    up::T
    right::T
    down::T
    left::T
end

################################################################################
########################### AbstractJoystickData ###############################

abstract type AbstractJoystickData end

################################################################################
############################# Joystick #########################################

mutable struct Joystick{D <: AbstractJoystickData} <: InputDevice
    const ptr::Ptr{SDL_Joystick}
    cache::D
end

IODevices.get_default_mapping(::Joystick) = GenericInputMapping()

is_connected(joy::Joystick) = is_connected(joy.ptr)
get_instance_id(joy::Joystick) = get_instance_id(joy.ptr)
get_guid(joy::Joystick) = get_guid(joy.ptr)

const connected_joysticks = Joystick[]

function update_connected_joysticks()
    @assert init_joysticks() == 0
    update_joysticks()
    filter!(is_connected, connected_joysticks)

    for idx in 0:num_joysticks()-1
        #these must be reevaluated on each iteration
        connected_ids = map(get_instance_id, connected_joysticks)
        connected_guids = map(get_guid, connected_joysticks)
        latest_id = (isempty(connected_ids) ? -1 : maximum(connected_ids))

        if get_instance_id(idx) > latest_id #newly connected
            guid = get_guid(idx)
            name = get_name(idx)
            if guid ∈ connected_guids #multiple identical devices unsupported
                @warn "Can't add $name, an identical device is already connected"
                continue
            end
            if haskey(supported_joysticks, guid)
                ptr = open(idx)
                cache = supported_joysticks[guid]()
                joy = Joystick(ptr, cache)
                push!(connected_joysticks, joy)
            else
                @warn "Can't add $name, device not supported"
                continue
            end

        end
    end

    return connected_joysticks
end

function IODevices.get_data!(joystick::Joystick{D}) where {D}
    sleep(0.01)
    update_joysticks()
    warn_not_connected(joystick.ptr)
    data = D(joystick)
    joystick.cache = data
    return data
end

################################################################################
########################### Thrustmaster T16000M ###############################

@kwdef struct T16000MAxes{T} <: AbstractAxisSet{4, T}
    stick_x::T
    stick_y::T
    stick_z::T
    throttle::T
end

@kwdef struct T16000MButtons{T} <: AbstractButtonSet{16, T}
    button_0::T; button_1::T; button_2::T; button_3::T;
    button_4::T; button_5::T; button_6::T; button_7::T;
    button_8::T; button_9::T; button_10::T; button_11::T;
    button_12::T; button_13::T; button_14::T; button_15::T;
end

@kwdef struct T16000MData <: AbstractJoystickData
    axes::T16000MAxes{Float64} = zeros(T16000MAxes{Float64})
    buttons::T16000MButtons{ButtonData} = zeros(T16000MButtons{ButtonData})
    hat::HatButtons{ButtonData} = zeros(HatButtons{ButtonData})
end

function T16000MData(joystick::Joystick{T16000MData})

    (; ptr, cache) = joystick

    axes_data = sacollect(T16000MAxes{Float64}, get_axis(ptr, i)/32768.0 for i in 0:3) |> rescale

    button_state = sacollect(T16000MButtons{Bool}, get_button(ptr, i) for i in 0:15)
    button_state_last = sacollect(T16000MButtons{Bool}, b.state for b in cache.buttons)
    button_change = get_change.(button_state, button_state_last) |> T16000MButtons
    button_data = sacollect(T16000MButtons, ButtonData(state, change) for (state, change) in zip(button_state, button_change))

    hat_state = HatButtons(get_hat(ptr, 0) .& SDL_HAT_4POS .!= 0)
    hat_state_last = sacollect(HatButtons{Bool}, b.state for b in cache.hat)
    hat_change = get_change.(hat_state, hat_state_last) |> HatButtons
    hat_data = sacollect(HatButtons, ButtonData(state, change) for (state, change) in zip(hat_state, hat_change))

    return T16000MData(axes_data, button_data, hat_data)

end

function rescale(data::T16000MAxes{T}) where {T<:AbstractFloat}
    (; stick_x, stick_y, stick_z, throttle) = data
    T16000MAxes{T}(; stick_x, stick_y, stick_z, throttle = 0.5*(1 - throttle))
end

const T16000M = Joystick{T16000MData}


################################################################################
############################# Thrustmaster TWCS ################################

@kwdef struct TWCSAxes{T} <: AbstractAxisSet{8, T}
    mini_stick_x::T
    mini_stick_y::T
    throttle::T
    right_pedal::T #only available when TFRP connected
    left_pedal::T #only available when TFRP connected
    rocker::T
    rudder::T #only available when TFRP connected
    antenna::T
end

@kwdef struct TWCSButtons{T} <: AbstractButtonSet{6, T}
    button_0::T #orange button below the hats on the right side
    button_1::T #orange button on the front, leftmost
    button_2::T #orange button on the front, rightmost
    button_3::T #vertical lever up
    button_4::T #vertical lever down
    button_5::T #mini stick push button
end

@kwdef struct TWCSData <: AbstractJoystickData
    axes::TWCSAxes{Float64} = zeros(TWCSAxes{Float64})
    buttons::TWCSButtons{ButtonData} = zeros(TWCSButtons{ButtonData})
    hat_top::HatButtons{ButtonData} = zeros(HatButtons{ButtonData})
    hat_middle::HatButtons{ButtonData} = zeros(HatButtons{ButtonData})
    hat_bottom::HatButtons{ButtonData} = zeros(HatButtons{ButtonData})
end

function TWCSData(joystick::Joystick{TWCSData})

    (; ptr, cache) = joystick

    axes = sacollect(TWCSAxes{Float64}, get_axis(ptr, i)/32768.0 for i in 0:7) |> rescale

    button_state = sacollect(TWCSButtons{Bool}, get_button(ptr, i) for i in 0:5)
    button_state_last = sacollect(TWCSButtons{Bool}, b.state for b in cache.buttons)
    button_change = get_change.(button_state, button_state_last) |> TWCSButtons
    buttons = sacollect(TWCSButtons, ButtonData(state, change) for (state, change) in zip(button_state, button_change))

    hat_top_state = HatButtons(get_hat(ptr, 0) .& SDL_HAT_4POS .!= 0)
    hat_top_state_last = sacollect(HatButtons{Bool}, b.state for b in cache.hat_top)
    hat_top_change = get_change.(hat_top_state, hat_top_state_last) |> HatButtons
    hat_top = sacollect(HatButtons, ButtonData(state, change) for (state, change) in zip(hat_top_state, hat_top_change))

    #buttons 6 to 9 correspond to the middle hat
    hat_middle_state = sacollect(HatButtons{Bool}, get_button(ptr, i) for i in 6:9)
    hat_middle_state_last = sacollect(HatButtons{Bool}, b.state for b in cache.hat_middle)
    hat_middle_change = get_change.(hat_middle_state, hat_middle_state_last) |> HatButtons
    hat_middle = sacollect(HatButtons, ButtonData(state, change) for (state, change) in zip(hat_middle_state, hat_middle_change))

    #buttons 10 to 13 correspond to the bottom hat
    hat_bottom_state = sacollect(HatButtons{Bool}, get_button(ptr, i) for i in 10:13)
    hat_bottom_state_last = sacollect(HatButtons{Bool}, b.state for b in cache.hat_bottom)
    hat_bottom_change = get_change.(hat_bottom_state, hat_bottom_state_last) |> HatButtons
    hat_bottom = sacollect(HatButtons, ButtonData(state, change) for (state, change) in zip(hat_bottom_state, hat_bottom_change))

    return TWCSData(; axes, buttons, hat_top, hat_middle, hat_bottom)

end

function rescale(data::TWCSAxes{T}) where {T<:AbstractFloat}
    (; mini_stick_x, mini_stick_y, throttle, right_pedal, left_pedal,
        rocker, rudder, antenna) = data
    TWCSAxes{T}(; mini_stick_x, mini_stick_y, throttle = 0.5*(1 - throttle),
                right_pedal, left_pedal, rocker, rudder, antenna = 0.5*(antenna + 1))
end

const TWCS = Joystick{TWCSData}



################################################################################
######################## VKBSim Gladiator NXT Evo ##############################

@kwdef struct GladiatorNXTEvoAxes{T} <: AbstractAxisSet{6, T}
    stick_x::T = 0.0
    stick_y::T = 0.0
    throttle::T = 0.0
    analog_hat_x::T = 0.0
    analog_hat_y::T = 0.0
    stick_z::T = 0.0
end

#when queried with GetJoystickButtons, Gladiator NXT Evo returns an array of
#132 values, from which the first 29 actually correspond to physical buttons. we
#only keep those here
@kwdef struct GladiatorNXTEvoButtons{T} <: AbstractButtonSet{29, T}
    fire_half::T = 0
    fire_full::T = 0
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


@kwdef struct GladiatorNXTEvoData <: AbstractJoystickData
    axes::GladiatorNXTEvoAxes{Float64} = zeros(GladiatorNXTEvoAxes{Float64})
    buttons::GladiatorNXTEvoButtons{ButtonData} = zeros(GladiatorNXTEvoButtons{ButtonData})
    hat::HatButtons{ButtonData} = zeros(HatButtons{ButtonData})
end

function GladiatorNXTEvoData(joystick::Joystick{GladiatorNXTEvoData})

    (; ptr, cache) = joystick

    axes_data = sacollect(GladiatorNXTEvoAxes{Float64}, get_axis(ptr, i)/32768.0 for i in 0:5) |> rescale

    button_state = sacollect(GladiatorNXTEvoButtons{Bool}, get_button(ptr, i) for i in 0:28)
    button_state_last = sacollect(GladiatorNXTEvoButtons{Bool}, b.state for b in cache.buttons)
    button_change = get_change.(button_state, button_state_last) |> GladiatorNXTEvoButtons
    button_data = sacollect(GladiatorNXTEvoButtons, ButtonData(state, change) for (state, change) in zip(button_state, button_change))

    hat_state = HatButtons(get_hat(ptr, 0) .& SDL_HAT_4POS .!= 0)
    hat_state_last = sacollect(HatButtons{Bool}, b.state for b in cache.hat)
    hat_change = get_change.(hat_state, hat_state_last) |> HatButtons
    hat_data = sacollect(HatButtons, ButtonData(state, change) for (state, change) in zip(hat_state, hat_change))

    return GladiatorNXTEvoData(axes_data, button_data, hat_data)

end

function rescale(data::GladiatorNXTEvoAxes{T}) where {T<:AbstractFloat}
    (; stick_x, stick_y, throttle, analog_hat_x, analog_hat_y, stick_z) = data
    GladiatorNXTEvoAxes{T}(; stick_x, stick_y, throttle = 0.5*(1 - throttle),
                            analog_hat_x, analog_hat_y, stick_z)
end

const GladiatorNXTEvo = Joystick{GladiatorNXTEvoData}


################################################################################
############################# Supported Joysticks ##############################

const supported_joysticks = Dict{SDL_JoystickGUID, Type}(
    SDL_JoystickGUID((0x03, 0x00, 0x00, 0x00, 0x4f, 0x04, 0x00, 0x00, 0x0a, 0xb1, 0x00, 0x00, 0x00, 0x05, 0x00, 0x00)) => T16000MData, #macOS
    SDL_JoystickGUID((0x03, 0x00, 0x00, 0x00, 0x4f, 0x04, 0x00, 0x00, 0x0a, 0xb1, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00)) => T16000MData, #Windows
    SDL_JoystickGUID((0x03, 0x00, 0x00, 0x00, 0x4f, 0x04, 0x00, 0x00, 0x87, 0xb6, 0x00, 0x00, 0x10, 0x01, 0x00, 0x00)) => TWCSData, #macOS
    SDL_JoystickGUID((0x03, 0x00, 0x00, 0x00, 0x4f, 0x04, 0x00, 0x00, 0x87, 0xb6, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00)) => TWCSData, #Windows
    SDL_JoystickGUID((0x03, 0x00, 0x00, 0x00, 0x1d, 0x23, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x4a, 0x21, 0x00, 0x00)) => GladiatorNXTEvoData, #macOS
    SDL_JoystickGUID((0x03, 0x00, 0x00, 0x00, 0x1d, 0x23, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00)) => GladiatorNXTEvoData, #Windows
)


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