module JNew

using StaticArrays
using UnPack
using SDL2_jll

using ..IODevices
using ..Types

# export JoystickData, Joystick
# export get_connected_joysticks
# export get_axis_value, exp_axis_curve
# export get_button_state, is_pressed
# export get_button_change, was_pressed, was_released

# export XBoxController, XBoxControllerData
# export T16000M, T16000MData
# export GladiatorNXTEvo, GladiatorNXTEvoData


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

abstract type AbstractAxisData{N} <: FieldVector{N, Float32} end

rescale(axes::AbstractAxisData) = axes

################################################################################
############################# AbstractButtonSet ################################

@kwdef struct ButtonData
    state::Bool = false
    change::ButtonChange = ButtonUnchanged
end

is_active(data::ButtonData) = data.state === true
was_pressed(data::ButtonData) = data.change === ButtonPressed
was_released(data::ButtonData) = data.change === ButtonReleased

abstract type AbstractButtonSetData{N} <: FieldVector{N, ButtonData} end


################################################################################
########################### AbstractJoystickData ###############################

abstract type AbstractJoystickData end

struct NoJoystickData <: AbstractJoystickData end


################################################################################
############################ SDL Initialization ################################

const SDL_INIT_JOYSTICK = Cuint(0x00000200)

const SDL_JoystickInstanceID = Int32

struct SDL_JoystickGUID
    data::NTuple{16, UInt8}
end

#initialize SDL2's joystick functionality
init_joysticks() = ccall((:SDL_Init, libsdl2), Cint, (UInt32,), SDL_INIT_JOYSTICK)

#update connected joysticks status and data
update_joysticks() = ccall((:SDL_JoystickUpdate, libsdl2), Cvoid, ())

num_joysticks() = ccall((:SDL_NumJoysticks, libsdl2), Cint, ())

function get_instance_id(device_index::Integer)
    ccall((:SDL_JoystickGetDeviceInstanceID, libsdl2), SDL_JoystickInstanceID, (Cint,), device_index)
end

function get_guid(device_index::Integer)
    ccall((:SDL_JoystickGetDeviceGUID, libsdl2), SDL_JoystickGUID, (Cint,), device_index)
end

function get_name(device_index::Integer)
    ccall((:SDL_JoystickNameForIndex, libsdl2), Ptr{Cchar}, (Cint,), device_index) |> unsafe_string
end


################################################################################
############################## SDL_Joystick ####################################

mutable struct SDL_Joystick end

function open(device_index::Integer)::Ptr{SDL_Joystick}
   ccall((:SDL_JoystickOpen, libsdl2), Ptr{SDL_Joystick}, (Cint,), device_index)
end

function get_instance_id(joystick::Ptr{SDL_Joystick})
    ccall((:SDL_JoystickInstanceID, libsdl2), SDL_JoystickInstanceID, (Ptr{SDL_Joystick},), joystick)
end

function get_guid(joystick::Ptr{SDL_Joystick})
    ccall((:SDL_JoystickGetGUID, libsdl2), SDL_JoystickGUID, (Ptr{SDL_Joystick},), joystick)
end

function get_name(joystick::Ptr{SDL_Joystick})
    ccall((:SDL_JoystickName, libsdl2), Ptr{Cchar}, (Ptr{SDL_Joystick},), joystick) |> unsafe_string
end

function close(joystick::Ptr{SDL_Joystick})
    ccall((:SDL_JoystickClose, libsdl2), Cvoid, (Ptr{SDL_Joystick},), joystick)
end

function is_connected(joystick::Ptr{SDL_Joystick})
    ccall((:SDL_JoystickGetAttached, libsdl2), Bool, (Ptr{SDL_Joystick},), joystick)
end

#observations:
#0) if a joystick is disconnected, update_joysticks() must be called for
#is_connected to return false. if we close it explicitly, it also returns false
#1) we unplug a joystick, call update_joysticks, plug it again and call
#update_joysticks. even if the second time its device_index is the same, its
#instance id will have increased. so an instance id is never reused even if it
#becomes available again, unless we restart the REPL


################################################################################
################################################################################


mutable struct Joystick{D <: AbstractJoystickData} <: InputDevice
    const ptr::Ptr{SDL_Joystick}
    data::D
end

get_instance_id(joy::Joystick) = get_instance_id(joy.ptr)
get_guid(joy::Joystick) = get_guid(joy.ptr)
is_connected(joy::Joystick) = is_connected(joy.ptr)

const connected_joysticks = Joystick[]

function update_connected_joysticks()
    @assert init_joysticks() == 0
    update_joysticks()
    filter!(is_connected, connected_joysticks)

    for device_index in 0:num_joysticks()-1
        #need to recompute these after adding any joystick
        connected_ids = map(get_instance_id, connected_joysticks)
        connected_guids = map(get_guid, connected_joysticks)
        latest_id = (isempty(connected_ids) ? -1 : maximum(connected_ids))

        if get_instance_id(device_index) > latest_id #newly connected
            if get_guid(device_index) ∈ connected_guids #already have an identical device
                @warn "Can't add $(get_name(device_index)), an identical device is already connected"
                continue
            end
            ptr = open(device_index)
            #check if it's supported and instantiate the appropriate Joystick subtype
            println("Adding $(get_name(ptr))")
            joy = Joystick(ptr, NoJoystickData())
            push!(connected_joysticks, joy)
        end
    end

    return connected_joysticks
end

function add_joystick()

end

#we could create three different IOMappings: TWCS, T16000, T16000M_TWCS. What we
#map when both are plugged is different

################################################################################
################################ T16000M #######################################

@kwdef struct T16000MAxisData <: AbstractAxisData{4}
    stick_x::Float32 = 0.0
    stick_y::Float32 = 0.0
    stick_z::Float32 = 0.0
    throttle::Float32 = 0.0
end

function rescale(data::T16000MAxisData)
    @unpack stick_x, stick_y, stick_z, throttle = data
    T16000MAxisData(; stick_x, stick_y, stick_z, throttle = 0.5*(1 - throttle))
end

@kwdef struct T16000MButtonSetData <: AbstractButtonSetData{20}
    button_0::ButtonData = ButtonData()
    button_1::ButtonData = ButtonData()
    button_2::ButtonData = ButtonData()
    button_3::ButtonData = ButtonData()
    button_4::ButtonData = ButtonData()
    button_5::ButtonData = ButtonData()
    button_6::ButtonData = ButtonData()
    button_7::ButtonData = ButtonData()
    button_8::ButtonData = ButtonData()
    button_9::ButtonData = ButtonData()
    button_10::ButtonData = ButtonData()
    button_11::ButtonData = ButtonData()
    button_12::ButtonData = ButtonData()
    button_13::ButtonData = ButtonData()
    button_14::ButtonData = ButtonData()
    button_15::ButtonData = ButtonData()
end

@kwdef struct T16000MHatData <: AbstractButtonSetData{4}
    up::ButtonData = ButtonData()
    right::ButtonData = ButtonData()
    down::ButtonData = ButtonData()
    left::ButtonData = ButtonData()
end

@kwdef struct T16000MData <: AbstractJoystickData
    axes::T16000MAxisData = T16000MAxisData()
    buttons::T16000MButtonSetData = T16000MButtonSetData()
    hat::T16000MHatData = T16000MHatData()
end

#NO: UN HAT TENGO QUE TRATARLO COMO BOTON DISCRETO. PARA UN TRIM ME INTERESA EL
#CAMBIO, NO EL VALOR INSTANTANEO. Voy a tener que definir las 4 posiciones
#posibles como botones, teniendo en cuenta que pueden ser ciertas mas de una a
#la vez! como methods auxiliares puedo definir is_center, was_center o lo que
#haga falta. pero vamos



# function IODevices.get_data!(joystick::Joystick{T}) where {T <: JoystickData{A, BS, BC}} where {A, BS, BC}

#     #from the axis values returned by GetJoystickAxes we only keep as many
#     #as defined by our AbstractButtonSet
#     axes = view(GetJoystickAxes(joystick.slot), 1:length(A)) |> A |> rescale

#     #from the button values returned by GetJoystickButtons we only keep as many
#     #as defined by our AbstractButtonSet
#     button_state = view(GetJoystickButtons(joystick.slot), 1:length(BS)) |> BS

#     button_state_last = joystick.cache.button_state

#     button_change = map(button_state, button_state_last) do current, last
#         (current && !last) && return ButtonPressed
#         (!current && last) && return ButtonReleased
#         return ButtonUnchanged
#     end |> BC

#     data = JoystickData(axes, button_state, button_change)

#     joystick.cache = data

#     return data
# end


# ################################################################################
# ############################ XBox Controller ###################################

# @kwdef struct XBoxControllerAxes <: AbstractAxisSet{6}
#     left_stick_x::Float32 = 0.0
#     left_stick_y::Float32 = 0.0
#     right_stick_x::Float32 = 0.0
#     right_stick_y::Float32 = 0.0
#     left_trigger::Float32 = 0.0
#     right_trigger::Float32 = 0.0
# end

# @kwdef struct XBoxControllerButtons{T} <: AbstractButtonSet{14, T}
#     button_A::T = 0
#     button_B::T = 0
#     button_X::T = 0
#     button_Y::T = 0
#     left_bumper::T = 0
#     right_bumper::T = 0
#     view::T = 0
#     menu::T = 0
#     left_stick_press::T = 0
#     right_stick_press::T = 0
#     dpad_up::T = 0
#     dpad_right::T = 0
#     dpad_down::T = 0
#     dpad_left::T = 0
# end

# const XBoxControllerData = JoystickData{XBoxControllerAxes,
#                                         XBoxControllerButtons{Bool},
#                                         XBoxControllerButtons{ButtonChange}}

# const XBoxController = Joystick{XBoxControllerData}

# get_XBoxControllerData() = JoystickData(XBoxControllerAxes(),
#                                         XBoxControllerButtons{Bool}(),
#                                         XBoxControllerButtons{ButtonChange}())

# get_XBoxController(slot::JoystickSlot) = Joystick(slot, get_XBoxControllerData())

# function rescale(axes::XBoxControllerAxes)
#     @unpack left_stick_x, left_stick_y, right_stick_x, right_stick_y,
#             left_trigger, right_trigger = axes
#     XBoxControllerAxes(; left_stick_x, left_stick_y, right_stick_x, right_stick_y,
#             left_trigger = 0.5*(1 + left_trigger),
#             right_trigger = 0.5*(1 + right_trigger))
# end


# ################################################################################
# ######################## VKBSim Gladiator NXT Evo ##############################

# @kwdef struct GladiatorNXTEvoAxes <: AbstractAxisSet{6}
#     stick_x::Float32 = 0.0
#     stick_y::Float32 = 0.0
#     throttle::Float32 = 0.0
#     analog_hat_x::Float32 = 0.0
#     analog_hat_y::Float32 = 0.0
#     stick_z::Float32 = 0.0
# end

# #when queried with GetJoystickButtons, the Gladiator NXT Evo returns an array of
# #132 values, from which the first 29 actually correspond to physical buttons. we
# #only keep those here
# @kwdef struct GladiatorNXTEvoButtons{T} <: AbstractButtonSet{29, T}
#     red_trigger_half::T = 0
#     red_trigger_full::T = 0
#     A2::T = 0
#     B1::T = 0
#     D1::T = 0
#     A3_up::T = 0
#     A3_right::T = 0
#     A3_down::T = 0
#     A3_left::T = 0
#     A3_press::T = 0
#     A4_up::T = 0
#     A4_right::T = 0
#     A4_down::T = 0
#     A4_left::T = 0
#     A4_press::T = 0
#     C1_up::T = 0
#     C1_right::T = 0
#     C1_down::T = 0
#     C1_left::T = 0
#     C1_press::T = 0
#     black_trigger_up::T = 0
#     black_trigger_down::T = 0
#     encoder_up::T = 0
#     encoder_down::T = 0
#     switch_up::T = 0
#     switch_down::T = 0
#     F1::T = 0
#     F2::T = 0
#     F3::T = 0
# end

# const GladiatorNXTEvoData = JoystickData{GladiatorNXTEvoAxes,
#                                         GladiatorNXTEvoButtons{Bool},
#                                         GladiatorNXTEvoButtons{ButtonChange}}

# const GladiatorNXTEvo = Joystick{GladiatorNXTEvoData}

# get_GladiatorNXTEvoData() = JoystickData(GladiatorNXTEvoAxes(),
#                                         GladiatorNXTEvoButtons{Bool}(),
#                                         GladiatorNXTEvoButtons{ButtonChange}())

# get_GladiatorNXTEvo(slot::JoystickSlot) = Joystick(slot, get_GladiatorNXTEvoData())

# function rescale(axes::GladiatorNXTEvoAxes)
#     @unpack stick_x, stick_y, throttle, analog_hat_x, analog_hat_y, stick_z = axes
#     GladiatorNXTEvoAxes(; stick_x, stick_y, throttle = 0.5*(1 - throttle),
#                         analog_hat_x, analog_hat_y, stick_z)
# end


# ################################################################################
# ############################# Initialization ###################################

# const connected_slots = Dict{JoystickSlot, Joystick}()

# #specify a Joystick constructor for each joystick model
# const supported_joysticks = Dict{String, Function}(
#     "Xbox Controller" => get_XBoxController, #XBox 360 Controller
#     "T.16000M" => get_T16000M, #Thrustmaster T16000M
#     "VKBsim Gladiator EVO  R" => get_GladiatorNXTEvo, #VKBSim Gladiator NXT Evo
# )

# function refresh_joysticks()

#     #calling PollEvents refreshes all the joystick slots, so JoystickPresent
#     #then returns their updated states. no need to explicitly handle the
#     #addition or removal of joysticks directly from a joystick_callback (which
#     #would also require a PollEvents call to work)
#     GLFW.PollEvents()
#     empty!(connected_slots)
#     for slot ∈ instances(JoystickSlot)
#         JoystickPresent(slot) && add_joystick(slot)
#     end

#     return connected_slots
# end

# get_connected_joysticks() = Tuple(values(refresh_joysticks()))

# function add_joystick(slot::JoystickSlot)

#     if !JoystickPresent(slot)
#         @warn("Failed to add joystick at slot $slot: not found")
#         return
#     end

#     joystick_model = GetJoystickName(slot) |> strip

#     if !(joystick_model ∈ keys(supported_joysticks))
#         @warn("Failed to add $joystick_model at slot $slot: not supported")
#         return
#     end

#     connected_slots[slot] = supported_joysticks[joystick_model](slot)

#     return connected_slots[slot]

# end

# function is_connected(joystick::Joystick)
#     slot = joystick.slot
#     #this joystick's slot must be listed in connected_slots, and the
#     #corresponding key must actually point to this joystick instance
#     return (slot ∈ keys(connected_slots) && connected_slots[slot] === joystick)
# end


# ################################################################################
# ############################ Misc Stuff ########################################

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



end #module