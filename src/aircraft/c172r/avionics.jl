using Flight.Engine.Systems
using Flight.Engine.IODevices
using Flight.Engine.Joysticks
using Flight.Engine.GUI
using Flight.Engine.Utils: Ranged

using Printf
using CImGui, CImGui.CSyntax, CImGui.CSyntax.CStatic

################################################################################
########################### MechanicalControls #################################

struct MechanicalControls <: AbstractAvionics end

#elevator↑ (stick forward) -> e↑ -> δe↑ -> trailing edge down -> Cm↓ -> pitch down
#aileron↑ (stick right) -> a↑ -> δa↑ -> left trailing edge down, right up -> Cl↓ -> roll right
#rudder↑ (right pedal forward) -> r↓ -> δr↓ -> rudder trailing edge right -> Cn↑ -> yaw right
#rudder↑ (right pedal forward) -> nose wheel steering right -> yaw right
#flaps↑ -> δf↑ -> flap trailing edge down -> CL↑
Base.@kwdef mutable struct MechanicalControlsU
    throttle::Ranged{Float64, 0, 1} = 0.0
    aileron_trim::Ranged{Float64, -1, 1} = 0.0
    aileron_offset::Ranged{Float64, -1, 1} = 0.0 #incremental command, for input devices
    elevator_trim::Ranged{Float64, -1, 1} = 0.0
    elevator_offset::Ranged{Float64, -1, 1} = 0.0 #incremental command, for input devices
    rudder_trim::Ranged{Float64, -1, 1} = 0.0
    rudder_offset::Ranged{Float64, -1, 1} = 0.0 #incremental command, for input devices
    brake_left::Ranged{Float64, 0, 1} = 0.0
    brake_right::Ranged{Float64, 0, 1} = 0.0
    flaps::Ranged{Float64, 0, 1} = 0.0
    mixture::Ranged{Float64, 0, 1} = 0.5
    eng_start::Bool = false
    eng_stop::Bool = false
end

Base.@kwdef struct MechanicalControlsY
    throttle::Float64 = 0.0
    aileron_trim::Float64 = 0.0
    aileron_offset::Float64 = 0.0
    elevator_trim::Float64 = 0.0
    elevator_offset::Float64 = 0.0
    rudder_trim::Float64 = 0.0
    rudder_offset::Float64 = 0.0
    brake_left::Float64 = 0.0
    brake_right::Float64 = 0.0
    flaps::Float64 = 0.0
    mixture::Float64 = 0.5
    eng_start::Bool = false
    eng_stop::Bool = false
end

Systems.init(::SystemU, ::MechanicalControls) = MechanicalControlsU()
Systems.init(::SystemY, ::MechanicalControls) = MechanicalControlsY()


########################### Update Methods #####################################

function Systems.f_ode!(avionics::System{MechanicalControls}, ::System{<:Airframe},
                ::KinematicData, ::AirData, ::System{<:AbstractTerrain})

    #MechanicalControls has no internal dynamics, just input-output feedthrough
    @unpack throttle, aileron_trim, aileron_offset, elevator_trim, elevator_offset,
            rudder_trim, rudder_offset, brake_left, brake_right, flaps, mixture,
            eng_start, eng_stop = avionics.u

    avionics.y = MechanicalControlsY(;
            throttle, aileron_trim, aileron_offset, elevator_trim, elevator_offset,
            rudder_trim, rudder_offset, brake_left, brake_right, flaps, mixture,
            eng_start, eng_stop)

end

#no digital components or state machines in MechanicalControls
@inline Systems.f_step!(::System{MechanicalControls}, ::System{<:Airframe}, ::KinematicSystem) = false
@inline Systems.f_disc!(::System{MechanicalControls}, ::System{<:Airframe}, ::KinematicSystem, Δt) = false


function map_controls!(airframe::System{<:Airframe}, avionics::System{MechanicalControls})

    @unpack throttle, aileron_trim, aileron_offset, elevator_trim, elevator_offset,
            rudder_trim, rudder_offset, brake_left, brake_right, flaps, mixture,
            eng_start, eng_stop = avionics.u

    @unpack aero, pwp, ldg = airframe

    pwp.u.engine.start = eng_start
    pwp.u.engine.stop = eng_stop
    pwp.u.engine.thr = throttle
    pwp.u.engine.mix = mixture
    ldg.u.nose.steering[] = (rudder_trim + rudder_offset) #rudder↑ (right pedal forward) -> nose wheel steering right
    ldg.u.left.braking[] = brake_left
    ldg.u.right.braking[] = brake_right
    aero.u.e = (elevator_trim + elevator_offset) #elevator↑ (stick forward) -> e↑ -> pitch down
    aero.u.a = (aileron_trim + aileron_offset) #aileron↑ (stick right) -> a↑ -> roll right
    aero.u.r = -(rudder_trim + rudder_offset) #rudder↑ (right pedal forward) -> r↓ -> yaw right
    aero.u.f = flaps #flaps↑ -> δf↑

    return nothing
end


############################ Joystick Mappings #################################

elevator_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
aileron_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)
rudder_curve(x) = exp_axis_curve(x, strength = 1.5, deadzone = 0.05)
brake_curve(x) = exp_axis_curve(x, strength = 1, deadzone = 0.05)

function IODevices.assign!(u::MechanicalControlsU,
            joystick::Joystick{XBoxControllerID}, ::DefaultMapping)

    u.aileron_offset = get_axis_value(joystick, :right_analog_x) |> aileron_curve
    u.elevator_offset = -get_axis_value(joystick, :right_analog_y) |> elevator_curve
    u.rudder_offset = get_axis_value(joystick, :left_analog_x) |> rudder_curve
    u.brake_left = get_axis_value(joystick, :left_trigger) |> brake_curve
    u.brake_right = get_axis_value(joystick, :right_trigger) |> brake_curve

    u.aileron_trim -= 0.01 * was_released(joystick, :dpad_left)
    u.aileron_trim += 0.01 * was_released(joystick, :dpad_right)
    u.elevator_trim -= 0.01 * was_released(joystick, :dpad_down)
    u.elevator_trim += 0.01 * was_released(joystick, :dpad_up)

    u.throttle += 0.1 * was_released(joystick, :button_Y)
    u.throttle -= 0.1 * was_released(joystick, :button_A)

    u.flaps += 0.3333 * was_released(joystick, :right_bumper)
    u.flaps -= 0.3333 * was_released(joystick, :left_bumper)

end


################################## GUI #########################################

struct InputState{T}
    enabled::Bool
    value::T
end

function Base.:(|)(target::Any, state::InputState{T}) where {T}
    state.enabled ? state.value : convert(T, target)
end

macro input_button(target, label, hue)
    return (quote
        # enable_flag = @cstatic check=false @c CImGui.Checkbox("Enable##"*($label), &check)
        # CImGui.SameLine()
        enable_flag = true
        CImGui.PushStyleColor(CImGui.ImGuiCol_Button, CImGui.HSV($hue, 0.6, 0.6))
        CImGui.PushStyleColor(CImGui.ImGuiCol_ButtonHovered, CImGui.HSV($hue, 0.7, 0.7))
        CImGui.PushStyleColor(CImGui.ImGuiCol_ButtonActive, CImGui.HSV($hue, 0.8, 0.8))
        CImGui.Button(($label))
        CImGui.PopStyleColor(3)
        enable_flag && ($(esc(target)) = CImGui.IsItemActive())
    end)
end

macro input_slider(target, label, lower_bound, upper_bound, default)
    return (quote
        enable_flag = @cstatic check=false @c CImGui.Checkbox($label, &check)
        CImGui.SameLine()
        slider_value = @cstatic f=Cfloat($default) @c CImGui.SliderFloat("##"*($label), &f, $lower_bound, $upper_bound)
        enable_flag && ($(esc(target)) = slider_value)
    end)
end

macro running_plot(source, label, lower_bound, upper_bound, default)
    window_height = 60
    return (quote
        @cstatic values=fill(Cfloat($default),90) values_offset=Cint(0) begin
            values[values_offset+1] = $(esc(source))
            values_offset = (values_offset+1) % length(values)
            CImGui.PlotLines(string($(esc(source)) |> Float32), values, length(values), values_offset,
                             $label, $lower_bound, $upper_bound, (0, $window_height))
        end
    end)
end


function GUI.draw!(sys::System{<:MechanicalControls}, label::String = "Cessna 172R Mechanical Controls")

    u = sys.u

    CImGui.Begin(label)

    CImGui.PushItemWidth(-60)

    @input_button(u.eng_start, "Engine Start", 0.4)
    @input_button(u.eng_stop, "Engine Stop", 0.0)
    @input_slider(u.throttle, "Throttle", 0, 1, 1)
    @input_slider(u.mixture, "Mixture", 0, 1, 0.5)
    @input_slider(u.brake_left, "Left Brake", 0, 1, 0.0)
    @input_slider(u.brake_right, "Right Brake", 0, 1, 0.0)
    @input_slider(u.aileron_offset, "Aileron Offset", -1, 1, 0.0)
    @input_slider(u.elevator_offset, "Elevator Offset", -1, 1, 0.0)
    @input_slider(u.rudder_offset, "Rudder Offset", -1, 1, 0.0)
    @input_slider(u.flaps, "Flaps", 0, 1, 0.0)
    @input_slider(u.aileron_trim, "Aileron Trim", -1, 1, 0.0)
    @input_slider(u.elevator_trim, "Elevator Trim", -1, 1, 0.0)
    @input_slider(u.rudder_trim, "Rudder Trim", -1, 1, 0.0)

    CImGui.PopItemWidth()

    CImGui.End()

    GUI.draw(sys, label)

end

function GUI.draw(sys::System{<:MechanicalControls}, label::String = "Cessna 172R Mechanical Controls")

    y = sys.y

    CImGui.Begin(label)

    CImGui.PushItemWidth(-60)

    @running_plot(y.throttle, "Throttle", 0, 1, 0.0)
    @running_plot(y.mixture, "Mixture", 0, 1, 0.5)
    @running_plot(y.brake_left, "Left Brake", 0, 1, 0.0)
    @running_plot(y.brake_right, "Right Brake", 0, 1, 0.0)
    @running_plot(y.aileron_offset, "Aileron Offset", -1, 1, 0.0)
    @running_plot(y.elevator_offset, "Elevator Offset", -1, 1, 0.0)
    @running_plot(y.rudder_offset, "Rudder Offset", -1, 1, 0.0)
    @running_plot(y.flaps, "Flaps", 0, 1, 0.0)
    @running_plot(y.aileron_trim, "Aileron Trim", -1, 1, 0.0)
    @running_plot(y.elevator_trim, "Elevator Trim", -1, 1, 0.0)
    @running_plot(y.rudder_trim, "Rudder Trim", -1, 1, 0.0)

    CImGui.PopItemWidth()

    CImGui.End()

end
