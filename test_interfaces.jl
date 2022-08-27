module TestIODevices

using Flight
using Flight.Essentials: PICompensator, PICompensatorU, PICompensatorY

struct TestMapping <: InputMapping end

function IODevices.assign!(u::PICompensatorU{1},
                            joystick::Joystick{XBoxControllerID},
                            ::TestMapping)

    u.input .= get_axis_value(joystick, :right_analog_y)
    u.reset .= was_released(joystick, :button_A)
    u.sat_enable .⊻= was_released(joystick, :button_Y)

end

Base.@kwdef struct PIDashboardState
    reset::Cuchar = false
    sat_enable::Cuchar = true
end

function GUI.draw_dashboard!(state::PIDashboardState, y::PICompensatorY{1})
    begin
        CImGui.Begin("PICompensator{1} Dashboard")
        CImGui.Text("Inputs:")
            @c CImGui.Checkbox("Reset", &(state.reset))
            @c CImGui.Checkbox("Enable Saturation", &(state.sat_enable))
        CImGui.Text("Outputs:")
            CImGui.Text("Output = $(data.out)")
            CImGui.Text("Saturation Status = $(data.sat_status)")
        CImGui.End()
    end
end

function IODevices.assign!(u::PICompensatorU{1}, dashboard::Dashboard, ::TestMapping)
    @unpack reset, sat_enable = dashboard.state
    u.reset .= reset
    u.sat_enable .= sat_enable
end

function test_input_pi()

    sys = PICompensator{1}(k_i = 0.5) |> System
    sim = Simulation(sys; t_end = 10, dt = 0.02)
    joystick = get_connected_joysticks()[1]
    joystick_interface = Sim.attach_io!(sim, joystick)
    dashboard = Dashboard{PIDashboardState}()
    dashboard_interface = Sim.attach_io!(sim, dashboard)

    @sync begin
        Threads.@spawn Sim.start!(joystick_interface)
        Threads.@spawn Sim.start!(dashboard_interface)
        Threads.@spawn Sim.run!(sim; rate = 1, verbose = true)
    end


    return sim
end

end #module