module TestIODevices

using UnPack
using CImGui
using CImGui.CSyntax

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

Base.@kwdef mutable struct PIDashboardState
    reset::Bool = false
    sat_enable::Bool = true
end

function GUI.draw_dashboard!(state::PIDashboardState, data::PICompensatorY{1})
    @show data
    # @show state
    begin
        CImGui.Begin("PICompensator{1} Dashboard")
        CImGui.Text("Inputs:")
            # reset_ref = Ref(true)
            # CImGui.Checkbox("Reset", reset_ref) #needs a Ref(Bool)
            # state.reset = reset_ref[]
            @c CImGui.Checkbox("Reset", &(state.reset))
            # @c CImGui.Checkbox("Enable Saturation", &(state.sat_enable))
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
    @show u.reset
end

function test_input_pi()

    sys = PICompensator{1}(k_i = 0.5) |> System
    sim = Simulation(sys; t_end = 30, dt = 1)
    # joystick = get_connected_joysticks()[1]
    # joystick_interface = Sim.attach_io!(sim, joystick; mapping = TestMapping())
    # dashboard = Dashboard{PIDashboardState}()
    # dashboard_interface = Sim.attach_io!(sim, dashboard; mapping = TestMapping())

    # @sync begin
        # errormonitor(Threads.@spawn IODevices.start!(joystick_interface))
        # # Threads.@spawn Sim.start!(dashboard_interface)
        # errormonitor(Threads.@spawn Sim.run!(sim; rate = 1, verbose = true))

    # end
    # return dashboard_interface, sim

    # @sync begin
        # Threads.@spawn IODevices._start!(joystick_interface)
        # Threads.@spawn IODevices._start!(dashboard_interface)
        Threads.@spawn Sim.run!(sim; rate = 1, verbose = true)
    # end

    # return joystick_interface, sim
end

end #module