module TestJoysticks

using UnPack
using StaticArrays

using Flight.FlightCore

export demo_joysticks

struct TestSystem <: SystemDefinition end

Systems.U(::TestSystem) = zeros(2)
Systems.Y(::TestSystem) = zeros(SVector{2,Float64})

@no_ode TestSystem
@no_disc TestSystem

Systems.f_step!(sys::System{<:TestSystem}) = (sys.y = SVector{2}(sys.u))

function GUI.draw!(sys::System{TestSystem})
    CImGui.Begin("Joystick Visualization")

    # Get window position and dimensions
    window_pos = CImGui.GetWindowPos()
    window_width = CImGui.GetWindowWidth()

    # Calculate center using the window's actual position on screen
    center_x = Float32(window_pos.x + window_width * 0.5)

    # Get the current position values (clamped between -1.0 and 1.0)
    x_pos = clamp(sys.y[1], -1.0, 1.0)
    y_pos = clamp(sys.y[2], -1.0, 1.0)

    # Draw horizontal bar
    CImGui.Text("X: $(round(x_pos, digits=3))")
    bar_width = Float32(200.0)
    bar_height = Float32(20.0)
    bar_center_x = center_x
    bar_left = Float32(bar_center_x - bar_width/2)
    cursor_y = Float32(CImGui.GetCursorScreenPos().y)

    # Background for horizontal bar
    draw_list = CImGui.GetWindowDrawList()
    CImGui.AddRectFilled(draw_list,
                        (bar_left, cursor_y),
                        (bar_left + bar_width, cursor_y + bar_height),
                        CImGui.ColorConvertFloat4ToU32(CImGui.HSV(0.0, 0.0, 0.2)))

    # Actual horizontal bar
    fill_width = Float32((x_pos + 1.0) * (bar_width / 2.0))
    CImGui.AddRectFilled(draw_list,
                        (bar_center_x, cursor_y),
                        (bar_center_x + fill_width, cursor_y + bar_height),
                        CImGui.ColorConvertFloat4ToU32(CImGui.HSV(0.4, 0.6, 0.6)))

    CImGui.Dummy((Float32(0.0), Float32(bar_height + 10)))

    # Draw vertical bar
    CImGui.Text("Y: $(round(y_pos, digits=3))")
    bar_height = Float32(200.0)
    bar_top = Float32(CImGui.GetCursorScreenPos().y)
    bar_width_v = Float32(10.0)  # Width of vertical bar

    # Background for vertical bar
    CImGui.AddRectFilled(draw_list,
                        (center_x - bar_width_v, bar_top),
                        (center_x + bar_width_v, bar_top + bar_height),
                        CImGui.ColorConvertFloat4ToU32(CImGui.HSV(0.0, 0.0, 0.2)))

    # Actual vertical bar
    fill_height = Float32((1.0 - y_pos) * (bar_height / 2.0))  # Y is inverted in screen coordinates
    bar_mid_y = Float32(bar_top + bar_height/2)

    CImGui.AddRectFilled(draw_list,
                        (center_x - bar_width_v, bar_mid_y),
                        (center_x + bar_width_v, bar_mid_y + fill_height),
                        CImGui.ColorConvertFloat4ToU32(CImGui.HSV(0.0, 0.7, 0.7)))

    CImGui.Dummy((Float32(0.0), Float32(bar_height + 20)))

    # Draw crosshair marker at the intersection
    marker_size = Float32(10.0)
    marker_x = Float32(bar_center_x + x_pos * (bar_width/2))
    marker_y = Float32(bar_top + bar_height/2 + y_pos * (bar_height/2))

    CImGui.AddCircleFilled(draw_list,
                         (marker_x, marker_y),
                         marker_size,
                         CImGui.ColorConvertFloat4ToU32(CImGui.HSV(0.13, 0.6, 0.9)))

    CImGui.End()
end


struct TestMapping <: IOMapping end

function Systems.assign_input!(sys::System{TestSystem}, data::XBoxController, ::TestMapping)
    sys.u[1] = get_axis_value(data, :stick_x)
    sys.u[2] = get_axis_value(data, :stick_y)
end


function demo_joysticks()

    sys = TestSystem() |> System
    sim = Simulation(sys; t_end = 30, dt = 0.02)

    for joystick in get_connected_joysticks()
        Sim.attach!(sim, joystick)
    end

    Sim.run_interactive!(sim)

end


end #module