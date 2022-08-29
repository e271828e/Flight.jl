module TestIODevices

using UnPack
using CImGui
using CImGui.CSyntax

using Flight
using Flight.Essentials: PICompensator, PICompensatorU, PICompensatorY

struct TestMapping <: InputMapping end

function IODevices.assign!(u::PICompensatorU{N},
                            joystick::Joystick{XBoxControllerID},
                            ::TestMapping) where {N}

    u.input .= get_axis_value(joystick, :right_analog_y)
    # u.reset .= was_released(joystick, :button_A)
    # u.sat_enable .⊻= was_released(joystick, :button_Y)
end


function test_input_pi()

    sys = PICompensator{2}(k_p = 0, k_i = 0.2) |> System
    sim = Simulation(sys; t_end = 120, dt = 0.05)
    joystick = get_connected_joysticks()[1]
    joystick_interface = Sim.attach_io!(sim, joystick; mapping = TestMapping())

    @sync begin
        Threads.@spawn IODevices.start!(joystick_interface)
        Threads.@spawn Sim.run_paced!(sim; rate = 5, verbose = true)
    end

end

end #module