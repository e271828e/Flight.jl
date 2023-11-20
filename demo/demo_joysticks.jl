module DemoJoysticks

using UnPack

using Flight.FlightCore.GUI
using Flight.FlightComponents.Control.Continuous: PIVector

export demo_joysticks

struct TestMapping <: InputMapping end

function IODevices.assign!(sys::System{<:PIVector{2}},
                            joystick::T16000M,
                            ::TestMapping) where {N}

    sys.u.setpoint[1] .= get_axis_value(joystick, :stick_x)
    sys.u.setpoint[2] .= get_axis_value(joystick, :stick_y)
end

function demo_joysticks()

    sys = PIVector{2}(k_p = 0, k_i = 0.2) |> System
    sim = Simulation(sys; t_end = 30, dt = 0.02)

    joy_interfaces = Vector{IODevices.Interface}()
    for joystick in get_connected_joysticks()
        push!(joy_interfaces, attach_io!(sim, joystick; mapping = TestMapping()))
    end

    @sync begin
        for interface in joy_interfaces
            Threads.@spawn IODevices.start!(interface)
        end

        # disable_gui!(sim)
        Threads.@spawn Sim.run_paced!(sim; pace = 1, verbose = true)
    end

end

end #module