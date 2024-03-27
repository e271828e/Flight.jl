module DemoUDP

using StructTypes, JSON3

using Flight
using Flight.FlightCore.Sim
using Flight.FlightComponents.Control.Continuous: PIVector, PIVectorInput

export demo_udp

function demo_udp()

    sys = PIVector{1}() |> System
    sim = Simulation(sys; t_end = 1, dt = 0.02)

    Sim.attach!(sim, Networking.DummyUDPSender())
    Sim.attach!(sim, Networking.DummyUDPReceiver())
    Sim.run_paced!(sim)

    #on the REPL, close the GUI prematurely. then run again until done. then run
    #again and get the isdone error. then reinit! and run again.

end

end #module