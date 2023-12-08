module TestGUI

using Test

using Flight.FlightCore.Systems
using Flight.FlightCore.GUI

using Flight.FlightComponents.Control

# using Flight.FlightAircraft.C172CAS
using Flight.FlightAircraft.C172RBase
# using Flight.FlightAircraft.C172RDirect

export test_gui

function test_gui()
    # target = PIVector{2}() |> System
    # target = C172CAS.PitchControl() |> System
    # target = Cessna172CAS() |> System
    # target = Cessna172RBase() |> System
    # target = Control.Discrete.PIDVector{3}() |> System
    # target = Control.Continuous.PIVector{3}() |> System
    target = Control.Discrete.LQRTracker{10, 2, 2}() |> System
    r = Renderer()
    GUI.init!(r)
    # GUI.run(r, GUI.draw!, target)
    GUI.run(r, GUI.draw, target)
end

end