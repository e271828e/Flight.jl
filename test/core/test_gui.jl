module TestGUI

using Test

using Flight.FlightCore.Systems
using Flight.FlightCore.GUI

using Flight.FlightAircraft.Control

export test_gui

function test_gui()
    # target = SimpleWorld(Cessna172Rv0(), SimpleEnvironment()) |> System
    target = PIContinuous{2}() |> System
    r = Renderer()
    GUI.init!(r)
    GUI.run(r, GUI.draw!, target)
end

end