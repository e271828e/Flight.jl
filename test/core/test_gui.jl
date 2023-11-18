module TestGUI

using Test

using Flight.FlightCore.Systems
using Flight.FlightCore.GUI

using Flight.FlightComponents.Control

# using Flight.FlightAircraft.C172FBWCAS
using Flight.FlightAircraft.C172RBase
# using Flight.FlightAircraft.C172RDirect

export test_gui

function test_gui()
    # target = PIContinuous{2}() |> System
    # target = C172FBWCAS.PitchControl() |> System
    # target = Cessna172FBWCAS() |> System
    target = Cessna172RBase() |> System
    r = Renderer()
    GUI.init!(r)
    GUI.run(r, GUI.draw!, target)
    # GUI.run(r, GUI.draw, target)
end

end