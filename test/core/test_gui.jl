module TestGUI

using Test

using Flight.FlightCore.Systems
using Flight.FlightCore.GUI

using Flight.FlightLib
# using Flight.FlightLib.Control
# using Flight.FlightLib.LandingGear

# using Flight.FlightAircraft.C172S
# using Flight.FlightAircraft.C172X

export test_gui

function test_gui()
    # target = Control.Discrete.PIDVector{2}() |> System
    # target = LandingGearUnit() |> System
    # target = Control.Continuous.PIVector{3}() |> System
    # target = C172CAS.PitchControl() |> System
    # target = Cessna172CAS() |> System
    # target = Cessna172SBase() |> System
    # target = Cessna172Sv0() |> System
    target = Atmosphere.ISAModel() |> System
    f_draw = let target = target
        () -> GUI.draw!(target)
        # return () -> GUI.draw!(target)
    end
    r = Renderer(; f_draw)
    GUI.render_loop(r)
    # GUI.run(r, GUI.draw!, target)
end

end