module TestGUI

using Test

using Flight.FlightCore.Systems
using Flight.FlightCore.GUI

using Flight.FlightLib.Control
using Flight.FlightLib.LandingGear

using Flight.FlightAircraft.C172R
using Flight.FlightAircraft.C172FBW
using Flight.FlightAircraft.C172RPA

export test_gui

function test_gui()
    # target = Control.Discrete.PIDVector{2}() |> System
    # target = LandingGearUnit() |> System
    # target = Control.Continuous.PIVector{3}() |> System
    # target = C172CAS.PitchControl() |> System
    # target = Cessna172CAS() |> System
    # target = Cessna172RBase() |> System
    # target = C172FBW.Actuation() |> System;
    target = Cessna172FBWv1() |> System
    # target = Cessna172R() |> System
    f_draw = let target = target
        () -> GUI.draw!(target)
        # return () -> GUI.draw!(target)
    end
    # target = C172FBWv1.LonControl() |> System
    # target = Cessna172FBW() |> System
    r = Renderer(; f_draw)
    GUI.render_loop(r)
    # GUI.run(r, GUI.draw!, target)
end

end