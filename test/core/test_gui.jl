module TestGUI

using Test
using Revise

using Flight.FlightCore

# using Flight.FlightLib
# using Flight.FlightAircraft


export test_gui

function test_gui(target::Model)
    f_draw = let target = target
        () -> GUI.draw!(target)
        # return () -> GUI.draw!(target)
    end
    r = Renderer(; f_draw)
    GUI.render_loop(r)
    # GUI.run(r, GUI.draw!, target)
end

end