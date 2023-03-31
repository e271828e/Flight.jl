module TestGUI

using Test

using Flight

export test_gui

# function test_gui()
#     ctl = C172RDirect.FeedthroughActuation() |> System
#     r = Renderer()
#     GUI.init!(r)
#     GUI.run(r, GUI.draw!, ctl)
# end

function test_gui()
    world = SimpleWorld() |> System
    r = Renderer()
    GUI.init!(r)
    GUI.run(r, GUI.draw!, world)
end

end