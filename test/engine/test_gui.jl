module TestGUI

using Test

using Flight

export test_gui

function test_gui()
    ctl = C172R.ReversibleControls() |> System
    r = Renderer()
    GUI.init!(r)
    GUI.run(r, GUI.draw!, ctl)
end

end