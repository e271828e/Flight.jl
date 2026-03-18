module TestGUI

using Test
using Revise

using FlightCore

export test_render_loop, test_gui

function test_gui()
    @testset verbose = true "GUI" begin
        @test_nowarn test_render_loop()
    end
end

struct TestSystem <: ModelDefinition end

Modeling.X(::TestSystem) = [0.0]

function GUI.draw!(mdl::Model{TestSystem})
    x = mdl.x
    CImGui.PushItemWidth(-50)
        x[1] = GUI.safe_slider("State", x[1], -10, 10, "%.3f")
    CImGui.PopItemWidth()
end

function test_render_loop(target = Model(TestSystem()), timeout::Real = 1.0, args...)
    f_draw = let target = target
        () -> GUI.draw!(target, args...)
    end
    r = Renderer(; f_draw)
    GUI.render_loop(r, timeout)
end

end