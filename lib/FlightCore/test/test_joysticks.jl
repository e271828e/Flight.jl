module TestJoysticks

using Test

using FlightCore

export test_joysticks

function test_joysticks()
    @testset verbose = true "Joysticks" begin
        test_t16000m()
    end
end

################################## Joystick ####################################

@kwdef struct TestSystem <: ModelDefinition end

@kwdef mutable struct TestSystemU
    input::Float64 = 0
end

@kwdef struct TestSystemY
    input::Float64 = 0
end

Modeling.U(::TestSystem) = TestSystemU()
Modeling.Y(::TestSystem) = TestSystemY()

@no_ode TestSystem
@no_step TestSystem

function Modeling.f_periodic!(::Unconditional, mdl::Model{<:TestSystem})
    mdl.y = TestSystemY(; input = mdl.u.input)
end

function GUI.draw(mdl::Model{TestSystem}, label::String = "TestSystem")

    (; input) = mdl.y

    BeginWindow(label)

        TextFormatted("input = $input")

    EndWindow()

end #function

function IODevices.assign_input!(mdl::Model{TestSystem}, ::IOMapping, data::Joysticks.T16000MData)
    mdl.u.input = data.axes.stick_x
end

function test_t16000m()

    @testset verbose = true "Joystick Input" begin

        mdl = TestSystem() |> Model
        sim = Simulation(mdl; t_end = 10.0)
        joystick = update_connected_joysticks()[1]
        Sim.attach!(sim, joystick)

        Sim.run!(sim; gui = true)

    end

end

end # module
