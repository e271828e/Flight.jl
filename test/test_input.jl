module TestInput

using Test
using BenchmarkTools
using UnPack
using ComponentArrays

using Flight
using Flight.Input: update!, assign!, AbstractJoystick

export TestInput

struct TestTarget end

function assign!(::TestTarget, xbox::XBoxController, ::DefaultInputMapping)
    return

end

function test_input()
    test_joysticks()
end

function test_joysticks()

    init_joysticks()
    for joy in values(connected_joysticks)
        test(joy)
    end

end

test(joy::AbstractJoystick) = println("No tests defined for $joy")

function test(joy::XBoxController)

    target = TestTarget()
    t0 = time()
    println("Testing $(typeof(joy))")
    while time() - t0 < 10
        update!(joy)
        assign!(target, joy)
        println()
        # show(joy)
        # joy.axes |> showfields
        sleep(1)
    end

end


end #module