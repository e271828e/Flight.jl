module TestInput

using Test
using BenchmarkTools
using UnPack
using ComponentArrays

using Flight
using Flight.Input: update!, assign!

export test_input

struct TestTarget end

function Input.assign!(::TestTarget, joy::Joystick{XBoxController}, ::DefaultInputMapping)
    println(joy)
end

function test_input()
    test_joysticks()
end

function test_joysticks()

    for joy in get_connected_joysticks()
        test(joy)
    end

end

test(joy::Joystick) = println("No tests defined for $(typeof(joy.id))")

function test(joy::Joystick{XBoxController})

    target = TestTarget()
    t0 = time()
    println("Testing $(typeof(joy.id))")
    while time() - t0 < 10
        update!(joy)
        assign!(target, joy)
        sleep(1)
    end

end


end #module