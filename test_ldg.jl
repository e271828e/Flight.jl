using StaticArrays
using BenchmarkTools
using Flight

# function test_basics()
wow_true = LandingGear.WoW(true)
wow_false = LandingGear.WoW(false)

roll = LandingGear.Rolling()
skid = LandingGear.Skidding(wet = LandingGear.StaticDynamic(0.2, 0.1))

dry = Terrain.Dry; wet = Terrain.Wet; icy = Terrain.Icy
LandingGear.get_μ(roll, dry)
LandingGear.get_μ(roll, wet)

LandingGear.get_μ(skid, dry)
LandingGear.get_μ(skid, wet)
LandingGear.get_μ(skid, icy)

LandingGear.get_μ(roll, dry, 0)
LandingGear.get_μ(roll, dry, 0.5)
LandingGear.get_μ(roll, dry, 1)

LandingGear.get_μ(skid, icy, 0.1)

contact = Contact()
LandingGear.ContactY()

get_x0(contact)
get_y0(contact)

contact_sys = System(contact)

@btime f_cont!($contact_sys, $wow_false)
@btime f_cont!($contact_sys, $wow_true, $wet, 1, @SVector [0, 0, 0])

damper = LandingGear.SimpleDamper()
LandingGear.get_damper_force(damper, -0.1, 0)

no_steering_sys = System(NoSteering())
set_steering_input(no_steering_sys, -1)
f_cont!(no_steering_sys)
get_steering_angle(no_steering_sys)

direct_steering_sys = System(DirectSteering())
set_steering_input(direct_steering_sys, -1)
f_cont!(direct_steering_sys)
try
    set_steering_input(direct_steering_sys, 2)
catch AssertionError
    print("OK")
end
get_steering_angle(direct_steering_sys)

no_braking_sys = System(NoBraking())
set_braking_input(no_braking_sys, 0)
f_cont!(no_braking_sys)
get_braking_coefficient(no_braking_sys)

direct_braking_sys = System(DirectBraking())
set_braking_input(direct_braking_sys, 0.1)
f_cont!(direct_braking_sys)
try
    set_braking_input(direct_braking_sys, -.5)
catch AssertionError
    print("OK")
end
get_braking_coefficient(direct_braking_sys)

# end