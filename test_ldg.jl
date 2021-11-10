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

@btime f_cont!(contact_sys, wow_false)
@btime f_cont!($contact_sys, $wow_true, $wet, 1, @SVector [0, 0, 0])
# @code_warntype f_cont!(contact_sys, wow_true, wet, 1, @SVector [0, 0, 0])



# end