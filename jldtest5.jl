using Flight
using StaticArrays
using Interpolations
using JLD

#for some reason, we need to load StaticArrays and Interpolations; other fields
#in Propeller and Piston.Engine, which are not visible either, don't seem to
#require their types to be explicitly loaded

prop = Propeller()
eng = Piston.Engine()

file = jldopen("jldtest.jld", "w")
write(file, "prop", prop)
write(file, "eng", eng)
close(file)

file = jldopen("jldtest.jld", "r")
prop_loaded = read(file, "prop")
eng_loaded = read(file, "eng")
close(file)
