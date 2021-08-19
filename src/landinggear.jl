module LandingGear

using StaticArrays
using ComponentArrays

using Flight.Airdata
using Flight.Dynamics
using Flight.Component
import Flight.Component: get_wr_Ob_b, get_h_Gc_b
import Flight.System: X, Y, U, f_cont!, f_disc!

export LandingGearLeg, LandingGearGroup

abstract type AbstractLandingGearLeg <: AbstractComponent end

Base.@kwdef struct LandingGearLeg <: AbstractLandingGearLeg
end

#################### AbstractSystem interface
X(::LandingGearLeg) = ComponentVector(state = 0.0)
Y(::LandingGearLeg) = ComponentVector(output = 0.0) #both are valid, but not having missing elements in the overall Y vector is slightly faster for some operations
U(::LandingGearLeg) = missing
f_cont!(y, ẋ, x, u, t, ldg::LandingGearLeg, trn = nothing) = (ẋ.state = 0.001x.state)
f_disc!(x, u, t, ldg::LandingGearLeg, trn = nothing) = false

#################### AbstractComponent interface
get_wr_Ob_b(y, comp::LandingGearLeg) = Wrench()
get_h_Gc_b(y, comp::LandingGearLeg) = SVector(0.0, 0, 0)


struct LandingGearGroup{C} <: AbstractComponentGroup{C} end

function LandingGearGroup(nt::NamedTuple{L, T}  where {L, T<:NTuple{N,AbstractLandingGearLeg} where {N}})
    LandingGearGroup{nt}()
end

end


#in System, define and extend f_branch!

# #individual Component
# f_branch!(y, dx, x, u, t, sys, args...) = f_branch!(Val(has_input(sys)), y, dx, x, u, t, args...)
# f_branch!(::Val{true}, y, dx, x, u, t, sys, args...) = f_cont!(y, dx, x, u, t, sys, args...)
# f_cont!(::HasInput, y, dx, x ,u, t, sys, args...) = f_cont!(y, dx, x, u, t, sys, args...)
# f_cont!(::HasNoInput, y, dx, x, u, t, sys, args...) = f_cont!(y, dx, x, t, sys, args...)

# #for a ComponentGroup
# f_cont!(MaybeInput(S), MaybeOutput(S), y, dx, x, u, t, sys, args...)
# f_cont!(::HasInput, ::HasOutput, y, dx, x ,u, t, sys, args...)
# #now, this method needs to consider the possibility for each component that it
# #may have or not Input or Output. so it must do
# for (label, component) in zip(keys(C), values(C))
#     if MaybeInput(typeof(component)) #need tocheck, because if it has no input, u[label] will not exist!
#         f_cont!(y_cmp, dx_cmp, x_cmp, u_cmp, t, cmp, args...)
#     else
#         f_cont!(y_cmp, dx_cmp, x_cmp, t, cmp, args...)
#     end
# end


#this method