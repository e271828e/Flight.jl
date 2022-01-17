
"""
#abstract type MaybeContinuousDynamics end struct
HasContinuousDynamics end struct HasNoContinuousDynamics end
MaybeContinuousDynamics(::Type{<:Any}) = NoContinuousDynamics
MaybeContinuousDynamics(::Type{<:EThruster}) = ContinuousDynamics
# MaybeContinuousDynamics(::Type{<:TunableISA}) = NoContinuousDynamics

#this is the first method that should be overridden by any System that
#declares having continuous dynamics. function f_cont!(s::System{C},
args...) where {C<:SystemDescriptor} f_cont!(s, MaybeContinuousDynamics(C),
args...) end function f_cont!(::S, ::ContinuousDynamics, args...) where
{S<:System} println("Warning: S declared itself as having
ContinuousDynamics but does not implement f_cont!") end f_cont!(::System,
::NoContinuousDynamics, args...) = nothing #its OK #this assumes that when some
type calls f_cont! with variable arguments, but the #second one is NOT of type
MaybeContinuousDynamics, the compiler will dispatch #to the first method first.
this method should have been overridden

#we could do the same with MaybeDiscreteDynamics


#what purpose does it serve? it avoids an unwanted and inadverted dispatch to
#the fallback f_cont! or f_disc due to wrong interface in its own
#implementation. obviously, the component is assumed to have actually
#implemented f_cont! and f_disc! and declared itself as having Continous or
#DiscreteDynamics. a component that neither declares ContinuousDynamics nor
#defines f_cont! will receive no warning

#we could go another step and do the same with init_x, which will also dispatch
#on MaybeContinuousDynamics, and init_d, on MaybeDiscreteDynamics. then
#MaybeInput, MaybeOutput for init_u and init_d

#however, all this does is force systems with continuous dyanmics to define the
#a MaybeContinuousDynamics method in addition to f_cont!. on exchange, systems
without continuous dynamics are not forced to implement f_cont!. not much to be
gained

"""