module ModelingNew

using Dates
using UnPack
using SciMLBase, OrdinaryDiffEq, DiffEqCallbacks
using ComponentArrays, RecursiveArrayTools

import Flight.Modeling: f_cont!, f_disc!, init_x, init_y, init_u, init_d
export SysDescNew, SysNew, SysGroupDescNew

############################# SysDescNew ###########################

abstract type SysDescNew end #anything from which we can build a System

init_x(::SysDescNew) = nothing
init_y(::SysDescNew) = nothing #sytems are not required to have outputs

############################# SysNew ############################

mutable struct SysNew{D <: SysDescNew, X <: Union{Nothing, AbstractVector{<:Float64}}, Y, P, S}
    ẋ::X #continuous state vector derivative
    x::X #continuous state vector (to be used as a buffer for f_cont! evaluation)
    y::Y #output state
    t::Base.RefValue{Float64} #this allows propagation of t updates down the subsystem hierarchy
    params::P
    subsystems::S
end

function SysNew(c::SysDescNew, ẋ = init_x(c), x = init_x(c), y = init_y(c), t = Ref(0.0))

    params = c #by default assign the system descriptor as System parameters
    subsystems = nothing
    SysNew{map(typeof, (c, x, y, params, subsystems))...}(
                                    ẋ, x, y, t, params, subsystems)
end

f_cont!(::SysNew, args...) = MethodError(f_cont!, args) |> throw
(f_disc!(::SysNew, args...)::Bool) = MethodError(f_disc!, args) |> throw


######################### SysGroupDescNew ###############################

abstract type SysGroupDescNew <: SysDescNew end

Base.keys(g::SysGroupDescNew) = propertynames(g)
Base.values(g::SysGroupDescNew) = map(λ -> getproperty(g, λ), keys(g))

function init_x(g::SysGroupDescNew)

    x_dict = Dict{Symbol, AbstractVector{<:Real}}()

    for label in keys(g)
        x_ss = init_x(getproperty(g, label))
        !isnothing(x_ss) ? x_dict[label] = x_ss : nothing
    end

    return !isempty(x_dict) ? ComponentVector(x_dict) : nothing

end
    # = NamedTuple{keys(g)}(init_x.(values(g))) |> ComponentVector
function init_y(g::SysGroupDescNew)

    y_dict = Dict{Symbol, Any}()

    for label in keys(g)
        y_ss = init_y(getproperty(g, label))
        !isnothing(y_ss) ? y_dict[label] = y_ss : nothing
    end

    return !isempty(y_dict) ? NamedTuple{Tuple(keys(y_dict))}(values(y_dict)) : nothing

end

function maybe_getproperty(x, label)
    #the x of a SysGroupDescNew will be either a ComponentVector or nothing.
    #if it's nothing it's because init_x returned nothing, and this is only
    #the case if all of its subsystems' init_x in turn returned nothing.
    #in this scenario, we can assign nothing to its subsystem's x
    #the same goes for y, but with a NamedTuple instead of a ComponentVector
    !isnothing(x) && (label in keys(x)) ? getproperty(x, label) : nothing
end

function SysNew(g::SysGroupDescNew, ẋ = init_x(g), x = init_x(g),
                    y = init_y(g), t = Ref(0.0))

    ss_labels = keys(g)
    ss_list = Vector{SysNew}()

    for label in ss_labels
        ẋ_ss = maybe_getproperty(ẋ, label) #if x is a ComponentVector, this returns a view
        x_ss = maybe_getproperty(x, label) #idem
        y_ss = maybe_getproperty(y, label)
        push!(ss_list, SysNew(getproperty(g, label), ẋ_ss, x_ss, y_ss, t))
        # push!(ss_list, SysNew(map((λ)->getproperty(λ, label), (g, ẋ, x, y, ))..., t))
    end

    desc = nothing
    subsystems = NamedTuple{ss_labels}(ss_list)

    SysNew{map(typeof, (g, x, y, desc, subsystems))...}(
                         ẋ, x, y, t, desc, subsystems)

end

end #module