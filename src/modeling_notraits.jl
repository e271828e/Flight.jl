module ModelingNoTraits

using Dates
using UnPack
using SciMLBase, OrdinaryDiffEq, DiffEqCallbacks
using ComponentArrays, RecursiveArrayTools
using DataStructures: OrderedDict

import Flight.Modeling: f_cont!, f_disc!, init_x, init_y, init_u, init_d
export SysDescNew, SysNew, SysGroupDescNew
export HasContinuousStates, HasNoContinuousStates
export HasOutputs, HasNoOutputs

############################# SysDescNew ###########################

abstract type SysDescNew end #anything from which we can build a System

# init_x(::T) where {T<:SysDescNew} = init_x(T)
init_x(::Type{T} where {T<:SysDescNew}) = nothing

# init_y(::T) where {T<:SysDescNew} = init_y(T)
init_y(::Type{T} where {T<:SysDescNew}) = nothing

############################# SysNew ############################

mutable struct SysNew{D <: SysDescNew, X <: Union{Nothing, AbstractVector{<:Float64}}, Y, P, S}
    ẋ::X #continuous state vector derivative
    x::X #continuous state vector (to be used as a buffer for f_cont! evaluation)
    y::Y #output state
    t::Base.RefValue{Float64} #this allows propagation of t updates down the subsystem hierarchy
    params::P
    subsystems::S
end

function SysNew(c::T, ẋ = init_x(T), x = init_x(T), y = init_y(T), t = Ref(0.0)) where {T<:SysDescNew}

    params = c #by default assign the system descriptor as System parameters
    subsystems = nothing
    SysNew{map(typeof, (c, x, y, params, subsystems))...}(
                                    ẋ, x, y, t, params, subsystems)
end

f_cont!(::SysNew, args...) = MethodError(f_cont!, args) |> throw
(f_disc!(::SysNew, args...)::Bool) = MethodError(f_disc!, args) |> throw


######################### SysGroupDescNew ###############################

abstract type SysGroupDescNew <: SysDescNew end

# Base.keys(g::SysGroupDescNew) = propertynames(g)
# Base.values(g::SysGroupDescNew) = map(λ -> getproperty(g, λ), keys(g))

function init_x(::Type{T}) where {T <: SysGroupDescNew}

    x_dict = OrderedDict{Symbol, AbstractVector{<:Real}}()

    for (label, type) in zip(fieldnames(T), T.types)
        x_ss = init_x(type)
        !isnothing(x_ss) ? x_dict[label] = x_ss : nothing
    end

    return !isempty(x_dict) ? ComponentVector(x_dict) : nothing

end

function init_y(::Type{T}) where {T <: SysGroupDescNew}

    y_dict = OrderedDict{Symbol, Any}()

    for (label, type) in zip(fieldnames(T), T.types)
        y_ss = init_y(type)
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

function SysNew(g::T, ẋ = init_x(T), x = init_x(T), y = init_y(T), t = Ref(0.0)) where {T<:SysGroupDescNew}

    ss_names = fieldnames(T)
    ss_list = Vector{SysNew}()

    for name in ss_names
        ẋ_ss = maybe_getproperty(ẋ, name) #if x is a ComponentVector, this returns a view
        x_ss = maybe_getproperty(x, name) #idem
        y_ss = maybe_getproperty(y, name)
        push!(ss_list, SysNew(getproperty(g, name), ẋ_ss, x_ss, y_ss, t))
        # push!(ss_list, SysNew(map((λ)->getproperty(λ, label), (g, ẋ, x, y, u,
        # d))..., t))
    end

    params = nothing
    subsystems = NamedTuple{ss_names}(ss_list)

    SysNew{map(typeof, (g, x, y, params, subsystems))...}(
                         ẋ, x, y, t, params, subsystems)

end


#NOT WORKING WITH REPL-DEFINED SYSGROUPDESC SUBTYPES
@inline @generated function f_cont!(sys::SysNew{T}, args...) where {T<:SysGroupDescNew}

    ss_labels = fieldnames(T)

    # Core.print("Generated function called")
    ex_main = Expr(:block)

    #call f_cont! on each subsystem
    ex_calls = Expr(:block)
    for label in ss_labels
        push!(ex_calls.args,
            :(f_cont!(sys.subsystems[$(QuoteNode(label))], args...)))
    end

    push!(ex_main.args, ex_calls)

    if !isnothing(init_y(T)) #find out if T has non-trivial output
        #retrieve the y from each subsystem and build a tuple from them
        # Core.println("Non-null output")
        ex_tuple = Expr(:tuple)
        for label in ss_labels
            push!(ex_tuple.args,
                :(sys.subsystems[$(QuoteNode(label))].y))
        end

        #build a NamedTuple from the subsystem's labels and the constructed tuple
        ex_y = Expr(:call, Expr(:curly, NamedTuple, ss_labels), ex_tuple)

        #assign the result to the parent system's y
        ex_assign_y = Expr(:(=), :(sys.y), ex_y)

        #push everything into the main block expression
        push!(ex_main.args, ex_assign_y)
    end

    push!(ex_main.args, :(return nothing))

    return ex_main

end


#WE NEED TO KNOW IN F_CONT! IF WE REALLY NEED TO ASSEMBLE Y. if init_y(T) is not
#nothing, we do

end #module