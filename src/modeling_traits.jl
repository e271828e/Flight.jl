module ModelingTraits

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

abstract type ContinuousStateTrait end
struct HasNoContinuousStates <: ContinuousStateTrait end
struct HasContinuousStates <: ContinuousStateTrait end

ContinuousStateTrait(::T) where {T<:SysDescNew} = ContinuousStateTrait(T)
ContinuousStateTrait(::Type{T}) where {T<:SysDescNew} = HasNoContinuousStates()

init_x(::T) where {T<:SysDescNew} = init_x(T)
init_x(::Type{T}) where {T<:SysDescNew} = init_x(T, ContinuousStateTrait(T))
init_x(::Type{T} where {T<:SysDescNew}, ::HasNoContinuousStates) = nothing
function init_x(::Type{T}, ::HasContinuousStates) where {T<:SysDescNew}
    error("$T was declared to have continuous states, so a method with signature
        init_x(::Type{$T})::AbstractVector{Float64} is required")
end

abstract type OutputTrait end
struct HasNoOutputs <: OutputTrait end
struct HasOutputs <: OutputTrait end

OutputTrait(::T) where {T<:SysDescNew} = OutputTrait(T)
OutputTrait(::Type{T}) where {T<:SysDescNew} = HasNoOutputs()

init_y(::T) where {T<:SysDescNew} = init_y(T)
init_y(::Type{T}) where {T<:SysDescNew} = init_y(T, OutputTrait(T))
init_y(::Type{T} where {T<:SysDescNew}, ::HasNoOutputs) = nothing
function init_x(::Type{T}, ::HasOutputs) where {T<:SysDescNew}
    error("$T was declared to have outputs, so a method with signature
        init_y(::Type{$T}) is required")
end

############################# SysNew ############################

mutable struct SysNew{D <: SysDescNew,
                    X <: Union{Nothing, AbstractFloat, AbstractVector{<:AbstractFloat}},
                    Y, P, S}
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

function ContinuousStateTrait(::Type{T}) where {T<:SysGroupDescNew}
    for field_type in T.types
        if ContinuousStateTrait(field_type) === HasContinuousStates()
            return HasContinuousStates()
        end
    end
    return HasNoContinuousStates()
end

function OutputTrait(::Type{T}) where {T<:SysGroupDescNew}
    for field_type in T.types
        if OutputTrait(field_type) === HasOutputs()
            return HasOutputs()
        end
    end
    return HasNoOutputs()
end

# @generated function ContinuousStateTrait(::T) where {T<:SysGroupDescNew}
#     ##this doesn't work, but it does when a Core.println statement is added!!
#     # for field_type in T.types
#     #     # if ContinuousStateTrait(field_type) === HasContinuousStates()
#     #     #     Core.print("Got states")
#     #     # end
#     # end
#     # return :(HasNoContinuousStates())
#     if any(isa.(ContinuousStateTrait.(T.types), HasContinuousStates))
#         return :(HasContinuousStates())
#     else
#         return :(HasNoContinuousStates())
#     end
# end

# @generated function OutputTrait(::T) where {T<:SysGroupDescNew}
#     # Core.println(T.types)
#     # has_outputs = false

#     # for field_type in T.types
#     #     if OutputTrait(field_type) === HasOutputs()
#     #         # Core.println(field_type)
#     #         has_outputs = true
#     #     end
#     # end
#     if any(isa.(OutputTrait.(T.types), HasOutputs))
#         return :(HasOutputs())
#     else
#         return :(HasNoOutputs())
#     end
# end

function init_x(::Type{T}, ::HasContinuousStates) where {T<:SysGroupDescNew}

    x_dict = OrderedDict{Symbol, Union{Real, AbstractVector{<:Real}}}()

    for (label, type) in zip(fieldnames(T), T.types)
        x_ss = init_x(type)
        !isnothing(x_ss) ? x_dict[label] = x_ss : nothing
    end

    @assert !isempty(x_dict) #T has continuous states, shouldn't end up empty

    return ComponentVector(x_dict)

end

function init_y(::Type{T}, ::HasOutputs) where {T<:SysGroupDescNew}

    y_dict = OrderedDict{Symbol, Any}()

    for (label, type) in zip(fieldnames(T), T.types)
        y_ss = init_y(type)
        !isnothing(y_ss) ? y_dict[label] = y_ss : nothing
    end

    @assert !isempty(y_dict) #T has outputs states, shouldn't end up empty

    return NamedTuple{Tuple(keys(y_dict))}(values(y_dict))

end

    #the x of a SysGroupDescNew will be either a ComponentVector or nothing.
    #if it's nothing it's because init_x returned nothing, and this is only
    #the case if all of its subsystems' init_x in turn returned nothing.
    #in this scenario, we can assign nothing to its subsystem's x
    #the same goes for y, but with a NamedTuple instead of a ComponentVector
function maybe_getproperty(x, label)
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
        # push!(ss_list, SysNew(map((λ)->getproperty(λ, label), (g, ẋ, x, y, ))..., t))
    end

    desc = nothing
    subsystems = NamedTuple{ss_names}(ss_list)

    SysNew{map(typeof, (g, x, y, desc, subsystems))...}(
                         ẋ, x, y, t, desc, subsystems)

end


@inline @generated function f_cont!(sys::SysNew{T}, args...) where {T<:SysGroupDescNew}

    ss_labels = fieldnames(T)
    ss_types = T.types

    # Core.print("Generated function called")
    ex_main = Expr(:block)

    #call f_cont! on each subsystem
    ex_calls = Expr(:block)
    for label in ss_labels
        push!(ex_calls.args,
            :(f_cont!(sys.subsystems[$(QuoteNode(label))], args...)))
    end

    push!(ex_main.args, ex_calls)

    if OutputTrait(T) === HasOutputs() #at least one child has output

        #tuple expression for the names of children with outputs
        ex_ss_labels = Expr(:tuple)
        #tuple expression for children's outputs
        ex_ss_outputs = Expr(:tuple)

        for (ss_label, ss_type) in zip(ss_labels, ss_types)
            if OutputTrait(ss_type) === HasOutputs() #does this child have output?
                #put its name and its output in the corresponding expressions
                push!(ex_ss_labels.args, :($(QuoteNode(ss_label))))
                push!(ex_ss_outputs.args, :(sys.subsystems[$(QuoteNode(ss_label))].y))
            end
        end

        #build a NamedTuple from the subsystem's labels and the constructed tuple
        ex_y = Expr(:call, Expr(:curly, NamedTuple, ex_ss_labels), ex_ss_outputs)

        #assign the result to the parent system's y
        ex_assign_y = Expr(:(=), :(sys.y), ex_y)

        #push everything into the main block expression
        push!(ex_main.args, ex_assign_y)
    end

    push!(ex_main.args, :(return nothing))

    return ex_main

end

#default implementation calls f_disc! on all Node subsystems with the same
#arguments provided to the parent Node's System, then ORs their outputs.
#can be overridden for specific SystemGroupDescriptor subtypes if needed
@inline @generated function (f_disc!(sys::SysNew{T}, args...)::Bool
    ) where {T<:SysGroupDescNew}

    # Core.print("Generated function called")
    ex = Expr(:block)
    push!(ex.args, :(x_mod = false))
    for label in fieldnames(T)
        #we need all f_disc! calls executed, so | must be used instead of ||
        push!(ex.args,
            :(x_mod = x_mod | f_disc!(sys.subsystems[$(QuoteNode(label))], args...)))
    end
    return ex

end

end #module