module Airframe #this module is really needed, do NOT merge it into aircraft

using LinearAlgebra
using StaticArrays, ComponentArrays
using UnPack

using Flight.Attitude
using Flight.Dynamics
using Flight.ModelingTools
import Flight.ModelingTools: System, init_x0, init_y0, init_u0, init_d0, f_cont!, f_disc!

using Flight.Plotting
import Flight.Plotting: plots

export AbstractAirframeComponent, AbstractAirframeNode, NullAirframeComponent
export get_wr_b, get_hr_b


abstract type AbstractAirframeComponent <: AbstractComponent end

function get_wr_b(::T) where {T<:System{<:AbstractAirframeComponent}}
    error("Method get_wr_b not implemented for type $T or incorrect call signature")
end
function get_hr_b(::T) where {T<:System{<:AbstractAirframeComponent}}
    error("Method hr_b not implemented for type $T or incorrect call signature")
end

######################### NullAirframeComponent ############################

struct NullAirframeComponent <: AbstractAirframeComponent end

@inline get_wr_b(::System{NullAirframeComponent}) = Wrench()
@inline get_hr_b(::System{NullAirframeComponent}) = zeros(SVector{3})

@inline f_cont!(::System{NullAirframeComponent}, args...) = nothing
@inline (f_disc!(::System{NullAirframeComponent}, args...)::Bool) = false

######################## AbstractAirframeNode #############################

#abstract type providing convenience methods for those
#AbstractAirframeComponents consisting only of AbstractAirframeComponent fields

abstract type AbstractAirframeNode <: AbstractAirframeComponent end

Base.keys(node::AbstractAirframeNode) = propertynames(node)
Base.values(node::AbstractAirframeNode) = map(λ -> getproperty(node, λ), keys(node))

init_x0(node::AbstractAirframeNode) = NamedTuple{keys(node)}(init_x0.(values(node))) |> ComponentVector
init_y0(node::AbstractAirframeNode) = NamedTuple{keys(node)}(init_y0.(values(node)))
init_u0(node::AbstractAirframeNode) = NamedTuple{keys(node)}(init_u0.(values(node)))
init_d0(node::AbstractAirframeNode) = NamedTuple{keys(node)}(init_d0.(values(node)))

function System(node::AbstractAirframeNode, ẋ = init_x0(node), x = init_x0(node),
                    y = init_y0(node), u = init_u0(node), d = init_d0(node), t = Ref(0.0))

    ss_list = Vector{System}()
    ss_labels = keys(node)
    # @show ss_labels
    for label in ss_labels
        push!(ss_list, System(map((λ)->getproperty(λ, label), (node, ẋ, x, y, u, d))..., t))
    end

    params = nothing
    subsystems = NamedTuple{ss_labels}(ss_list)

    System{map(typeof, (node, x, y, u, d, params, subsystems))...}(
                         ẋ, x, y, u, d, t, params, subsystems)

end

#default implementation calls f_cont! on all Node subsystems with the same
#arguments provided to the parent Node's System, then builds the NamedTuple
#can be overridden for specific AbstractAirframeNode subtypes if needed
@inline @generated function f_cont!(sys::System{C}, args...) where {C<:AbstractAirframeNode}

    # Core.print("Generated function called")
    ex_main = Expr(:block)

    #call f_cont! on each subsystem
    ex_calls = Expr(:block)
    ss_labels = fieldnames(C)
    for label in ss_labels
        push!(ex_calls.args,
            :(f_cont!(sys.subsystems[$(QuoteNode(label))], args...)))
    end

    #retrieve the y from each subsystem and build a tuple from them
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
    push!(ex_main.args, ex_calls)
    push!(ex_main.args, ex_assign_y)
    push!(ex_main.args, :(return nothing))

    return ex_main

end

#default implementation calls f_disc! on all Node subsystems with the same
#arguments provided to the parent Node's System, then ORs their outputs.
#can be overridden for specific AbstractAirframeNode subtypes if needed
@inline @generated function (f_disc!(sys::System{C}, args...)::Bool
    ) where {C<:AbstractAirframeNode}

    # Core.print("Generated function called")
    ex = Expr(:block)
    push!(ex.args, :(x_mod = false))
    for label in fieldnames(C)
        #we need all f_disc! calls executed, so | must be used instead of ||
        push!(ex.args,
            :(x_mod = x_mod | f_disc!(sys.subsystems[$(QuoteNode(label))], args...)))
    end
    return ex

end

#default implementation retrieves the wrench from each subsystem and sums them
@inline @generated function get_wr_b(sys::System{C}) where {C<:AbstractAirframeNode}

    # Core.print("Generated function called")
    ex = Expr(:block)
    push!(ex.args, :(wr = Wrench())) #allocate a zero wrench
    for label in fieldnames(C)
        push!(ex.args,
            :(wr += get_wr_b(sys.subsystems[$(QuoteNode(label))])))
    end
    return ex

end

#default implementation retrieves angular momentum from each subsystem and sums
#them
@inline @generated function get_hr_b(sys::System{C}) where {C<:AbstractAirframeNode}

    # Core.print("Generated function called")
    ex = Expr(:block)
    push!(ex.args, :(h = SVector(0., 0., 0.))) #allocate
    for label in fieldnames(C)
        push!(ex.args,
            :(h += get_hr_b(sys.subsystems[$(QuoteNode(label))])))
    end
    return ex

end

end #module