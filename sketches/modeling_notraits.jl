module ModelingNoTraits

using Dates
using UnPack
using SciMLBase, OrdinaryDiffEq, DiffEqCallbacks
using ComponentArrays, RecursiveArrayTools
using DataStructures: OrderedDict

import Flight.Modeling: f_cont!, f_disc!, init_x, init_y, init_u, init_d
export SysDescNew, SysNew, SysGroupDescNew


abstract type SysDescNew end #anything from which we can build a System

############################# SysNew ############################

mutable struct SysNew{  T <: SysDescNew,
                        X <: Union{Nothing, AbstractVector{Float64}},
                        Y, U, D, P, S}
    ẋ::X #continuous state vector derivative
    x::X #continuous state vector
    y::Y #output
    u::U #control input
    d::D #discrete state
    t::Base.RefValue{Float64} #this allows propagation of t updates down the subsystem hierarchy
    params::P
    subsystems::S
end

#this is dangerous, discourage defining init methods directly on instances
# init_x(::T) where {T<:SysDescNew} = init_x(T)
# init_y(::T) where {T<:SysDescNew} = init_y(T)
# init_u(::T) where {T<:SysDescNew} = init_u(T)
# init_d(::T) where {T<:SysDescNew} = init_d(T)

init_x(::Type{T} where {T<:SysDescNew}) = nothing
init_y(::Type{T} where {T<:SysDescNew}) = nothing
init_u(::Type{T} where {T<:SysDescNew}) = nothing
init_d(::Type{T} where {T<:SysDescNew}) = nothing

function SysNew(c::T, ẋ = init_x(T), x = init_x(T), y = init_y(T),
                u = init_u(T), d = init_d(T), t = Ref(0.0)) where {T<:SysDescNew}

    params = c #by default assign the system descriptor as System parameters
    subsystems = nothing
    SysNew{map(typeof, (c, x, y, u, d, params, subsystems))...}(
                                    ẋ, x, y, u, d, t, params, subsystems)
end

f_cont!(::SysNew, args...) = MethodError(f_cont!, args) |> throw
(f_disc!(::SysNew, args...)::Bool) = MethodError(f_disc!, args) |> throw


######################### SysGroupDescNew ###############################

abstract type SysGroupDescNew <: SysDescNew end

init_x(::Type{T}) where {T<:SysGroupDescNew} = init_cv(T, init_x)
init_y(::Type{T}) where {T<:SysGroupDescNew} = init_nt(T, init_y)
init_u(::Type{T}) where {T<:SysGroupDescNew} = init_nt(T, init_u)
init_d(::Type{T}) where {T<:SysGroupDescNew} = init_nt(T, init_d)

function init_cv(::Type{T}, f::Function) where {T<:SysGroupDescNew}

    dict = OrderedDict{Symbol, AbstractVector{Float64}}()

    for (ss_label, ss_type) in zip(fieldnames(T), T.types)
        ss_value = f(ss_type)
        !isnothing(ss_value) ? dict[ss_label] = ss_value : nothing
    end

    return !isempty(dict) ? ComponentVector(dict) : nothing

end

function init_nt(::Type{T}, f::Function) where {T<:SysGroupDescNew}
    dict = OrderedDict{Symbol, Any}()

    for (ss_label, ss_type) in zip(fieldnames(T), T.types)
        ss_value = f(ss_type)
        !isnothing(ss_value) ? dict[ss_label] = ss_value : nothing
    end

    return !isempty(dict) ? NamedTuple{Tuple(keys(dict))}(values(dict)) : nothing

end

#the x of a SysGroupDescNew will be either a ComponentVector or nothing. if it's
#nothing it's because init_x returned nothing, and this is only the case if all
#of its subsystems' init_x in turn returned nothing. in this scenario, we can
#assign nothing to its subsystem's x the same goes for y, but with a NamedTuple
#instead of a ComponentVector
function maybe_getproperty(x, label)
    !isnothing(x) && (label in keys(x)) ? getproperty(x, label) : nothing
end

function SysNew(g::T, ẋ = init_x(T), x = init_x(T), y = init_y(T),
                u = init_u(T), d = init_d(T), t = Ref(0.0)) where {T<:SysGroupDescNew}

    ss_names = fieldnames(T)
    ss_list = Vector{SysNew}()

    for name in ss_names
        ẋ_ss = maybe_getproperty(ẋ, name) #if x is a ComponentVector, this returns a view
        x_ss = maybe_getproperty(x, name) #idem
        y_ss = maybe_getproperty(y, name)
        u_ss = maybe_getproperty(u, name)
        d_ss = maybe_getproperty(d, name)
        push!(ss_list, SysNew(getproperty(g, name), ẋ_ss, x_ss, y_ss, u_ss, d_ss, t))
        # push!(ss_list, SysNew(map((λ)->getproperty(λ, label), (g, ẋ, x, y, ))..., t))
    end

    params = nothing
    subsystems = NamedTuple{ss_names}(ss_list)

    SysNew{map(typeof, (g, x, y, u, d, params, subsystems))...}(
                         ẋ, x, y, u, d, t, params, subsystems)

end

@inline @generated function (f_cont!(sys::SysNew{T, X, Y}, args...)
    where {T<:SysGroupDescNew, X, Y <: Union{Nothing, NamedTuple{L, M}}} where {L, M})

    # Core.println("Generated function called")
    ex_main = Expr(:block)

    #call f_cont! on every subsystem
    ex_calls = Expr(:block)
    for label in fieldnames(T)
        push!(ex_calls.args,
            :(f_cont!(sys.subsystems[$(QuoteNode(label))], args...)))
    end

    push!(ex_main.args, ex_calls)

    #when Y is not Nothing, L will hold the labels of those subsystems that have
    #output. therefore, we know we can retrieve the y fields of those subsystems
    #and assemble them into a NamedTuple, which will have the same type as Y
    if Y !== Nothing

        Core.println("Assembling outputs")
        #tuple expression for the names of children with outputs
        ex_ss_outputs = Expr(:tuple) #tuple expression for children's outputs

        for label in L
            push!(ex_ss_outputs.args,
                :(sys.subsystems[$(QuoteNode(label))].y))
        end

        #build a NamedTuple from the subsystem's labels and the constructed tuple
        ex_y = Expr(:call, Expr(:curly, NamedTuple, L), ex_ss_outputs)

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