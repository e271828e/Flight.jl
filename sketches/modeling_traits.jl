module ModelingTraits

using Dates
using UnPack
using SciMLBase, OrdinaryDiffEq, DiffEqCallbacks
using ComponentArrays, RecursiveArrayTools
using DataStructures: OrderedDict

import Flight.Modeling: f_cont!, f_disc!, init_x, init_y, init_u, init_d
export SysDescNew, SysNew, SysGroupDescNew
export HasContinuousState, HasNoContinuousState
export HasOutput, HasNoOutput


abstract type SysDescNew end #anything from which we can build a System

############################# System Traits ###########################

abstract type ContinuousStateTrait end
struct HasNoContinuousState <: ContinuousStateTrait end
struct HasContinuousState <: ContinuousStateTrait end

ContinuousStateTrait(::T) where {T<:SysDescNew} = ContinuousStateTrait(T)
ContinuousStateTrait(::Type{T}) where {T<:SysDescNew} = HasNoContinuousState()

init_x(::T) where {T<:SysDescNew} = init_x(T)
init_x(::Type{T}) where {T<:SysDescNew} = init_x(T, ContinuousStateTrait(T))
init_x(::Type{T} where {T<:SysDescNew}, ::HasNoContinuousState) = nothing
function init_x(::Type{T}, ::HasContinuousState) where {T<:SysDescNew}
    error("$T was declared to have continuous state, so a method with signature
        init_x(::Type{$T})::AbstractVector{Float64} is required")
end

abstract type OutputTrait end
struct HasNoOutput <: OutputTrait end
struct HasOutput <: OutputTrait end

OutputTrait(::T) where {T<:SysDescNew} = OutputTrait(T)
OutputTrait(::Type{T}) where {T<:SysDescNew} = HasNoOutput()

init_y(::T) where {T<:SysDescNew} = init_y(T)
init_y(::Type{T}) where {T<:SysDescNew} = init_y(T, OutputTrait(T))
init_y(::Type{T} where {T<:SysDescNew}, ::HasNoOutput) = nothing
function init_y(::Type{T}, ::HasOutput) where {T<:SysDescNew}
    error("$T was declared to have output, so a method with signature
        init_y(::Type{$T}) is required")
end

abstract type InputTrait end
struct HasNoInput <: InputTrait end
struct HasInput <: InputTrait end

InputTrait(::T) where {T<:SysDescNew} = InputTrait(T)
InputTrait(::Type{T}) where {T<:SysDescNew} = HasNoInput()

init_u(::T) where {T<:SysDescNew} = init_u(T)
init_u(::Type{T}) where {T<:SysDescNew} = init_u(T, InputTrait(T))
init_u(::Type{T} where {T<:SysDescNew}, ::HasNoInput) = nothing
function init_u(::Type{T}, ::HasInput) where {T<:SysDescNew}
    error("$T was declared to have input, so a method with signature
        init_u(::Type{$T}) is required")
end

abstract type DiscreteStateTrait end
struct HasNoDiscreteState <: DiscreteStateTrait end
struct HasDiscreteState <: DiscreteStateTrait end

DiscreteStateTrait(::T) where {T<:SysDescNew} = DiscreteStateTrait(T)
DiscreteStateTrait(::Type{T}) where {T<:SysDescNew} = HasNoDiscreteState()

init_d(::T) where {T<:SysDescNew} = init_d(T)
init_d(::Type{T}) where {T<:SysDescNew} = init_d(T, DiscreteStateTrait(T))
init_d(::Type{T} where {T<:SysDescNew}, ::HasNoDiscreteState) = nothing
function init_d(::Type{T}, ::HasDiscreteState) where {T<:SysDescNew}
    error("$T was declared to have input, so a method with signature
        init_d(::Type{$T}) is required")
end

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


function ContinuousStateTrait(::Type{T}) where {T<:SysGroupDescNew}
    for field_type in T.types
        if ContinuousStateTrait(field_type) === HasContinuousState()
            return HasContinuousState()
        end
    end
    return HasNoContinuousState()
end

function OutputTrait(::Type{T}) where {T<:SysGroupDescNew}
    for field_type in T.types
        if OutputTrait(field_type) === HasOutput()
            return HasOutput()
        end
    end
    return HasNoOutput()
end

function InputTrait(::Type{T}) where {T<:SysGroupDescNew}
    for field_type in T.types
        if InputTrait(field_type) === HasInput()
            return HasInput()
        end
    end
    return HasNoInput()
end

function DiscreteStateTrait(::Type{T}) where {T<:SysGroupDescNew}
    for field_type in T.types
        if DiscreteStateTrait(field_type) === HasDiscreteState()
            return HasDiscreteState()
        end
    end
    return HasNoDiscreteState()
end

init_x(::Type{T}, ::HasContinuousState) where {T<:SysGroupDescNew} = init_cv(T, init_x)
init_y(::Type{T}, ::HasOutput) where {T<:SysGroupDescNew} = init_nt(T, init_y)
init_u(::Type{T}, ::HasInput) where {T<:SysGroupDescNew} = init_nt(T, init_u)
init_d(::Type{T}, ::HasDiscreteState) where {T<:SysGroupDescNew} = init_nt(T, init_d)

function init_nt(::Type{T}, f::Function) where {T<:SysGroupDescNew}
    dict = OrderedDict{Symbol, Any}()

    for (label, type) in zip(fieldnames(T), T.types)
        ss_value = f(type)
        !isnothing(ss_value) ? dict[label] = ss_value : nothing
    end

    return NamedTuple{Tuple(keys(dict))}(values(dict))
end

function init_cv(::Type{T}, f::Function) where {T<:SysGroupDescNew}

    dict = OrderedDict{Symbol, AbstractVector{Float64}}()

    for (label, type) in zip(fieldnames(T), T.types)
        ss_value = f(type)
        !isnothing(ss_value) ? dict[label] = ss_value : nothing
    end

    return ComponentVector(dict)

end

# function init_x(::Type{T}, ::HasContinuousState) where {T<:SysGroupDescNew}

#     dict = OrderedDict{Symbol, AbstractVector{Float64}}()

#     for (label, type) in zip(fieldnames(T), T.types)
#         x_ss = init_x(type)
#         !isnothing(x_ss) ? dict[label] = x_ss : nothing
#     end

#     return ComponentVector(dict)

# end

# function init_y(::Type{T}, ::HasOutput) where {T<:SysGroupDescNew}

#     dict = OrderedDict{Symbol, Any}()

#     for (label, type) in zip(fieldnames(T), T.types)
#         y_ss = init_y(type)
#         !isnothing(y_ss) ? dict[label] = y_ss : nothing
#     end

#     return NamedTuple{Tuple(keys(dict))}(values(dict))

# end

# function init_u(::Type{T}, ::HasInput) where {T<:SysGroupDescNew}

#     dict = OrderedDict{Symbol, Any}()

#     for (label, type) in zip(fieldnames(T), T.types)
#         u_ss = init_u(type)
#         !isnothing(u_ss) ? dict[label] = u_ss : nothing
#     end

#     return NamedTuple{Tuple(keys(dict))}(values(dict))

# end

# function init_d(::Type{T}, ::HasDiscreteState) where {T<:SysGroupDescNew}

#     dict = OrderedDict{Symbol, Any}()

#     for (label, type) in zip(fieldnames(T), T.types)
#         u_ss = init_d(type)
#         !isnothing(u_ss) ? dict[label] = u_ss : nothing
#     end

#     @assert !isempty(dict) #T has outputs states, shouldn't end up empty

#     return NamedTuple{Tuple(keys(dict))}(values(dict))

# end
# @generated function ContinuousStateTrait(::T) where {T<:SysGroupDescNew}
#     ##this doesn't work, but it does when a Core.println statement is added!!
#     # for field_type in T.types
#     #     # if ContinuousStateTrait(field_type) === HasContinuousState()
#     #     #     Core.print("Got states")
#     #     # end
#     # end
#     # return :(HasNoContinuousState())
#     if any(isa.(ContinuousStateTrait.(T.types), HasContinuousState))
#         return :(HasContinuousState())
#     else
#         return :(HasNoContinuousState())
#     end
# end

# @generated function OutputTrait(::T) where {T<:SysGroupDescNew}
#     # Core.println(T.types)
#     # has_outputs = false

#     # for field_type in T.types
#     #     if OutputTrait(field_type) === HasOutput()
#     #         # Core.println(field_type)
#     #         has_outputs = true
#     #     end
#     # end
#     if any(isa.(OutputTrait.(T.types), HasOutput))
#         return :(HasOutput())
#     else
#         return :(HasNoOutput())
#     end
# end


    #the x of a SysGroupDescNew will be either a ComponentVector or nothing.
    #if it's nothing it's because init_x returned nothing, and this is only
    #the case if all of its subsystems' init_x in turn returned nothing.
    #in this scenario, we can assign nothing to its subsystem's x
    #the same goes for y, but with a NamedTuple instead of a ComponentVector
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


@inline @generated function f_cont!(sys::SysNew{T, X, Y}, args...) where {T<:SysGroupDescNew, X, Y <: Union{Nothing, NamedTuple{L, M}}} where {L, M}
# @inline @generated function f_cont!(sys::SysNew{T, X, Y}, args...) where {T<:SysGroupDescNew, X <: Union{Nothing, AbstractVector{Float64}}, Y <: Union{Nothing, NamedTuple{L, M}}} where {L, M}

    #when Y is not Nothing, L will hold the labels of those subsystems that have
    #output. therefore, we know we can retrieve the y fields of those subsystems
    #and assemble them into a NamedTuple, which will have the same type as Y

    Core.println("Generated function called")
    ex_main = Expr(:block)

    #call f_cont! on every subsystem
    ex_calls = Expr(:block)
    for label in fieldnames(T)
        push!(ex_calls.args,
            :(f_cont!(sys.subsystems[$(QuoteNode(label))], args...)))
    end

    push!(ex_main.args, ex_calls)

    if Y !== Nothing #then it is a NamedTuple type with keys L

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