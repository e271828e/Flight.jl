module System

using ComponentArrays

export SystemDescriptor, SystemDescriptorGroup, HybridSystem, X, U, f_cont!, f_disc!
# export plotlog

abstract type SystemDescriptor end #anything that can go in a HybridSystem

no_extend_error(f::Function, ::Type{S}) where {S} = error("Function $f not implemented for subtype $S or incorrect call signature")

X(::C) where {C<:SystemDescriptor} = no_extend_error(X, C)
D(::C) where {C<:SystemDescriptor} = nothing #systems are not required to have discrete states
U(::C) where {C<:SystemDescriptor} = nothing #sytems are not required to have control inputs

#need the C type parameter for dispatch, the rest for type stability
struct HybridSystem{C<:SystemDescriptor, X, D, U, P, S}
    ẋ::X #continuous state vector derivative
    x::X #continuous state vector (to be used as a buffer for f_cont! evaluation)
    d::D #discrete state vector
    u::U #control inputs
    t::Base.RefValue{Float64} #this allows propagation of t updates down the subsystem hierarchy
    params::P
    subsystems::S
end

f_cont!(::S, args...) where {S<:HybridSystem} = no_extend_error(f_cont!, S)

#this method is free to modify the system's discrete state, control inputs and
#continuous state. if it modifies the latter, it must return true, false
#otherwise. it is dangerous to provide a default fallback for f_disc!, because
#if the intended f_disc! implementation for the System has the wrong interface,
#the dispatch will revert to the fallback, which may not be obvious at all. it
#is safer to force each concrete System that does not require an actual f_disc!
#to implement a trivial f_disc! that returns false
(f_disc!(::S, args...)::Bool) where {S<:HybridSystem} = no_extend_error(f_cont!, S)


######################### SystemDescriptorGroup ##############################

struct SystemDescriptorGroup{T<:SystemDescriptor,N,L} <: SystemDescriptor
    components::NamedTuple{L, M} where {L, M <: NTuple{N, T} where {N, T<:SystemDescriptor}}
    function SystemDescriptorGroup(nt::NamedTuple{L, M}) where {L, M<:NTuple{N, T}} where {N, T<:SystemDescriptor}
        new{T,N,L}(nt)
    end
end
# (G::Type{<:SystemDescriptorGroup})(;kwargs...) = G((; kwargs...))
SystemDescriptorGroup(;kwargs...) = SystemDescriptorGroup((; kwargs...))
Base.length(::SystemDescriptorGroup{T,N,L}) where {T,N,L} = N
Base.getindex(g::SystemDescriptorGroup, i) = getindex(getfield(g,:components), i)
Base.getproperty(g::SystemDescriptorGroup, i::Symbol) = getproperty(getfield(g,:components), i)
Base.keys(::SystemDescriptorGroup{T,N,L}) where {T,N,L} = L
Base.values(g::SystemDescriptorGroup) = values(getfield(g,:components))

X(g::SystemDescriptorGroup{T,N,L}) where {T,N,L} = ComponentVector(NamedTuple{L}(X.(values(g))))
D(g::SystemDescriptorGroup{T,N,L}) where {T,N,L} = NamedTuple{L}(D.(values(g)))
U(g::SystemDescriptorGroup{T,N,L}) where {T,N,L} = NamedTuple{L}(U.(values(g)))

function HybridSystem(g::SystemDescriptorGroup{T,N,L},
    ẋ = X(g), x = X(g), d = D(g), u = U(g), t = Ref(0.0)) where {T,N,L}
    #having L allows us to know the length of g and therefore the number of
    #expressions we need to generate
    s_list = Vector{HybridSystem}()
    for label in L
        s_cmp = HybridSystem(map((λ)->getproperty(λ, label), (g, ẋ, x, d, u))..., t)
        push!(s_list, s_cmp)
    end
    params = nothing #everything is already stored in the subsystem's parameters
    subsystems = NamedTuple{L}(s_list)
    HybridSystem{map(typeof, (g, x, d, u, params, subsystems))...}(ẋ, x, d, u, t, params, subsystems)
end

@inline @generated function f_cont!(sys::HybridSystem{C}, args...) where {C<:SystemDescriptorGroup{T,N,L}} where {T,N,L}
    ex_tuple = Expr(:tuple) #construct a tuple around all the individual function calls
    for label in L
        label = QuoteNode(label)
        ex_ss = quote
            f_cont!(sys.subsystems[$label], args...)
        end
        push!(ex_tuple.args, ex_ss)
    end
    # Core.println(ex_tuple)
    #construct a NamedTuple from the labels L and the constructed tuple
    ex_all = Expr(:call, Expr(:curly, NamedTuple, L), ex_tuple)
    return ex_all
end

@inline @generated function (f_disc!(sys::HybridSystem{C}, args...)::Bool) where {C<:SystemDescriptorGroup{T,N,L}} where {T,N,L}
    ex = Expr(:block)
    push!(ex.args, :(x_mod = false))
    for label in L
        label = QuoteNode(label)
        ex_ss = quote
            x_mod = x_mod || f_disc!(sys.subsystems[$label], args...)
        end
        push!(ex.args, ex_ss)
    end
    return ex
end

# #replace this with the appropriate overloads, Plot recipes, whatever
# plotlog(log, sys::AbstractSystem) = extend_error(S)




end