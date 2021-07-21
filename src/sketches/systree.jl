using Base: Float64, @kwdef
using Flight.LBV

@define_node XSubsystemA (a1 = LBVLeaf{3}, a2 = LBVLeaf{2})
@define_node XSubsystemB (b1 = LBVLeaf{1}, b2 = LBVLeaf{4})
@define_node XNodeSystem (a = XSubsystemA, b = XSubsystemB)

Base.@kwdef struct SubsystemAParams
    kA::Float64 = 145.0
end

Base.@kwdef struct SubsystemBParams
    kB::Float64 = 239.0
end

Base.@kwdef struct NodeSystemParams
    a::SubsystemAParams = SubsystemAParams()
    b::SubsystemBParams = SubsystemBParams()
    kS::Float64 = π
end

struct SubsystemA
    x::XSubsystemA{Float64}
    kA::Float64
    function SubsystemA(x::XSubsystemA{Float64} = XSubsystemA(),
                        params::SubsystemAParams = SubsystemAParams())
        new(x, params.kA)
    end
end

struct SubsystemB
    x::XSubsystemB{Float64}
    kB::Float64
    function SubsystemB(x::XSubsystemB{Float64} = XSubsystemB(),
                        params::SubsystemBParams = SubsystemBParams())
        new(x, params.kB)
    end
end

struct NodeSystem
    x::XNodeSystem{Float64}
    a::SubsystemA
    b::SubsystemB
    kS::Float64
    function NodeSystem(x::XNodeSystem{Float64} = XNodeSystem(),
                        params::NodeSystemParams = NodeSystemParams())
        a = SubsystemA(x.a, params.a)
        b = SubsystemB(x.b, params.b)
        kS = params.kS
        new(x, a, b, kS)
    end
end
#we must not accept externally defined SubsystemA or SubsystemB instances as
#inputs to NodeSystem. instead, these need to be instantiated internally by
#NodeSystem's constructor from their respective params descriptors, passing
#NodeSystem's x.a and x.b blocks as state vectors to the SubsystemA and
#SubsystemB constructors. if we wanted to do this a posteriori, systems would
#have to be mutable, which would probably lead to type instability. this is
#understandable. for example, NodeSystem.x is declared as XNodeSystem{Float64},
#which still leaves the last type parameter in LBV undefined. if we make
#NodeSystem mutable, we could assign to x values with different type parameters
#during runtime

function demo()
    #and now we can do...
    params_b = SubsystemBParams(π/3)
    #thanks to the macro-generated kw constructor we can set only some parameters:
    params_node = NodeSystemParams(b = params_b, kS = 31)
    #preallocate an external state vector
    x = XNodeSystem(zeros(length(XNodeSystem)))
    @show x
    #and pass it to the constructor
    sys = NodeSystem(x, params_node)
    @show sys
    #now, subsystem a holds block x.a, so we can modify x (AKA sys.x) from a
    sys.a.x.a1 .= rand(3)
    @show x
    #and this can be verified to be type stable:
    @code_warntype sys.a.x .+ 1
end

#now, the question is... for SubsystemA and SubsystemB, which are Leaf systems,
#do we really need to pass a Params type? or could we simply pass the actual
#parameters as individual arguments? well, we could, but since these parameters
#will be part of a subsystem parameter hierarchy, why not store them in an
#ad-hoc struct defined specifically for that Subsystem type.