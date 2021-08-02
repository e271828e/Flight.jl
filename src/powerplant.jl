# module Powerplant

using LinearAlgebra
using StaticArrays
using ComponentArrays
using DifferentialEquations
using Unitful
using UnPack

using Flight.Dynamics

#=TODO:
2) add saving callback to system constructor
3)

=#

#export XTemplate, UTemplate, XType, UType, f!, get_outputs #this is done by each module defining
#a systemdescriptor

abstract type SystemDescriptor end

struct System{D}

    integrator::OrdinaryDiffEq.ODEIntegrator #in fact everything is stored within the integrator!

    function System(desc::D, x₀::X, u₀::U, method, kwargs...) where {D<:SystemDescriptor, X, U}
        @assert X <: XType(D)
        @assert U <: UType(D)
        params = (u = u₀, d = desc)
        problem = ODEProblem{true}(f!, x₀, (0, Inf), params)
        integrator = init(problem, method; save_everystep = false, kwargs...)
        new{D}(integrator)
    end

# if construction of the integrator needs knowledge of a particular function for
#that concrete System, we can always define a method that dispatches on the
#System type parameter and returns the function that it needs. for example, to
#provide the saving function for the saving callback, we can do
#get_output_function(System{ElectricPwp})

        #need to add saving callbacks!!!

        ### REMEMBER TO call U_MODIFIED!!!!!!!!!!!!!!

        # makes sense to disable save_everystep.
        #many of the state vector variables themselves (q_el, q_lb...) aren't
        #particularly intuitive. we are already defining a saving callback that
        #outputs all the needed variables, including the values of the input
        #vector at that specific instant. the saving callback must call the
        #output function for the system. however, this can still be overridden
        #upon construction

end

System(desc::D; x = XTemplate(D), u = UTemplate(D), method = Tsit5(), kwargs...) where {D} =
    System(desc, x, u, method, kwargs...)

Base.getproperty(sys::System, s::Symbol) = getproperty(sys, Val(s))

#forward everything to the integrator...
Base.getproperty(sys::System, ::Val{S}) where {S} = getproperty(getfield(sys, :integrator), S)
Base.getproperty(sys::System, ::Val{:integrator}) = getfield(sys, :integrator)
Base.getproperty(sys::System, ::Val{:descriptor}) = sys.integrator.p.descriptor

#...except for x and u (because DiffEqs calls the state u, instead of x)
Base.getproperty(sys::System, ::Val{:u}) = sys.integrator.p.u #input vector
Base.getproperty(sys::System, ::Val{:x}) = sys.integrator.u #state vector

#need to forward everything to integrator except x and u

#this is the signature required by the integrator. it will call the specific
#method defined by each System subtype
DifferentialEquations.step!(sys::System) = step!(sys.integrator)
DifferentialEquations.step!(sys::System, dt::Real) = step!(sys.integrator, dt, stop_at_tdt = true)


################## ELECTRIC POWERPLANT ######################


@enum TurnSense begin
    CW = 1
    CCW = -1
end

Base.@kwdef struct SimpleProp
    sense::TurnSense = CW
    Jxx::Float64 = 1.0
    kF::Float64 = 1
    kM::Float64 = 0.01
end

h_Gc_xc(prop::SimpleProp, ω::Real) = prop.Jxx * ω

function Dynamics.Wrench(prop::SimpleProp, ω::Real) # airdata should also be an input
    @unpack sense, kF, kM = prop
    F_ext_Os_s = kF * ω^2 * SVector(1,0,0)
    M_ext_Os_s = -Int(sense) * kM * ω^2 * SVector(1,0,0)
    Wrench(F = F_ext_Os_s, M = M_ext_Os_s)
end

Base.@kwdef struct ElectricMotor #defaults from Hacker Motors Q150-4M-V2
    i₀::Float64 = 6.78
    R::Float64 = 0.004
    kV::Float64 = 13.09 #rad/s/V
    Vb::Float64 = 48
    J::Float64 = 0.003 #kg*m^2 #ballpark figure, assuming a cylinder
end

function torque(eng::ElectricMotor, throttle::Real, ω::Real)
    @unpack i₀, R, kV, Vb = eng
    return ((throttle * Vb - ω / kV) / R - i₀) / kV
end

Base.@kwdef struct ElectricPwp <: SystemDescriptor
    frame::FrameTransform = FrameTransform()
    engine::ElectricMotor = ElectricMotor()
    propeller::SimpleProp = SimpleProp()
    gear_ratio::Float64 = 1.0 #propRPM / engRPM, typically <=1
end

#the parametric type definition is useful to allow a ComponentVector holding
#either a view or an actual array as a method argument. this is in turn required
#if we want to pass block views from a parent component array to a method
#operating on it
const XElectricPwpTemplate = ComponentVector(ω = 0.0)
const XElectricPwpAxes = typeof(getaxes(XElectricPwpTemplate))
const XElectricPwp = ComponentVector{Float64, D, XElectricPwpAxes} where {D <: AbstractVector{Float64}}
# XElectricPwpType() = similar(XElectricPwpTemplate)

const UElectricPwpTemplate = ComponentVector(throttle = 0.0)
const UElectricPwpAxes = typeof(getaxes(UElectricPwpTemplate))
const UElectricPwp = ComponentVector{Float64, D, UElectricPwpAxes} where {D <: AbstractVector{Float64}}

XType(::Type{ElectricPwp}) = XElectricPwp
UType(::Type{ElectricPwp}) = UElectricPwp
XTemplate(::Type{ElectricPwp}) = XElectricPwpTemplate
UTemplate(::Type{ElectricPwp}) = UElectricPwpTemplate
# get_update_function(::Type{ElectricPwp}) = update!
# get_output_function(::Type{ElectricPwp}) = output

#this is the signature required by the integrator. it will call the specific
#method defined by each System subtype

# update!(ẋ, x, p, t) = f!(ẋ, x, p.u, p.d, t)
f!(ẋ, x, p, t) = f!(ẋ, x, p.u, t, p.d)

# saving callback method signature
h(x, t, integrator) = h(x, integrator.p.u, t, integrator.p.d)

function f!(ẋ::XElectricPwp, x::XElectricPwp, u::UElectricPwp, t, desc::ElectricPwp)
    ẋ.ω = 0.1
    # println(u, desc)
end
function h(x::XElectricPwp, u::UElectricPwp, t, desc::ElectricPwp) #output function dispatch
    println(x,u,t)
end

#TRY UNITFUL OUT, EVALUATE PERFORMANCE.
#conclusions:
#1) performance of arrays with components of different dimensions seems to be
#   (understandably) bad

# https://discourse.julialang.org/t/building-large-modular-models-in-julia/20755/7