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
    log::SavedValues

    function System(desc::D, x₀, u₀, method, kwargs...) where {D<:SystemDescriptor}

        @assert typeof(x₀) <: x_type(D)
        @assert typeof(u₀) <: u_type(D)

        params = (u = u₀, d = desc)
        log = SavedValues(Float64, save_type(D))
        scb = SavingCallback(f_output, log)

        problem = ODEProblem{true}(f_update!, x₀, (0, Inf), params)
        integrator = init(problem, method; callback = scb, save_everystep = false, kwargs...)
        new{D}(integrator, log)
    end


        # makes sense to disable save_everystep.
        #many of the state vector variables themselves (q_el, q_lb...) aren't
        #particularly intuitive. we are already defining a saving callback that
        #outputs all the needed variables, including the values of the input
        #vector at that specific instant. the saving callback must call the
        #output function for the system. however, this can still be overridden
        #upon construction

end

System(desc::D; x = x₀(D), u = u₀(D), method = Tsit5(), kwargs...) where {D} =
    System(desc, x, u, method, kwargs...)

Base.getproperty(sys::System, s::Symbol) = getproperty(sys, Val(s))

Base.getproperty(sys::System, ::Val{:descriptor}) = sys.integrator.p.descriptor
Base.getproperty(sys::System, ::Val{:integrator}) = getfield(sys, :integrator)
Base.getproperty(sys::System, ::Val{:log}) = getfield(sys, :log)

#forward everything else to the integrator...
Base.getproperty(sys::System, ::Val{S}) where {S} = getproperty(getfield(sys, :integrator), S)

#...except for x and u (because DiffEqs calls the state u, instead of x)
Base.getproperty(sys::System, ::Val{:u}) = sys.integrator.p.u #input vector
Base.getproperty(sys::System, ::Val{:x}) = sys.integrator.u #state vector

f_update!(ẋ, x, p, t) = f!(ẋ, x, p.u, t, p.d)
f_output(x, t, integrator) = h(x, integrator.p.u, t, integrator.p.d)

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
    kF::Float64 = 0.1
    kM::Float64 = 0.01
    J::Float64 = 1.0
end


function Dynamics.Wrench(prop::SimpleProp, ω::Real) # airdata should also be an input
    @unpack kF, kM = prop
    F_ext_Os_s = kF * ω^2 * SVector(1,0,0)
    M_ext_Os_s = -sign(ω) * kM * ω^2 * SVector(1,0,0)
    Wrench(F = F_ext_Os_s, M = M_ext_Os_s)
end

Base.@kwdef struct ElectricMotor #defaults from Hacker Motors Q150-4M-V2
    i₀::Float64 = 6.78
    R::Float64 = 0.004
    kV::Float64 = 13.09 #rad/s/V
    Vb::Float64 = 48
    J::Float64 = 0.003 #kg*m^2 #ballpark figure, assuming a cylinder
    s::TurnSense = CW
end

function torque(eng::ElectricMotor, throttle::Real, ω::Real)
    @unpack i₀, R, kV, Vb, s = eng
    V = i₀ * R + throttle * Vb
    return Int(s) * ((V - Int(s)*ω/kV) / R - i₀) / kV
end

Base.@kwdef struct Gearbox
    n::Float64 = 1.0 #gear ratio
    η::Float64 = 1.0 #efficiency
end

Base.@kwdef struct ElectricPwp <: SystemDescriptor
    frame::FrameTransform = FrameTransform()
    engine::ElectricMotor = ElectricMotor()
    propeller::SimpleProp = SimpleProp()
    gearbox::Gearbox = Gearbox()
end

#the parametric type definition is useful to allow a ComponentVector holding
#either a view or an actual array as a method argument. this is in turn required
#if we want to pass block views from a parent component array to a method
#operating on it
const XElectricPwpTemplate = ComponentVector(ω_shaft = 0.0)
const XElectricPwpAxes = typeof(getaxes(XElectricPwpTemplate))
const XElectricPwp = ComponentVector{Float64, D, XElectricPwpAxes} where {D <: AbstractVector{Float64}}

const UElectricPwpTemplate = ComponentVector(throttle = 0.0)
const UElectricPwpAxes = typeof(getaxes(UElectricPwpTemplate))
const UElectricPwp = ComponentVector{Float64, D, UElectricPwpAxes} where {D <: AbstractVector{Float64}}

x_type(::Type{ElectricPwp}) = XElectricPwp
u_type(::Type{ElectricPwp}) = UElectricPwp
x₀(::Type{ElectricPwp}) = XElectricPwpTemplate
u₀(::Type{ElectricPwp}) = UElectricPwpTemplate
save_type(::Type{ElectricPwp}) = YElectricPwp

#a descriptor only has to provide its x template. then, system must construct

struct YElectricPwp
    throttle::Float64
    ω_shaft_dot::Float64
    ω_shaft::Float64
    ω_prop::Float64
    wr_Oc_c::Wrench
    wr_Ob_b::Wrench
    h_Gc_b::SVector{3, Float64}
end

function f!(ẋ::XElectricPwp, x::XElectricPwp, u::UElectricPwp, t, desc::ElectricPwp)
    out = h(x, u, t, desc)
    ẋ.ω_shaft = out.ω_shaft_dot
end

function h(x::XElectricPwp, u::UElectricPwp, t, desc::ElectricPwp) #output function dispatch

    @unpack frame, engine, propeller, gearbox = desc
    @unpack n, η = gearbox

    ω_shaft = x.ω_shaft
    ω_prop = ω_shaft / n

    wr_Oc_c = Wrench(propeller, ω_prop)
    wr_Ob_b = frame * wr_Oc_c

    M_eng_shaft = torque(engine, u.throttle, ω_shaft)
    M_air_prop = wr_Oc_c.M[1]

    ω_shaft_dot = (M_eng_shaft + M_air_prop/(η*n)) / (engine.J + propeller.J/(η*n^2))

    h_Gc_c = SVector(engine.J * ω_shaft + propeller.J * ω_prop, 0, 0)
    h_Gc_b = frame.q_bc * h_Gc_c

    YElectricPwp(u.throttle, ω_shaft_dot, ω_shaft, ω_prop, wr_Oc_c, wr_Ob_b, h_Gc_b)

end

#TRY UNITFUL OUT, EVALUATE PERFORMANCE.
#conclusions:
#1) performance of arrays with components of different dimensions seems to be
#   (understandably) bad

# https://discourse.julialang.org/t/building-large-modular-models-in-julia/20755/7