module Powerplant

using LinearAlgebra
using StaticArrays
using ComponentArrays
using UnPack

using Flight.Dynamics
using Flight.System
import Flight.System: init_x, init_u

export SimpleProp, ElectricMotor, ElectricThruster, ElectricPowerplant
export ElectricThrusterOutput
export init_x, init_u, f_update!, f_output, step!
export ElectricThrusterSystem, ElectricPowerplantSystem

@enum TurnSense begin
    CW = 1
    CCW = -1
end

Base.@kwdef struct SimpleProp
    kF::Float64 = 0.1
    kM::Float64 = 0.01
    J::Float64 = 1.0
end

Base.@kwdef struct Gearbox
    n::Float64 = 1.0 #gear ratio
    η::Float64 = 1.0 #efficiency
end

Base.@kwdef struct ElectricMotor #defaults from Hacker Motors Q150-4M-V2
    i₀::Float64 = 6.78
    R::Float64 = 0.004
    kV::Float64 = 13.09 #rad/s/V
    Vb::Float64 = 48
    J::Float64 = 0.003 #kg*m^2 #ballpark figure, assuming a cylinder
    s::TurnSense = CW
end

Base.@kwdef struct ElectricThruster <: System.Descriptor
    frame::FrameTransform = FrameTransform()
    motor::ElectricMotor = ElectricMotor()
    gearbox::Gearbox = Gearbox()
    propeller::SimpleProp = SimpleProp()
end

function wrench(prop::SimpleProp, ω::Real) # airdata should also be an input
    @unpack kF, kM = prop
    F_ext_Os_s = kF * ω^2 * SVector(1,0,0)
    M_ext_Os_s = -sign(ω) * kM * ω^2 * SVector(1,0,0)
    Wrench(F = F_ext_Os_s, M = M_ext_Os_s)
end

function torque(eng::ElectricMotor, throttle::Real, ω::Real)
    @unpack i₀, R, kV, Vb, s = eng
    V = i₀ * R + throttle * Vb
    return Int(s) * ((V - Int(s)*ω/kV) / R - i₀) / kV
end

struct ElectricThrusterOutput
    throttle::Float64
    ω_shaft_dot::Float64
    ω_shaft::Float64
    ω_prop::Float64
    wr_Oc_c::Wrench
    wr_Ob_b::Wrench
    h_Gc_b::SVector{3, Float64}
end

init_x(::ElectricThruster) = ComponentVector(ω_shaft = 0.0)
init_u(::ElectricThruster) = ComponentVector(throttle = 0.0)

#interface required by DifferentialEquations
f_update!(ẋ, x, p, t) = f_update!(ẋ, x, p.u, t, p.d)
f_output(x, t, integrator) = f_output(x, integrator.p.u, t, integrator.p.d)

function f_update!(ẋ, x, u, t, desc::ElectricThruster)
    out = f_output(x, u, t, desc)
    ẋ.ω_shaft = out.ω_shaft_dot
end

function f_output(x, u, t, desc::ElectricThruster)

    @unpack frame, motor, propeller, gearbox = desc
    @unpack n, η = gearbox

    ω_shaft = x.ω_shaft
    ω_prop = ω_shaft / n

    wr_Oc_c = wrench(propeller, ω_prop)
    wr_Ob_b = frame * wr_Oc_c

    M_eng_shaft = torque(motor, u.throttle, ω_shaft)
    M_air_prop = wr_Oc_c.M[1]

    ω_shaft_dot = (M_eng_shaft + M_air_prop/(η*n)) / (motor.J + propeller.J/(η*n^2))

    h_Gc_c = SVector(motor.J * ω_shaft + propeller.J * ω_prop, 0, 0)
    h_Gc_b = frame.q_bc * h_Gc_c

    ElectricThrusterOutput(u.throttle, ω_shaft_dot, ω_shaft, ω_prop, wr_Oc_c, wr_Ob_b, h_Gc_b)

end

ElectricThrusterSystem(d = ElectricThruster(); kwargs...) =
    System.Continuous(d, init_x(d), init_u(d), f_update!, f_output; kwargs...)


struct ElectricPowerplant{N} <: System.Descriptor
    labels::NTuple{N, Symbol}
    thrusters::NTuple{N, ElectricThruster}
    function ElectricPowerplant(nt::NamedTuple) #Dicts are not ordered, so they won't do
        N = length(nt)
        @assert eltype(nt) == ElectricThruster
        new{N}(keys(nt), values(nt))
    end
end
ElectricPowerplant() = ElectricPowerplant((main = ElectricThruster(),))

function init_x(p::ElectricPowerplant)
    #we know that x is non-empty for all NamedTuple components, but this allows
    #us to generalize to the case of an heterogeneous powerplant where some
    #elements have a state and others don't
    # labels = p.labels
    # blocks = init_x.(p.thrusters)
    # with_x = collect(blocks.!=nothing)
    # nt = NamedTuple{names[with_x]}(blocks[with_x])
    ComponentVector(NamedTuple{p.labels}(init_x.(p.thrusters)))
end

function init_u(p::ElectricPowerplant)
    # names = keys(p.thrusters)
    # blocks = init_u.(values(p.thrusters))
    # with_u = collect(blocks.!=nothing)
    # nt = NamedTuple{names[with_u]}(blocks[with_u])
    # ComponentVector(nt)
    ComponentVector(NamedTuple{p.labels}(init_u.(p.thrusters)))
end

function f_update!(ẋ, x, u, t, pwp::ElectricPowerplant)
    for (name, thruster) in zip(pwp.labels, pwp.thrusters)
        f_update!(getproperty(ẋ, name), getproperty(x, name), getproperty(u, name), t, thruster)
    end
end

function f_output(x, u, t, pwp::ElectricPowerplant{N}) where {N}
    # v = Vector{ElectricThrusterOutput}(undef, N)
    v = ElectricThrusterOutput[]
    for (name, thruster) in zip(pwp.labels, pwp.thrusters)
        push!(v, f_output(getproperty(x, name), getproperty(u, name), t, thruster))
        # v[1] = f_output(getproperty(x, name), getproperty(u, name), t, thruster)
    end
    return NamedTuple{pwp.labels}(tuple(v...))
    # return v
end

ElectricPowerplantSystem(d = ElectricPowerplant(); kwargs...) =
    System.Continuous(d, init_x(d), init_u(d), f_update!, f_output; kwargs...)

end #module