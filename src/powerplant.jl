module Powerplant

using LinearAlgebra
using StaticArrays
using ComponentArrays
using UnPack

using Flight.Dynamics
using Flight.AirData
using Flight.System
import Flight.System: x_template, u_template, d_template, f_output!
import Flight.System: y_type

export SimpleProp, Gearbox, ElectricMotor
export EThruster, EThrusterX, EThrusterU, EThrusterY, EThrusterD, EThrusterSys

export EPowerplant, EPowerplantSys

export Group

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

Base.@kwdef struct EThruster <: AbstractComponent
    frame::FrameTransform = FrameTransform()
    motor::ElectricMotor = ElectricMotor()
    gearbox::Gearbox = Gearbox()
    propeller::SimpleProp = SimpleProp()
end

function wrench(prop::SimpleProp, ω::Real, ::AirDataSensed) #air data just for interface demo
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

Base.@kwdef struct EThrusterY
    throttle::Float64 = 0.0
    ω_shaft::Float64 = 0.0
    ω_prop::Float64 = 0.0
    wr_Oc_c::Wrench = Wrench()
    wr_Ob_b::Wrench = Wrench()
    h_Gc_b::SVector{3, Float64} = SVector{3}(0,0,0)
end

Base.@kwdef struct EThrusterD #external data sources (other than control inputs)
    air::AirDataSensed = AirDataSensed()
end

x_template(::Type{EThruster}) = ComponentVector(ω_shaft = 0.0)
u_template(::Type{EThruster}) = ComponentVector(throttle = 0.0)
d_template(::Type{EThruster}) = EThrusterD() #only called by the System constructor when the component is simulated alone
y_type(::Type{EThruster}) = EThrusterY

const EThrusterXAxes = typeof(getaxes(x_template(EThruster)))
const EThrusterUAxes = typeof(getaxes(u_template(EThruster)))
const EThrusterX{D} = ComponentVector{Float64, D, EThrusterXAxes} where {D<:AbstractVector{Float64}}
const EThrusterU{D} = ComponentVector{Float64, D, EThrusterUAxes} where {D<:AbstractVector{Float64}}

function f_output!(ẋ::EThrusterX, x::EThrusterX, u::EThrusterU, ::Real, data::EThrusterD, thr::EThruster)

    @unpack frame, motor, propeller, gearbox = thr
    @unpack n, η = gearbox

    throttle = u.throttle
    ω_shaft = x.ω_shaft
    ω_prop = ω_shaft / n

    wr_Oc_c = wrench(propeller, ω_prop, data.air)
    wr_Ob_b = frame * wr_Oc_c

    M_eng_shaft = torque(motor, throttle, ω_shaft)
    M_air_prop = wr_Oc_c.M[1]

    ω_shaft_dot = (M_eng_shaft + M_air_prop/(η*n)) / (motor.J + propeller.J/(η*n^2))

    h_Gc_c = SVector(motor.J * ω_shaft + propeller.J * ω_prop, 0, 0)
    h_Gc_b = frame.q_bc * h_Gc_c

    #update out.ẋ
    ẋ.ω_shaft = ω_shaft_dot

    EThrusterY(throttle, ω_shaft, ω_prop, wr_Oc_c, wr_Ob_b, h_Gc_b)

end

EThrusterSys(thr::EThruster = EThruster(); kwargs...) =
    System.Continuous(thr; kwargs...)

end #module