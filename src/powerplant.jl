module Powerplant

using LinearAlgebra
using StaticArrays
using ComponentArrays
using UnPack

using Flight.Airframe
using Flight.Airdata

# import ComponentArrays: ComponentVector
import Flight.System: X, Y, U, D, f_cont!, f_disc!
import Flight.Airframe: get_wr_Ob_b, get_h_Gc_b

export SimpleProp, Gearbox, ElectricMotor, CW, CCW
export EThruster, EThrusterX, EThrusterU, EThrusterY, EThrusterD, EThrusterSys

@enum TurnSense begin
    CW = 1
    CCW = -1
end

Base.@kwdef struct SimpleProp
    kF::Float64 = 0.1
    kM::Float64 = 0.1
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

################ Electric Thruster ###################

abstract type AbstractThruster <: AbstractComponent end

Base.@kwdef struct EThruster <: AbstractThruster
    frame::ComponentFrame = ComponentFrame()
    motor::ElectricMotor = ElectricMotor()
    gearbox::Gearbox = Gearbox()
    propeller::SimpleProp = SimpleProp()
end

const EThrusterXTemplate = ComponentVector(ω_shaft = 0.0)
const EThrusterUTemplate = ComponentVector(throttle = 0.0)
const EThrusterYTemplate = ComponentVector(
    throttle = 0.0, ω_shaft = 0.0, ω_prop = 0.0,
    wr_Oc_c = ComponentVector(Wrench()), wr_Ob_b = ComponentVector(Wrench()),
    h_Gc_b = zeros(3))

const EThrusterX{D} = ComponentVector{Float64, D, typeof(getaxes(EThrusterXTemplate))} where {D<:AbstractVector{Float64}}
const EThrusterU{D} = ComponentVector{Float64, D, typeof(getaxes(EThrusterUTemplate))} where {D<:AbstractVector{Float64}}
const EThrusterY{D} = ComponentVector{Float64, D, typeof(getaxes(EThrusterYTemplate))} where {D<:AbstractVector{Float64}}
const EThrusterD{D} = AirDataY{D} where {D}

#AbstractSystem interface
X(::EThruster) = copy(EThrusterXTemplate)
U(::EThruster) = copy(EThrusterUTemplate)
Y(::EThruster) = copy(EThrusterYTemplate)
D(::EThruster) = Y(AirData())

function wrench(prop::SimpleProp, ω::Real, ::EThrusterD) #air data just for interface demo
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

get_wr_Ob_b(y::EThrusterY, ::EThruster) = y.wr_Ob_b
get_h_Gc_b(y::EThrusterY, ::EThruster) = y.h_Gc_b

function f_cont!(y::EThrusterY, ẋ::EThrusterX, x::EThrusterX, u::EThrusterU, ::Real, data::EThrusterD, thr::EThruster)

    @unpack frame, motor, propeller, gearbox = thr
    @unpack n, η = gearbox

    ω_shaft = x.ω_shaft
    ω_prop = ω_shaft / n
    throttle = u.throttle
    M_eng_shaft = torque(motor, throttle, ω_shaft)

    wr_Oc_c = wrench(propeller, ω_prop, data)
    wr_Ob_b = frame * wr_Oc_c
    M_air_prop = wr_Oc_c.M[1]

    h_Gc_c = SVector(motor.J * ω_shaft + propeller.J * ω_prop, 0, 0)
    h_Gc_b = frame.q_bc * h_Gc_c

    ω_shaft_dot = (M_eng_shaft + M_air_prop/(η*n)) / (motor.J + propeller.J/(η*n^2))
    ẋ.ω_shaft = ω_shaft_dot

    @pack! y = throttle, ω_shaft, ω_prop, wr_Oc_c, wr_Ob_b, h_Gc_b

    return nothing

end


EThrusterSys(thr::EThruster = EThruster(); kwargs...) =
    System.Continuous(thr; kwargs...)

end #module