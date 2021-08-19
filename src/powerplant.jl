module Powerplant

using LinearAlgebra
using StaticArrays
using ComponentArrays
using Unitful
using UnPack

using Flight.Component
import Flight.Component: get_wr_Ob_b, get_h_Gc_b
using Flight.Airdata

import Flight.System: X, Y, U, f_cont!, f_disc!

export SimpleProp, Gearbox, ElectricMotor, Battery, CW, CCW
export EThruster, PowerplantGroup

@enum TurnSense begin
    CW = 1
    CCW = -1
end

Base.@kwdef struct SimpleProp
    kF::Float64 = 2e-3
    kM::Float64 = 5e-5
    J::Float64 = 1.0
end

function wrench(prop::SimpleProp, ω::Real, air::AirDataY) #air data just for interface demo
    @unpack kF, kM = prop
    F_ext_Os_s = kF * ω^2 * SVector(1,0,0)
    M_ext_Os_s = -tanh(ω/1.0) * kM * ω^2 * SVector(1,0,0) #choose ω_ref = 1.0
    Wrench(F = F_ext_Os_s, M = M_ext_Os_s)
end

Base.@kwdef struct Gearbox
    n::Float64 = 1.0 #gear ratio
    η::Float64 = 1.0 #efficiency
end

Base.@kwdef struct ElectricMotor #defaults from Hacker Motors Q150-4M-V2
    i₀::Float64 = 6.78
    R::Float64 = 0.004
    kV::Float64 = 125u"rpm*V" |> upreferred |> ustrip #rad/s/V, using the same value for kM
    J::Float64 = 0.003 #kg*m^2 #ballpark figure, assuming a cylinder
    α::TurnSense = CW
end

back_emf(m::ElectricMotor, ω::Real) = Int(m.α) * ω / m.kV
torque(m::ElectricMotor, i::Real, ω::Real) = (Int(m.α) * i - tanh(ω) * m.i₀) / m.kV
R(m::ElectricMotor) = m.R

Base.@kwdef struct Battery
    n_cells::Int64 = 14 #number of cells in series
    V_cell::Float64 = 4.2 #no-load cell voltage at C=Cmax (V)
    R_cell::Float64 = 5e-3 #internal cell resistance (Ω)
    Cmax::Float64 = 50000u"mA*hr" |> upreferred |> ustrip #capacity (Coulomb)
end

voltage_curve(b::Battery, charge_ratio::Real) = 1 #c = charge_ratio
voltage_open(b::Battery, charge_ratio::Real) = b.n_cells * b.V_cell * voltage_curve(b, charge_ratio)
R(b::Battery) = b.n_cells * b.R_cell
ċ(b::Battery, i::Real) = -i/b.Cmax

################ Electric Thruster ###################

abstract type AbstractThruster <: AbstractComponent end

Base.@kwdef struct EThruster <: AbstractThruster
    frame::Frame = Frame()
    battery::Battery = Battery()
    motor::ElectricMotor = ElectricMotor()
    gearbox::Gearbox = Gearbox()
    propeller::SimpleProp = SimpleProp()
end

const EThrusterXTemplate = ComponentVector(ω_shaft = 0.0, c_bat = 1.0)
const EThrusterUTemplate = ComponentVector(throttle = 0.0)
const EThrusterYTemplate = ComponentVector(
    throttle = 0.0, ω_shaft = 0.0, ω_prop = 0.0, i = 0.0, c_bat = 1.0,
    wr_Oc_c = ComponentVector(Wrench()), wr_Ob_b = ComponentVector(Wrench()),
    h_Gc_b = zeros(3))

#if we wanted to dispatch on the specific ComponentVector subtype...
# const EThrusterX{D} = ComponentVector{Float64, D, typeof(getaxes(EThrusterXTemplate))} where {D<:AbstractVector{Float64}}
# const EThrusterU{D} = ComponentVector{Float64, D, typeof(getaxes(EThrusterUTemplate))} where {D<:AbstractVector{Float64}}
# const EThrusterY{D} = ComponentVector{Float64, D, typeof(getaxes(EThrusterYTemplate))} where {D<:AbstractVector{Float64}}

#AbstractSystem interface
X(::EThruster) = copy(EThrusterXTemplate)
U(::EThruster) = copy(EThrusterUTemplate)
Y(::EThruster) = copy(EThrusterYTemplate)

get_wr_Ob_b(y, ::EThruster) = y.wr_Ob_b
get_h_Gc_b(y, ::EThruster) = y.h_Gc_b

f_disc!(x, u, t, thr::EThruster) = false

function f_cont!(y, ẋ, x, u, t, thr::EThruster, air::AirDataY = Y(AirData()))

    @unpack frame, battery, motor, propeller, gearbox = thr
    @unpack n, η = gearbox
    @unpack ω_shaft, c_bat = x

    # ω_shaft = x.ω_shaft
    # c_bat = x.c_bat
    throttle = u.throttle

    ω_prop = ω_shaft / n
    wr_Oc_c = wrench(propeller, ω_prop, air)
    wr_Ob_b = frame * wr_Oc_c

    i = (throttle * voltage_open(battery, c_bat) - back_emf(motor, ω_shaft)) /
        (R(battery) + R(motor))

    M_eng_shaft = torque(motor, i, ω_shaft)
    M_air_prop = wr_Oc_c.M[1]

    h_Gc_c = SVector(motor.J * ω_shaft + propeller.J * ω_prop, 0, 0)
    h_Gc_b = frame.q_bc * h_Gc_c

    ω_shaft_dot = (M_eng_shaft + M_air_prop/(η*n)) / (motor.J + propeller.J/(η*n^2))
    ẋ.ω_shaft = ω_shaft_dot
    ẋ.c_bat = ċ(battery, i)

    @pack! y = throttle, ω_shaft, ω_prop, i, c_bat, wr_Oc_c, wr_Ob_b, h_Gc_b

    return nothing

end

struct PowerplantGroup{C} <: AbstractComponentGroup{C} end

function PowerplantGroup(nt::NamedTuple{L, T}  where {L, T<:NTuple{N,AbstractThruster} where {N}})
    PowerplantGroup{nt}()
end
#= #interestingly, this does not work:
PowerplantGroup(nt::NamedTuple{L, NTuple{N, T}  where {L,N,T<:NTuple{N,
Powerplant.AbstractThruster}}) = PowerplantGroup{nt}()

#the reason is probably that NamedTuple, unlike Tuple (and therefore NTuple) is
#NOT covariant. that is:
#(EThruster(), EThruster()) isa NTuple{N, AbstractThruster} where {N} = true
#however:
#(a=EThruster(), b=NThruster) isa NamedTuple{L, NTuple{N, AbstractThruster} where {N} = false
=#
end #module