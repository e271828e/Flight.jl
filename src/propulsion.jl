module Propulsion

using LinearAlgebra
using StaticArrays
using ComponentArrays
using StructArrays
using RecursiveArrayTools
using Unitful
using UnPack
using Plots

using Flight.Airdata
using Flight.Airframe
using Flight.System
import Flight.Airframe: get_wr_b, get_hr_b
import Flight.System: HybridSystem, get_x0, get_d0, get_u0, f_cont!, f_disc!

export SimpleProp, Gearbox, ElectricMotor, Battery, CW, CCW
export EThruster, EThrusterU, EThrusterD, EThrusterY

# abstract type AbstractThruster <: SystemDescriptor end

@enum TurnSense begin
    CW = 1
    CCW = -1
end

################ EThruster Component ###################

Base.@kwdef struct SimpleProp
    kF::Float64 = 2e-3
    kM::Float64 = 5e-5
    J::Float64 = 0.5
end

function get_wrench(prop::SimpleProp, ω::Real, air::AirY) #air data just for interface demo
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

Base.@kwdef struct EThruster <: AbstractAirframeComponent
    frame::FrameSpec = FrameSpec()
    battery::Battery = Battery()
    motor::ElectricMotor = ElectricMotor()
    gearbox::Gearbox = Gearbox()
    propeller::SimpleProp = SimpleProp()
end

const EThrusterXTemplate = ComponentVector(ω_shaft = 0.0, c_bat = 1.0)
const EThrusterX{T, D} = ComponentVector{T, D, typeof(getaxes(EThrusterXTemplate))} where {T,D}

#disallow default values to avoid subtle bugs when failing to change a
#constructor call in user code after changing the struct definition
Base.@kwdef struct EThrusterY <: AbstractY{EThruster}
    throttle::Float64
    ω_shaft::Float64
    ω_prop::Float64
    i::Float64
    c_bat::Float64
    wr_c::Wrench
    wr_b::Wrench
    hr_b::SVector{3,Float64}
end

Base.@kwdef mutable struct EThrusterU <: AbstractU{EThruster}
    throttle::Float64 = 0.0
end

#required to make EThruster compatible with HybridSystem
Base.@kwdef mutable struct EThrusterD <: AbstractD{EThruster} end

get_x0(::EThruster) = copy(EThrusterXTemplate)
get_d0(::EThruster) = EThrusterD()
get_u0(::EThruster) = EThrusterU()


################ EThruster HybridSystem ###################

function HybridSystem(thr::EThruster, ẋ = get_x0(thr), x = get_x0(thr),
                        d = get_d0(thr), u = get_u0(thr), t = Ref(0.0))
    params = thr #params is the component itself
    subsystems = nothing #no subsystems to define
    HybridSystem{map(typeof, (thr, x, d, u, params, subsystems))...}(ẋ, x, d, u, t, params, subsystems)
end

get_wr_b(y::EThrusterY) = y.wr_b
get_hr_b(y::EThrusterY) = y.hr_b

f_disc!(sys::HybridSystem{EThruster}) = false

function f_cont!(sys::HybridSystem{EThruster}, air::AirY)

    @unpack ẋ, x, u, params = sys #no need for subsystems
    @unpack frame, battery, motor, propeller, gearbox = params
    @unpack n, η = gearbox
    @unpack ω_shaft, c_bat = x

    throttle = u.throttle

    ω_prop = ω_shaft / n
    wr_c = get_wrench(propeller, ω_prop, air)
    wr_b = frame(wr_c)

    i = (throttle * voltage_open(battery, c_bat) - back_emf(motor, ω_shaft)) /
        (R(battery) + R(motor))

    M_eng_shaft = torque(motor, i, ω_shaft)
    M_air_prop = wr_c.M[1]

    hr_c = SVector(motor.J * ω_shaft + propeller.J * ω_prop, 0, 0)
    hr_b = frame.q_bc(hr_c)

    ω_shaft_dot = (M_eng_shaft + M_air_prop/(η*n)) / (motor.J + propeller.J/(η*n^2))
    ẋ.ω_shaft = ω_shaft_dot
    ẋ.c_bat = ċ(battery, i)

    return EThrusterY(throttle, ω_shaft, ω_prop, i, c_bat, wr_c, wr_b, hr_b)

end

"""
function System.rplot(t::AbstractVector{<:Real}, y::AbstractVector{EThrusterY}, args...)

    y_sa = StructArray(y)
    @unpack throttle, ω_shaft, i, c_bat, wr_c, wr_b, hr_b = y_sa


    #throttle,
    #battery charge: same plot, two subplots
    #wr_b: two subplots, one for F and another for M, all components

    throttle_plot = plot(t, throttle, title = "Throttle", xlabel = "u")
    # plot(t, ω_shaft)

    display(throttle_plot)

    hr_b = convert(Array, VectorOfArray(hr_b))

    #my use case is not really "i have a weird custom type that i want plotted
    #in a specific way" (this would be a User Recipe or a Type Recipe), but more
    #like "i have multiple custom types. each of which consists of fields, and
    #for each of these fields i want to generate a separate plot which will be
    #one in a set of predefined layouts". this is more of a Plot Recipeo

    #would like to have recipes to plot
    #a) A n-element vector in a single subplot, provide legend and colors
    #b) A n-element vector split in multiple subplots, vertical or horizontal
    #c) Multiple scalars in multiple subplots, vertical or horizontal

    #no, porque un Plot recipe supone que ya tiene como argumentos los datos
    #x,y,z de una serie. y en el caso de un n-element vector no es una serie,
    #son 3 realmente. y una User Recipe se puede usar para definir un cierto
    #tipo de layout. defino ese layout como un nuevo data type con @userplot
    #o a mano. o sea que si! SI QUE SON USER RECIPES lo que necesito.
    #porque son varias series.

    #ojo: A MENOS QUE dentro del plotattributes[:y] que recibe un Plot recipe si
    #que llegue un array de 3 filas. no. un Plot recibe recibe una serie, y una
    #serie es 1-D. http://docs.juliaplots.org/latest/input_data/#columns-are-series

end

"""





# function PropulsionGroup(nt::NamedTuple{L, T}  where {L, T<:NTuple{N,AbstractThruster} where {N}})
#     PropulsionGroup{nt}()
# end
#= #interestingly, this does not work:
PropulsionGroup(nt::NamedTuple{L, NTuple{N, T}  where {L,N,T<:NTuple{N,
Propulsion.AbstractThruster}}) = PropulsionGroup{nt}()

#the reason is that NamedTuple, unlike Tuple (and therefore NTuple) is
#NOT covariant. that is:
#(EThruster(), EThruster()) isa NTuple{N, AbstractThruster} where {N} == true
#however:
#(a=EThruster(), b=NThruster) isa NamedTuple{L, NTuple{N, AbstractThruster} where {N} == false
=#
end #module