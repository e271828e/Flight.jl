module Propulsion

using LinearAlgebra
using StaticArrays, ComponentArrays, StructArrays, RecursiveArrayTools
using Unitful
using UnPack
using Plots

using Flight.Airdata
using Flight.Dynamics
using Flight.Airframe
using Flight.ModelingTools
import Flight.Airframe: get_wr_b, get_hr_b
import Flight.ModelingTools: System, get_x0, get_y0, get_u0, get_d0, f_cont!, f_disc!

using Flight.Plotting
import Flight.Plotting: plots

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

function get_wrench(prop::SimpleProp, ω::Real, air::AirData) #air data just for interface demo
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

#maybe disallow default values to avoid subtle bugs when failing to change a
#constructor call in user code after changing the struct definition
Base.@kwdef struct EThrusterY
    throttle::Float64 = 0
    ω_shaft::Float64 = 0
    ω_prop::Float64 = 0
    i::Float64 = 0
    c_bat::Float64 = 0
    wr_c::Wrench = Wrench()
    wr_b::Wrench = Wrench()
    hr_b::SVector{3,Float64} = zeros(SVector{3})
end

Base.@kwdef mutable struct EThrusterU
    throttle::Float64 = 0.0
end

Base.@kwdef mutable struct EThrusterD end

get_x0(::EThruster) = copy(EThrusterXTemplate)
get_d0(::EThruster) = EThrusterD()
get_u0(::EThruster) = EThrusterU()
get_y0(::EThruster) = EThrusterY()


################ EThruster System ###################

function System(thr::EThruster, ẋ = get_x0(thr), x = get_x0(thr),
                        y = get_y0(thr), u = get_u0(thr), d = get_d0(thr), t = Ref(0.0))
    params = thr #params is the component itself
    subsystems = nothing #no subsystems to define
    System{map(typeof, (thr, x, y, u, d, params, subsystems))...}(ẋ, x, y, u, d, t, params, subsystems)
end

get_wr_b(sys::System{EThruster}) = sys.y.wr_b
get_hr_b(sys::System{EThruster}) = sys.y.hr_b

f_disc!(sys::System{EThruster}) = false

function f_cont!(sys::System{EThruster}, air::AirData)

    @unpack ẋ, x, y, u, params = sys #no need for subsystems
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

    sys.y = EThrusterY(throttle, ω_shaft, ω_prop, i, c_bat, wr_c, wr_b, hr_b)

    return nothing

end

function plots(t, data::AbstractVector{<:EThrusterY}; mode, save_path, kwargs...)

    @unpack wr_c, wr_b = StructArray(data)
    # sa_wr_c = StructArray(wr_c)
    # sa_wr_b = StructArray(wr_b)

    # plt_F_c = thplot(t, sa_wr_c.F;
    #     ylabel = L"$F \ (N)$",
    #     plot_title = "Thruster Force [Thruster Frame]",
    #     th_split = :h,
    #     kwargs...)

    # plt_M_c = thplot(t, sa_wr_c.M;
    #     ylabel = L"$M \ (Nm)$",
    #     plot_title = "Thruster Moment [Thruster Frame]",
    #     th_split = :h,
    #     kwargs...)

    # savefig(plt_F_c, joinpath(save_path, "F_Oc_c.png"))
    # savefig(plt_M_c, joinpath(save_path, "M_Oc_c.png"))

    pd = Dict{String, Plots.Plot}()

    pd["01_wr_Oc_c"] = thplot(t, wr_c;
        plot_title = "Thruster Wrench [Thruster Frame]",
        wr_source = "thr", wr_frame = "c",
        kwargs...)

    pd["02_wr_Ob_b"] = thplot(t, wr_b;
        plot_title = "Thruster Wrench [Airframe]",
        wr_source = "thr", wr_frame = "b",
        kwargs...)

    save_plots(pd; save_path)

end

##################### Hand of God ##########################

# Base.@kwdef struct HandOfGod <: AbstractAirframeComponent
#     frame::FrameSpec = FrameSpec()
#     enabled::Bool = false
# end

# Base.@kwdef struct HandOfGodY
#     wr_c::Wrench = Wrench()
#     wr_b::Wrench = Wrench()
#     hr_b::SVector{3,Float64} = zeros(SVector{3})
# end

# Base.@kwdef mutable struct HandOfGodU
#     F::MVector{3,Float64} = 0.0
#     M::MVector{3,Float64} = 0.0
# end

# #required to make EThruster compatible with System
# Base.@kwdef mutable struct EThrusterD end

# get_x0(::EThruster) = copy(EThrusterXTemplate)
# get_d0(::EThruster) = EThrusterD()
# get_u0(::EThruster) = EThrusterU()
# get_y0(::EThruster) = EThrusterY()

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