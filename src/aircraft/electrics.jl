module Electrics

using LinearAlgebra
using StaticArrays, ComponentArrays, StructArrays, RecursiveArrayTools
using Unitful
using UnPack
using Plots

using Flight.Utils
using Flight.Systems
using Flight.Plotting

using Flight.Air
using Flight.Kinematics
using Flight.RigidBody

import Flight.Systems: init, f_cont!, f_disc!
import Flight.Plotting: make_plots
import Flight.RigidBody: MassTrait, WrenchTrait, AngularMomentumTrait, get_wr_b, get_hr_b

export EThruster

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

function get_wrench(prop::SimpleProp, ω::Real, air::AirflowData) #air data just for interface demo
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

Base.@kwdef struct EThruster <: SystemDescriptor
    frame::FrameTransform = FrameTransform()
    battery::Battery = Battery()
    motor::ElectricMotor = ElectricMotor()
    gearbox::Gearbox = Gearbox()
    propeller::SimpleProp = SimpleProp()
end

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

init(::SystemX, ::EThruster) = init_x(ω_shaft = 0.0, c_bat = 1.0)
init(::SystemU, ::EThruster) = EThrusterU()
init(::SystemY, ::EThruster) = EThrusterY()

################ EThruster System ###################

get_wr_b(sys::System{EThruster}) = sys.y.wr_b
get_hr_b(sys::System{EThruster}) = sys.y.hr_b

f_disc!(sys::System{EThruster}) = false

function f_cont!(sys::System{EThruster}, kin::KinematicData, air::AirflowData)

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
    hr_b = frame.q(hr_c)

    ω_shaft_dot = (M_eng_shaft + M_air_prop/(η*n)) / (motor.J + propeller.J/(η*n^2))
    ẋ.ω_shaft = ω_shaft_dot
    ẋ.c_bat = ċ(battery, i)

    sys.y = EThrusterY(throttle, ω_shaft, ω_prop, i, c_bat, wr_c, wr_b, hr_b)

    return nothing

end


############################### Electrics ######################################


function make_plots(th::TimeHistory{<:EThrusterY}; kwargs...)

    pd = OrderedDict{Symbol, Plots.Plot}()

    pd[:wr_Oc_c] = plot(th.wr_c;
        plot_title = "Thruster Wrench [Thruster Frame]",
        wr_source = "thr", wr_frame = "c",
        kwargs...)

    pd[:wr_Ob_b] = plot(th.wr_b;
        plot_title = "Thruster Wrench [Vehicle Axes]",
        wr_source = "thr", wr_frame = "b",
        kwargs...)

    return pd

end



##################### Hand of God ##########################

# Base.@kwdef struct HandOfGod <: SystemDescriptor
#     frame::FrameTransform = FrameTransform()
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

#
end #module