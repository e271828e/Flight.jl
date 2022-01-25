module Aerodynamics

# using LinearAlgebra
using StaticArrays, ComponentArrays
using UnPack

using Flight.Modeling
using Flight.Plotting
using Flight.Attitude
using Flight.Kinematics
using Flight.Dynamics
using Flight.Airdata
using Flight.Components

import Flight.Modeling: init_x, init_y, init_u, init_d, f_cont!, f_disc!
import Flight.Plotting: plots
import Flight.Components: WrenchTrait, AngularMomentumTrait, get_wr_b

export AbstractAerodynamics, SimpleDrag


abstract type AbstractAerodynamics <: SystemDescriptor end

WrenchTrait(::System{<:AbstractAerodynamics}) = HasWrench()
AngularMomentumTrait(::System{<:AbstractAerodynamics}) = HasNoAngularMomentum()

#################### SimpleDrag ######################

Base.@kwdef struct SimpleDrag <: AbstractAerodynamics
    c_d::Float64 = 0.5
    A::Float64 = 1.0
end

Base.@kwdef struct SimpleDragY
    D::Float64 = 0.0
    wr_b::Wrench = Wrench()
end

init_y(::Type{SimpleDrag}) = SimpleDragY()

#in case of type instability use this one instead (or maybe a function barrier)
# function f_cont!(sys::System{SimpleDrag}, air::AirData, ::KinData, ::Any, ::Any)
function f_cont!(sys::System{SimpleDrag}, air::AirData, args...)

    @unpack α_b, β_b = air
    q_bw = get_wind_axes(α_b, β_b)

    D = 0.5 * air.q * sys.params.c_d * sys.params.A
    F_Ob_w = SVector{3,Float64}(-D, 0.0, 0.0)
    F_Ob_b = q_bw(F_Ob_w)

    #now the wrench is trivially computed
    wr_b = Wrench(F = F_Ob_b)

    sys.y = SimpleDragY(D, wr_b)

end

f_disc!(::System{SimpleDrag}) = false

Components.get_wr_b(sys::System{SimpleDrag}) = sys.y.wr_b

function plots(t, data::AbstractVector{SimpleDragY}; mode, save_path, kwargs...)

    @unpack D, wr_b = StructArray(data)

    pd = Dict{String, Plots.Plot}()

    pd["01_D"] = thplot(t, D;
        plot_title = "Aerodynamic Drag",
        label = "",
        ylabel = L"$D (N)$",
        kwargs...)

    pd["02_wr_Ob_b"] = thplot(t, wr_b;
        plot_title = "Aerodynamic Wrench [Airframe]",
        wr_source = "aero", wr_frame = "b",
        kwargs...)

    save_plots(pd; save_path)

end

#################### TestAerodynamics ################

#this subtype showcases the solutions to two potential scenarios:
#1) what to do if the aerodynamic data is expressed in a different frame (a)
#   than that chosen for the airframe (b), and we actually care about the
#   velocity lever arm between them

# answer: having v_wOb_b in AirData, any SystemDescriptor to which AirData is
# passed is free to transform v_wOb_b into its own frame. this may be
# particularly useful for an Aerodynamics component. if its frame is f(Oa, εa):

# v_wOb_b = v_eOb_b - v_ew_b
# v_eOa_b = v_eOb_b + ω_eb_b × r_ObOa_b
# v_wOa_b = v_eOa_b - v_ew_b = v_eOb_b + ω_eb_b × r_ObOa_b - v_ew_b
# v_wOa_b = v_wOb_b + ω_eb_b × r_ObOa_b
# v_wOa_a = q_ba'(v_wOa_b)

#generally, v_wOa ≈ v_wOb, so the whole conversion won't be necessary; at most,
#we will need to reproject v_wOb_b into v_wOb_a to comply with the axes used by
#the aerodynamics database, and from v_wOb_a we can compute airflow angles for
#those axes

#2) what to do if there is dependence of α_dot and β_dot in the aerodynamic
#   coefficients

# answer: we make the AerodynamicsSystem stateful, with filtered versions of α
# and β as states.

# we implement a first order filter as:
# ẋ.α = 1/τ * (α_in - x.α)
# ẋ.β = 1/τ * (β_in - x.β)

# now we have filtered versions of α_dot and β_dot available to use in the
# aerodynamics data.

# note that having strictly α_dot and β_dot as inputs (not their filtered
# counterparts) would make the complete world model implicit: α_dot and β_dot
# depend on v̇_eOb_b and v̇_ew_n, which in turn depend on aerodynamic forces and
# moments, which in turn depend on α_dot and β_dot. not only that, but this
# formulation would mean that any discontinuous changes in wind velocity would
# cause these derivatives to go to infinity. so... having filtered versions of
# these variables is the way to go.

#the remaining question is: what do we input as airflow angles to the
#aerodynamic database? the plain α_c, β_c inputs or the filtered states α, β?
#the first option seems better, because they are direct feedthrough, so unlike
#the filtered states, they do not introduce delays. all aerodynamic
#contributions not involving α_dot and β_dot are then delay-free. however, we
#should assume that in any case τ will have been chosen so that the bandwidth of
#the airflow angle filter extends well beyond that of aircraft dynamics



Base.@kwdef struct TestAerodynamics <: AbstractAerodynamics
    frame::FrameTransform = FrameTransform()
    τ::Float64 = 0.1 #time constant for first-order airflow angle filters
end

Base.@kwdef struct TestAerodynamicsY
    α_dot::Float64 = 0
    α_in::Float64 = 0
    α::Float64 = 0
    β_dot::Float64 = 0
    β_in::Float64 = 0
    β::Float64 = 0
    wr_c::Wrench = Wrench()
    wr_b::Wrench = Wrench()
end

init_x(::TestAerodynamics) = ComponentVector(α = 0.0, β = 0.0)
init_y(::TestAerodynamics) = TestAerodynamicsY()

function f_cont!(sys::System{TestAerodynamics}, air::AirData, kin::KinData, ::Any, ::Any)

    @unpack ẋ, x, y, params = sys
    @unpack α, β = x

    t_bc = params.frame

    # v_wOc_b = v_wOb_b #simply use the airframe origin velocity
    v_wOc_b = air.v_wOb_b + kin.ω_eb_b × t_bc.r
    v_wOc_c = t_bc.q'(v_wOc_b)
    α_in, β_in = get_airflow_angles(v_wOc_c) #airflow angles in component frame

    α_dot = 1/τ * (α_in - x.α)
    β_dot = 1/τ * (β_in - x.β)

    #here we would actually compute aerodynamic forces and moments using either
    #α_c, β_c or α_c_in and β_c_in
    wr_c = Wrench()
    wr_b = t_bc(wr_c)

    ẋ.α = α_dot
    ẋ.β = β_dot

    sys.y = TestAerodynamicsY(α_dot, α_in, α, β_dot, β_in, β, wr_c, wr_b)

end

end #module