module Aerodynamics

# using LinearAlgebra
# using StaticArrays
# using UnPack


using Flight.Attitude
using Flight.Kinematics
using Flight.Dynamics
using Flight.AirData
using Flight.Airframe
import Flight.Airframe: get_wr_b, get_hr_b

using Flight.Plotting
import Flight.Plotting: plots


abstract type AbstractAerodynamics <: AbstractAirframeComponent end

#compute airflow angles at frame c from the c-frame aerodynamic velocity
function get_airflow_angles(v_wOc_c::AbstractVector{<:Real})
    return(
        α = atan(v_wOc_c[3], v_wOc_c[1]),
        β = atan(v_wOc_c[2], √(v_wOc_c[1]^2 + v_wOc_c[3]^2)))
end

#AbstractAerodynamics subtypes don't produce angular momentum
get_hr_b(::System{<:AbstractAerodynamics}) = zeros(SVector{3})

#################### SimpleDrag ######################

#a System{<:AbstractAerodynamics} must at least store its wr_c and wr_b in its
#y, so that it can be retrieved when get_wr_b is called. we always assume that
#f_cont! will be called before get_wr_b

#as an AirframeComponent, this needs to implement get_wr_b, get_hr_b
Base.@kwdef struct SimpleDrag <: AbstractAerodynamics
    frame::FrameSpec = FrameSpec()
    c_d::Float64 = 0.5
    A::Float64 = 1.0
end

Base.@kwdef struct SimpleDragY
    α_c::Float64 = 0.0
    β_c::Float64 = 0.0
    wr_c::Wrench = Wrench()
    wr_b::Wrench = Wrench()
end

#don't care about kinematics, control surfaces or terrain, so we drop them with args...
function f_cont!(sys::System{SimpleDrag}, air::AirData, args...)

    error("Implement this")
    v_wOc_b = v_wOb_b #use the airframe origin velocity
    v_wOc_c = frame.q_bc'(v_wOc_b)
    α_c, β_c = get_airflow_angles(v_wOc_c)

    #compute drag in wind axes, back to component axes

    #OJO: en general, para un Aerodynamics mas sofisticado lo primero que tendre
    #que hacer es reproyectar v_wOb_b en ejes c. despues, calcular α_c y β_c.
    #con eso, wr_c. despues, transformar de vuelta
end

f_disc!(::System{SimpleDrag}) = false



# having v_wOb_b in AirData, any AbstractComponent to which AirData is passed is
# free to transform v_wOb_b into its own frame. this may be particularly useful
# for an Aerodynamics component. if its frame is f(Oa, εa):

# v_wOb_b = v_eOb_b - v_ew_b
# v_eOa_b = v_eOb_b + ω_eb_b × r_ObOa_b
# v_wOa_b = v_eOa_b - v_ew_b = v_eOb_b + ω_eb_b × r_ObOa_b - v_ew_b
# v_wOa_b = v_wOb_b + ω_eb_b × r_ObOa_b
# v_wOa_a = q_ba'(v_wOa_b)

#generally, v_wOa ≈ v_wOb, so the whole conversion won't be necessary; at most,
#we may need to reproject v_wOb_b into v_wOb_a to comply with the axes used by
#the aerodynamics database


    #TURN THIS INTO A RECIPE FOR AIRFLOW ANGLES

    # pd["05_α_β"] = thplot(t, rad2deg.(hcat(α_b, β_b));
    #     plot_title = "Airflow Angles [Airframe]",
    #     label = "",
    #     title = ["Angle of Attack" "Sideslip Angle"],
    #     ylabel = [L"$\alpha_b \ (deg)$" L"$\beta_b \ (deg)$"],
    #     th_split = :h, link = :none,
    #     kwargs...)



end #module