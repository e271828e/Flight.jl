module Friction

using StaticArrays
using LinearAlgebra
using UnPack
using Plots

using Flight.Utils
using Flight.Plotting
using Flight.Systems

import Flight.Systems: init, f_cont!, f_disc!
import Flight.Plotting: make_plots

export get_μ

########################### Parameters #############################

struct Parameters
    μ_s::Float64 #static friction coefficient
    μ_d::Float64 #dynamic friction coefficient (μ_d < μ_s)
    v_s::Float64 #static friction upper velocity threshold
    v_d::Float64 #dynamic friction lower velocity threshold (v_d > v_s)

    function Parameters(; μ_s = 0, μ_d = 0, v_s = 0, v_d = 1)
        @assert μ_s >= μ_d
        @assert v_s < v_d
        new(μ_s, μ_d, v_s, v_d)
    end
end

function get_μ(fr::Parameters, v::Real)
    @unpack v_s, v_d, μ_s, μ_d = fr
    κ_sd = clamp((norm(v) - v_s) / (v_d - v_s), 0, 1)
    return κ_sd * μ_d + (1 - κ_sd) * μ_s
end


######################### Regulator###################################

struct Regulator{N} <: SystemDescriptor
    k_p::SVector{N,Float64} #proportional gain
    k_i::SVector{N,Float64} #integral gain
    k_l::SVector{N,Float64} #integrator leak factor
end

Regulator{N}(; k_p = 5.0, k_i = 400.0, k_l = 0.2) where {N} = Regulator{N}(k_p, k_i, k_l)

function Regulator{N}(k_p::Real, k_i::Real, k_l::Real) where {N}
    Regulator{N}(fill(k_p, N), fill(k_i, N), fill(k_l, N))
end

Base.@kwdef struct RegulatorU{N}
    reset::MVector{N,Bool} = zeros(Bool, N)
end

Base.@kwdef struct RegulatorY{N}
    reset::SVector{N,Bool} = zeros(SVector{N, Bool}) #reset input
    v::SVector{N,Float64} = zeros(SVector{N}) #constraint velocity
    s::SVector{N,Float64} = zeros(SVector{N})  #constraint velocity integral
    α_p::SVector{N,Float64} = zeros(SVector{N}) #proportional constraint force scale factor
    α_i::SVector{N,Float64} = zeros(SVector{N}) #integral constraint force scale factor
    α_raw::SVector{N,Float64} = zeros(SVector{N}) #total scale factor, raw
    α::SVector{N,Float64} = zeros(SVector{N}) #total scale factor, clipped
    sat::SVector{N,Bool} = zeros(SVector{N, Bool}) #scale factor saturation flag
end

init(::Regulator{N}, ::SystemX) where {N} = zeros(N) #v regulator integrator states
init(::Regulator{N}, ::SystemY) where {N} = RegulatorY{N}()
init(::Regulator{N}, ::SystemU) where {N} = RegulatorU{N}()

f_cont!(sys::System{<:Regulator{1}}, v_in::Real) = f_cont!(sys, SVector{1, Float64}(v_in))

function f_cont!(sys::System{<:Regulator{N}}, v_in::AbstractVector{<:Real}) where {N}

    @unpack k_p, k_i, k_l = sys.params

    v = SVector{N, Float64}(v_in)
    s = SVector{N, Float64}(sys.x)
    reset = SVector{N, Bool}(sys.u.reset)

    α_p = -k_p .* v
    α_i = -k_i .* s
    α_raw = α_p + α_i #raw μ scaling
    α = clamp.(α_raw, -1, 1) #clipped μ scaling
    sat = abs.(α_raw) .> abs.(α) #saturated?
    sys.ẋ .= (v - k_l .* s) .* .!sat .* .!reset #if not, integrator accumulates

    sys.y = RegulatorY(; reset, v, s, α_p, α_i, α_raw, α, sat)

end

function f_disc!(sys::System{<:Regulator{N}}) where {N}

    x = sys.x
    reset = sys.u.reset

    #not vectorized to avoid allocations
    x_mod = false
    for i in 1:N
        x_tmp = x[i]
        x[i] *= !reset[i]
        x_mod = x_mod || (x[i] != x_tmp)
    end
    return x_mod

end

############################## Plotting ########################################

function make_plots(th::TimeHistory{<:RegulatorY}; kwargs...)

    pd = OrderedDict{Symbol, Plots.Plot}()

    splt_v = plot(th.v; title = "Velocity",
        ylabel = L"$v \ (m/s)$", kwargs...)

    splt_s = plot(th.s; title = "Velocity Integral",
        ylabel = L"$s \ (m)$", kwargs...)

    splt_α_p = plot(th.α_p; title = "Proportional Term",
        ylabel = L"$\alpha_p$", kwargs...)

    splt_α_i = plot(th.α_i; title = "Integral Term",
        ylabel = L"$\alpha_i$", kwargs...)

    splt_α_raw = plot(th.α_raw; title = "Raw Output",
        ylabel = L"$\alpha_{raw}$", kwargs...)

    splt_α = plot(th.α; title = "Clipped Output",
        ylabel = L"$\alpha$", kwargs...)

    splt_sat = plot(th.sat; title = "Saturation",
        ylabel = L"$S$", kwargs...)

    pd[:vs] = plot(splt_v, splt_s, splt_sat;
        plot_title = "Contact Point Kinematics",
        layout = (1,3),
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    pd[:pi] = plot(splt_α_p, splt_α_i, splt_sat;
        plot_title = "Proportional and Integral Terms",
        layout = (1,3),
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    pd[:output] = plot(splt_α_raw, splt_α, splt_sat;
        plot_title = "Regulator Output",
        layout = (1,3),
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    return pd

end


end #module