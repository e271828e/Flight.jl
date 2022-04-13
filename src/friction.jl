module Friction

using StaticArrays
using LinearAlgebra
using UnPack

using Flight.Modeling
using Flight.Plotting

import Flight.Modeling: init, f_cont!, f_disc!
import Flight.Plotting: plots

export FrictionParameters, FrictionRegulator, get_μ

########################### FrictionParameters #############################

struct FrictionParameters
    μ_s::Float64 #static friction coefficient
    μ_d::Float64 #dynamic friction coefficient (μ_d < μ_s)
    v_s::Float64 #static friction upper velocity threshold
    v_d::Float64 #dynamic friction lower velocity threshold (v_d > v_s)

    function FrictionParameters(; μ_s, μ_d, v_s, v_d)
        @assert μ_s > μ_d
        @assert v_s < v_d
        new(μ_s, μ_d, v_s, v_d)
    end
end

function get_μ(fr::FrictionParameters, v::Real)
    @unpack v_s, v_d, μ_s, μ_d = fr
    κ_sd = clamp((norm(v) - v_s) / (v_d - v_s), 0, 1)
    return κ_sd * μ_d + (1 - κ_sd) * μ_s
end


######################### FrictionRegulator###################################

Base.@kwdef struct FrictionRegulator{N} <: SystemDescriptor
    k_p::SVector{N,Float64} #proportional gain
    k_i::SVector{N,Float64} #integral gain
    k_l::SVector{N,Float64} #integrator leak factor
end

# FrictionRegulator{N}(; k_p, k_i, k_l) where {N} = FrictionRegulator{N}(k_p, k_i, k_l)

function FrictionRegulator{N}(k_p::Real, k_i::Real, k_l::Real) where {N}
    FrictionRegulator{N}(fill(k_p, N), fill(k_i, N), fill(k_l, N))
end

Base.@kwdef struct FrictionRegulatorU{N}
    reset::MVector{N,Bool} = zeros(Bool, N)
end

Base.@kwdef struct FrictionRegulatorY{N}
    reset::SVector{N,Bool} = zeros(SVector{N, Bool}) #reset input
    v::SVector{N,Float64} = zeros(SVector{N}) #constraint velocity
    s::SVector{N,Float64} = zeros(SVector{N})  #constraint velocity integral
    α_p::SVector{N,Float64} = zeros(SVector{N}) #proportional constraint force scale factor
    α_i::SVector{N,Float64} = zeros(SVector{N}) #integral constraint force scale factor
    α_raw::SVector{N,Float64} = zeros(SVector{N}) #total scale factor, raw
    α::SVector{N,Float64} = zeros(SVector{N}) #total scale factor, clipped
    sat::SVector{N,Bool} = zeros(SVector{N, Bool}) #scale factor saturation flag
end

init(::FrictionRegulator{N}, ::SystemX) where {N} = zeros(N) #v regulator integrator states
init(::FrictionRegulator{N}, ::SystemY) where {N} = FrictionRegulatorY{N}()
init(::FrictionRegulator{N}, ::SystemU) where {N} = FrictionRegulatorU{N}()

function f_cont!(sys::System{<:FrictionRegulator{N}}, v_in::AbstractVector{<:Real}) where {N}

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

    sys.y = FrictionRegulatorY(; reset, v, s, α_p, α_i, α_raw, α, sat)

end

function f_disc!(sys::System{<:FrictionRegulator{N}}) where {N}

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

function plots(t, data::AbstractVector{<:FrictionRegulatorY}; mode, save_path, kwargs...)

    @unpack v, s, α_p, α_i, α_raw, α, sat = StructArray(data)

    pd = Dict{String, Plots.Plot}()

    splt_v = thplot(t, v; title = "Velocity",
        ylabel = L"$v \ (m/s)$", kwargs...)
    splt_s = thplot(t, s; title = "Velocity Integral",
        ylabel = L"$s \ (m)$", kwargs...)
    splt_α_p = thplot(t, α_p; title = "Proportional Term",
        ylabel = L"$\alpha_p$", label = "", kwargs...)
    splt_α_i = thplot(t, α_i; title = "Integral Term",
        ylabel = L"$\alpha_i$", label = "", kwargs...)
    splt_α_raw = thplot(t, α_raw; title = "Raw Output",
        ylabel = L"$\alpha_{raw}$", label = "", kwargs...)
    splt_α = thplot(t, α; title = "Clipped Output",
        ylabel = L"$\alpha$", label = "", kwargs...)
    splt_sat = thplot(t, sat; title = "Saturation",
        ylabel = L"$S$", label = "", kwargs...)

    pd["01_vs"] = plot(splt_v, splt_s;
        plot_title = "Contact Point Kinematics",
        layout = (1,2),
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    pd["02_pi"] = plot(splt_α_p, splt_α_i, splt_sat;
        plot_title = "Proportional and Integral Terms",
        layout = (1,3),
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    pd["03_pi"] = plot(splt_α_raw, splt_α, splt_sat;
        plot_title = "Regulator Output",
        layout = (1,3),
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    save_plots(pd; save_path)

end

end #module