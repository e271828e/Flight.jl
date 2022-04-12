module Friction

using StaticArrays
using UnPack

using Flight.Modeling

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

struct FrictionRegulator{N} <: SystemDescriptor
    k_p::SVector{N,Float64} #proportional gain
    k_i::SVector{N,Float64} #integral gain
    k_l::SVector{N,Float64} #integrator leak factor
end

FrictionRegulator{N}(; k_p, k_i, k_l = 0.0) where {N} = FrictionRegulator{N}(k_p, k_i, k_l)

function FrictionRegulator{N}(k_p::Real, k_i::Real, k_l::Real) where {N}
    FrictionRegulator{N}(fill(k_p, N), fill(k_i, N), fill(k_l, N))
end

Base.@kwdef struct FrictionRegulatorY{N}
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

function f_cont!(sys::System{<:FrictionRegulator{N}}, v_in::Union{Real, AbstractVector{<:Real}}) where {N}

    @unpack k_p, k_i, k_l = sys.params

    v = SVector{N, Float64}(v_in)
    s = SVector{N, Float64}(sys.x)

    α_p = -k_p .* v
    α_i = -k_i .* s
    α_raw = α_p + α_i #raw μ scaling
    α = clamp.(α_raw, -1, 1) #clipped μ scaling
    sat = abs.(α_raw) .> abs.(α) #saturated?
    sys.ẋ .= (v - k_l .* s) .* .!sat #if not, integrator accumulates

    sys.y = FrictionRegulatorY(; v, s, α_p, α_i, α_raw, α, sat)

end

function f_disc!(sys::System{FrictionRegulator}, reset::Bool)

    if !reset
        return false
    else
        x_mod = (sys.x != 0.0 ? true : false)
        sys.x .= 0.0 #reset integrator states
        return x_mod
    end

end


end #module