module Stochastic

using ComponentArrays
using StaticArrays
using UnPack
using Random

using Flight.Engine.Systems

export DiscreteGaussianWN, SampledGaussianWN, OrnsteinUhlenbeck, GaussianDoubleIntegrator


################################################################################
########################### GaussianProcess ##################################

#An N-dimensional discrete-time Gaussian stochastic process
abstract type GaussianProcess{N} <: Component end

#randomize state
function Random.randn!(rng::AbstractRNG, sys::System{<:GaussianProcess};
                       σ::Union{Real, AbstractVector{<:Real}})
    !isnothing(sys.x) && (randn!(rng, sys.x); sys.x .*= σ)
end

function Systems.f_disc!(sys::System{<:GaussianProcess}, Δt::Real,
                         rng::AbstractRNG)
    randn!(rng, sys.u) #generate a N(0,1) sample and apply it to the System's input
    _f_disc!(sys, Δt)
end

function Systems.f_disc!(sys::System{<:GaussianProcess}, Δt::Real,
                         u::Union{Real, AbstractVector{<:Real}})
    sys.u .= u #apply a directly provided N(0,1) sample to the System's input
    _f_disc!(sys, Δt)
end


################################################################################
############################ SampledGaussianWN #################################

"""
Sampled continuous Gaussian white noise process

We consider a continuous white noise process with PSD = N0/2. Before being
sampled at a frequency f_s = 1/Δt, this process must be band-limited to the
Nyquist frequency W = f_s/2 = 1/(2Δt). The result is a discrete noise process
with variance σ² = N0 * W = (2PSD) * 1/(2Δt) = PSD / Δt [Kay, Chapter 17.8]
"""
Base.@kwdef struct SampledGaussianWN{N} <: GaussianProcess{N}
    PSD::SVector{N,Float64} = ones(SVector{N})
end

Systems.init(::SystemU, cmp::SampledGaussianWN{N}) where {N} = zeros(N)
Systems.init(::SystemY, cmp::SampledGaussianWN{N}) where {N} = zeros(SVector{N,Float64})

function _f_disc!(sys::System{<:SampledGaussianWN{N}}, Δt::Real) where {N}

    σ = .√(sys.params.PSD / Δt)
    u = SVector{N,Float64}(sys.u)
    sys.y = SVector{N,Float64}(σ .* u)

    return false #no x

end


################################################################################
########################### DiscreteGaussianWN #################################

"""
Discrete Gaussian white noise process
"""
Base.@kwdef struct DiscreteGaussianWN{N} <: GaussianProcess{N}
    σ::SVector{N,Float64} = ones(SVector{N})
end

Systems.init(::SystemU, cmp::DiscreteGaussianWN{N}) where {N} = zeros(N)
Systems.init(::SystemY, cmp::DiscreteGaussianWN{N}) where {N} = zeros(SVector{N,Float64})

function _f_disc!(sys::System{<:DiscreteGaussianWN{N}}, ::Real) where {N}

    u = SVector{N,Float64}(sys.u)
    sys.y = SVector{N,Float64}(sys.params.σ .* u)

    return false #no x

end


################################################################################
############################# OrnsteinUhlenbeck ################################

"""
An exact discretization of the Ornstein-Uhlenbeck process:
dx = -1/T_c * x * dt + k_w * dW

T_c is a time constant, W is the Wiener process and k_w is a noise power
constant, which can be interpreted as the square root PSD of the white noise
process k_w * dW/dt (dW/dt is unit-PSD continuous white noise)
"""
Base.@kwdef struct OrnsteinUhlenbeck{N} <: GaussianProcess{N}
    T_c::SVector{N,Float64} = ones(SVector{N}) #time constant
    k_w::SVector{N,Float64} = ones(SVector{N}) #noise PSD square root
end

σ²(sys::System{<:OrnsteinUhlenbeck}) = (sys.params.k_w.^2 .* sys.params.T_c/2)
σ(sys::System{<:OrnsteinUhlenbeck}) = sqrt.(σ²(sys))

Systems.init(::SystemU, cmp::OrnsteinUhlenbeck{N}) where {N} = zeros(N)
Systems.init(::SystemX, cmp::OrnsteinUhlenbeck{N}) where {N} = zeros(N)
Systems.init(::SystemY, cmp::OrnsteinUhlenbeck{N}) where {N} = zeros(SVector{N,Float64})

function _f_disc!(sys::System{<:OrnsteinUhlenbeck{N}}, Δt::Real) where {N}

    @unpack x, u, params = sys
    @unpack T_c, k_w = params

    α = exp.(-Δt ./ T_c)
    β = .√(σ²(sys) .* (1 .- α.^2))

    x .= α .* x .+ β .* u

    sys.y = SVector{N, Float64}(x)

    return true #x modified

end


################################################################################
############################# GaussianDoubleIntegrator ######################################

"""
Gaussian stochastic double integrator with embedded velocity-acceleration and
position-acceleration feedback.
"""
struct GaussianDoubleIntegrator{N} <: GaussianProcess{N}
    k_u::SVector{N,Float64} #noise gain
    k_av::SVector{N,Float64} #acceleration-velocity feedback (<0 stabilizes)
    k_ap::SVector{N,Float64} #acceleration-position feedback (<0 stabilizes)
end

function GaussianDoubleIntegrator{N}(; k_u::Real = 1.0, k_av::Real = 0., k_ap::Real = 0.) where {N}
    GaussianDoubleIntegrator{N}(map(x-> fill(x,N), (k_u, k_av, k_ap))...)
end

Base.@kwdef struct GaussianDoubleIntegratorY{N}
    a::SVector{N,Float64} = zeros(SVector{N})
    v::SVector{N,Float64} = zeros(SVector{N})
    p::SVector{N,Float64} = zeros(SVector{N})
end

Systems.init(::SystemU, cmp::GaussianDoubleIntegrator{N}) where {N} = zeros(N)
Systems.init(::SystemY, cmp::GaussianDoubleIntegrator{N}) where {N} = GaussianDoubleIntegratorY{N}()

function Systems.init(::SystemX, cmp::GaussianDoubleIntegrator{N}) where {N}
    ComponentVector(v = zeros(N), p = zeros(N))
end

function _f_disc!(sys::System{<:GaussianDoubleIntegrator{N}}, Δt::Real) where {N}

    @unpack x, u, params = sys
    @unpack k_u, k_av, k_ap = params

    (v, p, u) = map(SVector{N,Float64}, (x.v, x.p, u))

    a = k_av .* v .+ k_ap .* p .+ k_u .* u
    x.v += Δt .* a #broadcasted assignment .= allocates
    x.p += Δt .* v #broadcasted assignment .= allocates

    sys.y = GaussianDoubleIntegratorY(; a, v = SVector{N}(x.v), p = SVector{N}(x.p))

    return true #x modified

end

end #module