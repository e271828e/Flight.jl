module Stochastic

using ComponentArrays
using StaticArrays
using UnPack
using Random
using DataStructures

using Flight.Engine.Systems

export StochasticProcess
export DiscreteGWN, SampledGWN, SampledGRW, SampledOU, DoubleIntegrator


################################################################################
########################### StochasticProcess ##################################

abstract type StochasticProcess <: Component end

function Random.randn!(sys::System{<:StochasticProcess}, args...)
    randn!(Random.default_rng(), sys, args...)
end

function Random.randn!(rng::AbstractRNG,
                       sys::System{<:StochasticProcess, X};
                       x̄::Union{Real, AbstractVector{<:Real}},
                       σ::Union{Real, AbstractVector{<:Real}}) where {X <: AbstractVector{Float64}}
    randn!(rng, sys.x)
    sys.x .*= σ
    sys.x .+= x̄
end

#the System's input is driven by standard Gaussian white noise
function Systems.f_disc!(sys::System{<:StochasticProcess, X, Y, U},
                         Δt::Real, rng::AbstractRNG) where {X, Y, U <: AbstractVector{Float64}}
    randn!(rng, sys.u) #generate a N(0,1) sample and apply it to the System's input
    f_disc!(sys, Δt)
end

function Systems.f_disc!(sys::System{<:StochasticProcess, X, Y, U}, Δt::Real,
                         u::Union{Real, AbstractVector{<:Real}}) where {X, Y, U <: AbstractVector{Float64}}
    sys.u .= u #apply a directly provided N(0,1) sample to the System's input
    f_disc!(sys, Δt)
end

################################################################################
########################### DiscreteGWN #################################

"""
Discrete Gaussian white noise process
"""
Base.@kwdef struct DiscreteGWN{N} <: StochasticProcess
    σ::SVector{N,Float64} = ones(SVector{N})
end

Systems.init(::SystemU, cmp::DiscreteGWN{N}) where {N} = zeros(N)
Systems.init(::SystemY, cmp::DiscreteGWN{N}) where {N} = zeros(SVector{N,Float64})

@inline σ²(sys::System{<:DiscreteGWN}) = σ(sys).^2
@inline σ(sys::System{<:DiscreteGWN}) = sys.params.σ

function Systems.f_disc!(sys::System{<:DiscreteGWN}, ::Real)
    u = SVector{N,Float64}(sys.u)
    sys.y = σ(sys) .* u
    return false
end

################################################################################
############################ SampledGWN #################################

"""
Sampled continuous Gaussian white noise process

We consider a continuous white noise process with PSD = N0/2. Before being
sampled at a frequency f_s = 1/Δt, this process must be band-limited to the
Nyquist frequency W = f_s/2 = 1/(2Δt). The result is a discrete noise process
with variance σ² = N0 * W = (2PSD) * 1/(2Δt) = PSD / Δt [Kay, Chapter 17.8]
"""
Base.@kwdef struct SampledGWN{N} <: StochasticProcess
    PSD::SVector{N,Float64} = ones(SVector{N})
end

Systems.init(::SystemU, cmp::SampledGWN{N}) where {N} = zeros(N)
Systems.init(::SystemY, cmp::SampledGWN{N}) where {N} = zeros(SVector{N,Float64})

@inline σ²(sys::System{<:SampledGWN}, Δt::Real) = SVector(sys.params.PSD ./ Δt)
@inline σ(sys::System{<:SampledGWN}, Δt::Real) = .√(σ²(sys, Δt))

function Systems.f_disc!(sys::System{<:SampledGWN{N}}, Δt::Real) where {N}
    u = SVector{N,Float64}(sys.u)
    sys.y = σ(sys, Δt) .* u
    return false #no x
end


################################################################################
############################# SampledGRW ################################

"""
Sampled continuous Gaussian random walk

This is just a Wiener process scaled with a non-unit noise PSD.
"""
struct SampledGRW{N} <: StochasticProcess
    k_w::SVector{N,Float64} #noise PSD square root
end

SampledGRW{N}(; k_w::Real = 1.0) where {N} = SampledGRW{N}(fill(k_w, N))

Systems.init(::SystemU, cmp::SampledGRW{N}) where {N} = zeros(N)
Systems.init(::SystemX, cmp::SampledGRW{N}) where {N} = zeros(N)
Systems.init(::SystemY, cmp::SampledGRW{N}) where {N} = zeros(SVector{N,Float64})

function Systems.f_disc!(sys::System{<:SampledGRW{N}}, Δt::Real) where {N}

    @unpack x, u, params = sys
    @unpack k_w = params

    x .= x .+ Δt .* k_w .* u

    sys.y = SVector{N, Float64}(x)

    return true #x modified

end

################################################################################
############################# SampledOU ################################

"""
An exact discretization of the Ornstein-Uhlenbeck process:
dx = -1/T_c * x * dt + k_w * dW

T_c is a time constant, W is the Wiener process and k_w is a noise power
constant, which can be interpreted as the square root PSD of the white noise
process k_w * dW/dt (dW/dt is unit-PSD continuous white noise)
"""
struct SampledOU{N} <: StochasticProcess
    T_c::SVector{N,Float64} #time constant
    k_w::SVector{N,Float64} #noise PSD square root
end

function SampledOU{N}(; T_c::Real = 1.0, k_w::Real = 1.0) where {N}
    SampledOU{N}(map(x-> fill(x,N), (T_c, k_w))...)
end

#stationary variance and standard deviation
@inline σ²(sys::System{<:SampledOU}) = (sys.params.k_w.^2 .* sys.params.T_c/2)
@inline σ(sys::System{<:SampledOU}) = sqrt.(σ²(sys))

Systems.init(::SystemU, cmp::SampledOU{N}) where {N} = zeros(N)
Systems.init(::SystemX, cmp::SampledOU{N}) where {N} = zeros(N)
Systems.init(::SystemY, cmp::SampledOU{N}) where {N} = zeros(SVector{N,Float64})

function Systems.f_disc!(sys::System{<:SampledOU{N}}, Δt::Real) where {N}

    @unpack x, u, params = sys
    @unpack T_c, k_w = params

    α = exp.(-Δt ./ T_c)
    β = .√(σ²(sys) .* (1 .- α.^2))

    x .= α .* x .+ β .* u

    sys.y = SVector{N, Float64}(x)

    return true #x modified

end


################################################################################
############################# DoubleIntegrator ######################################

"""
Gaussian stochastic double integrator with embedded velocity-acceleration and
position-acceleration feedback.
"""
struct DoubleIntegrator{N} <: StochasticProcess
    k_u::SVector{N,Float64} #noise gain
    k_av::SVector{N,Float64} #velocity feedback gain (>0 stabilizes)
    k_ap::SVector{N,Float64} #position feedback gain (>0 stabilizes)
end

function DoubleIntegrator{N}(;
                k_u::Real = 1.0, k_av::Real = 0., k_ap::Real = 0.) where {N}
    DoubleIntegrator{N}(map(x-> fill(x,N), (k_u, k_av, k_ap, σ0_v, σ0_p))...)
end

Base.@kwdef struct DoubleIntegratorY{N}
    a::SVector{N,Float64} = zeros(SVector{N})
    v::SVector{N,Float64} = zeros(SVector{N})
    p::SVector{N,Float64} = zeros(SVector{N})
end

Systems.init(::SystemU, cmp::DoubleIntegrator{N}) where {N} = zeros(N)
Systems.init(::SystemY, cmp::DoubleIntegrator{N}) where {N} = DoubleIntegratorY{N}()
function Systems.init(::SystemX, cmp::DoubleIntegrator{N}) where {N}
    ComponentVector(v = zeros(N), p = zeros(N))
end

function Systems.f_disc!(sys::System{<:DoubleIntegrator{N}}, Δt::Real) where {N}

    @unpack x, u, params = sys
    @unpack k_u, k_av, k_ap = params

    (v, p, u) = map(SVector{N,Float64}, (x.v, x.p, u))

    a = -k_av .* v .- k_ap .* p .+ k_u .* u
    x.v += Δt .* a #broadcasted assignment .= allocates
    x.p += Δt .* v #broadcasted assignment .= allocates

    sys.y = DoubleIntegratorY(; a, v = SVector{N}(x.v), p = SVector{N}(x.p))

    return true #x modified

end

end #module