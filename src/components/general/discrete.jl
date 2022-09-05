module Discrete

using ComponentArrays, StaticArrays, UnPack, LinearAlgebra

using Flight.Engine.Systems

################################################################################
############################# OrnsteinUhlenbeck ################################

"""
An exact discretization of the Ornstein-Uhlenbeck process:
dx = -1/T_c * x * dt + k_w * dW

Where T_c is a time constant, W is the Wiener process and k_w is a noise
power constant, which can be interpreted as the square root PSD of the white
noise process k_w * dW/dt (dW/dt is unit-PSD continuous white noise)

Input must be driven by a rng with standard normal distribution.
"""

Base.@kwdef struct OrnsteinUhlenbeck{N} <: Component
    T_c::SVector{N,Float64} = ones(SVector{N}) #time constant
    k_w::SVector{N,Float64} = ones(SVector{N}) #noise PSD square root
end

function stationary_var(sys::System{<:OrnsteinUhlenbeck})
    @unpack T_c, k_w = sys.params
    return k_w.^2 .* T_c/2
end

Systems.init(::SystemU, cmp::OrnsteinUhlenbeck{N}) where {N} = zeros(N)
Systems.init(::SystemX, cmp::OrnsteinUhlenbeck{N}) where {N} = zeros(N)
Systems.init(::SystemY, cmp::OrnsteinUhlenbeck{N}) where {N} = zeros(SVector{N,Float64})

function Systems.f_disc!(sys::System{<:OrnsteinUhlenbeck{N}}, Δt::Real) where {N}

    @unpack x, u, params = sys
    @unpack T_c, k_w = params

    σ²x = stationary_var(sys)
    α = exp.(-Δt ./ T_c)
    β = sqrt.(σ²x .* (1 .- α.^2))

    x .= α .* x .+ β .* u

    sys.y = SVector{N, Float64}(x)

    return true #x modified

end

################################################################################
############################# SecondOrder ######################################

"""
Discrete autoregressive second order model with velocity-acceleration and
position-acceleration feedback. When the System's u is fed by a rng, it works as
a discrete-time stochastic process, with k_u controlling the σ of each noise
component
"""

struct SecondOrder{N} <: Component
    k_u::SVector{N,Float64} #input gain
    k_av::SVector{N,Float64} #acceleration-velocity feedback (<0 stabilizes)
    k_ap::SVector{N,Float64} #acceleration-position feedback (<0 stabilizes)
end

function SecondOrder{N}(; k_u::Real = 1.0, k_av::Real = 0., k_ap::Real = 0.) where {N}
    SecondOrder{N}(map(x-> fill(x,N), (k_u, k_av, k_ap))...)
end

Base.@kwdef struct SecondOrderY{N}
    a::SVector{N,Float64} = zeros(SVector{N})
    v::SVector{N,Float64} = zeros(SVector{N})
    p::SVector{N,Float64} = zeros(SVector{N})
end

function Systems.init(::SystemX, cmp::SecondOrder{N}) where {N}
    ComponentVector(v = zeros(N), p = zeros(N))
end

Systems.init(::SystemU, cmp::SecondOrder{N}) where {N} = zeros(N)
Systems.init(::SystemY, cmp::SecondOrder{N}) where {N} = SecondOrderY{N}()

function Systems.f_disc!(sys::System{<:SecondOrder{N}}, Δt::Real) where {N}

    @unpack x, u, params = sys
    @unpack k_u, k_av, k_ap = params

    (v, p, u) = map(SVector{N,Float64}, (x.v, x.p, u))

    a = k_av .* v .+ k_ap .* p .+ k_u .* u
    x.v += Δt .* a #broadcasted assignment .= allocates
    x.p += Δt .* v #broadcasted assignment .= allocates

    sys.y = SecondOrderY(; a, v = SVector{N}(x.v), p = SVector{N}(x.p))

    return true #x modified

end



end #module