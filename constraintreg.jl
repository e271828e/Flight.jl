using StaticArrays

# Base.@kwdef struct FrictionRegulator{T <: Union{Float64, SVector}}
#     k_p::T = 5.0 #proportional gain
#     k_i::T = 400.0 #integral gain
#     k_l::T = 0.0 #integrator leak factor

#     FrictionRegulator(k_p::Real, k_i, k_l) = new{Float64}(k_p, k_i, k_l)

#     function FrictionRegulator(k_p::AbstractVector{<:Real}, k_i, k_l)
#         N = length(k_p)
#         new{SVector{N,Float64}}(k_p, k_i, k_l)
#     end
# end


struct FrictionRegulator{N}
    k_p::SVector{N,Float64} #proportional gain
    k_i::SVector{N,Float64} #integral gain
    k_l::SVector{N,Float64} #integrator leak factor
end

FrictionRegulator{N}(; k_p, k_i, k_l) where {N} = FrictionRegulator{N}(k_p, k_i, k_l)

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

init(::FrictionRegulator{N}, ::Float64) where {N} = zeros(N)
init(::FrictionRegulator{N}) where {N} = FrictionRegulatorY{N}()