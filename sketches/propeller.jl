using UnPack
using Roots
using Plots
using Flight
using Trapz

const π² = π^2

abstract type AbstractDistribution end
get_value(::T, args...) where {T<:AbstractDistribution} =
    error("method get_value not implemented for AbstractDistribution $T")
(d::AbstractDistribution)(args...) = get_value(d, args...)

struct ConstantDistribution <: AbstractDistribution
    a::Float64
end
get_value(d::ConstantDistribution, ::Real) = d.a

Base.@kwdef struct EllipticDistribution <: AbstractDistribution
    a::Float64 = 0.075
end
get_value(d::EllipticDistribution, ζ::Real) = d.a*√(1 - ζ^2)

abstract type AbstractCoefficient end
get_value(::T, args...) where {T<:AbstractCoefficient} =
    error("method get_value not implemented for AbstractCoefficient $T")
(c::AbstractCoefficient)(ζ; α) = get_value(c, ζ, α)

struct DefaultCL <: AbstractCoefficient end
get_value(::DefaultCL, ::Real, α::Real) = (α < 0.25 ? 2π*α : π/2 * cos(α)/cos(0.25))

struct DefaultCLSlope <: AbstractCoefficient end
get_value(::DefaultCLSlope, ::Real, α::Real) = (α < 0.25 ? 2π : -π/2 * sin(α)/cos(0.25))

struct DefaultCD <: AbstractCoefficient end
function get_value(::DefaultCD, ::Real, α::Real)
    if α < 0.25
        0.006 + 0.224α^2
    elseif α < 0.3
        -1.0234 + 16.6944α^2
    else
        π/2 * sin(α)/cos(0.25)
    end
end


abstract type AbstractPropellerFamily <: SystemDescriptor end
# a propeller family is characterized by a series of scale-independent parameters
# given as functions of non-dimensional radial parameter ζ ∈ [ζ_h, 1]

struct FixedPitchPropellerFamily{ C, P, A <: AbstractDistribution,
                                  L, D, S <: AbstractCoefficient} <: AbstractPropellerFamily
    k::Int #number of blades
    ζ_h::Float64 #hub to blade diameter ratio
    c̃::C #chord to diameter ratio distribution c_d(ζ), ĉb = k*c_d
    p̃::P #chord-line-pitch to diameter ratio distribution k_c(ζ)
    α_0::A #airfoil zero-lift angle of attack distribution α_0(ζ)
    cL::L #airfoil lift coefficient cl(ζ, α) (α measured wrt zero-lift line)
    cD::D #airfoil drag coefficient cd(ζ, α)
    cL_α::S #airfoil lift coefficient slope cl_α(ζ, α)
end

function FixedPitchPropellerFamily(; k = 2, ζ_h = 0.1,
        c̃ = EllipticDistribution(0.075), p̃ = ConstantDistribution(0.5),
        α_0 = ConstantDistribution(deg2rad(-2.1)),
        cL = DefaultCL(), cD = DefaultCD(), cL_α = DefaultCLSlope())
    FixedPitchPropellerFamily{map(typeof, (c̃, p̃, α_0, cL, cD, cL_α))...}(k, ζ_h, c̃, p̃, α_0, cL, cD, cL_α)
end

function f_ind(k, β_t, c̃, cL, ζ, β, ε_inf, ε_i)
    α = β - ε_inf - ε_i
    k*c̃(ζ) / (8ζ) * cL(ζ; α) - acos((exp(-k*(1-ζ)/(2sin(β_t))))) * tan(ε_i) * sin(ε_inf + ε_i)
end

function compute_coefficients_CW(pf::FixedPitchPropellerFamily, J::Real, n_ζ = 101)

    @unpack k, ζ_h, c̃, p̃, α_0, cL, cD, cL_α = pf

    β_t = atan(p̃(1) / π) - α_0(1)

    ζ = range(ζ_h, 1, length = n_ζ)
    dC_Fx = similar(ζ); dC_Fz_α = similar(ζ)
    dC_Mx = similar(ζ); dC_Mz_α = similar(ζ)
    dC_P = similar(ζ)

    f_ζ = let k = k, β_t = β_t, c̃ = c̃, cL = cL
        (ζ, β, ε_inf, ε_i) -> f_ind(k, β_t, c̃, cL, ζ, β, ε_inf, ε_i)
    end

    for (i, ζ) in enumerate(ζ)

        β = atan(p̃(ζ) / (π*ζ)) - α_0(ζ)
        ε_inf = atan(J / (π*ζ))

        f = let ζ = ζ, β = β, ε_inf = ε_inf
            ε_i -> f_ζ(ζ, β, ε_inf, ε_i)
        end

        ε_i = find_zero(f, 1) #1 appears to be a suitable initial guess
        # @show ε_i

        ε = ε_inf + ε_i
        α = β - ε

        kc̃ = k * c̃(ζ)
        ζ² = ζ^2; ζ³ = ζ^3
        cos_ε = cos(ε); sin_ε = sin(ε)
        cos²ε_i = cos(ε_i)^2; cos²ε_inf = cos(ε_inf)^2
        tan_ε_inf = tan(ε_inf); tan²ε_inf = tan_ε_inf^2

        let cL = cL(ζ; α), cD = cD(ζ; α), cL_α = cL_α(ζ; α)

            dC_Fx[i] = π²/4 * ζ² * kc̃ * cos²ε_i / cos²ε_inf * (cL * cos_ε - cD * sin_ε)

            #the aerodynamic torque on a CW prop is negative along x (hence the
            #sign change wrt Phillips)
            dC_Mx[i] = -π²/8 * ζ³ * kc̃ * cos²ε_i / cos²ε_inf * (cD * cos_ε + cL * sin_ε)

            #a positive angle of attack yields a negative force along the
            #propeller's z axis (hence the sign change wrt Phillips)
            dC_Fz_α[i] = -π²/8 * ζ² * kc̃ * cos²ε_i * (2tan_ε_inf * (cD * cos_ε + cL * sin_ε) - tan²ε_inf * (cL * cos_ε - (cL_α + cD) * sin_ε))

            dC_Mz_α[i] = -π²/16 * ζ³ * kc̃ * cos²ε_i * (2tan_ε_inf * (cL * cos_ε - cD * sin_ε) + tan²ε_inf * ((cL_α + cD) * cos_ε + cL * sin_ε))

            #the aerodynamic power on any prop is negative since it opposes its
            #angular velocity
            dC_P[i] = 2π*dC_Mx[i]

        end
    end

    C_Fx = trapz(ζ, dC_Fx)
    C_Mx = trapz(ζ, dC_Mx)
    C_P = trapz(ζ, dC_P)
    C_Fz_α = trapz(ζ, dC_Fz_α)
    C_Mz_α = trapz(ζ, dC_Mz_α)
    C_Fy_β = C_Fz_α
    C_My_β = C_Mz_α

    η_p = - J * C_Fx / C_P
    # @show dC_Fx
    # @show dC_Mx
    # @show dC_P

    @show C_Fx
    @show C_Mx

    @show C_Fz_α
    @show C_Mz_α

    @show C_P
    @show η_p



end

    #the number of blades affects the total blade section circulation at a given
    #radial distance in the Goldstein's condition, and therefore the more blades
    #we have, the larger the induced velocity will be, and the less is to be
    #gained from adding further blades. in fact, while the traction coefficient
    #increases with the number of blades, the propulsive efficiency decreases

function compute_coefficient_tables_CW()
    #when building the tables, start increasing J from 0 until C_Fx goes below a
    #minimum threshold, for example 1e-3. that last data point, we don't include
    #in the tables
end

struct PropellerPerformanceData end
#we will need a type parameter for the type of Interpolation we are using (the
#same for all coefficients)

abstract type TurnSense end
(::Type{T})(s::TurnSense) where {T} = convert(T, s)

struct CW <: TurnSense end
Base.convert(::Type{T}, ::CW) where {T<:Real} = T(1)

struct CCW <: TurnSense end
Base.convert(::Type{T}, ::CCW) where {T<:Real} = T(-1)

struct FixedPitchPropeller{S <: TurnSense} <: SystemDescriptor
    data::PropellerPerformanceData
    d::Float64 #diameter
end

#maybe don't carry all the family type, instead simply create a

#when using the propeller in a PistonPowerplant, we must require that the Sense type
#parameters of the propeller and the engine are the same

    # family::FixedPropellerFamily
    # d::Float64 #propeller diameter
    # s::TurnSense

    #where do we put the performance tables? do we compute them within the
    #propeller family constructor? if we do so, we only need to instantiate the
    #family once, then use within multiple sizes and instnaces of actual propellers



# function test_propeller_analysis()

#     prop = FixedPitchPropellerFamily()
#     compute_coefficients(prop, 0, 0, 0)
# end
