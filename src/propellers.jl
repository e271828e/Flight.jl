module Propellers

using UnPack
using Roots
using Plots
using Trapz
using StaticArrays
using StructArrays
using Interpolations
using LinearAlgebra

import Interpolations: knots, bounds

using Flight.Modeling, Flight.Misc
using Flight.Kinematics, Flight.Dynamics, Flight.Airdata

import Flight.Modeling: init_x, init_y, init_u, init_d, f_cont!, f_disc!
import Flight.Dynamics: MassTrait, WrenchTrait, AngularMomentumTrait, get_hr_b, get_wr_b
import Flight.Plotting: plots


export PropCoefficients, PropBlade, Propeller, FixedPitch, VariablePitch


const π² = π^2

################################################################################
############################## Propeller Blade #################################

abstract type AbstractFunction{N} end

get_value(f::AbstractFunction{N}, ::Vararg{Any, M}) where {N, M} =
    error("A $(typeof(f)) was called with $M arguments, requires $N")

get_value(f::AbstractFunction{N}, ::Vararg{Any, N}) where {N} =
    error("Method get_value not implemented for $(typeof(f))")

(d::AbstractFunction)(args...) = get_value(d, args...)

struct ConstantFunction <: AbstractFunction{1}
    a::Float64
end
get_value(d::ConstantFunction, ::Real) = d.a

Base.@kwdef struct EllipticFunction <: AbstractFunction{1}
    a::Float64 = 0.075
end
get_value(d::EllipticFunction, ζ::Real) = d.a*√(1 - ζ^2)


abstract type AbstractAirfoil end

struct DefaultAirfoil <: AbstractAirfoil end

α_0(::DefaultAirfoil) = deg2rad(-2.1)

function cL(a::DefaultAirfoil, α::Real, M::Real = 0.0)
    if M <= 0.8
        (α < 0.25 ? 2π*α : π/2 * cos(α)/cos(0.25)) / √(1 - M^2)
    elseif M >= 1.2
        # 4α/√(M^2 - 1)
        # 4α/√(M^2 - 1)
        (α < 0.25 ? 4*α : cos(α)/cos(0.25)) / √(M^2 - 1)
    else
        cL(a, α, 0.8) + (cL(a, α, 1.2) - cL(a, α, 0.8)) / 0.4 * (M - 0.8) #recursive
    end
end

function cL_α(a::DefaultAirfoil, α::Real, M::Real = 0.0)
    if M <= 0.8
        (α < 0.25 ? 2π : -π/2 * sin(α)/cos(0.25)) / √(1 - M^2)
    elseif M >= 1.2
        (α < 0.25 ? 4 : -sin(α)/cos(0.25)) / √(M^2 - 1)
    else
        cL_α(a, α, 0.8) + (cL_α(a, 1.2) - cL_α(a, α, 0.8)) / 0.4 * (M - 0.8) #recursive
    end
end

function cD(::DefaultAirfoil, α::Real, M::Real = 0.0)

    # @assert M <= 1.2 "cD data only valid for M <= for 1.2"

    #incompressible
    if α < 0.25
        cD_inc = 0.006 + 0.224α^2
    elseif α < 0.3
        cD_inc = -1.0234 + 16.6944α^2
    else
        cD_inc = π/2 * sin(α)/cos(0.25)
    end

    #drag divergence correction
    if M <= 0.8
        κ_dd = 1.0
    elseif M <= 0.95
        κ_dd = 1.0 + 160000*(M-0.8)^4/27
    elseif M <= 1.0
        κ_dd = 6.0 - 800*(1 - M)^2
    else #not valid beyond M = 1.2
        κ_dd = 6 - 5(M - 1)
    end

    return κ_dd * cD_inc

end


#ζ ∈ [ζ_h, 1]
Base.@kwdef struct PropBlade{C, P, A <: AbstractAirfoil}
    ζ_h::Float64 = 0.1 #hub diameter to blade diameter ratio
    c̃::C = EllipticFunction(0.075) #chord to diameter ratio c_d(ζ)
    p̃::P = ConstantFunction(0.9)#chord-line-pitch to diameter ratio k_c(ζ)
    airfoil::A = DefaultAirfoil()
end

#pitch angle measured relative to the airfoil chord line
get_βc(b::PropBlade, ζ::Real, Δβ::Real) = atan(b.p̃(ζ) / (π*ζ)) + Δβ

#pitch angle measured relative to the airfoil zero-lift line (aerodynamic)
get_βa(b::PropBlade, ζ::Real, Δβ::Real) = get_βc(b, ζ, Δβ) - α_0(b.airfoil)

################################################################################
############################ PropCoefficients ##################################

Base.@kwdef struct PropCoefficients{T}
    C_Fx::T
    C_Mx::T
    C_Fz_α::T
    C_Mz_α::T
    C_P::T
    η_p::T
end

#compute coefficients for n_blades propeller blades with the given properties at a
#specific advance ratio and blade pitch offset setting. CW is assumed.
function PropCoefficients(blade::PropBlade, n_blades::Int; J::Real, M_t::Real, Δβ::Real, n_ζ = 201)

    airfoil = blade.airfoil

    ζ = range(blade.ζ_h, 1, length = n_ζ)
    βa_t = get_βa(blade, 1.0, Δβ) #add the pitch angle offset to the blade's own

    dC_Fx = similar(ζ); dC_Fz_α = similar(ζ)
    dC_Mx = similar(ζ); dC_Mz_α = similar(ζ)

    ε_i = 1.0 #1 appears to be a suitable initial guess

    for (i, ζ) in enumerate(ζ)

        ε_inf = atan(J / (π*ζ))
        βa = get_βa(blade, ζ, Δβ)
        c̃ = blade.c̃(ζ)

        f = let n_blades = n_blades, c̃ = c̃, airfoil = airfoil, βa_t = βa_t, J = J, M_t = M_t, βa = βa, ε_inf = ε_inf, ζ = ζ
            ε_i -> induced_angle_eq(; n_blades, c̃, airfoil, βa_t, J, M_t, βa, ε_inf, ζ, ε_i)
        end

        ε_i = find_zero(f, ε_i) #start at the solution for the previous radial location
        ε = ε_inf + ε_i
        α = βa - ε
        M = M_section(; J, M_t, ζ, ε_i)
        @assert (α < π/2 && α > -π/3) "α out of bounds after solving for ε_i"

        kc̃ = n_blades * c̃
        ζ² = ζ^2; ζ³ = ζ^3
        cos_ε = cos(ε); sin_ε = sin(ε)
        cos²ε_i = cos(ε_i)^2; cos²ε_inf = cos(ε_inf)^2
        tan_ε_inf = tan(ε_inf); tan²ε_inf = tan_ε_inf^2

        #the aerodynamic torque on a CW prop is negative along x (hence the sign
        #change wrt Phillips). a positive angle of attack yields a negative
        #force along the propeller's z axis (hence the sign change wrt Phillips)
        let cL = cL(airfoil, α, M), cD = cD(airfoil, α, M), cL_α = cL_α(airfoil, α, M)

            dC_Fx[i] = π²/4 * ζ² * kc̃ * cos²ε_i / cos²ε_inf * (cL * cos_ε - cD * sin_ε)
            dC_Mx[i] = -π²/8 * ζ³ * kc̃ * cos²ε_i / cos²ε_inf * (cD * cos_ε + cL * sin_ε)
            dC_Fz_α[i] = -π²/8 * ζ² * kc̃ * cos²ε_i * (2tan_ε_inf * (cD * cos_ε + cL * sin_ε) - tan²ε_inf * (cL * cos_ε - (cL_α + cD) * sin_ε))
            dC_Mz_α[i] = -π²/16 * ζ³ * kc̃ * cos²ε_i * (2tan_ε_inf * (cL * cos_ε - cD * sin_ε) + tan²ε_inf * ((cL_α + cD) * cos_ε + cL * sin_ε))

            # @show α, M, cL, cL_α, cD, dC_Fx[i], dC_Mx[i]

        end
    end

    C_Fx = trapz(ζ, dC_Fx)
    C_Mx = trapz(ζ, dC_Mx)
    C_Fz_α = trapz(ζ, dC_Fz_α)
    C_Mz_α = trapz(ζ, dC_Mz_α)
    C_P = 2π * C_Mx #all of these are for a CW propelller, we'll deal with the signs later
    η_p = (C_Fx > 0 ? -J * C_Fx / C_P : 0.0)

    PropCoefficients(C_Fx, C_Mx, C_Fz_α, C_Mz_α, C_P, η_p)

end

function M_section(; J, M_t, ζ, ε_i)
    J² = J^2
    ζ² = ζ^2
    return M_t * √((π² * ζ² + J²) / (π² + J²)) * cos(ε_i)
end

function induced_angle_eq(; n_blades, c̃, airfoil, βa_t, J, M_t, βa, ε_inf, ζ, ε_i)
    α = βa - ε_inf - ε_i
    M = M_section(; J, M_t, ζ, ε_i)
    return n_blades*c̃ / (8ζ) * cL(airfoil, α, M) - acos((exp(-n_blades*(1-ζ)/(2sin(βa_t))))) * tan(ε_i) * sin(ε_inf + ε_i)
end

################################################################################
################################ PitchControl ###################################

abstract type PitchControl end

struct FixedPitch <: PitchControl end

Base.@kwdef struct VariablePitch <: PitchControl
    bounds::NTuple{2, Float64} = (0.0, deg2rad(20))
end

################################################################################
################################ PropDataset ###################################

struct PropDataset{P <: PitchControl, T <: Interpolations.Extrapolation}
    _data::PropCoefficients{T}
end
PropDataset{P}(data::PropCoefficients{T}) where {P, T} = PropDataset{P, T}(data)

Base.getproperty(dataset::PropDataset, s::Symbol) = getproperty(dataset, Val(s))
@generated function Base.getproperty(dataset::PropDataset, ::Val{S}) where {S}
    if S === :data
        return :(getfield(dataset, :_data))
    elseif S ∈ fieldnames(PropCoefficients)
        return :(getfield(getfield(dataset, :_data), $(QuoteNode(S))))
    else
        error("PropDataset has no property $S")
    end
end

Interpolations.knots(dataset::PropDataset) = Interpolations.knots(dataset.C_Fx)
Interpolations.bounds(dataset::PropDataset) = Interpolations.bounds(dataset.C_Fx.itp)

#the number of blades affects the total blade section circulation at a given
#radial distance in the Goldstein's condition, and therefore the more blades
#we have, the larger the induced velocity will be, and the less is to be
#gained from adding further blades. in fact, while the traction coefficient
#increases with the number of blades, the propulsive efficiency decreases

#compute coefficient dataset for a fixed pitch propeller with n_blades blades of the
#specified geometry
function PropDataset(::FixedPitch, blade::PropBlade, n_blades::Int; opts...)

    n_ζ = get(opts, :n_ζ, 201)

    n_J = get(opts, :n_J, 41)
    J_max = get(opts, :J_max, 1.5)
    J_range = range(0, J_max, length = n_J)

    n_M_t = get(opts, :n_M_t, 41)
    M_t_max = get(opts, :M_t_max, 1.5)
    M_t_range = range(0, M_t_max, length = n_M_t)

    data = Array{PropCoefficients{Float64}}(undef, length(J_range), length(M_t_range))

    for (j, M_t) in enumerate(M_t_range)
        for (i, J) in enumerate(J_range)
            data[i,j] = PropCoefficients(blade, n_blades; J, M_t, Δβ=0, n_ζ)
        end
    end

    data_sa = data |> StructArray |> StructArrays.components

    interps = [ LinearInterpolation((J_range, M_t_range), c, extrapolation_bc = Flat())  for c in data_sa]
    data = PropCoefficients(interps...)

    neg_ratio = sum((data_sa.C_Fx) .< 0) / length(data_sa.C_Fx)
    if neg_ratio > 0.3
        println("Warning: $(neg_ratio * 100)% of dataset points have negative traction")
    end

    PropDataset{FixedPitch}(data)

end

#compute coefficient dataset for a variable pitch propeller with n_blades blades of the
#specified geometry
function PropDataset(p::VariablePitch, blade::PropBlade, n_blades::Int; opts...)

    n_ζ = get(opts, :n_ζ, 201)

    n_J = get(opts, :n_J, 41)
    J_max = get(opts, :J_max, 1.5)
    J_range = range(0, J_max, length = n_J)

    n_M_t = get(opts, :n_M_t, 41)
    M_t_max = get(opts, :M_t_max, 1.5)
    M_t_range = range(0, M_t_max, length = n_M_t)

    n_Δβ = get(opts, :n_Δβ, 5)
    Δβ_range = range(p.bounds[1], p.bounds[2], length = n_Δβ)

    @show J_range, M_t_range, Δβ_range

    data = Array{PropCoefficients{Float64}}(undef, (length(J_range), length(M_t_range), length(Δβ_range)) )

    for (k, Δβ) in enumerate(Δβ_range)
        for (j, M_t) in enumerate(M_t_range)
            for (i, J) in enumerate(J_range)
                data[i,j,k] = PropCoefficients(blade, n_blades; J, M_t, Δβ, n_ζ)
            end
        end
    end

    data_sa = data |> StructArray |> StructArrays.components

    interps = [ LinearInterpolation((J_range, M_t_range, Δβ_range), c, extrapolation_bc = Flat())
        for c in data_sa]
    data = PropCoefficients(interps...)

    neg_ratio = sum((data_sa.C_Fx) .< 0) / length(data_sa.C_Fx)
    if neg_ratio > 0.3
        println("Warning: $(neg_ratio * 100)% of dataset points have negative traction")
    end

    PropDataset{VariablePitch}(data)

end

function PropCoefficients(dataset::PropDataset{FixedPitch}, J::Real, M_tip::Real)

    @unpack C_Fx, C_Mx, C_Fz_α, C_Mz_α, C_P, η_p = dataset

    PropCoefficients(   C_Fx = C_Fx(J, M_tip), C_Mx = C_Mx(J, M_tip), C_Fz_α = C_Fz_α(J, M_tip),
                        C_Mz_α = C_Mz_α(J, M_tip), C_P = C_P(J, M_tip), η_p = η_p(J, M_tip))
end

function PropCoefficients(dataset::PropDataset{VariablePitch}, J::Real, M_tip::Real, Δβ::Real)

    @unpack C_Fx, C_Mx, C_Fz_α, C_Mz_α, C_P, η_p = dataset

    PropCoefficients(  C_Fx = C_Fx(J, M_tip, Δβ), C_Mx = C_Mx(J, M_tip, Δβ),
                        C_Fz_α = C_Fz_α(J, M_tip, Δβ), C_Mz_α = C_Mz_α(J, M_tip, Δβ),
                        C_P = C_P(J, M_tip, Δβ), η_p = η_p(J, M_tip, Δβ))
end

################################################################################
################################ Propeller #####################################

@enum TurnSense begin
    CW = 1
    CCW = -1
end

abstract type AbstractPropeller <: SystemDescriptor end

MassTrait(::System{<:AbstractPropeller}) = HasNoMass()
WrenchTrait(::System{<:AbstractPropeller}) = GetsExternalWrench()
AngularMomentumTrait(::System{<:AbstractPropeller}) = HasAngularMomentum()

struct Propeller{P <: PitchControl, B <: PropBlade,  D <: PropDataset{P}} <: AbstractPropeller
    pitch::P
    blade::B
    n_blades::Int #number of blades
    d::Float64 #diameter
    Ixx::Float64 #axial moment of inertia, labelled Ixx to avoid confusion with advance ratio
    sense::TurnSense
    t_bp::FrameTransform
    dataset::D
end

function Propeller(pitch = FixedPitch(), blade = PropBlade(), n_blades = 2,
                   d = 2.0, Ixx = 0.3, sense = CW, t_bp = FrameTransform(); dataset_opts...)

    dataset = PropDataset(pitch, blade, n_blades; dataset_opts...)

    Propeller{typeof(pitch), typeof(blade), typeof(dataset)}(
        pitch, blade, n_blades, d, Ixx, sense, t_bp, dataset)
end


Base.@kwdef struct PropellerY
    v_wOp_p::SVector{3,Float64} = zeros(SVector{3}) #local aerodynamic velocity, propeller axes
    ω::Float64 = 0 #angular velocity
    J::Float64 = 0 #advance ratio
    wr_p::Wrench = Wrench() #resulting aerodynamic Wrench, propeller frame
    wr_b::Wrench = Wrench() #resulting aerodynamic Wrench, airframe
    P::Float64 = 0.0 #power produced by the propeller
    η_p::Float64 = 0.0 #propulsive efficiency
end

Base.@kwdef mutable struct VariablePitchU
    pitch_setting::Bounded{Float64, 0, 1} = 0.0 #elevator control input (+ pitch down)
end
init_u(::Type{<:Propeller{FixedPitch}}) = nothing
init_u(::Type{<:Propeller{VariablePitch}}) = VariablePitchU()
init_y(::Type{<:Propeller}) = PropellerY()


function PropCoefficients(sys::System{<:Propeller{FixedPitch}}, J::Real)
    PropCoefficients(sys.params.dataset, J)
end

function PropCoefficients(sys::System{<:Propeller{VariablePitch}}, J::Real)

    Δβ = linear_scaling(sys.u.pitch_setting, sys.params.pitch.bounds)
    # @show Δβ
    PropCoefficients(sys.params.dataset, J, Δβ)
end

function f_cont!(sys::System{<:Propeller}, kin::KinData, air::AirData, ω::Real)

    @unpack d, t_bp, sense, dataset = sys.params

    v_wOp_b = air.v_wOb_b + kin.vel.ω_eb_b × t_bp.r
    v_wOp_p = t_bp.q'(v_wOp_b)

    #compute advance ratio. here we use the velocity vector magnitude rather
    #than its axial component. J must be positive, so we need the abs for CCW
    #propellers. also, we must prevent division by zero in a non-rotating
    #propeller
    abs_ω_min = 1.0
    v_J = norm(v_wOp_p)
    ω_J = max(abs(ω), abs_ω_min)
    J = 2π * v_J / (ω_J * d)

    coeffs = PropCoefficients(sys, J)

    @unpack C_Fx, C_Mx, C_Fz_α, C_Mz_α, C_P, η_p = coeffs
    C_Fy_β = C_Fz_α #by y/z symmetry
    C_My_β = C_Mz_α #by y/z symmetry

    α_p, β_p = get_airflow_angles(v_wOp_p)

    #datasets are computed for a CW propeller. by symmetry considerations, one
    #can reason that, for a CCW propeller, force coefficients remain the same,
    #while moment coefficients must change sign
    C_F = SVector{3,Float64}(C_Fx, C_Fy_β * β_p, C_Fz_α * α_p)
    C_M = Int(sense) * SVector{3,Float64}(C_Mx, C_My_β * β_p, C_Mz_α * α_p)

    ρ = air.ρ
    f = ω/2π; f² = f^2; f³ = f * f²
    d⁴ = d^4; d⁵ = d * d⁴

    F_Op_p = ρ * f² * d⁴ * C_F
    M_Op_p = ρ * f² * d⁵ * C_M
    P      = Int(sense) * ρ * abs(f³) * d⁵ * C_P

    wr_p = Wrench(F_Op_p, M_Op_p)
    wr_b = t_bp(wr_p)

    sys.y = PropellerY(; v_wOp_p, ω, J, wr_p, wr_b, P, η_p)

end

f_disc!(::System{<:Propeller}, args...) = false

get_wr_b(sys::System{<:Propeller}) = sys.y.wr_b

end #module