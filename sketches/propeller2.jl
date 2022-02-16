module Propulsion2

using UnPack
using Roots
using Plots
using Trapz
using StaticArrays
using StructArrays
using Interpolations
using LinearAlgebra

using Flight.Modeling, Flight.Misc
using Flight.Kinematics, Flight.Dynamics, Flight.Airdata

import Flight.Modeling: init_x, init_y, init_u, init_d, f_cont!, f_disc!
import Flight.Dynamics: MassTrait, WrenchTrait, AngularMomentumTrait, get_hr_b, get_wr_b
import Flight.Plotting: plots

export PropCoefficients, PropBlade, PropFamily, FixedPitch, VariablePitch

###############################################
# en f_cont! de pwp no hay que olvidar asignar u = omega_shaft / n



# error("Write some fucking tests")

const π² = π^2


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



abstract type AbstractCoefficient <: AbstractFunction{2} end

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

(c::AbstractCoefficient)(ζ; α) = get_value(c, ζ, α)

#ζ ∈ [ζ_h, 1]
Base.@kwdef struct PropBlade{C, P, A <: AbstractFunction{1}, L, D, S <: AbstractCoefficient}
    ζ_h::Float64 = 0.1 #hub diameter to blade diameter ratio
    c̃::C = EllipticFunction(0.075) #chord to diameter ratio c_d(ζ)
    p̃::P = ConstantFunction(0.9)#chord-line-pitch to diameter ratio k_c(ζ)
    α_0::A = ConstantFunction(deg2rad(-2.1)) #airfoil zero-lift angle of attack α_0(ζ)
    cL::L = DefaultCL() #airfoil lift coefficient cl(ζ, α) (α measured wrt zero-lift line)
    cD::D = DefaultCD() #airfoil drag coefficient cd(ζ, α)
    cL_α::S = DefaultCLSlope() #airfoil lift coefficient slope cl_α(ζ, α)
end

#############################################################################

abstract type PitchTrait end

struct FixedPitch <: PitchTrait end

Base.@kwdef struct VariablePitch <: PitchTrait
    bounds::NTuple{2, Float64} = (0.0, deg2rad(20))
end

#############################################################################

Base.@kwdef struct PropCoefficients{T}
    C_Fx::T
    C_Mx::T
    C_Fz_α::T
    C_Mz_α::T
    C_P::T
    η_p::T
end

#compute coefficients for k propeller blades with the given properties at a
#specific advance ratio and blade pitch offset setting
function PropCoefficients(k::Real, blade::PropBlade, J::Real, Δβc::Real; n_ζ = 201)

    @unpack ζ_h, c̃, p̃, α_0, cL, cD, cL_α = blade

    ζ = range(ζ_h, 1, length = n_ζ)
    β = atan.(p̃.(ζ) ./ (π*ζ)) - α_0.(ζ) .+ Δβc #aerodynamic pitch angle β = βc - α_0 + Δβc
    β_t = β[end] #aerodynamic pitch angle at the tip

    dC_Fx = similar(ζ); dC_Fz_α = similar(ζ)
    dC_Mx = similar(ζ); dC_Mz_α = similar(ζ)

    ε_i = 1.0 #1 appears to be a suitable initial guess

    for (i, (ζ, β)) in enumerate(zip(ζ, β))

        ε_inf = atan(J / (π*ζ))

        f = let k = k, β_t = β_t, c̃ = c̃, cL = cL, ζ = ζ, β = β, ε_inf = ε_inf
            ε_i -> induced_angle_eq(k, β_t, c̃, cL, ζ, β, ε_inf, ε_i)
        end

        ε_i = find_zero(f, ε_i) #start at the solution for the previous radial location

        ε = ε_inf + ε_i
        α = β - ε
        # @assert (α < π/2 && α > -π/3) "α out of bounds after solving for ε_i"
        # @show β |> rad2deg
        # @show ε_i |> rad2deg
        # @show α |> rad2deg

        kc̃ = k * c̃(ζ)
        ζ² = ζ^2; ζ³ = ζ^3
        cos_ε = cos(ε); sin_ε = sin(ε)
        cos²ε_i = cos(ε_i)^2; cos²ε_inf = cos(ε_inf)^2
        tan_ε_inf = tan(ε_inf); tan²ε_inf = tan_ε_inf^2

        #the aerodynamic torque on a CW prop is negative along x (hence the sign
        #change wrt Phillips). a positive angle of attack yields a negative
        #force along the propeller's z axis (hence the sign change wrt Phillips)
        let cL = cL(ζ; α), cD = cD(ζ; α), cL_α = cL_α(ζ; α)

            dC_Fx[i] = π²/4 * ζ² * kc̃ * cos²ε_i / cos²ε_inf * (cL * cos_ε - cD * sin_ε)
            dC_Mx[i] = -π²/8 * ζ³ * kc̃ * cos²ε_i / cos²ε_inf * (cD * cos_ε + cL * sin_ε)
            dC_Fz_α[i] = -π²/8 * ζ² * kc̃ * cos²ε_i * (2tan_ε_inf * (cD * cos_ε + cL * sin_ε) - tan²ε_inf * (cL * cos_ε - (cL_α + cD) * sin_ε))
            dC_Mz_α[i] = -π²/16 * ζ³ * kc̃ * cos²ε_i * (2tan_ε_inf * (cL * cos_ε - cD * sin_ε) + tan²ε_inf * ((cL_α + cD) * cos_ε + cL * sin_ε))

        end
    end

    C_Fx = trapz(ζ, dC_Fx)
    C_Mx = trapz(ζ, dC_Mx)
    C_Fz_α = trapz(ζ, dC_Fz_α)
    C_Mz_α = trapz(ζ, dC_Mz_α)
    C_P = 2π * C_Mx
    η_p = (C_Fx > 0 ? -J * C_Fx / C_P : 0.0)

    PropCoefficients(C_Fx, C_Mx, C_Fz_α, C_Mz_α, C_P, η_p)

end

function induced_angle_eq(k, β_t, c̃, cL, ζ, β, ε_inf, ε_i)
    α = β - ε_inf - ε_i
    k*c̃(ζ) / (8ζ) * cL(ζ; α) - acos((exp(-k*(1-ζ)/(2sin(β_t))))) * tan(ε_i) * sin(ε_inf + ε_i)
end

# function PropCoefficients(k::Int, blade::PropBlade, J_bounds, Δβc_bounds; opts...)
# end
#compute coefficient dataset for a fixed pitch propeller with k blades of the
#specified geometry
function PropCoefficients(::FixedPitch, k::Int, blade::PropBlade; opts...)

    n_ζ = get(opts, :n_ζ, 201)
    n_J = get(opts, :n_J, 101)
    J_max = get(opts, :J_max, 1.5)

    J_range = range(0, J_max, length = n_J)
    data = Vector{PropCoefficients}(undef, length(J_range))

    for (i, J) in enumerate(J_range)
        data[i] = PropCoefficients(k, blade, J, 0; n_ζ)
    end

    data_sa = data |> StructArray |> StructArrays.components

    neg_ratio = sum((data_sa.C_Fx) .< 0) / length(data_sa.C_Fx)
    if neg_ratio > 0.3
        println("Warning: $(neg_ratio * 100)% of dataset points have negative traction")
    end

    #flat extrapolation by default
    interps = [ LinearInterpolation(J_range, c, extrapolation_bc = Flat())  for c in data_sa]

    PropCoefficients(interps...)

end

function PropCoefficients(p::VariablePitch, k::Int, blade::PropBlade; opts...)

    n_ζ = get(opts, :n_ζ, 201)
    n_J = get(opts, :n_J, 101)
    n_Δβc = get(opts, :n_Δβc, 21)
    J_max = get(opts, :J_max, 1.5)

    J_range = range(0, J_max, length = n_J)
    Δβc_range = range(p.bounds[1], p.bounds[2], length = n_Δβc)
    data = Array{PropCoefficients{Float64}}(undef, (length(J_range), length(Δβc_range)) )

    for (j, Δβc) in enumerate(Δβc_range)
        for (i, J) in enumerate(J_range)
            data[i,j] = PropCoefficients(k, blade, J, Δβc; n_ζ)
        end
    end

    data_sa = data |> StructArray |> StructArrays.components

    neg_ratio = sum((data_sa.C_Fx) .< 0) / length(data_sa.C_Fx)
    if neg_ratio > 0.3
        println("Warning: $(neg_ratio * 100)% of dataset points have negative traction")
    end

    interps = [ LinearInterpolation((J_range, Δβc_range), c, extrapolation_bc = Flat()) #flat extrapolation by default
        for c in data_sa]

    PropCoefficients(interps...)

end

function PropCoefficients(data::PropCoefficients{<:Interpolations.Extrapolation}, J::Real)

    @unpack C_Fx, C_Mx, C_Fz_α, C_Mz_α, C_P, η_p = data

    PropCoefficients(   C_Fx = C_Fx(J), C_Mx = C_Mx(J), C_Fz_α = C_Fz_α(J),
                        C_Mz_α = C_Mz_α(J), C_P = C_P(J), η_p = η_p(J))
end

function PropCoefficients(data::PropCoefficients{<:Interpolations.Extrapolation}, J::Real, Δβc::Real)

    @unpack C_Fx, C_Mx, C_Fz_α, C_Mz_α, C_P, η_p = data

    PropCoefficients(  C_Fx = C_Fx(J, Δβc), C_Mx = C_Mx(J, Δβc),
                        C_Fz_α = C_Fz_α(J, Δβc), C_Mz_α = C_Mz_α(J, Δβc),
                        C_P = C_P(J, Δβc), η_p = η_p(J, Δβc))
end

############################################################################

#the number of blades affects the total blade section circulation at a given
#radial distance in the Goldstein's condition, and therefore the more blades
#we have, the larger the induced velocity will be, and the less is to be
#gained from adding further blades. in fact, while the traction coefficient
#increases with the number of blades, the propulsive efficiency decreases


struct PropFamily{P <: PitchTrait, B <: PropBlade, D <: PropCoefficients{<:Interpolations.Extrapolation}}
    k::Int #number of blades
    pitch::P
    blade::B
    data::D

    function PropFamily(; k = 2, blade = PropBlade(), pitch = FixedPitch(), dataset_opts...)
        data = PropCoefficients(pitch, k, blade; dataset_opts...)
        new{typeof(pitch), typeof(blade), typeof(data)}(k, pitch, blade, data)
    end

end

PropCoefficients(pf::PropFamily{FixedPitch}, J::Real) = PropCoefficients(pf.data, J)
PropCoefficients(pf::PropFamily{VariablePitch}, J::Real, Δβc::Real) = PropCoefficients(pf.data, J, Δβc)

##########################################################################

@enum TurnSense begin
    CW = 1
    CCW = -1
end

abstract type AbstractPropeller <: SystemDescriptor end

MassTrait(::System{<:AbstractPropeller}) = HasNoMass()
WrenchTrait(::System{<:AbstractPropeller}) = GetsExternalWrench()
AngularMomentumTrait(::System{<:AbstractPropeller}) = HasAngularMomentum()

struct Propeller{F <: PropFamily} <: AbstractPropeller
    d::Float64 #diameter
    Ixx::Float64 #axial moment of inertia, labelled Ixx to avoid confusion with advance ratio
    t_bp::FrameTransform
    sense::TurnSense
    family::F
end

function Propeller(; d = 2.0, Ixx = 0.3, t_bp = FrameTransform(), sense = CW, family = PropFamily())
    Propeller(d, Ixx, t_bp, sense, family)
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
init_u(::Type{<:Propeller{<:PropFamily{FixedPitch}}}) = nothing
init_u(::Type{<:Propeller{<:PropFamily{VariablePitch}}}) = VariablePitchU()
init_y(::Type{<:Propeller}) = PropellerY()


function PropCoefficients(sys::System{Propeller{PropFamily{FixedPitch}}}, J::Real, ::Real)
    PropCoefficients(sys.params.family, J)
end

function PropCoefficients(sys::System{<:Propeller{<:PropFamily{VariablePitch}}}, J::Real)

    Δβc = linear_scaling(sys.u.pitch_setting, sys.params.family.pitch.bounds)
    # @show Δβc
    PropCoefficients(sys.params.family, J, Δβc)
end

function f_cont!(sys::System{<:Propeller}, kin::KinData, air::AirData, ω::Real)

    @unpack d, t_bp, sense, family = sys.params

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

    #the dataset is computed for a CW propeller. by symmetry considerations, one
    #can reason that, for a CCW propeller, force coefficients remain the same,
    #while moment coefficients must change sign
    C_F = SVector{3,Float64}(C_Fx, C_Fy_β * β_p, C_Fz_α * α_p)
    C_M = Int(sense) * SVector{3,Float64}(C_Mx, C_My_β * β_p, C_Mz_α * α_p)

    ρ = air.ρ
    f = ω/2π; f² = f^2; f³ = f * f²
    d⁴ = d^4; d⁵ = d * d⁴

    F_Op_p = ρ * f² * d⁴ * C_F
    M_Op_p = ρ * f² * d⁵ * C_M
    P      = ρ * f³ * d⁵ * C_P

    wr_p = Wrench(F_Op_p, M_Op_p)
    wr_b = t_bp(wr_p)

    sys.y = PropellerY(; v_wOp_p, ω, J, wr_p, wr_b, P, η_p)

end

f_disc!(::System{<:Propeller}, args...) = false

get_wr_b(sys::System{<:Propeller}) = sys.y.wr_b

end #module