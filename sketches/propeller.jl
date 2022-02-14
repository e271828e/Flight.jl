module Propulsion2

using UnPack
using Roots
using Plots
using Trapz
using StaticArrays
using StructArrays
using Interpolations
using LinearAlgebra

using Flight.Modeling
using Flight.Kinematics, Flight.Dynamics, Flight.Airdata

import Flight.Modeling: init_x, init_y, init_u, init_d, f_cont!, f_disc!
import Flight.Dynamics: MassTrait, WrenchTrait, AngularMomentumTrait, get_hr_b, get_wr_b
import Flight.Plotting: plots

export FixedPitchFamily, FixedPitchPropeller


    ###############################################
    ###############################################
    ###############################################
    ###############################################
    ###############################################
    # en f_cont! de pwp no hay que olvidar asignar u = omega_shaft / n


"""
    como implementar VariablePitchPropeller
    0) VariablePitchFamily constara de una FixedPitchFamily mas un campo
       Δβc_range, que nos dice el incremento de βc (positivo o negativo) que se
       le puede aplicar
    1) VariablePitchDataset ira en funcion de dos variables, J y Δβ_c. pero ojo,
       ya no podemos simplemente cortar en un J para el que  C_Fx sea menor que
       un umbral, porque eso es para un Δβc determinado, para otro puede ser
       necesario seguir
    2) ahora VariablePitchPropeller sera igual que FixedPitchPropeller, solo que
       en dataset tendra VariablePitchDataset
    3) definimos una U que sea Bounded{0,1}, donde 0 corresponde a Δβc_range[1]
       y 1 a Δβc_range[2]
    4) en f_cont! extraemos u, y hacemos PropellerCoefficients(J, Δβc)
    5) a partir de ahi, todo es similar a FixedPitch
"""

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



abstract type AbstractPropellerFamily <: SystemDescriptor end
# a propeller family is characterized by a series of scale-independent parameters
# given as functions of non-dimensional radial parameter ζ ∈ [ζ_h, 1]

#the number of blades affects the total blade section circulation at a given
#radial distance in the Goldstein's condition, and therefore the more blades
#we have, the larger the induced velocity will be, and the less is to be
#gained from adding further blades. in fact, while the traction coefficient
#increases with the number of blades, the propulsive efficiency decreases

struct FixedPitchFamily{C, P, A <: AbstractFunction{1},
                        L, D, S <: AbstractCoefficient} <: AbstractPropellerFamily
    k::Int #number of blades
    ζ_h::Float64 #hub to blade diameter ratio
    c̃::C #chord to diameter ratio c_d(ζ)
    p̃::P #chord-line-pitch to diameter ratio k_c(ζ)
    α_0::A #airfoil zero-lift angle of attack α_0(ζ)
    cL::L #airfoil lift coefficient cl(ζ, α) (α measured wrt zero-lift line)
    cD::D #airfoil drag coefficient cd(ζ, α)
    cL_α::S #airfoil lift coefficient slope cl_α(ζ, α)
end

function FixedPitchFamily(; k = 2,
                            ζ_h = 0.1,
                            c̃ = EllipticFunction(0.075),
                            p̃ = ConstantFunction(0.9),
                            α_0 = ConstantFunction(deg2rad(-2.1)),
                            cL = DefaultCL(),
                            cD = DefaultCD(),
                            cL_α = DefaultCLSlope())

    FixedPitchFamily{map(typeof, (c̃, p̃, α_0, cL, cD, cL_α))...}(k, ζ_h, c̃, p̃, α_0, cL, cD, cL_α)

end



Base.@kwdef struct PropellerCoefficients #CW propeller
    C_Fx::Float64 = 0.0
    C_Mx::Float64 = 0.0 #normally, C_Mx < 0 (negative along x_p for a CW propeller)
    C_Fz_α::Float64 = 0.0
    C_Mz_α::Float64 = 0.0
    C_P::Float64 = 0.0 #normally, C_P < 0 (the propeller does negative work along the shaft, ω * Mx < 0)
    η_p::Float64 = 0.0
end

function PropellerCoefficients(pf::FixedPitchFamily, J::Real, n_ζ = 201)

    @unpack k, ζ_h, c̃, p̃, α_0, cL, cD, cL_α = pf

    β_t = atan(p̃(1) / π) - α_0(1)

    ζ = range(ζ_h, 1, length = n_ζ)
    dC_Fx = similar(ζ); dC_Fz_α = similar(ζ)
    dC_Mx = similar(ζ); dC_Mz_α = similar(ζ)

    ε_i = 1.0 #1 appears to be a suitable initial guess

    for (i, ζ) in enumerate(ζ)

        β = atan(p̃(ζ) / (π*ζ)) - α_0(ζ)
        ε_inf = atan(J / (π*ζ))

        f = let k = k, β_t = β_t, c̃ = c̃, cL = cL, ζ = ζ, β = β, ε_inf = ε_inf
            ε_i -> induced_angle_eq(k, β_t, c̃, cL, ζ, β, ε_inf, ε_i)
        end

        ε_i = find_zero(f, ε_i) #start at the solution for the previous radial location

        ε = ε_inf + ε_i
        α = β - ε

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
    η_p = (C_Fx > 0 ? -J * C_Fx / C_P : 0)

    PropellerCoefficients(; C_Fx, C_Mx, C_Fz_α, C_Mz_α, C_P, η_p)

end

function induced_angle_eq(k, β_t, c̃, cL, ζ, β, ε_inf, ε_i)
    α = β - ε_inf - ε_i
    k*c̃(ζ) / (8ζ) * cL(ζ; α) - acos((exp(-k*(1-ζ)/(2sin(β_t))))) * tan(ε_i) * sin(ε_inf + ε_i)
end

struct FixedPitchDataset{T <: Interpolations.Extrapolation}
    C_Fx::T
    C_Mx::T
    C_Fz_α::T
    C_Mz_α::T
    C_P::T
    η_p::T
end

function FixedPitchDataset(pf::FixedPitchFamily; n_ζ = 201, ΔJ = 0.01, J_max = 1.5, C_Fx_min = -0.02)

    data = Vector{PropellerCoefficients}()

    for J in range(0, J_max; step = ΔJ)
        coefs = PropellerCoefficients(pf, J, n_ζ)
        coefs.C_Fx < C_Fx_min ? break : nothing
        push!(data, coefs)
    end

    J_valid = range(0; step = ΔJ, length = length(data))

    interps = [ LinearInterpolation(J_valid, c, extrapolation_bc = Flat()) #flat extrapolation by default
        for c in data |> StructArray |> StructArrays.components]

    FixedPitchDataset(interps...)

end

function PropellerCoefficients(dataset::FixedPitchDataset, J::Real)

    @unpack C_Fx, C_Mx, C_Fz_α, C_Mz_α, C_P, η_p = dataset

    PropellerCoefficients(  C_Fx = C_Fx(J), C_Mx = C_Mx(J),
                            C_Fz_α = C_Fz_α(J), C_Mz_α = C_Mz_α(J),
                            C_P = C_P(J), η_p = η_p(J)
        )

end



@enum TurnSense begin
    CW = 1
    CCW = -1
end

abstract type AbstractPropeller <: SystemDescriptor end

MassTrait(::System{<:AbstractPropeller}) = HasNoMass()
WrenchTrait(::System{<:AbstractPropeller}) = GetsExternalWrench()
AngularMomentumTrait(::System{<:AbstractPropeller}) = HasAngularMomentum()



struct FixedPitchPropeller{D <: FixedPitchDataset} <: AbstractPropeller
    d::Float64 #diameter
    Ixx::Float64 #axial moment of inertia, labelled Ixx to avoid confusion with advance ratio
    t_bp::FrameTransform
    sense::TurnSense
    dataset::D
end

function FixedPitchPropeller(; d = 2.0, Ixx = 0.3, t_bp = FrameTransform(), sense = CW, family = FixedPitchFamily())
    FixedPitchPropeller(d, Ixx, t_bp, sense, FixedPitchDataset(family))
end

Base.@kwdef struct FixedPitchPropellerY
    v_wOp_p::SVector{3,Float64} = zeros(SVector{3}) #local aerodynamic velocity, propeller axes
    ω::Float64 = 0 #angular velocity
    J::Float64 = 0 #advance ratio
    α_p::Float64 = 0 #propeller angle of attack
    β_p::Float64 = 0 #propeller angle of sideslip
    coeffs::PropellerCoefficients = PropellerCoefficients()
    wr_p::Wrench = Wrench() #resulting aerodynamic Wrench, propeller frame
    wr_b::Wrench = Wrench() #resulting aerodynamic Wrench, airframe
    P::Float64 = 0.0 #power produced by the propeller
    η_p::Float64 = 0.0 #propulsive efficiency
end

init_y(::Type{<:FixedPitchPropeller}) = FixedPitchPropellerY()

function f_cont!(sys::System{<:FixedPitchPropeller}, kin::KinData, air::AirData, ω::Real)

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

    coeffs = PropellerCoefficients(dataset, J)

    @unpack C_Fx, C_Mx, C_Fz_α, C_Mz_α, C_P, η_p = coeffs
    C_Fy_β = C_Fz_α #by y/z symmetry
    C_My_β = C_Mz_α #by y/z symmetry

    α_p, β_p = get_airflow_angles(v_wOp_p)

    #the dataset is computed for a CW propeller. by symmetry considerations, one
    #can reason that, for a CCW propeller, force coefficients remain the same,
    #while moment coefficients must change sign
    C_F = SVector{3,Float64}(C_Fx, C_Fy_β * β_p, C_Fz_α * α_p)
    C_M = Int(sense) * SVector{3,Float64}(C_Fx, C_My_β * β_p, C_Mz_α * α_p)

    ρ = air.ρ
    f = ω/2π; f² = f^2; f³ = f * f²
    d⁴ = d^4; d⁵ = d * d⁴

    F_Op_p = ρ * f² * d⁴ * C_F
    M_Op_p = ρ * f² * d⁵ * C_M
    P      = ρ * f³ * d⁵ * C_P

    wr_p = Wrench(F_Op_p, M_Op_p)
    wr_b = t_bp(wr_p)

    sys.y = FixedPitchPropellerY(; v_wOp_p, ω, J, α_p, β_p, coeffs, wr_p, wr_b, P, η_p)

end

f_disc!(::System{<:FixedPitchPropeller}, args...) = false

get_wr_b(sys::System{<:FixedPitchPropeller}) = sys.y.wr_b

end #module