module Propellers

using LinearAlgebra, StaticArrays, StructArrays, Interpolations, UnPack
using Roots: find_zero
using Trapz: trapz
using Plots

import Interpolations: knots, bounds

using Flight.Utils
using Flight.Systems
using Flight.Plotting
using Flight.Kinematics, Flight.Dynamics, Flight.Air

import Flight.Systems: init, f_cont!, f_disc!
import Flight.Plotting: make_plots
import Flight.Dynamics: MassTrait, WrenchTrait, AngularMomentumTrait, get_hr_b, get_wr_b

export FixedPitch, VariablePitch, Propeller

const π² = π^2

############################## Propeller Blade #################################
################################################################################

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
        cL_α(a, α, 0.8) + (cL_α(a, α, 1.2) - cL_α(a, α, 0.8)) / 0.4 * (M - 0.8) #recursive
    end
end

function cD(::DefaultAirfoil, α::Real, M::Real = 0.0)

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

Base.@kwdef struct Blade{C, P, A <: AbstractAirfoil}
    ζ_h::Float64 = 0.2 #hub diameter to blade diameter ratio, ζ ∈ [ζ_h, 1]
    c̃::C = EllipticFunction(0.075) #chord to diameter ratio c_d(ζ)
    p̃::P = ConstantFunction(0.8)#chord-line-pitch to diameter ratio k_c(ζ)
    airfoil::A = DefaultAirfoil()
end

#pitch angle measured relative to the airfoil chord line
get_βc(b::Blade, ζ::Real, Δβ::Real) = atan(b.p̃(ζ) / (π*ζ)) + Δβ

#pitch angle measured relative to the airfoil zero-lift line (aerodynamic)
get_βa(b::Blade, ζ::Real, Δβ::Real) = get_βc(b, ζ, Δβ) - α_0(b.airfoil)


################################################################################
############################ Coefficients ##################################

#note: the number of blades affects the total blade section circulation at a
#given radial distance in the Goldstein's condition, and therefore the more
#blades we have, the larger the induced velocity will be, and the less is to be
#gained from adding further blades. in fact, while the traction coefficient
#increases with the number of blades, the propulsive efficiency decreases

Base.@kwdef struct Coefficients{T}
    C_Fx::T
    C_Mx::T
    C_Fz_α::T
    C_Mz_α::T
    C_P::T
    η_p::T
end

#compute coefficients for n_blades propeller blades with the given properties at a
#specific advance ratio and blade pitch offset setting. CW is assumed.
function Coefficients(blade::Blade, n_blades::Int; J::Real, M_tip::Real, Δβ::Real, n_ζ = 201)

    airfoil = blade.airfoil

    ζ = range(blade.ζ_h, 1, length = n_ζ)
    βa_t = get_βa(blade, 1.0, Δβ) #add the pitch angle offset to the blade's own

    dC_Fx = similar(ζ); dC_Fz_α = similar(ζ)
    dC_Mx = similar(ζ); dC_Mz_α = similar(ζ)

    ε_i = 1 #appears to be a suitable initial guess

    for (i, ζ) in enumerate(ζ)

        ε_inf = atan(J / (π*ζ))
        βa = get_βa(blade, ζ, Δβ)
        c̃ = blade.c̃(ζ)

        βa > π/2 ?  println("Warning: Aerodynamic AoA at $ζ is $(rad2deg(βa))° " *
            "for pitch offset $(rad2deg(Δβ))°") : nothing

        f = let n_blades = n_blades, c̃ = c̃, airfoil = airfoil, βa_t = βa_t, J = J, M_tip = M_tip, βa = βa, ε_inf = ε_inf, ζ = ζ
            ε_i -> induced_angle_eq(; n_blades, c̃, airfoil, βa_t, J, M_tip, βa, ε_inf, ζ, ε_i)
        end

        ε_i = find_zero(f, ε_i) #start at the solution for the previous radial location
        ε = ε_inf + ε_i
        α = βa - ε
        M = M_section(; J, M_tip, ζ, ε_i)


        @assert (α < π/2 && α > -π/3) "α = $α out of bounds after solving for ε_i "*
                "with J = $J, M_tip = $M_tip, Δβ = $Δβ, βa = $βa, ε_i = $ε_i"

        kc̃ = n_blades * c̃
        ζ² = ζ^2; ζ³ = ζ^3
        cos_ε = cos(ε); sin_ε = sin(ε)
        cos²ε_i = cos(ε_i)^2; cos²ε_inf = cos(ε_inf)^2
        tan_ε_inf = tan(ε_inf); tan²ε_inf = tan_ε_inf^2

        let cL = cL(airfoil, α, M), cD = cD(airfoil, α, M), cL_α = cL_α(airfoil, α, M)

            dC_Fx[i] = π²/4 * ζ² * kc̃ * cos²ε_i / cos²ε_inf * (cL * cos_ε - cD * sin_ε)
            dC_Mx[i] = -π²/8 * ζ³ * kc̃ * cos²ε_i / cos²ε_inf * (cD * cos_ε + cL * sin_ε)
            dC_Fz_α[i] = -π²/8 * ζ² * kc̃ * cos²ε_i * (2tan_ε_inf * (cD * cos_ε + cL * sin_ε) - tan²ε_inf * (cL * cos_ε - (cL_α + cD) * sin_ε))
            dC_Mz_α[i] = -π²/16 * ζ³ * kc̃ * cos²ε_i * (2tan_ε_inf * (cL * cos_ε - cD * sin_ε) + tan²ε_inf * ((cL_α + cD) * cos_ε + cL * sin_ε))

        end
    end

    #these are for a CW propeller, we'll deal with CCW signs somewhere else
    C_Fx = trapz(ζ, dC_Fx)
    C_Mx = trapz(ζ, dC_Mx)
    C_Fz_α = trapz(ζ, dC_Fz_α)
    C_Mz_α = trapz(ζ, dC_Mz_α)
    C_P = 2π * C_Mx
    η_p = (C_Fx > 0 ? -J * C_Fx / C_P : 0.0)

    Coefficients(C_Fx, C_Mx, C_Fz_α, C_Mz_α, C_P, η_p)

end

function M_section(; J, M_tip, ζ, ε_i)
    J² = J^2
    ζ² = ζ^2
    return M_tip * √((π² * ζ² + J²) / (π² + J²)) * cos(ε_i)
end

function induced_angle_eq(; n_blades, c̃, airfoil, βa_t, J, M_tip, βa, ε_inf, ζ, ε_i)
    α = βa - ε_inf - ε_i
    M = M_section(; J, M_tip, ζ, ε_i)
    return n_blades*c̃ / (8ζ) * cL(airfoil, α, M) - acos((exp(-n_blades*(1-ζ)/(2sin(βa_t))))) * tan(ε_i) * sin(ε_inf + ε_i)
end

################################################################################
################################ PitchControl ###################################

abstract type PitchControl end

struct FixedPitch <: PitchControl end

Base.@kwdef struct VariablePitch <: PitchControl
    bounds::NTuple{2, Float64} = (deg2rad(0), deg2rad(10))
end

################################################################################
################################ Dataset ###################################

struct Dataset{T <: Interpolations.Extrapolation}
    _data::Coefficients{T}
end

Base.getproperty(dataset::Dataset, s::Symbol) = getproperty(dataset, Val(s))
@generated function Base.getproperty(dataset::Dataset, ::Val{S}) where {S}
    if S === :data
        return :(getfield(dataset, :_data))
    elseif S ∈ fieldnames(Coefficients)
        return :(getfield(getfield(dataset, :_data), $(QuoteNode(S))))
    else
        error("Dataset has no property $S")
    end
end

Interpolations.knots(dataset::Dataset) = Interpolations.knots(dataset.C_Fx)
Interpolations.bounds(dataset::Dataset) = Interpolations.bounds(dataset.C_Fx.itp)

#compute coefficient dataset for a propeller with n_blades of the specified
#geometry and specified pitch control
function Dataset(p::PitchControl, blade::Blade, n_blades::Int; opts...)

    n_ζ = get(opts, :n_ζ, 101)

    J_range = range(0, get(opts, :J_max, 1.5), length = get(opts, :n_J, 21))
    J_mode = BSpline(Linear())
    J_scaling = J_range

    M_tip_range = range(0, get(opts, :M_tip_max, 1.5), length = get(opts, :n_M_tip, 21))
    M_tip_mode = BSpline(Linear())
    M_tip_scaling = M_tip_range

    if p isa VariablePitch
        Δβ_range = range(p.bounds[1], p.bounds[2], length = get(opts, :n_Δβ, 10))
        Δβ_mode = BSpline(Linear())
        Δβ_scaling = Δβ_range
    else
        Δβ_range = range(0, 0, length = 1)
        Δβ_mode = NoInterp()
        Δβ_scaling = 1:1
    end

    data = [Coefficients(blade, n_blades; J, M_tip, Δβ, n_ζ)
            for (J, M_tip, Δβ) in Iterators.product(J_range, M_tip_range, Δβ_range)] |> StructArray |> StructArrays.components

    interps = [extrapolate(scale(interpolate(coef, (J_mode, M_tip_mode, Δβ_mode)), J_scaling, M_tip_scaling, Δβ_scaling), (Flat(), Flat(), Flat()))
            for coef in data]

    Dataset(Coefficients(interps...))

end

function Coefficients(dataset::Dataset, J::Real, M_tip::Real, Δβ::Real)

    @unpack C_Fx, C_Mx, C_Fz_α, C_Mz_α, C_P, η_p = dataset

    Coefficients(  C_Fx = C_Fx(J, M_tip, Δβ), C_Mx = C_Mx(J, M_tip, Δβ),
                        C_Fz_α = C_Fz_α(J, M_tip, Δβ), C_Mz_α = C_Mz_α(J, M_tip, Δβ),
                        C_P = C_P(J, M_tip, Δβ), η_p = η_p(J, M_tip, Δβ))
end

(d::Dataset)(J::Real, M_tip::Real, Δβ::Real) = Coefficients(d, J, M_tip, Δβ)

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

struct Propeller{P <: PitchControl, B <: Blade,  D <: Dataset} <: AbstractPropeller
    pitch::P
    blade::B
    dataset::D
    n_blades::Int
    sense::TurnSense
    d::Float64 #diameter
    J_xx::Float64 #axial moment of inertia, J_xx label avoids confusion with advance ratio
    t_bp::FrameTransform
end

function Propeller(; pitch = FixedPitch(), blade = Blade(), n_blades = 2,
                   d = 2.0, J_xx = 0.3, sense = CW, t_bp = FrameTransform(), dataset_opts = ())

    dataset = Dataset(pitch, blade, n_blades; dataset_opts...)

    Propeller{typeof(pitch), typeof(blade), typeof(dataset)}(
        pitch, blade, dataset, n_blades, sense, d, J_xx, t_bp)
end


Base.@kwdef struct PropellerY
    v_wOp_p::SVector{3,Float64} = zeros(SVector{3}) #local aerodynamic velocity, propeller axes
    ω::Float64 = 0 #angular velocity
    J::Float64 = 0 #advance ratio
    M_tip::Float64 = 0 #blade tip Mach number
    Δβ::Float64 = 0 #blade pitch offset
    wr_p::Wrench = Wrench() #resulting aerodynamic Wrench, propeller frame
    wr_b::Wrench = Wrench() #resulting aerodynamic Wrench, vehicle frame
    hr_p::SVector{3,Float64} = zeros(SVector{3})#angular momentum, propeller frame
    hr_b::SVector{3,Float64} = zeros(SVector{3}) #angular momentum, vehicle frame
    P::Float64 = 0.0 #power produced by the propeller
    η_p::Float64 = 0.0 #propulsive efficiency
end

init(::Propeller{FixedPitch}, ::SystemU) = nothing
init(::Propeller{VariablePitch}, ::SystemU) = Ref(Ranged(0.0, 0, 1))
init(::Propeller, ::SystemY) = PropellerY()

get_Δβ(sys::System{<:Propeller{FixedPitch}}) = 0.0
get_Δβ(sys::System{<:Propeller{VariablePitch}}) = linear_scaling(sys.u[], sys.params.pitch.bounds)

function f_cont!(sys::System{<:Propeller}, kin::KinData, air::AirflowData, ω::Real)

    @unpack d, J_xx, t_bp, sense, dataset = sys.params
    #remove this, it may happen due to friction overshoot at low RPMs
    # @assert sign(ω) * Int(sys.params.sense) >= 0 "Propeller turning in the wrong sense"

    v_wOp_b = air.v_wOb_b + kin.vel.ω_eb_b × t_bp.r
    v_wOp_p = t_bp.q'(v_wOp_b)

    #compute advance ratio. here we use the velocity vector magnitude rather
    #than its axial component. J must be positive, so we need the abs for CCW
    #propellers. also, prevent division by zero in a non-rotating propeller
    abs_ω_min = 1.0
    v_J = norm(v_wOp_p)
    ω_J = max(abs(ω), abs_ω_min)
    J = 2π * v_J / (ω_J * d)

    M_tip = abs(ω)*(d/2) / air.a

    Δβ = get_Δβ(sys)
    coeffs = Coefficients(dataset, J, M_tip, Δβ)

    @unpack C_Fx, C_Mx, C_Fz_α, C_Mz_α, C_P, η_p = coeffs
    C_Fy_β = C_Fz_α #by y/z symmetry
    C_My_β = C_Mz_α #by y/z symmetry

    α_p, β_p = get_airflow_angles(v_wOp_p)

    C_F = SVector{3,Float64}(C_Fx, C_Fy_β * β_p, C_Fz_α * α_p)
    C_M = Int(sense) * SVector{3,Float64}(C_Mx, C_My_β * β_p, C_Mz_α * α_p)

    ρ = air.ρ
    f = ω/2π; f² = f^2; f³ = f * f²
    d⁴ = d^4; d⁵ = d * d⁴

    F_Op_p = ρ * f² * d⁴ * C_F
    M_Op_p = ρ * f² * d⁵ * C_M
    P      = ρ * abs(f³) * d⁵ * C_P

    wr_p = Wrench(F_Op_p, M_Op_p)
    wr_b = t_bp(wr_p)

    hr_p = SVector(J_xx * ω, 0, 0)
    hr_b = t_bp.q(hr_p)

    sys.y = PropellerY(; v_wOp_p, ω, J, M_tip, Δβ, wr_p, wr_b, hr_p, hr_b, P, η_p)

end

f_disc!(::System{<:Propeller}, args...) = false

get_wr_b(sys::System{<:Propeller}) = sys.y.wr_b
get_hr_b(sys::System{<:Propeller}) = sys.y.hr_b


################################################################################
####################### ConstantSpeedPropeller #################################

# abstract type AbstractGovernor <: SystemDescriptor end

# struct ConstantSpeedPropeller{P <: Propeller{VariablePitch}, G <: AbstractGovernor} <: AbstractPropeller
#     propeller::P
#     governor::G
# end

# struct DummyGovernor end
# @inline f_cont!(::System{DummyGovernor}, args...) = nothing
# @inline (f_disc!(::System{NullSystemDescriptor}, args...)::Bool) = false

################################################################################
############################ Plots #############################################

#called as: plot_airfoil(airfoil; Plotting.defaults...)
function plot_airfoil(airfoil::Propellers.AbstractAirfoil; plot_settings...)

    α = range(-π/6, π/3, length = 100)
    M = range(0, 1.5, length = 6)
    iter = Iterators.product(α, M)

    cL_data = [cL(airfoil, α, M) for (α, M) in iter]
    cD_data = [cD(airfoil, α, M) for (α, M) in iter]
    cL_α_data = [cL_α(airfoil, α, M) for (α, M) in iter]

    label = latexstring.("M = ".*string.(M'))
    titles = (L"c_L", L"c_{L, \alpha}", L"c_D")
    x_label = L"\alpha \ (rad)"

    p = Vector{Plots.Plot}()
    for (i, data) in enumerate((cL_data, cD_data, cL_α_data))
        push!(p, plot(α, data; label = label, title = titles[i], xlabel = x_label, plot_settings...))
    end

    return p

end

function plot_J_Δβ(dataset::Propellers.Dataset, M_tip::Real = 0.0; plot_settings...)

    J = knots(dataset).iterators[1] |> collect
    Δβ_bounds = bounds(dataset)[3]
    Δβ = range(Δβ_bounds[1], Δβ_bounds[2], length = 5)

    data = [dataset(J, M_tip, Δβ) for (J, Δβ) in Iterators.product(J, Δβ)]
    data = data |> StructArray |> StructArrays.components

    @unpack C_Fx, C_Mx, C_Fz_α, C_Mz_α, C_P, η_p = data

    label = latexstring.("\$ \\Delta \\beta = " .* string.(rad2deg.(Δβ')) .* "\\degree \$")
    label_pos = [:bottomleft, :topleft, :bottomleft, :bottomleft, :topright, :topleft]
    x_label = L"J"
    y_labels = [L"C_{Fx}", L"C_{Mx}", L"C_{Fz, \alpha}", L"C_{Mz, \alpha}", L"M_{tip}", L"\eta_p"]
    titles = ["Traction Coefficient ", "Torque Coefficient", "Off-Axis Force Coefficient Derivative",
                "Off-Axis Moment Coefficient Derivative", "Power Coefficient", "Propulsive Efficiency"] .*
                " (Blade Tip Mach Number = $M_tip)"

    p = Vector{Plots.Plot}()
    for (i, c) in enumerate((C_Fx, C_Mx, C_Fz_α, C_Mz_α, C_P, η_p))
        push!(p, plot( J, c;
                        title = titles[i], label = label, legend = label_pos[i],
                        xlabel = x_label, ylabel = y_labels[i], plot_settings...))
    end

    return p

end

function plot_M_J(dataset::Propellers.Dataset, Δβ::Real = 0.0; plot_settings...)

    M_tip = knots(dataset).iterators[2] |> collect
    J_bounds = bounds(dataset)[1]
    J = range(J_bounds[1], J_bounds[2], length = 5)

    data = [dataset(J, M_tip, Δβ) for (M_tip, J) in Iterators.product(M_tip, J)]
    data = data |> StructArray |> StructArrays.components

    @unpack C_Fx, C_Mx, C_Fz_α, C_Mz_α, C_P, η_p = data

    label = latexstring.("J = ".*string.(J'))
    label_pos = [:bottomleft, :topleft, :bottomleft, :bottomleft, :topright, :topleft]
    x_label = L"M_{tip}"
    y_labels = [L"C_{Fx}", L"C_{Mx}", L"C_{Fz, \alpha}", L"C_{Mz, \alpha}", L"M_{tip}", L"\eta_p"]
    titles = ["Traction Coefficient ", "Torque Coefficient", "Off-Axis Force Coefficient Derivative",
                "Off-Axis Moment Coefficient Derivative", "Power Coefficient", "Propulsive Efficiency"] .*
                " (Blade Pitch Offset = $(rad2deg(Δβ))°)"

    p = Vector{Plots.Plot}()
    for (i, c) in enumerate((C_Fx, C_Mx, C_Fz_α, C_Mz_α, C_P, η_p))
        push!(p, plot(M_tip, c;
                        title = titles[i], label = label, legend = label_pos[i],
                        xlabel = x_label, ylabel = y_labels[i],
                        plot_settings...))
    end

    return p

end

function plot_J_M(dataset::Propellers.Dataset, Δβ::Real = 0.0; plot_settings...)

    J = knots(dataset).iterators[1] |> collect
    M_tip_bounds = bounds(dataset)[2]
    M_tip = range(M_tip_bounds[1], M_tip_bounds[2], length = 5)

    data = [dataset(J, M_tip, Δβ) for (J, M_tip) in Iterators.product(J, M_tip)]
    data = data |> StructArray |> StructArrays.components

    @unpack C_Fx, C_Mx, C_Fz_α, C_Mz_α, C_P, η_p = data

    label = latexstring.("M_{tip} = ".*string.(M_tip'))
    label_pos = [:bottomleft, :topleft, :bottomleft, :bottomleft, :topright, :topleft]
    x_label = L"J"
    y_labels = [L"C_{Fx}", L"C_{Mx}", L"C_{Fz, \alpha}", L"C_{Mz, \alpha}", L"M_{tip}", L"\eta_p"]
    titles = ["Traction Coefficient ", "Torque Coefficient", "Off-Axis Force Coefficient Derivative",
                "Off-Axis Moment Coefficient Derivative", "Power Coefficient", "Propulsive Efficiency"] .*
                " (Blade Pitch Offset = $(rad2deg(Δβ))°)"

    p = Vector{Plots.Plot}()
    for (i, c) in enumerate((C_Fx, C_Mx, C_Fz_α, C_Mz_α, C_P, η_p))
        push!(p, plot(J, c;
                        title = titles[i], label = label, legend = label_pos[i],
                        xlabel = x_label, ylabel = y_labels[i],
                        plot_settings...))
    end

    return p

end

end #module