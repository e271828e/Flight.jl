module Propellers

using LinearAlgebra, StaticArrays, StructArrays, Interpolations, UnPack
using Logging
using HDF5
using Roots: find_zero
using Trapz: trapz

using Flight.FlightCore
using Flight.FlightCore.Utils

using Flight.FlightPhysics

export FixedPitch, VariablePitch, Propeller

const π² = π^2

radpersec2RPM(ω) = ω/(π/30)

############################## Propeller Blade #################################
################################################################################

abstract type AbstractDistribution{N} end

get_value(f::AbstractDistribution{N}, ::Vararg{Any, M}) where {N, M} =
    error("A $(typeof(f)) was called with $M arguments, requires $N")

get_value(f::AbstractDistribution{N}, ::Vararg{Any, N}) where {N} =
    error("Method get_value not implemented for $(typeof(f))")

(f::AbstractDistribution)(args...) = get_value(f, args...)

struct ConstantDistribution <: AbstractDistribution{1}
    a::Float64
end
get_value(f::ConstantDistribution, ::Real) = f.a

@kwdef struct EllipticDistribution <: AbstractDistribution{1}
    a::Float64 = 0.075
end
get_value(f::EllipticDistribution, ζ::Real) = f.a*√(1 - ζ^2)


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

@kwdef struct Blade{C, P, A <: AbstractAirfoil}
    ζ_h::Float64 = 0.2 #hub diameter to blade diameter ratio, ζ ∈ [ζ_h, 1]
    c̃::C = EllipticDistribution(0.075) #chord to diameter ratio c_d(ζ)
    p̃::P = ConstantDistribution(0.8)#chord-line-pitch to diameter ratio k_c(ζ)
    airfoil::A = DefaultAirfoil()
end

#pitch angle measured relative to the airfoil chord line
get_βc(b::Blade, ζ::Real, Δβ::Real) = atan(b.p̃(ζ) / (π*ζ)) + Δβ

#pitch angle measured relative to the airfoil zero-lift line (aerodynamic)
get_βa(b::Blade, ζ::Real, Δβ::Real) = get_βc(b, ζ, Δβ) - α_0(b.airfoil)


################################################################################
############################ Coefficients ##################################

@kwdef struct Coefficients{T}
    C_Fx::T
    C_Mx::T
    C_Fz_α::T
    C_Mz_α::T
    C_P::T
    η_p::T
end

function Base.NamedTuple(coefs::Coefficients)
    names = fieldnames(typeof(coefs))
    fields = map(n -> getfield(coefs, n), names)
    NamedTuple{names}(fields)
end

#compute propeller coefficients for n blades with the specified properties, at a
#given advance ratio, blade tip Mach number and blade pitch offset setting. CW
#is assumed.
function Coefficients(blade::Blade, n_blades::Int, J::Real, Mt::Real, Δβ::Real, n_ζ::Integer = 101)

    @assert J >= 0 "Advance ratio must be 0 or positive"
    @assert Mt >= 0 "Blade tip Mach number must be 0 or positive"
    @assert n_ζ > 50 "At least 51 radial discretization points must be used"

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

        (βa > π/2) &&  @warn("Aerodynamic AoA at $ζ is $(rad2deg(βa))° " *
            "for pitch offset $(rad2deg(Δβ))°")

        f = let n_blades = n_blades, c̃ = c̃, airfoil = airfoil, βa_t = βa_t, J = J, Mt = Mt, βa = βa, ε_inf = ε_inf, ζ = ζ
            ε_i -> induced_angle_eq(; n_blades, c̃, airfoil, βa_t, J, Mt, βa, ε_inf, ζ, ε_i)
        end

        ε_i = find_zero(f, ε_i) #start at the solution for the previous radial location
        ε = ε_inf + ε_i
        α = βa - ε
        M = M_section(; J, Mt, ζ, ε_i)


        @assert (α < π/2 && α > -π/3) "α = $α out of bounds after solving for ε_i "*
                "with J = $J, Mt = $Mt, Δβ = $Δβ, βa = $βa, ε_i = $ε_i"

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

    #these are for a CW propeller, we'll deal with the CCW case through symmetry
    #inside the propeller update function
    C_Fx = trapz(ζ, dC_Fx)
    C_Mx = trapz(ζ, dC_Mx)
    C_Fz_α = trapz(ζ, dC_Fz_α)
    C_Mz_α = trapz(ζ, dC_Mz_α)
    C_P = 2π * C_Mx
    η_p = (C_Fx > 0 ? -J * C_Fx / C_P : 0.0)

    Coefficients(C_Fx, C_Mx, C_Fz_α, C_Mz_α, C_P, η_p)

end

function M_section(; J, Mt, ζ, ε_i)
    J² = J^2
    ζ² = ζ^2
    return Mt * √((π² * ζ² + J²) / (π² + J²)) * cos(ε_i)
end

function induced_angle_eq(; n_blades, c̃, airfoil, βa_t, J, Mt, βa, ε_inf, ζ, ε_i)
    α = βa - ε_inf - ε_i
    M = M_section(; J, Mt, ζ, ε_i)
    return n_blades*c̃ / (8ζ) * cL(airfoil, α, M) - acos((exp(-n_blades*(1-ζ)/(2sin(βa_t))))) * tan(ε_i) * sin(ε_inf + ε_i)
end



################################################################################
################################## Lookup ######################################

struct Lookup{T <: Interpolations.Extrapolation}
    interps::Coefficients{T}
    data::Coefficients{Array{Float64,3}}
    J_bounds::NTuple{2, Float64}
    Mt_bounds::NTuple{2, Float64}
    Δβ_bounds::NTuple{2, Float64}
end

Base.getproperty(lookup::Lookup, s::Symbol) = getproperty(lookup, Val(s))
@generated function Base.getproperty(lookup::Lookup, ::Val{S}) where {S}
    if S ∈ fieldnames(Lookup)
        return :(getfield(lookup, $(QuoteNode(S))))
    elseif S ∈ fieldnames(Coefficients)
        return :(getfield(getfield(lookup, :interps), $(QuoteNode(S))))
    else
        error("Lookup has no property $S")
    end
end

function Lookup(blade::Blade = Blade(), n_blades::Int = 2;
                J_range::AbstractRange{Float64} = range(0, 1.5, length = 21),
                Mt_range::AbstractRange{Float64} = range(0, 1.5, length = 21),
                Δβ_range::AbstractRange{Float64} = range(0, 0, length = 1),
                n_ζ::Integer = 101)

    data_points = [Coefficients(blade, n_blades, J, Mt, Δβ, n_ζ)
            for (J, Mt, Δβ) in Iterators.product(J_range, Mt_range, Δβ_range)]

    J_bounds = (J_range[1], J_range[end])
    Mt_bounds = (Mt_range[1], Mt_range[end])
    Δβ_bounds = (Δβ_range[1], Δβ_range[end])

    data = Coefficients(StructArrays.components(StructArray(data_points))...)

    Lookup(data, J_bounds, Mt_bounds, Δβ_bounds)

end

function Lookup(data::Coefficients{Array{Float64, 3}},
                J_bounds::NTuple{2, Float64},
                Mt_bounds::NTuple{2, Float64},
                Δβ_bounds::NTuple{2, Float64})

    J_length, Mt_length, Δβ_length = size(data.C_Fx)
    J_range = range(J_bounds..., length = J_length)
    Mt_range = range(Mt_bounds..., length = Mt_length)
    Δβ_range = range(Δβ_bounds..., length = Δβ_length)

    (J, Mt, Δβ) =  map((J_range, Mt_range, Δβ_range)) do range
        (mode, scaling) = length(range) > 1 ? (BSpline(Linear()), range) : (NoInterp(), 1:1)
        return (mode = mode, scaling = scaling)
    end

    interps = [extrapolate(scale(interpolate(coef,
                                            (J.mode, Mt.mode, Δβ.mode)),
                                 J.scaling, Mt.scaling, Δβ.scaling),
                           (Flat(), Flat(), Flat())) for coef in NamedTuple(data)]

    Lookup(Coefficients(interps...), data, J_bounds, Mt_bounds, Δβ_bounds)

end

(d::Lookup)(J::Real, Mt::Real, Δβ::Real) = Coefficients(d, J, Mt, Δβ)

function Coefficients(lookup::Lookup, J::Real, Mt::Real, Δβ::Real)

    @unpack C_Fx, C_Mx, C_Fz_α, C_Mz_α, C_P, η_p = lookup

    Coefficients(   C_Fx = C_Fx(J, Mt, Δβ),
                    C_Mx = C_Mx(J, Mt, Δβ),
                    C_Fz_α = C_Fz_α(J, Mt, Δβ),
                    C_Mz_α = C_Mz_α(J, Mt, Δβ),
                    C_P = C_P(J, Mt, Δβ),
                    η_p = η_p(J, Mt, Δβ))
end

function save_lookup(lookup::Lookup, fname::String)

    fid = h5open(fname, "w")

    foreach(pairs(NamedTuple(lookup.data))) do (name, data)
        fid[string(name)] = data
    end

    fid["J_start"], fid["J_end"] = lookup.J_bounds
    fid["Mt_start"], fid["Mt_end"] = lookup.Mt_bounds
    fid["Δβ_start"], fid["Δβ_end"] = lookup.Δβ_bounds

    close(fid)
end

function load_lookup(fname::String)

    fid = h5open(fname, "r")

    coef_data = map(fieldnames(Coefficients)) do name
        read(fid, string(name))
    end
    data = Coefficients(coef_data...)

    J_bounds = (read(fid["J_start"]), read(fid["J_end"]))
    Mt_bounds = (read(fid["Mt_start"]), read(fid["Mt_end"]))
    Δβ_bounds = (read(fid["Δβ_start"]), read(fid["Δβ_end"]))

    close(fid)

    return Lookup(data, J_bounds, Mt_bounds, Δβ_bounds)

end

################################################################################
################################ Propeller #####################################

@enum TurnSense begin
    CW = 1
    CCW = -1
end

abstract type PitchStyle end
struct FixedPitch <: PitchStyle end
struct VariablePitch <: PitchStyle end

abstract type AbstractPropeller <: SystemDefinition end

struct Propeller{P <: PitchStyle, L <: Lookup} <: AbstractPropeller
    pitch::P
    lookup::L
    sense::TurnSense
    d::Float64 #diameter
    J_xx::Float64 #axial moment of inertia, J_xx label avoids confusion with advance ratio
    t_bp::FrameTransform #vehicle frame to propeller frame
end

function Propeller(lookup::Lookup = Lookup();
                   sense ::TurnSense = CW, d::Real = 2.0, J_xx::Real = 0.3,
                   t_bp::FrameTransform = FrameTransform())

    Δβ_bounds = lookup.Δβ_bounds
    pitch =  Δβ_bounds[1] == Δβ_bounds[end] ? FixedPitch() : VariablePitch()
    Propeller{typeof(pitch), typeof(lookup)}(pitch, lookup, sense, d, J_xx, t_bp)
end

@kwdef struct PropellerY
    v_wOp_p::SVector{3,Float64} = zeros(SVector{3}) #local aerodynamic velocity, propeller axes
    ω::Float64 = 0 #angular velocity
    J::Float64 = 0 #advance ratio
    Mt::Float64 = 0 #blade tip Mach number
    Δβ::Float64 = 0 #blade pitch offset
    wr_p::Wrench = Wrench() #resulting aerodynamic Wrench, propeller frame
    wr_b::Wrench = Wrench() #resulting aerodynamic Wrench, vehicle frame
    hr_p::SVector{3,Float64} = zeros(SVector{3})#angular momentum, propeller frame
    hr_b::SVector{3,Float64} = zeros(SVector{3}) #angular momentum, vehicle frame
    P::Float64 = 0.0 #power produced by the propeller
    η_p::Float64 = 0.0 #propulsive efficiency
end

Systems.U(::Propeller{FixedPitch}) = nothing
Systems.U(::Propeller{VariablePitch}) = Ref(Ranged(0.0, 0., 1.))
Systems.Y(::Propeller) = PropellerY()

function get_Δβ(sys::System{<:Propeller{FixedPitch}})
    Δβ_bounds = sys.constants.lookup.Δβ_bounds
    @assert Δβ_bounds[1] == Δβ_bounds[2]
    return Δβ_bounds[1]
end

function get_Δβ(sys::System{<:Propeller{VariablePitch}})
    Δβ_bounds = sys.constants.lookup.Δβ_bounds
    @assert Δβ_bounds[2] > Δβ_bounds[1]
    return linear_scaling(sys.u[], Δβ_bounds)
end

function Systems.f_ode!(sys::System{<:Propeller}, kin::KinData, air::AirData, ω::Real)

    @unpack d, J_xx, t_bp, sense, lookup = sys.constants
    #remove this, it may happen due to friction overshoot at low RPMs
    # @assert sign(ω) * Int(sys.constants.sense) >= 0 "Propeller turning in the wrong sense"

    v_wOp_b = air.v_wOb_b + kin.ω_eb_b × t_bp.r
    v_wOp_p = t_bp.q'(v_wOp_b)

    #compute advance ratio. here we use the velocity vector magnitude rather
    #than its axial component. J must be positive, so we need the abs for CCW
    #propellers
    abs_ω_min = 1.0 #prevent division by zero in a non-rotating propeller
    v_J = norm(v_wOp_p)
    ω_J = max(abs(ω), abs_ω_min)
    J = 2π * v_J / (ω_J * d)

    Mt = abs(ω)*(d/2) / air.a

    Δβ = get_Δβ(sys)
    coeffs = Coefficients(lookup, J, Mt, Δβ)

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

    sys.y = PropellerY(; v_wOp_p, ω, J, Mt, Δβ, wr_p, wr_b, hr_p, hr_b, P, η_p)

end

Dynamics.get_mp_b(::System{<:Propeller}) = MassProperties()
Dynamics.get_wr_b(sys::System{<:Propeller}) = sys.y.wr_b
Dynamics.get_hr_b(sys::System{<:Propeller}) = sys.y.hr_b


################################################################################
####################### ConstantSpeedPropeller #################################

# abstract type AbstractGovernor <: SystemDefinition end

# struct ConstantSpeedPropeller{P <: Propeller{VariablePitch}, G <: AbstractGovernor} <: AbstractPropeller
#     propeller::P
#     governor::G
# end

# struct DummyGovernor end

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

function plot_J_Δβ(lookup::Propellers.Lookup, Mt::Real = 0.0; plot_settings...)

    @unpack J_bounds, Mt_bounds, Δβ_bounds = lookup._data

    @assert Mt_bounds[1] <= Mt <= Mt_bounds[2]

    J = range(J_bounds..., length = 50)
    Δβ = range(Δβ_bounds..., length = 5)

    data = [lookup(J, Mt, Δβ) for (J, Δβ) in Iterators.product(J, Δβ)]
    data = data |> StructArray |> StructArrays.components

    @unpack C_Fx, C_Mx, C_Fz_α, C_Mz_α, C_P, η_p = data

    label = latexstring.("\$ \\Delta \\beta = " .* string.(rad2deg.(Δβ')) .* "\\degree \$")
    label_pos = [:bottomleft, :topleft, :bottomleft, :bottomleft, :topright, :topleft]
    x_label = L"J"
    y_labels = [L"C_{Fx}", L"C_{Mx}", L"C_{Fz, \alpha}", L"C_{Mz, \alpha}", L"M_{tip}", L"\eta_p"]
    titles = ["Traction Coefficient ", "Torque Coefficient", "Off-Axis Force Coefficient Derivative",
                "Off-Axis Moment Coefficient Derivative", "Power Coefficient", "Propulsive Efficiency"] .*
                " (Blade Tip Mach Number = $Mt)"

    p = Vector{Plots.Plot}()
    for (i, c) in enumerate((C_Fx, C_Mx, C_Fz_α, C_Mz_α, C_P, η_p))
        push!(p, plot( J, c;
                        title = titles[i], label = label, legend = label_pos[i],
                        xlabel = x_label, ylabel = y_labels[i], plot_settings...))
    end

    return p

end

function plot_M_J(lookup::Propellers.Lookup, Δβ::Real = 0.0; plot_settings...)

    @unpack J_bounds, Mt_bounds, Δβ_bounds = lookup._data

    @assert Δβ_bounds[1] <= Δβ <= Δβ_bounds[2]

    Mt = range(Mt_bounds..., length = 50)
    J = range(lookup._data.J_bounds..., length = 5)

    data = [lookup(J, Mt, Δβ) for (Mt, J) in Iterators.product(Mt, J)]
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
        push!(p, plot(Mt, c;
                        title = titles[i], label = label, legend = label_pos[i],
                        xlabel = x_label, ylabel = y_labels[i],
                        plot_settings...))
    end

    return p

end

function plot_J_M(lookup::Propellers.Lookup, Δβ::Real = 0.0; plot_settings...)

    @unpack J_bounds, Mt_bounds, Δβ_bounds = lookup._data

    @assert Δβ_bounds[1] <= Δβ <= Δβ_bounds[2]

    J = range(J_bounds..., length = 100)
    Mt = range(Mt_bounds..., length = 5)

    data = [lookup(J, Mt, Δβ) for (J, Mt) in Iterators.product(J, Mt)]
    data = data |> StructArray |> StructArrays.components

    @unpack C_Fx, C_Mx, C_Fz_α, C_Mz_α, C_P, η_p = data

    label = latexstring.("M_{tip} = ".*string.(Mt'))
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

function GUI.draw(sys::System{<:Propeller}, p_open::Ref{Bool} = Ref(true),
                    window_label::String = "Propeller")

    @unpack ω, J, Mt, Δβ, wr_p, hr_p, P, η_p = sys.y

    CImGui.Begin(window_label, p_open)

        CImGui.Text(@sprintf("Angular Rate (Propeller/Body): %.7f RPM", radpersec2RPM(ω)))
        CImGui.Text(@sprintf("Advance Ratio: %.7f", J))
        CImGui.Text(@sprintf("Blade Tip Mach Number: %.7f", Mt))
        CImGui.Text(@sprintf("Blade Pitch Offset: %.7f deg", rad2deg(Δβ)))
        CImGui.Text(@sprintf("Axial Angular Momentum: %.7f kg*(m^2)/s", hr_p[1]))
        CImGui.Text(@sprintf("Power: %.7f kW", 1e-3*P))
        CImGui.Text(@sprintf("Propulsive Efficiency: %.7f", η_p))
        GUI.draw(wr_p.F, "Aerodynamic Force (Op) [Propeller]", "N")
        GUI.draw(wr_p.M, "Aerodynamic Torque (Op) [Propeller]", "N*m")
        GUI.draw(hr_p, "Axial Angular Momentum [Propeller]", "N*m")

    CImGui.End()

end


end #module