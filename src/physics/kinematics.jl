module Kinematics

using StaticArrays, StructArrays, ComponentArrays, LinearAlgebra, UnPack

using Flight.FlightCore

using ..Attitude
using ..Geodesy

export AbstractKinematicDescriptor, ECEF, LTF, NED
export KinInit, KinData, KinSystem

const v_min_χγ = 0.1 #minimum speed for valid χ, γ


############################## Initialization ##################################
################################################################################

#user-friendly conditions for kinematic state initialization
struct Initializer
    q_nb::RQuat
    Ob::Geographic{NVector, Ellipsoidal}
    ω_lb_b::SVector{3, Float64}
    v_eOb_n::SVector{3, Float64}
    Δx::Float64
    Δy::Float64
end

const KinInit = Initializer

function Initializer(;
    q_nb::Abstract3DRotation = RQuat(), loc::Abstract2DLocation = LatLon(),
    h::Altitude = HOrth(), ω_lb_b::AbstractVector{<:Real} = zeros(SVector{3}),
    v_eOb_n::AbstractVector{<:Real} = zeros(SVector{3}), Δx::Real = 0.0, Δy::Real = 0.0)

    Ob = Geographic(loc, h)
    Initializer(q_nb, Ob, ω_lb_b, v_eOb_n, Δx, Δy)
end


######################### AbstractKinematicDescriptor ##########################
################################################################################

abstract type AbstractKinematicDescriptor <: SystemDefinition end
const KinSystem = System{<:AbstractKinematicDescriptor}

const XVelTemplate = ComponentVector(ω_eb_b = zeros(3), v_eOb_b = zeros(3))
Systems.U(::AbstractKinematicDescriptor) = zero(XVelTemplate)


################################# KinematicsY ##################################
################################################################################

@kwdef struct KinData
    e_nb::REuler = REuler()
    q_nb::RQuat = RQuat()
    q_eb::RQuat = RQuat()
    q_en::RQuat = RQuat()
    ϕ_λ::LatLon = LatLon()
    n_e::NVector = NVector()
    h_e::Altitude{Ellipsoidal} = HEllip(0.0)
    h_o::Altitude{Orthometric} = HOrth(0.0)
    Δxy::SVector{2,Float64} = zeros(SVector{2})
    r_eOb_e::SVector{3,Float64} = zeros(SVector{3})
    ω_lb_b::SVector{3,Float64} = zeros(SVector{3})
    ω_eb_b::SVector{3,Float64} = zeros(SVector{3})
    ω_ie_b::SVector{3,Float64} = zeros(SVector{3})
    ω_ib_b::SVector{3,Float64} = zeros(SVector{3})
    v_eOb_b::SVector{3,Float64} = zeros(SVector{3})
    v_eOb_n::SVector{3,Float64} = zeros(SVector{3})
    v_gnd::Float64 = 0.0 #ground speed
    χ_gnd::Float64 = 0.0 #course angle
    γ_gnd::Float64 = 0.0 #flight path angle
end

function KinData(ic::KinInit)

    @unpack q_nb, Ob, ω_lb_b, v_eOb_n, Δx, Δy = ic

    e_nb = REuler(q_nb)
    q_en = ltf(Ob)
    q_eb = q_en ∘ q_nb

    n_e = NVector(Ob)
    ϕ_λ = LatLon(Ob)
    h_e = HEllip(Ob)
    h_o = HOrth(h_e, n_e)
    Δxy = SVector(Δx, Δy)
    r_eOb_e = Cartesian(Ob)

    ω_el_n = get_ω_el_n(v_eOb_n, Ob)
    ω_el_b = q_nb'(ω_el_n)
    ω_eb_b = ω_el_b + ω_lb_b

    ω_ie_e = SVector(0., 0., ω_ie)
    ω_ie_b = q_eb'(ω_ie_e)
    ω_ib_b = ω_ie_b + ω_eb_b

    v_eOb_b = q_nb'(v_eOb_n)

    v_gnd = norm(v_eOb_n)
    χ_gnd = v_gnd > v_min_χγ ? azimuth(v_eOb_n) : 0.0
    γ_gnd = v_gnd > v_min_χγ ? inclination(v_eOb_n) : 0.0

    KinData(;   e_nb, q_nb, q_eb, q_en, ϕ_λ, n_e, h_e, h_o, Δxy, r_eOb_e,
                ω_lb_b, ω_eb_b, ω_ie_b, ω_ib_b,
                v_eOb_b, v_eOb_n, v_gnd, χ_gnd, γ_gnd)

end


struct KinematicsY{S}
    data::KinData #general, implementation-agnostic kinematic data
    impl::S #implementation-specific outputs
end

KinData(y::KinematicsY) = y.data
KinData(sys::KinSystem) = KinData(sys.y)

Base.getproperty(y::KinematicsY, s::Symbol) = getproperty(y, Val(s))

@generated function Base.getproperty(y::KinematicsY{T}, ::Val{S}) where {T, S}
    if S === :data || S === :impl
        return :(getfield(y, $(QuoteNode(S))))
    elseif S ∈ fieldnames(KinData)
        return :(getfield(getfield(y, :data), $(QuoteNode(S))))
    elseif S ∈ fieldnames(T)
        return :(getfield(getfield(y, :impl), $(QuoteNode(S))))
    else
        error("$(typeof(y)) has no property $S")
    end
end

function normalize_block!(x, ε) #returns true if norm was corrected
    norm_x = norm(x)
    (abs(norm_x - 1.0) > ε) && (x ./= norm_x)
    return nothing
end

function Geodesy.gravity(kin_data::KinData)
    gravity(Geographic(kin_data.n_e, kin_data.h_e))
end

function Geodesy.G_n(kin_data::KinData)
    G_n(Geographic(kin_data.n_e, kin_data.h_e))
end


########################### LTF-based Kinematics #########################
##########################################################################

#fast, singularity-free (all-attitude, all-latitude) kinematic mechanization,
#appropriate for simulation. local tangent frame is the wander-azimuth frame

struct LTF <: AbstractKinematicDescriptor end

@kwdef struct LTFData
    q_lb::RQuat = RQuat()
    q_el::RQuat = RQuat()
    ω_el_l::SVector{3,Float64} = zeros(SVector{3})
end

Systems.X(::LTF) = ComponentVector(
    q_lb = zeros(4), q_el = zeros(4), Δx = 0.0, Δy = 0.0, h_e = 0.0)

Systems.Y(::LTF) = KinematicsY(KinData(), LTFData())

function Systems.init!(sys::System{LTF}, ic::Initializer = Initializer())

    @unpack x, u = sys
    @unpack q_nb, Ob, ω_lb_b, v_eOb_n, Δx, Δy = ic

    ω_el_n = get_ω_el_n(v_eOb_n, Ob)
    ω_el_b = q_nb'(ω_el_n)
    ω_eb_b = ω_el_b + ω_lb_b
    v_eOb_b = q_nb'(v_eOb_n)
    h_e = HEllip(Ob)

    q_lb = q_nb #arbitrarily initializes wander angle ψ_nl to 1

    u.ω_eb_b = ω_eb_b
    u.v_eOb_b = v_eOb_b

    x.q_lb = q_lb[:]
    x.q_el = ltf(Ob)[:]
    x.Δx = Δx
    x.Δy = Δy
    x.h_e = h_e

    f_ode!(sys) #update ẋ and y, not strictly necessary

end


function Systems.f_ode!(sys::System{LTF})

    @unpack ẋ, x, u = sys

    q_lb = RQuat(x.q_lb, normalization = false)
    q_el = RQuat(x.q_el, normalization = false)
    ω_eb_b = SVector{3}(u.ω_eb_b)
    v_eOb_b = SVector{3}(u.v_eOb_b)
    h_e = HEllip(x.h_e[1])
    Δxy = SVector(x.Δx, x.Δy)

    ψ_nl = get_ψ_nl(q_el)
    q_nl = Rz(ψ_nl)
    q_nb = q_nl ∘ q_lb
    q_eb = q_el ∘ q_lb
    q_en = q_eb ∘ q_nb'
    e_nb = REuler(q_nb)

    n_e = NVector(q_el)
    ϕ_λ = LatLon(n_e)
    h_o = HOrth(h_e, n_e)

    v_eOb_n = q_nb(v_eOb_b)
    Ob = Geographic(n_e, h_e)
    r_eOb_e = Cartesian(Ob)
    ω_el_n = get_ω_el_n(v_eOb_n, Ob)

    ω_el_l = q_nl'(ω_el_n)
    ω_el_b = q_lb'(ω_el_l)
    ω_lb_b = ω_eb_b - ω_el_b

    ω_ie_e = SVector(0., 0., ω_ie)
    ω_ie_b = q_eb'(ω_ie_e)
    ω_ib_b = ω_ie_b + ω_eb_b

    v_gnd = norm(v_eOb_n)
    χ_gnd = v_gnd > v_min_χγ ? azimuth(v_eOb_n) : 0.0
    γ_gnd = v_gnd > v_min_χγ ? inclination(v_eOb_n) : 0.0

    ẋ.q_lb = Attitude.dt(q_lb, ω_lb_b)
    ẋ.q_el = Attitude.dt(q_el, ω_el_l)

    ẋ.Δx = v_eOb_n[1]
    ẋ.Δy = v_eOb_n[2]
    ẋ.h_e = -v_eOb_n[3]

    sys.y = KinematicsY(
        KinData( e_nb, q_nb, q_eb, q_en, ϕ_λ, n_e, h_e, h_o, Δxy, r_eOb_e,
                 ω_lb_b, ω_eb_b, ω_ie_b, ω_ib_b,
                 v_eOb_b, v_eOb_n, v_gnd, χ_gnd, γ_gnd),
        LTFData(; q_lb, q_el, ω_el_l)
    )

end


function Systems.f_step!(sys::System{LTF}, ε = 1e-8)
    normalize_block!(sys.x.q_lb, ε)
    normalize_block!(sys.x.q_el, ε)
end


@inline function get_ω_el_n(v_eOb_n::AbstractVector{<:Real}, Ob::Geographic)

    (R_N, R_E) = radii(Ob)
    h_e = HEllip(Ob)

    return SVector{3}(
        v_eOb_n[2] / (R_E + Float64(h_e)),
        -v_eOb_n[1] / (R_N + Float64(h_e)),
        0.0)

end

########################## ECEF-based Kinematics #########################
##########################################################################

#fast, singularity-free (all-attitude, all-latitude) kinematic mechanization,
#appropriate for simulation.

struct ECEF <: AbstractKinematicDescriptor end

@kwdef struct ECEFData
    q_en::RQuat = RQuat()
    ω_el_n::SVector{3,Float64} = zeros(SVector{3})
end

Systems.X(::ECEF) = ComponentVector(
    q_eb = zeros(4), n_e = zeros(3), Δx = 0.0, Δy = 0.0, h_e = 0.0)

Systems.Y(::ECEF) = KinematicsY(KinData(), ECEFData())

function Systems.init!(sys::System{ECEF}, ic::Initializer = Initializer())

    @unpack x, u = sys
    @unpack q_nb, Ob, ω_lb_b, v_eOb_n, Δx, Δy = ic

    n_e = NVector(Ob)
    h_e = HEllip(Ob)

    q_en = ltf(n_e)
    q_eb = q_en ∘ q_nb

    ω_el_n = get_ω_el_n(v_eOb_n, Ob)
    ω_el_b = q_nb'(ω_el_n)
    ω_eb_b = ω_el_b + ω_lb_b
    v_eOb_b = q_nb'(v_eOb_n)

    u.ω_eb_b = ω_eb_b
    u.v_eOb_b = v_eOb_b

    x.q_eb = q_eb[:]
    x.n_e = n_e[:]
    x.Δx = Δx
    x.Δy = Δy
    x.h_e = h_e

    f_ode!(sys) #update ẋ and y, not strictly necessary

end

function Systems.f_ode!(sys::System{ECEF})

    @unpack ẋ, x, u = sys

    q_eb = RQuat(x.q_eb, normalization = false)
    n_e = NVector(x.n_e, normalization = false)
    ω_eb_b = SVector{3}(u.ω_eb_b)
    v_eOb_b = SVector{3}(u.v_eOb_b)
    Δxy = SVector(x.Δx, x.Δy)
    h_e = HEllip(x.h_e[1])
    h_o = HOrth(h_e, n_e)
    ϕ_λ = LatLon(n_e)

    q_en = ltf(n_e)
    q_nb = q_en' ∘ q_eb
    e_nb = REuler(q_nb)

    Ob = Geographic(n_e, h_e)
    r_eOb_e = Cartesian(Ob)
    v_eOb_n = q_nb(v_eOb_b)
    ω_el_n = get_ω_el_n(v_eOb_n, Ob)
    ω_el_b = q_nb'(ω_el_n)
    ω_lb_b = ω_eb_b - ω_el_b

    ω_ie_e = SVector(0., 0., ω_ie)
    ω_ie_b = q_eb'(ω_ie_e)
    ω_ib_b = ω_ie_b + ω_eb_b

    v_gnd = norm(v_eOb_n)
    χ_gnd = v_gnd > v_min_χγ ? azimuth(v_eOb_n) : 0.0
    γ_gnd = v_gnd > v_min_χγ ? inclination(v_eOb_n) : 0.0

    ẋ.q_eb = Attitude.dt(q_eb, ω_eb_b)
    ẋ.n_e = q_en(ω_el_n × SVector{3,Float64}(0,0,-1))
    ẋ.Δx = v_eOb_n[1]
    ẋ.Δy = v_eOb_n[2]
    ẋ.h_e = -v_eOb_n[3]

    sys.y = KinematicsY(
        KinData(e_nb, q_nb, q_eb, q_en, ϕ_λ, n_e, h_e, h_o, Δxy, r_eOb_e,
                ω_lb_b, ω_eb_b, ω_ie_b, ω_ib_b,
                v_eOb_b, v_eOb_n, v_gnd, χ_gnd, γ_gnd),
        ECEFData(; q_en, ω_el_n)
    )

end


function Systems.f_step!(sys::System{ECEF}, ε = 1e-8)
    normalize_block!(sys.x.q_eb, ε)
    normalize_block!(sys.x.n_e, ε)
end


################################ NED Kinematics ################################
################################################################################

#non singularity-free kinematic mechanization. useful mostly for aircraft
#control analysis and design

struct NED <: AbstractKinematicDescriptor end

@kwdef struct NEDData
    ω_nb_b::SVector{3,Float64} = zeros(SVector{3})
    ω_en_n::SVector{3,Float64} = zeros(SVector{3})
end

Systems.X(::NED) = ComponentVector(ψ_nb = 0.0, θ_nb = 0.0, φ_nb = 0.0,
                                ϕ = 0.0, λ = 0.0, Δx = 0.0, Δy = 0.0, h_e = 0.0)

Systems.Y(::NED) = KinematicsY(KinData(), NEDData())


function Systems.init!(sys::System{NED}, ic::Initializer = Initializer())

    @unpack x, u = sys
    @unpack q_nb, Ob, ω_lb_b, v_eOb_n, Δx, Δy = ic

    ω_el_n = get_ω_el_n(v_eOb_n, Ob)
    ω_el_b = q_nb'(ω_el_n)
    ω_eb_b = ω_el_b + ω_lb_b
    v_eOb_b = q_nb'(v_eOb_n)

    e_nb = REuler(q_nb)
    ϕ_λ = LatLon(Ob)
    h_e = HEllip(Ob)

    u.ω_eb_b = ω_eb_b
    u.v_eOb_b = v_eOb_b

    x.ψ_nb = e_nb.ψ
    x.θ_nb = e_nb.θ
    x.φ_nb = e_nb.φ
    x.ϕ = ϕ_λ.ϕ
    x.λ = ϕ_λ.λ
    x.Δx = Δx
    x.Δy = Δy
    x.h_e = h_e

    f_ode!(sys) #update ẋ and y, not strictly necessary

end

function Systems.f_ode!(sys::System{NED})

    @unpack ẋ, x, u = sys

    e_nb = REuler(x.ψ_nb, x.θ_nb, x.φ_nb)
    ϕ_λ = LatLon(x.ϕ, x.λ)
    h_e = HEllip(x.h_e[1])
    Δxy = SVector(x.Δx, x.Δy)
    ω_eb_b = SVector{3}(u.ω_eb_b)
    v_eOb_b = SVector{3}(u.v_eOb_b)

    n_e = NVector(ϕ_λ)
    h_o = HOrth(h_e, n_e)

    q_nb = RQuat(e_nb)
    q_en = ltf(n_e)
    q_eb = q_en ∘ q_nb

    v_eOb_n = q_nb(v_eOb_b)
    Ob = Geographic(n_e, h_e)
    r_eOb_e = Cartesian(Ob)

    ω_en_n = get_ω_en_n(v_eOb_n, Ob)
    ω_en_b = q_nb'(ω_en_n)
    ω_nb_b = ω_eb_b - ω_en_b

    ω_el_n = get_ω_el_n(v_eOb_n, Ob)
    ω_el_b = q_nb'(ω_el_n)
    ω_lb_b = ω_eb_b - ω_el_b

    ω_ie_e = SVector{3}(0, 0, ω_ie)
    ω_ie_b = q_eb'(ω_ie_e)
    ω_ib_b = ω_ie_b + ω_eb_b

    v_gnd = norm(v_eOb_n)
    χ_gnd = azimuth(v_eOb_n)
    γ_gnd = inclination(v_eOb_n)

    ė_nb = Attitude.dt(e_nb, ω_nb_b)
    ϕ_λ_dot = Geodesy.dt(ϕ_λ, ω_en_n)

    ẋ.ψ_nb = ė_nb.ψ
    ẋ.θ_nb = ė_nb.θ
    ẋ.φ_nb = ė_nb.φ
    ẋ.ϕ = ϕ_λ_dot.ϕ
    ẋ.λ = ϕ_λ_dot.λ
    ẋ.Δx = v_eOb_n[1]
    ẋ.Δy = v_eOb_n[2]
    ẋ.h_e = -v_eOb_n[3]

    sys.y = KinematicsY(
        KinData( e_nb, q_nb, q_eb, q_en, ϕ_λ, n_e, h_e, h_o, Δxy, r_eOb_e,
                 ω_lb_b, ω_eb_b, ω_ie_b, ω_ib_b,
                 v_eOb_b, v_eOb_n, v_gnd, χ_gnd, γ_gnd),
        NEDData(; ω_nb_b, ω_en_n)
    )

end

Systems.f_step!(::System{NED}) = nothing


@inline function get_ω_en_n(v_eOb_n::AbstractVector{<:Real}, Ob::Geographic)

    (R_N, R_E) = radii(Ob)
    h_e = HEllip(Ob)
    ϕ = LatLon(Ob).ϕ

    return SVector{3}(
        v_eOb_n[2] / (R_E + Float64(h_e)),
        -v_eOb_n[1] / (R_N + Float64(h_e)),
        -v_eOb_n[2] * tan(ϕ) / (R_E + Float64(h_e))
        )
end


################################# Plotting #####################################

# #@userplot allows defining a custom plot for a specific dataset without having
# #to create a custom type for dispatch. we just wrap the data in the userplot
# #type generated by Plots, and is received inside the recipe in its field "args"
@userplot Trajectory3D
@recipe function f(t3d::Trajectory3D)

    # https://daschw.github.io/recipes/#series_recipes

    path = t3d.args
    n = length(path)

    x_ext, y_ext, z_ext = map(extrema, (path.x, path.y, path.z))
    x_mid, y_mid, z_mid = map(v -> 0.5sum(v), (x_ext, y_ext, z_ext))
    x_span, y_span, z_span = map(v -> v[2] - v[1], (x_ext, y_ext, z_ext))
    span = max(x_span, y_span, z_span)

    x_bounds = (x_mid - 0.5span, x_mid + 0.5span)
    y_bounds = (y_mid - 0.5span, y_mid + 0.5span)
    z_bounds = (z_mid - 0.5span, z_mid + 0.5span)

    path_xp = StructArray((x = fill(x_bounds[1], n), y = path.y, z = path.z))
    path_yp = StructArray((x = path.x, y = fill(y_bounds[1], n), z = path.z))
    path_zp = StructArray((x = path.x, y = path.y, z = fill(z_bounds[1], n)))

    #--> sets default values, which can be overridden by each series using :=
    linewidth --> 3
    markersize --> 8
    xguide --> L"$\Delta x\ (m)$"
    yguide --> L"$\Delta y\ (m)$"
    zguide --> L"$h\ (m)$"
    legend --> false

    xlims --> x_bounds
    ylims --> y_bounds
    zlims --> z_bounds

    yflip --> true #needed for local NED frame

    #marker projection lines
    linecolor --> :lightgray
    for i in (1,n)
        @series begin
            linestyle := :dashdotdot
            linewidth := 2
            [path.x[i], path_xp.x[1]], [path.y[i], path.y[i]], [path.z[i], path.z[i]]
        end
        @series begin
            linestyle := :dashdotdot
            linewidth := 2
            [path.x[i], path.x[i]], [path.y[i], path_yp.y[1]], [path.z[i], path.z[i]]
        end
        @series begin
            linestyle := :dashdotdot
            linewidth := 2
            [path.x[i], path.x[i]], [path.y[i], path.y[i]], [path.z[i], path_zp.z[1]]
        end
    end

    #3D path, path projections and start-end markers
    for (p, c) in zip((path_xp, path_yp, path_zp, path),
                      (:darkgray, :darkgray, :darkgray, :blue))

        @series begin
            linestyle := :solid
            seriestype := :path3d
            linecolor := c
            p.x, p.y, p.z
        end
        @series begin
            seriestype := :scatter3d
            markercolor := :green
            [p.x[1]], [p.y[1]], [p.z[1]]
        end
        @series begin
            seriestype := :scatter3d
            markercolor := :red
            [p.x[end]], [p.y[end]], [p.z[end]]
        end

    end

    return nothing

end

function Plotting.make_plots(ts::TimeSeries{<:KinematicsY}; kwargs...)

    return OrderedDict(
        :data => make_plots(ts.data; kwargs...),
    )

end

function Plotting.make_plots(ts::TimeSeries{<:KinData}; kwargs...)

    pd = OrderedDict{Symbol, Plots.Plot}()

    plot_level = get(kwargs, :plot_level, :full)

    #example of capturing the plot_level keyword to control which plots are generated
    if plot_level == :simplified
        return pd #nothing also works
    end

    pd[:e_nb] = plot(
        ts.q_nb; #will automatically be converted to REuler for plotting
        plot_title = "Attitude (Vehicle/NED)",
        rot_ref = "n", rot_target = "b",
        kwargs...)

    #will be automatically converted to LatLon for plotting
    #remove the title added by the LatLon TH recipe
    subplot_latlon = plot(ts.n_e; title = "", ts_split = :v, kwargs...)

    #remove the title added by the Altitude TH recipe
    subplot_h = plot(ts.h_e; title = "", kwargs...)
                plot!(ts.h_o; title = "", kwargs...)

    subplot_xy = plot(
        ts.Δxy;
        label = [L"$\int v_{eb}^{x_n} dt$" L"$\int v_{eb}^{y_n} dt$"],
        ylabel = [L"$\Delta x\ (m)$" L"$\Delta y \ (m)$"],
        ts_split = :h, link = :none, kwargs...)

    pd[:Ob_geo] = plot(
        subplot_latlon, subplot_h;
        layout = grid(1, 2, widths = [0.67, 0.33]),
        plot_title = "Position (WGS84)",
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    pd[:Ob_xyh] = plot(
        subplot_xy, subplot_h;
        layout = grid(1, 2, widths = [0.67, 0.33]),
        plot_title = "Position (Local Cartesian)",
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    #when a plot is assembled from multiple subplots, the plot_titlefontsize
    #attribute no longer works, and it is titlefontisze what determines the font
    #size of the overall figure title (which normally is used for subplots).

    ts_Δx, ts_Δy = get_components(ts.Δxy)
    path = StructArray((x = ts_Δx._data, y = ts_Δy._data, z = Float64.(ts.h_e._data)))

    pd[:Ob_t3d] = plot(
        Trajectory3D(path);
        plot_title = "Trajectory (Local Cartesian, Ellipsoidal Altitude)",
        titlefontsize = 20,
        camera = (30, 15),
        kwargs...,
        size = (1920, 1920)
        )

    pd[:ω_lb_b] = plot(
        ts.ω_lb_b;
        plot_title = "Angular Velocity (Vehicle/LTF) [Vehicle Axes]",
        label = ["Roll Rate" "Pitch Rate" "Yaw Rate"],
        ylabel = [L"$p \ (rad/s)$" L"$q \ (rad/s)$" L"$r \ (rad/s)$"],
        ts_split = :h,
        kwargs...)

    pd[:v_eOb_n] = plot(
        ts.v_eOb_n;
        plot_title = "Velocity (Vehicle/ECEF) [NED Axes]",
        label = ["North" "East" "Down"],
        ylabel = [L"$v_{eb}^{N} \ (m/s)$" L"$v_{eb}^{E} \ (m/s)$" L"$v_{eb}^{D} \ (m/s)$"],
        ts_split = :h, link = :none,
        kwargs...)

    subplot_v_gnd = plot(ts.v_gnd; title = "Ground Speed",
        ylabel = L"$v_{gnd} \ (m/s)$", label = "", kwargs...)
    subplot_χ = plot(ts._t, rad2deg.(ts.χ_gnd._data); title = "Course Angle",
        ylabel = L"$\chi_{gnd} \ (deg)$", label = "", kwargs...)
    subplot_γ = plot(ts._t, rad2deg.(ts.γ_gnd._data); title = "Flight Path Angle",
        ylabel = L"$\gamma_{cv} \ (deg)$", label = "", kwargs...)

    pd[:vχγ] = plot(subplot_v_gnd, subplot_χ, subplot_γ;
        plot_title = "Velocity (Vehicle/ECEF) [NED Axes]",
        layout = (1,3),
        link = :none,
        kwargs...)

    pd[:v_eOb_b] = plot(
        ts.v_eOb_b;
        plot_title = "Velocity (Vehicle/ECEF) [Vehicle Axes]",
        ylabel = [L"$v_{eb}^{x_b} \ (m/s)$" L"$v_{eb}^{y_b} \ (m/s)$" L"$v_{eb}^{z_b} \ (m/s)$"],
        ts_split = :h,
        kwargs...)

    # pd[:ω_el_n] = plot(
    #     ts.ω_el_n;
    #     plot_title = "Local Tangent Frame Transport Rate (LTF/ECEF) [NED Axes]",
    #     ylabel = L"$\omega_{el}^{l} \ (rad/s)$",
    #     ts_split = :h,
    #     kwargs...)

    return pd

end


################################################################################
################################# GUI ##########################################

function GUI.draw(sys::KinSystem, p_open::Ref{Bool} = Ref(true),
                    label::String = "Kinematics")

    @unpack e_nb, ϕ_λ, h_e, h_o, Δxy, ω_lb_b, ω_eb_b, ω_ib_b, v_eOb_b, v_eOb_n  = sys.y.data

    CImGui.Begin(label, p_open)

    if CImGui.TreeNode("Angular Velocity (Body / LTF) [Body]")

        CImGui.Text(@sprintf("Yaw Rate: %.7f deg/s", rad2deg(ω_lb_b[1])))
        CImGui.Text(@sprintf("Pitch Rate: %.7f deg/s", rad2deg(ω_lb_b[2])))
        CImGui.Text(@sprintf("Roll Rate: %.7f deg/s", rad2deg(ω_lb_b[3])))

        CImGui.TreePop()
    end

    GUI.draw(rad2deg.(ω_eb_b - ω_lb_b), "Transport Rate (LTF / ECEF) [Body]", "deg/s")
    GUI.draw(rad2deg.(ω_ib_b),"Angular Velocity (Body / ECI) [Body]", "deg/s")
    GUI.draw(v_eOb_n, "Velocity (O / ECEF) [NED]", "m/s")
    GUI.draw(v_eOb_b, "Velocity (O / ECEF) [Body]", "m/s")

    if CImGui.TreeNode("Attitude (Body / NED)")

        @unpack ψ, θ, φ = e_nb
        CImGui.Text(@sprintf("Heading: %.7f deg", rad2deg(ψ)))
        CImGui.Text(@sprintf("Inclination: %.7f deg", rad2deg(θ)))
        CImGui.Text(@sprintf("Bank: %.7f deg", rad2deg(φ)))

        CImGui.TreePop()
    end

    if CImGui.TreeNode("Position (O / ECEF)")

        @unpack ϕ, λ = ϕ_λ
        CImGui.Text(@sprintf("Latitude: %.7f deg", rad2deg(ϕ)))
        CImGui.Text(@sprintf("Longitude: %.7f deg", rad2deg(λ)))
        CImGui.Text(@sprintf("Northward Increment: %.7f m", Δxy[1]))
        CImGui.Text(@sprintf("Eastward Increment: %.7f m", Δxy[2]))
        CImGui.Text(@sprintf("Altitude (Ellipsoidal): %.7f m", Float64(h_e)))
        CImGui.Text(@sprintf("Altitude (Orthometric): %.7f m", Float64(h_o)))

        CImGui.TreePop()
    end

    CImGui.End()

end


end #module