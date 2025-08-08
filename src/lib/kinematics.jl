module Kinematics

using StaticArrays, StructArrays, ComponentArrays, LinearAlgebra, UnPack
using Plots, LaTeXStrings, DataStructures

using Flight.FlightCore

using ..Attitude
using ..Geodesy

export AbstractKinematicDescriptor, ECEF, WA, NED
export KinInit, KinData, KinSystem

const v_min_χγ = 0.1 #minimum speed for valid χ, γ


############################## Initialization ##################################
################################################################################

#user-friendly kinematic conditions for body frame initialization
struct Initializer
    q_nb::RQuat #attitude with respect to NED frame
    n_e::NVector #2D location, n-vector
    h_e::Altitude{Ellipsoidal} #ellipsoidal altitude
    ω_wb_b::SVector{3, Float64} #angular velocity with respect to local level frame, body coordinates
    v_eb_n::SVector{3, Float64} #Earth-relative velocity, NED coordinates
    Δx::Float64 #Northward velocity integral
    Δy::Float64 #Eastward velocity integral
end

const KinInit = Initializer

function Initializer(;
    q_nb::Abstract3DRotation = RQuat(), location::Abstract2DLocation = LatLon(),
    h::Altitude = HOrth(), ω_wb_b::AbstractVector{<:Real} = zeros(SVector{3}),
    v_eb_n::AbstractVector{<:Real} = zeros(SVector{3}), Δx::Real = 0.0, Δy::Real = 0.0)

    n_e = NVector(location)
    h_e = HEllip(h, n_e)

    Initializer(q_nb, n_e, h_e, ω_wb_b, v_eb_n, Δx, Δy)
end


################################### KinData ####################################
################################################################################

struct KinData
    e_nb::REuler #body frame attitude with respect to NED frame, Euler angles
    q_nb::RQuat #body frame attitude with respect to NED frame, quaternion
    q_eb::RQuat #body frame attitude with respect to ECEF frame, quaternion
    q_en::RQuat #NED frame attitude with respect to ECEF frame, quaternion
    ϕ_λ::LatLon #2D location, latitude / longitude
    n_e::NVector #2D location, n-vector
    h_e::Altitude{Ellipsoidal} #ellipsoidal altitude
    h_o::Altitude{Orthometric} #orthometric altitude
    Δxy::SVector{2,Float64} #horizontal velocity integral
    r_eb_e::SVector{3,Float64} #Cartesian ECEF position
    ω_wb_b::SVector{3,Float64} #angular velocity with respect to local level frame, body coordinates
    ω_eb_b::SVector{3,Float64} #angular velocity with respect to ECEF frame, body coordinates
    v_eb_b::SVector{3,Float64} #Earth-relative velocity, body coordinates
    v_eb_n::SVector{3,Float64} #Earth-relative velocity, NED coordinates
    v_gnd::Float64 #Earth-relative speed
    χ_gnd::Float64 #Earth-relative course angle
    γ_gnd::Float64 #Earth-relative flight path angle
end

function KinData(ic::KinInit = KinInit())

    @unpack q_nb, n_e, h_e, ω_wb_b, v_eb_n, Δx, Δy = ic

    Ob = Geographic(n_e, h_e)
    e_nb = REuler(q_nb)
    q_en = ltf(n_e)
    q_eb = q_en ∘ q_nb

    ϕ_λ = LatLon(n_e)
    h_o = HOrth(h_e, n_e)
    Δxy = SVector(Δx, Δy)
    r_eb_e = Cartesian(Ob)

    ω_ew_n = get_ω_ew_n(v_eb_n, Ob)
    ω_ew_b = q_nb'(ω_ew_n)
    ω_eb_b = ω_ew_b + ω_wb_b

    v_eb_b = q_nb'(v_eb_n)

    v_gnd = norm(v_eb_n)
    χ_gnd = v_gnd > v_min_χγ ? azimuth(v_eb_n) : 0.0
    γ_gnd = v_gnd > v_min_χγ ? inclination(v_eb_n) : 0.0

    KinData(   e_nb, q_nb, q_eb, q_en, ϕ_λ, n_e, h_e, h_o, Δxy, r_eb_e,
                ω_wb_b, ω_eb_b, v_eb_b, v_eb_n, v_gnd, χ_gnd, γ_gnd)

end

Base.getproperty(data::KinData, s::Symbol) = getproperty(data, Val(s))

@generated function Base.getproperty(data::KinData, ::Val{S}) where {S}
    if S ∈ fieldnames(KinData)
        return :(getfield(data, $(QuoteNode(S))))
    elseif S === :psi || S === :ψ
        return :(getfield(data, :e_nb).ψ)
    elseif S === :theta || S === :θ
        return :(getfield(data, :e_nb).θ)
    elseif S === :phi || S === :φ
        return :(getfield(data, :e_nb).φ)
    elseif S === :lat || S === :ϕ
        return :(getfield(data, :ϕ_λ).ϕ)
    elseif S === :lon || S === :λ
        return :(getfield(data, :ϕ_λ).λ)
    else
        error("$(typeof(data)) has no property $S")
    end
end


function normalize_block!(x, ε)
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


######################### AbstractKinematicDescriptor ##########################
################################################################################

abstract type AbstractKinematicDescriptor <: ModelDefinition end
@no_disc AbstractKinematicDescriptor

const XVelTemplate = ComponentVector(ω_eb_b = zeros(3), v_eb_b = zeros(3))

Modeling.U(::AbstractKinematicDescriptor) = zero(XVelTemplate)
Modeling.Y(::AbstractKinematicDescriptor) = KinData()

const KinSystem = Model{<:AbstractKinematicDescriptor}

KinData(mdl::KinSystem) = mdl.y

########################### WA-based Kinematics #########################
##########################################################################

#fast, singularity-free (all-attitude, all-latitude) kinematic mechanization,
#appropriate for simulation. WA = wander-azimuth frame

struct WA <: AbstractKinematicDescriptor end

Modeling.X(::WA) = ComponentVector(
    q_wb = zeros(4), q_ew = zeros(4), Δx = 0.0, Δy = 0.0, h_e = 0.0)

function Modeling.init!(mdl::Model{WA}, ic::Initializer = Initializer())

    @unpack x, u = mdl
    @unpack q_nb, n_e, h_e, ω_wb_b, v_eb_n, Δx, Δy = ic

    Ob = Geographic(n_e, h_e)
    ω_ew_n = get_ω_ew_n(v_eb_n, Ob)
    ω_ew_b = q_nb'(ω_ew_n)
    ω_eb_b = ω_ew_b + ω_wb_b
    v_eb_b = q_nb'(v_eb_n)
    h_e = HEllip(Ob)

    q_wb = q_nb #arbitrarily initializes wander angle ψ_nw to 1

    u.ω_eb_b = ω_eb_b
    u.v_eb_b = v_eb_b

    x.q_wb = q_wb[:]
    x.q_ew = ltf(n_e)[:]
    x.Δx = Δx
    x.Δy = Δy
    x.h_e = h_e

    f_ode!(mdl) #update ẋ and y, not strictly necessary

end


function Modeling.f_ode!(mdl::Model{WA})

    @unpack ẋ, x, u = mdl

    q_wb = RQuat(x.q_wb, normalization = false)
    q_ew = RQuat(x.q_ew, normalization = false)
    ω_eb_b = SVector{3}(u.ω_eb_b)
    v_eb_b = SVector{3}(u.v_eb_b)
    h_e = HEllip(x.h_e[1])
    Δxy = SVector(x.Δx, x.Δy)

    ψ_nw = get_ψ_nw(q_ew)
    q_nw = Rz(ψ_nw)
    q_nb = q_nw ∘ q_wb
    q_eb = q_ew ∘ q_wb
    q_en = q_eb ∘ q_nb'
    e_nb = REuler(q_nb)

    n_e = NVector(q_ew)
    ϕ_λ = LatLon(n_e)
    h_o = HOrth(h_e, n_e)

    v_eb_n = q_nb(v_eb_b)
    Ob = Geographic(n_e, h_e)
    r_eb_e = Cartesian(Ob)
    ω_ew_n = get_ω_ew_n(v_eb_n, Ob)

    ω_ew_w = q_nw'(ω_ew_n)
    ω_ew_b = q_wb'(ω_ew_w)
    ω_wb_b = ω_eb_b - ω_ew_b

    v_gnd = norm(v_eb_n)
    χ_gnd = v_gnd > v_min_χγ ? azimuth(v_eb_n) : 0.0
    γ_gnd = v_gnd > v_min_χγ ? inclination(v_eb_n) : 0.0

    ẋ.q_wb = Attitude.dt(q_wb, ω_wb_b)
    ẋ.q_ew = Attitude.dt(q_ew, ω_ew_w)

    ẋ.Δx = v_eb_n[1]
    ẋ.Δy = v_eb_n[2]
    ẋ.h_e = -v_eb_n[3]

    mdl.y = KinData( e_nb, q_nb, q_eb, q_en, ϕ_λ, n_e, h_e, h_o, Δxy, r_eb_e,
                 ω_wb_b, ω_eb_b, v_eb_b, v_eb_n, v_gnd, χ_gnd, γ_gnd)

end


function Modeling.f_step!(mdl::Model{WA}, ε = 1e-8)
    normalize_block!(mdl.x.q_wb, ε)
    normalize_block!(mdl.x.q_ew, ε)
end


@inline function get_ω_ew_n(v_eb_n::AbstractVector{<:Real}, Ob::Geographic)

    (R_N, R_E) = radii(Ob)
    h_e = HEllip(Ob)

    return SVector{3}(
        v_eb_n[2] / (R_E + Float64(h_e)),
        -v_eb_n[1] / (R_N + Float64(h_e)),
        0.0)

end

########################## ECEF-based Kinematics #########################
##########################################################################

#fast, singularity-free (all-attitude, all-latitude) kinematic mechanization,
#appropriate for simulation.

struct ECEF <: AbstractKinematicDescriptor end

Modeling.X(::ECEF) = ComponentVector(
    q_eb = zeros(4), n_e = zeros(3), Δx = 0.0, Δy = 0.0, h_e = 0.0)

function Modeling.init!(mdl::Model{ECEF}, ic::Initializer = Initializer())

    @unpack x, u = mdl
    @unpack q_nb, n_e, h_e, ω_wb_b, v_eb_n, Δx, Δy = ic

    Ob = Geographic(n_e, h_e)

    q_en = ltf(n_e)
    q_eb = q_en ∘ q_nb

    ω_ew_n = get_ω_ew_n(v_eb_n, Ob)
    ω_ew_b = q_nb'(ω_ew_n)
    ω_eb_b = ω_ew_b + ω_wb_b
    v_eb_b = q_nb'(v_eb_n)

    u.ω_eb_b = ω_eb_b
    u.v_eb_b = v_eb_b

    x.q_eb = q_eb[:]
    x.n_e = n_e[:]
    x.Δx = Δx
    x.Δy = Δy
    x.h_e = h_e

    f_ode!(mdl) #update ẋ and y, not strictly necessary

end

function Modeling.f_ode!(mdl::Model{ECEF})

    @unpack ẋ, x, u = mdl

    q_eb = RQuat(x.q_eb, normalization = false)
    n_e = NVector(x.n_e, normalization = false)
    ω_eb_b = SVector{3}(u.ω_eb_b)
    v_eb_b = SVector{3}(u.v_eb_b)
    Δxy = SVector(x.Δx, x.Δy)
    h_e = HEllip(x.h_e[1])
    h_o = HOrth(h_e, n_e)
    ϕ_λ = LatLon(n_e)

    q_en = ltf(n_e)
    q_nb = q_en' ∘ q_eb
    e_nb = REuler(q_nb)

    Ob = Geographic(n_e, h_e)
    r_eb_e = Cartesian(Ob)
    v_eb_n = q_nb(v_eb_b)
    ω_ew_n = get_ω_ew_n(v_eb_n, Ob)
    ω_ew_b = q_nb'(ω_ew_n)
    ω_wb_b = ω_eb_b - ω_ew_b

    v_gnd = norm(v_eb_n)
    χ_gnd = v_gnd > v_min_χγ ? azimuth(v_eb_n) : 0.0
    γ_gnd = v_gnd > v_min_χγ ? inclination(v_eb_n) : 0.0

    ẋ.q_eb = Attitude.dt(q_eb, ω_eb_b)
    ẋ.n_e = q_en(ω_ew_n × SVector{3,Float64}(0,0,-1))
    ẋ.Δx = v_eb_n[1]
    ẋ.Δy = v_eb_n[2]
    ẋ.h_e = -v_eb_n[3]

    mdl.y = KinData(e_nb, q_nb, q_eb, q_en, ϕ_λ, n_e, h_e, h_o, Δxy, r_eb_e,
                ω_wb_b, ω_eb_b, v_eb_b, v_eb_n, v_gnd, χ_gnd, γ_gnd)

end

function Modeling.f_step!(mdl::Model{ECEF}, ε = 1e-8)
    normalize_block!(mdl.x.q_eb, ε)
    normalize_block!(mdl.x.n_e, ε)
end


################################ NED Kinematics ################################
################################################################################

#non singularity-free kinematic mechanization. useful mostly for aircraft
#control analysis and design

struct NED <: AbstractKinematicDescriptor end

Modeling.X(::NED) = ComponentVector(ψ_nb = 0.0, θ_nb = 0.0, φ_nb = 0.0,
                                ϕ = 0.0, λ = 0.0, Δx = 0.0, Δy = 0.0, h_e = 0.0)


function Modeling.init!(mdl::Model{NED}, ic::Initializer = Initializer())

    @unpack x, u = mdl
    @unpack q_nb, n_e, h_e, ω_wb_b, v_eb_n, Δx, Δy = ic

    Ob = Geographic(n_e, h_e)
    ω_ew_n = get_ω_ew_n(v_eb_n, Ob)
    ω_ew_b = q_nb'(ω_ew_n)
    ω_eb_b = ω_ew_b + ω_wb_b
    v_eb_b = q_nb'(v_eb_n)

    e_nb = REuler(q_nb)
    ϕ_λ = LatLon(n_e)

    u.ω_eb_b = ω_eb_b
    u.v_eb_b = v_eb_b

    x.ψ_nb = e_nb.ψ
    x.θ_nb = e_nb.θ
    x.φ_nb = e_nb.φ
    x.ϕ = ϕ_λ.ϕ
    x.λ = ϕ_λ.λ
    x.Δx = Δx
    x.Δy = Δy
    x.h_e = h_e

    f_ode!(mdl) #update ẋ and y, not strictly necessary

end

function Modeling.f_ode!(mdl::Model{NED})

    @unpack ẋ, x, u = mdl

    e_nb = REuler(x.ψ_nb, x.θ_nb, x.φ_nb)
    ϕ_λ = LatLon(x.ϕ, x.λ)
    h_e = HEllip(x.h_e[1])
    Δxy = SVector(x.Δx, x.Δy)
    ω_eb_b = SVector{3}(u.ω_eb_b)
    v_eb_b = SVector{3}(u.v_eb_b)

    n_e = NVector(ϕ_λ)
    h_o = HOrth(h_e, n_e)

    q_nb = RQuat(e_nb)
    q_en = ltf(n_e)
    q_eb = q_en ∘ q_nb

    v_eb_n = q_nb(v_eb_b)
    Ob = Geographic(n_e, h_e)
    r_eb_e = Cartesian(Ob)

    ω_en_n = get_ω_en_n(v_eb_n, Ob)
    ω_en_b = q_nb'(ω_en_n)
    ω_nb_b = ω_eb_b - ω_en_b

    ω_ew_n = get_ω_ew_n(v_eb_n, Ob)
    ω_ew_b = q_nb'(ω_ew_n)
    ω_wb_b = ω_eb_b - ω_ew_b

    v_gnd = norm(v_eb_n)
    χ_gnd = azimuth(v_eb_n)
    γ_gnd = inclination(v_eb_n)

    ė_nb = Attitude.dt(e_nb, ω_nb_b)
    ϕ_λ_dot = Geodesy.dt(ϕ_λ, ω_en_n)

    ẋ.ψ_nb = ė_nb.ψ
    ẋ.θ_nb = ė_nb.θ
    ẋ.φ_nb = ė_nb.φ
    ẋ.ϕ = ϕ_λ_dot.ϕ
    ẋ.λ = ϕ_λ_dot.λ
    ẋ.Δx = v_eb_n[1]
    ẋ.Δy = v_eb_n[2]
    ẋ.h_e = -v_eb_n[3]

    mdl.y = KinData( e_nb, q_nb, q_eb, q_en, ϕ_λ, n_e, h_e, h_o, Δxy, r_eb_e,
                 ω_wb_b, ω_eb_b, v_eb_b, v_eb_n, v_gnd, χ_gnd, γ_gnd)

end

Modeling.f_step!(::Model{NED}) = nothing


@inline function get_ω_en_n(v_eb_n::AbstractVector{<:Real}, Ob::Geographic)

    (R_N, R_E) = radii(Ob)
    h_e = HEllip(Ob)
    ϕ = LatLon(Ob).ϕ

    return SVector{3}(
        v_eb_n[2] / (R_E + Float64(h_e)),
        -v_eb_n[1] / (R_N + Float64(h_e)),
        -v_eb_n[2] * tan(ϕ) / (R_E + Float64(h_e))
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

    pd[:ω_wb_b] = plot(
        ts.ω_wb_b;
        plot_title = "Angular Velocity (Vehicle/WA) [Vehicle Axes]",
        label = ["Roll Rate" "Pitch Rate" "Yaw Rate"],
        ylabel = [L"$p \ (rad/s)$" L"$q \ (rad/s)$" L"$r \ (rad/s)$"],
        ts_split = :h,
        kwargs...)

    pd[:v_eb_n] = plot(
        ts.v_eb_n;
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

    pd[:v_eb_b] = plot(
        ts.v_eb_b;
        plot_title = "Velocity (Vehicle/ECEF) [Vehicle Axes]",
        ylabel = [L"$v_{eb}^{x_b} \ (m/s)$" L"$v_{eb}^{y_b} \ (m/s)$" L"$v_{eb}^{z_b} \ (m/s)$"],
        ts_split = :h,
        kwargs...)

    # pd[:ω_ew_n] = plot(
    #     ts.ω_ew_n;
    #     plot_title = "Wander-Azimuth Transport Rate (WA/ECEF) [NED Axes]",
    #     ylabel = L"$\omega_{el}^{l} \ (rad/s)$",
    #     ts_split = :h,
    #     kwargs...)

    return pd

end


################################################################################
################################# GUI ##########################################


GUI.draw(dyn::Model{<:AbstractKinematicDescriptor}) = GUI.draw(KinData(dyn))

function GUI.draw(data::KinData, p_open::Ref{Bool} = Ref(true),
                    label::String = "Kinematic Data")

    @unpack e_nb, ϕ_λ, h_e, h_o, Δxy, ω_wb_b, ω_eb_b, v_eb_b, v_eb_n  = data

    CImGui.Begin(label, p_open)

    if CImGui.TreeNode("Angular Velocity (Body / WA) [Body]")

        CImGui.Text(@sprintf("Yaw Rate: %.7f deg/s", rad2deg(ω_wb_b[1])))
        CImGui.Text(@sprintf("Pitch Rate: %.7f deg/s", rad2deg(ω_wb_b[2])))
        CImGui.Text(@sprintf("Roll Rate: %.7f deg/s", rad2deg(ω_wb_b[3])))

        CImGui.TreePop()
    end

    GUI.draw(rad2deg.(ω_eb_b - ω_wb_b), "Transport Rate (WA / ECEF) [Body]", "deg/s")
    GUI.draw(v_eb_n, "Velocity (Ob / ECEF) [NED]", "m/s")
    GUI.draw(v_eb_b, "Velocity (Ob / ECEF) [Body]", "m/s")

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


################################################################################
############################### XPlane12Control ######################################

function Network.XPlanePose(kin_data::KinData)

    @unpack ϕ_λ, e_nb, h_o = kin_data

    ϕ = rad2deg(ϕ_λ.ϕ)
    λ = rad2deg(ϕ_λ.λ)
    h = Float64(h_o)

    ψ = rad2deg(e_nb.ψ)
    θ = rad2deg(e_nb.θ)
    φ = rad2deg(e_nb.φ)

    Network.XPlanePose(ϕ, λ, h, ψ, θ, φ)

end

end #module