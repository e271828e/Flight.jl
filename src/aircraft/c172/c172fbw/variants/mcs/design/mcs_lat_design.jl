module MCSLatDesign

using Flight
using Flight.FlightCore

using Flight.FlightPhysics
using Flight.FlightComponents
using Flight.FlightComponents.Control.Continuous: LinearizedSS
using Flight.FlightComponents.Control.Discrete: PIDParams, LQRTrackerParams
using Flight.FlightComponents.Control.PIDOpt: Settings, Metrics, optimize_PID, build_PID, check_results

using Flight.FlightAircraft.C172
using Flight.FlightAircraft.C172FBW
using Flight.FlightAircraft.C172FBWMCS

using HDF5
using Logging
using UnPack
using ControlSystems
using RobustAndOptimalControl
using StaticArrays
using StructArrays
using ComponentArrays
using LinearAlgebra


function generate_lookups(
    EAS_range::AbstractRange{Float64} = range(25, 50, length = 3),
    h_range::AbstractRange{Float64} = range(10, 1000, length = 2);
    channel::Symbol = :lat,
    folder::String = dirname(@__DIR__)) #save to parent folder by default

    if channel === :lon
        f_opt = design_lon
    elseif channel === :lat
        f_opt = design_lat
    else
        error("Valid values for channel keyword: :lon, :lat")
    end

    results = map(Iterators.product(EAS_range, h_range)) do (EAS, h)

        println("Designing controllers for EAS = $EAS, h = $h")

        #all other design point parameters at default
        flaps = EAS < 30 ? 1.0 : 0.0
        design_point = C172.TrimParameters(; Ob = Geographic(LatLon(), HOrth(h)), EAS, flaps)

        results = f_opt(design_point)
        return results

    end |> StructArray |> StructArrays.components

    return results

    filenames = joinpath.(dirname(@__DIR__), "data", string.(keys(results)) .* "_lookup.h5")

    bounds = (EAS = (EAS_range[1], EAS_range[end]), h = (h_range[1], h_range[end]))

    foreach(values(results), filenames) do results, fname

        save_lookup(data, bounds, joinpath(folder, fname))
    end

end

#save lookup. assumes uniformly scaled grid.

#this would work for a Union of Array{<:LQRTrackerParams} and PID Params

function save_lookup(data::Array{<:LQRTrackerParams, N}) where {N}
# function save_lookup(data::Array{<:LQRTrackerParams, N},
#                     bounds::NTuple{N, Tuple{Real,Real}},
#                     fname::String) where {N}

    data_nt = StructArrays.components(StructArray(data))

    fname = joinpath.(dirname(@__DIR__), "data", "test_lookup2.h5")

    fid = h5open(fname, "w")
    foreach(keys(data_nt), values(data_nt)) do k, v
        fid[string(k)] = stack(v)
    end
    close(fid)



    ############################## READBACK ####################################
    fid = h5open(fname, "r")

    #read fieldnames as ordered in LQRTrackerParams and into an instance
    params_stacked = LQRTrackerParams(map(name -> read(fid[string(name)]), fieldnames(LQRTrackerParams))...)

    close(fid)

    @unpack x_trim, u_trim, z_trim, C_fbk, C_fwd, C_int = params_stacked

    #from each 1+N dimensional stacked array (x_trim, u_trim, z_trim) read from
    #the HDF5, generate a N-dimensional array of SVectors of appropriate size
    x_trim_data, u_trim_data, z_trim_data = map((x_trim, u_trim, z_trim)) do a
        map(SVector{size(a)[1]}, eachslice(a; dims = Tuple(2:ndims(a))))
    end

    #from each 2+N dimensional stacked array (C_fbk, C_fwd, C_int) read from the
    #HDF5, generate a N-dimensional array of SMatrices of appropriate size
    C_fbk_data, C_fwd_data, C_int_data = map((C_fbk, C_fwd, C_int)) do a
        map(SMatrix{size(a)[1], size(a)[2]}, eachslice(a; dims = Tuple(3:ndims(a))))
    end

    #sizes of interpolating dimensions
    interp_sizes = size(x_trim_data)
    return interp_sizes


    # return LQRTrackerParams(; x_trim, u_trim, z_trim, C_fbk, C_fwd, C_int)

    #now we generate the interpolators. we can also do it for an arbitrary
    #number of interpolation dimensions.

end

function Control.Discrete.LQRTrackerParams(a::Array{<:LQRTrackerParams, N}) where {N}
    println("Hi")

end


function design_lat(design_point::C172.TrimParameters = C172.TrimParameters())

    ac = Cessna172FBWBase(NED()) |> System #linearization requires NED kinematics

    P_lss_lat = Control.Continuous.LinearizedSS(ac, design_point; model = :lat);

    x_labels_lat = keys(P_lss_lat.x0) |> collect
    y_labels_lat = keys(P_lss_lat.y0) |> collect
    u_labels_lat = keys(P_lss_lat.u0) |> collect

    #remove heading from state vector to avoid quasi-static equilibrium
    x_labels = deleteat!(x_labels_lat, findfirst(isequal(:ψ), x_labels_lat))
    y_labels = deleteat!(y_labels_lat, findfirst(isequal(:ψ), y_labels_lat))
    u_labels = u_labels_lat

    #ensure consistency in component selection and ordering between our design model
    #and MCS avionics implementation for state and control vectors
    @assert tuple(x_labels...) === propertynames(C172FBWMCS.XLat())
    @assert tuple(u_labels...) === propertynames(C172FBWMCS.ULat())

    #extract design model
    P_lss = Control.Continuous.submodel(P_lss_lat; x = x_labels, u = u_labels, y = y_labels)
    P_nss = named_ss(P_lss)


    ############################### φ + β ######################################

    z_labels = [:φ, :β]
    @assert tuple(z_labels...) === propertynames(C172FBWMCS.ZLatPhiBeta())

    x_trim = P_lss.x0
    u_trim = P_lss.u0
    z_trim = P_lss.y0[z_labels]

    n_x = length(P_lss.x0)
    n_u = length(P_lss.u0)
    n_z = length(z_labels)

    F = P_lss.A
    G = P_lss.B
    Hx = P_lss.C[z_labels, :]
    Hu = P_lss.D[z_labels, :]

    Hx_int = Hx[:β, :]'
    Hu_int = Hu[:β, :]'
    n_int, _ = size(Hx_int)

    F_aug = [F zeros(n_x, n_int); Hx_int zeros(n_int, n_int)]
    G_aug = [G; Hu_int]
    Hx_aug = [Hx zeros(n_z, n_int)]
    Hu_aug = Hu

    P_aug = ss(F_aug, G_aug, Hx_aug, Hu_aug)

    @unpack v_x, v_y = P_lss.x0
    v_norm = norm([v_x, v_y])

    #weight matrices
    Q = ComponentVector(p = 0, r = 0.1, φ = 0.15, v_x = 0/v_norm, v_y = 0.1/v_norm, β_filt = 0, ail_v = 0, ail_p = 0, rud_v = 0, rud_p = 0, ξ_β = 0.001) |> diagm
    R = C172FBWMCS.ULat(aileron_cmd = 0.1, rudder_cmd = 0.05) |> diagm

    #compute gain matrix
    C_aug = lqr(P_aug, Q, R)

    #extract system state and integrator blocks from the feedback matrix
    C_x = C_aug[:, 1:n_x]
    C_ξ = C_aug[:, n_x+1:end]

    #construct
    A = [F G; Hx Hu]
    B = inv(A)
    B_12 = B[1:n_x, n_x+1:end]
    B_22 = B[n_x+1:end, n_x+1:end]

    C_fbk = C_x
    C_fwd = B_22 + C_x * B_12
    C_int = ComponentMatrix(zeros(n_u, n_z), Axis(u_labels), Axis(z_labels))
    C_int[:, :β] .= C_ξ

    #convert everything to plain arrays
    params_φβ = LQRTrackerParams(;
        C_fbk = Matrix(C_fbk), C_fwd = Matrix(C_fwd), C_int = Matrix(C_int),
        x_trim = Vector(x_trim), u_trim = Vector(u_trim), z_trim = Vector(z_trim))

    #some useful signal labels
    u_labels_fbk = Symbol.(string.(u_labels) .* "_fbk")
    u_labels_fwd = Symbol.(string.(u_labels) .* "_fwd")
    u_labels_sum = Symbol.(string.(u_labels) .* "_sum")
    u_labels_int_u = Symbol.(string.(u_labels) .* "_int_u")
    u_labels_int = Symbol.(string.(u_labels) .* "_int")
    u_labels_ξ = Symbol.(string.(u_labels) .* "_ξ")

    z_labels_sp = Symbol.(string.(z_labels) .* "_sp")
    z_labels_sp1 = Symbol.(string.(z_labels) .* "_sp1")
    z_labels_sp2 = Symbol.(string.(z_labels) .* "_sp2")
    z_labels_err = Symbol.(string.(z_labels) .* "_err")
    z_labels_sum = Symbol.(string.(z_labels) .* "_sum")
    z_labels_sp_fwd = Symbol.(string.(z_labels) .* "_sp_fwd")
    z_labels_sp_sum = Symbol.(string.(z_labels) .* "_sp_sum")


    C_fbk_ss = named_ss(ss(C_fbk), u = x_labels, y = u_labels_fbk)

    C_fwd_ss = named_ss(ss(C_fwd), u = z_labels_sp_fwd, y = u_labels_fwd)

    C_int_ss = named_ss(ss(C_int), u = z_labels_err, y = u_labels_int_u)

    int_ss = named_ss(ss(tf(1, [1,0])) .* I(2),
                        x = u_labels_ξ,
                        u = u_labels_int_u,
                        y = u_labels_int);

    P_nss = named_ss(P_lss);

    φ_err_sum = sumblock("φ_err = φ_sum - φ_sp_sum")
    β_err_sum = sumblock("β_err = β_sum - β_sp_sum")

    aileron_cmd_sum = sumblock("aileron_cmd_sum = aileron_cmd_fwd - aileron_cmd_fbk - aileron_cmd_int")
    rudder_cmd_sum = sumblock("rudder_cmd_sum = rudder_cmd_fwd - rudder_cmd_fbk - rudder_cmd_int")

    φ_sp_splitter = splitter(:φ_sp, 2)
    β_sp_splitter = splitter(:β_sp, 2)

    connections = vcat(
        Pair.(x_labels, x_labels),
        Pair.(z_labels, z_labels_sum),
        Pair.(z_labels_sp1, z_labels_sp_sum),
        Pair.(z_labels_sp2, z_labels_sp_fwd),
        Pair.(z_labels_err, z_labels_err),
        Pair.(u_labels_sum, u_labels),
        Pair.(u_labels_fwd, u_labels_fwd),
        Pair.(u_labels_fbk, u_labels_fbk),
        Pair.(u_labels_int, u_labels_int),
        Pair.(u_labels_int_u, u_labels_int_u),
        )

    Logging.disable_logging(Logging.Warn)
    P_nss_φβ = connect([P_nss, int_ss, C_fwd_ss, C_fbk_ss, C_int_ss,
                        φ_err_sum, β_err_sum,
                        aileron_cmd_sum, rudder_cmd_sum,
                        φ_sp_splitter, β_sp_splitter], connections;
                        w1 = z_labels_sp, z1 = P_nss.y)
    Logging.disable_logging(Logging.LogLevel(typemin(Int32)))


    ############################### p + β ######################################

    P_φ2p = P_nss_φβ[:p, :φ_sp];

    p2φ_int = tf(1, [1, 0]) |> ss
    P_p2φ_opt = series(p2φ_int, ss(P_φ2p))

    t_sim_p2φ = 10
    lower_bounds = PIDParams(; k_p = 0.1, k_i = 0.0, k_d = 0.0, τ_f = 0.01)
    upper_bounds = PIDParams(; k_p = 10.0, k_i = 35.0, k_d = 1.5, τ_f = 0.01)
    settings = Settings(; t_sim = t_sim_p2φ, lower_bounds, upper_bounds)
    weights = Metrics(; Ms = 1, ∫e = 10, ef = 2, ∫u = 0.1, up = 0.00)
    params_0 = PIDParams(; k_p = 1.5, k_i = 3, k_d = 0.1, τ_f = 0.01)

    p2φ_results = optimize_PID(P_p2φ_opt; params_0, settings, weights, global_search = false)
    params_p2φ = p2φ_results.params
    if !check_results(p2φ_results, Metrics(; Ms = 1.3, ∫e = 0.1, ef = 0.02, ∫u = Inf, up = Inf))
        println("Warning: Checks failed for roll rate control PID, design point $(design_point)")
        println(p2φ_results.metrics)
    end

    p2φ_PID = build_PID(p2φ_results.params)
    C_p2φ = named_ss(series(p2φ_int, p2φ_PID), :C_p2φ; u = :p_err, y = :φ_sp)

    p2φ_sum = sumblock("p_err = p_sp - p")
    P_nss_pβ = connect([p2φ_sum, C_p2φ, P_nss_φβ], [:p_err=>:p_err, :p=>:p, :φ_sp=>:φ_sp], w1 = [:p_sp, :β_sp], z1 = P_nss_φβ.y)

    # @unpack k_p, k_i, k_d, T_i, T_d = p2φ_results.params
    # @unpack k_p, k_i, k_d, T_i, T_d = p2φ_results.params
    # @show k_p, k_i, k_d, T_i, T_d
    # @show p2φ_results.metrics
    # @show p2φ_results.exit_flag

    return (φβ = params_φβ, p2φ = params_p2φ)


end


end #module