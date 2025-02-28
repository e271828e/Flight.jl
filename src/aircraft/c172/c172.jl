module C172

using LinearAlgebra, StaticArrays, ComponentArrays, UnPack, HDF5, Interpolations
using Logging
using Reexport
using NLopt

using Flight.FlightCore
using Flight.FlightLib


################################################################################
################################ Airframe #####################################

struct Airframe <: SystemDefinition end

# This component represents the platform's structure, together with any
# components rigidly attached to it, such as powerplant or landing gear, but not
# payload or fuel contents. Its mass corresponds roughly to the aircraft's
# Standard Empty Weight

#Airframe mass properties computed in the vehicle reference frame b
const mp_b_afm = let
    #define the airframe as a RigidBodyDistribution
    afm_c = RigidBodyDistribution(767.0, SA[820.0 0 0; 0 1164.0 0; 0 0 1702.0])
    #a RigidBodyDistribution is always specified in a reference frame c with
    #origin at its center of mass. now, define the transform from vehicle
    #reference frame b to reference frame c (pure translation)
    t_bc = FrameTransform(r = SVector{3}(0.056, 0, 0.582))
    #translate the airframe's mass properties to frame b
    MassProperties(afm_c, t_bc)
end

#the airframe itself receives no external actions. these are considered to act
#upon the vehicle's aerodynamics, power plant and landing gear. the same goes
#for rotational angular momentum.
Dynamics.get_mp_b(::System{Airframe}) = mp_b_afm
Dynamics.get_wr_b(::System{Airframe}) = Wrench()
Dynamics.get_hr_b(::System{Airframe}) = zeros(SVector{3})


################################################################################
############################ Aerodynamics ######################################

function generate_aero_data(fname = joinpath(@__DIR__, "aero.h5"))

    h5open(fname, "w") do fid

        ############################## C_D ##################################

        create_group(fid, "C_D")
        C_D = fid["C_D"]

        C_D["zero"] = 0.027

        create_group(C_D, "δe")
        C_D["δe"]["δe"] = [-1.0 0.0 1.0] |> vec
        C_D["δe"]["data"] = [ 0.06 0 0.06] |> vec

        create_group(C_D, "β")
        C_D["β"]["β"] = [-1.0 0.0 1.0] |> vec
        C_D["β"]["data"] = [ 0.17 0 0.17] |> vec

        create_group(C_D, "ge")
        C_D["ge"]["Δh_nd"] = [ 0.0000 0.1000 0.1500 0.2000 0.3000 0.4000 0.5000 0.6000 0.7000 0.8000 0.9000 1.0000 1.1000 ] |> vec
        C_D["ge"]["data"] = [ 0.4800 0.5150 0.6290 0.7090 0.8150 0.8820 0.9280 0.9620 0.9880 1.0000 1.0000 1.0000 1.0000 ] |> vec

        create_group(C_D, "δf")
        C_D["δf"]["δf"] = deg2rad.([ 0.0000	10.0000 20.0000 30.0000 ]) |> vec
        C_D["δf"]["data"] = [ 0.0000 0.0070 0.0120 0.0180 ] |> vec

        create_group(C_D, "α_δf")
        C_D["α_δf"]["α"] = [ -0.0873 -0.0698 -0.0524 -0.0349 -0.0175 0.0000	0.0175	0.0349	0.0524	0.0698	0.0873 0.1047	0.1222	0.1396	0.1571	0.1745	0.1920	0.2094	0.2269	0.2443	0.2618	0.2793	0.2967	0.3142	0.3316	0.3491] |> vec
        C_D["α_δf"]["δf"] = deg2rad.([ 0.0000	10.0000	20.0000	30.0000 ]) |> vec
        C_D["α_δf"]["data"] = [ 0.0041	0.0000	0.0005	0.0014
                                0.0013	0.0004	0.0025	0.0041
                                0.0001	0.0023	0.0059	0.0084
                                0.0003	0.0057	0.0108	0.0141
                                0.0020	0.0105	0.0172	0.0212
                                0.0052	0.0168	0.0251	0.0299
                                0.0099	0.0248	0.0346	0.0402
                                0.0162	0.0342	0.0457	0.0521
                                0.0240	0.0452	0.0583	0.0655
                                0.0334	0.0577	0.0724	0.0804
                                0.0442	0.0718	0.0881	0.0968
                                0.0566	0.0874	0.1053	0.1148
                                0.0706	0.1045	0.1240	0.1343
                                0.0860	0.1232	0.1442	0.1554
                                0.0962	0.1353	0.1573	0.1690
                                0.1069	0.1479	0.1708	0.1830
                                0.1180	0.1610	0.1849	0.1975
                                0.1298	0.1746	0.1995	0.2126
                                0.1424	0.1892	0.2151	0.2286
                                0.1565	0.2054	0.2323	0.2464
                                0.1727	0.2240	0.2521	0.2667
                                0.1782	0.2302	0.2587	0.2735
                                0.1716	0.2227	0.2507	0.2653
                                0.1618	0.2115	0.2388	0.2531
                                0.1475	0.1951	0.2214	0.2351
                                0.1097	0.1512	0.1744	0.1866
        ]

        ############################## C_Y ##################################

        create_group(fid, "C_Y")
        C_Y = fid["C_Y"]

        C_Y["δr"] = 0.1870
        C_Y["δa"] = 0.0

        create_group(C_Y, "β_δf")
        C_Y["β_δf"]["β"] = [-0.3490 0 0.3490] |> vec
        C_Y["β_δf"]["δf"] = deg2rad.([0 30]) |> vec
        C_Y["β_δf"]["data"] = [
                            0.1370	0.1060
                            0.0000	0.0000
                            -0.1370	-0.1060
        ]
        create_group(C_Y, "p")
        C_Y["p"]["α"] = [0.0 0.094] |> vec
        C_Y["p"]["δf"] = deg2rad.([0 30]) |> vec
        C_Y["p"]["data"] = [
                            -0.0750	-0.1610
                            -0.1450	-0.2310
        ]
        create_group(C_Y, "r")
        C_Y["r"]["α"] = [0.0 0.094] |> vec
        C_Y["r"]["δf"] = deg2rad.([0 30]) |> vec
        C_Y["r"]["data"] = [
                            0.2140	0.1620
                            0.2670	0.2150
        ]


        ############################### C_L #################################

        create_group(fid, "C_L")
        C_L = fid["C_L"]

        C_L["δe"] = 0.4300
        C_L["q"] = 3.900
        C_L["α_dot"] = 1.700

        create_group(C_L, "ge")
        C_L["ge"]["Δh_nd"] = [ 0.0000 0.1000 0.1500 0.2000 0.3000 0.4000 0.5000 0.6000 0.7000 0.8000 0.9000 1.0000 1.1000 ] |> vec
        C_L["ge"]["data"] = [ 1.2030 1.1270 1.0900 1.0730 1.0460 1.0550 1.0190 1.0130 1.0080 1.0060 1.0030 1.0020 1.0000 ] |> vec

        create_group(C_L, "α")
        C_L["α"]["α"] = [ -0.0900 0.0000	0.0900	0.1000	0.1200	0.1400	0.1600	0.1700	0.1900	0.2100	0.2400	0.2600	0.2800	0.3000	0.3200	0.3400	0.3600	] |> vec
        C_L["α"]["stall"] = [0.0 1.0] |> vec
        C_L["α"]["data"] = [-0.2200	-0.2200
                           	0.2500	0.2500
                           	0.7300	0.7300
                           	0.8300	0.7800
                           	0.9200	0.7900
                           	1.0200	0.8100
                           	1.0800	0.8200
                           	1.1300	0.8300
                           	1.1900	0.8500
                           	1.2500	0.8600
                           	1.3500	0.8800
                           	1.4400	0.9000
                           	1.4700	0.9200
                           	1.4300	0.9500
                           	1.3800	0.9900
                           	1.3000	1.0500
                           	1.1500	1.1500
        ]

        create_group(C_L, "δf")
        C_L["δf"]["δf"] = deg2rad.([ 0.0000	10.0000 20.0000 30.0000 ]) |> vec
        C_L["δf"]["data"] = [ 0.0000 0.2 0.3 0.35] |> vec


        ############################### C_l #################################

        create_group(fid, "C_l")
        C_l = fid["C_l"]

        C_l["δa"] = 0.229
        C_l["δr"] = 0.0147
        C_l["β"] = -0.09226
        C_l["p"] = -0.4840

        create_group(C_l, "r")
        C_l["r"]["α"] = [0.0 0.094] |> vec
        C_l["r"]["δf"] = deg2rad.([0 30]) |> vec
        C_l["r"]["data"] = [
                            0.0798	0.1246
                            0.1869	0.2317
        ]

        ############################# C_m ###################################

        create_group(fid, "C_m")
        C_m = fid["C_m"]

        C_m["zero"] = 0.100
        C_m["δe"] = -1.1220
        C_m["α"] = -1.8000
        C_m["q"] = -12.400
        C_m["α_dot"] = -7.2700

        create_group(C_m, "δf")
        C_m["δf"]["δf"] = deg2rad.([0 10 20 30]) |> vec
        C_m["δf"]["data"] = [ 0.0000 -0.0654 -0.0981 -0.1140 ] |> vec

        ############################# C_n ###################################

        create_group(fid, "C_n")
        C_n = fid["C_n"]

        C_n["δr"] = -0.0430
        C_n["δa"] = -0.0053
        C_n["β"] = 0.05874
        C_n["p"] = -0.0278
        C_n["r"] = -0.0937

    end

end

function generate_aero_lookup(fname = joinpath(@__DIR__, "aero.h5"))

    fid = h5open(fname, "r")

    gr_C_D = fid["C_D"]
    gr_C_Y = fid["C_Y"]
    gr_C_L = fid["C_L"]
    gr_C_l = fid["C_l"]
    gr_C_m = fid["C_m"]
    gr_C_n = fid["C_n"]

    C_D = (
        z = gr_C_D["zero"] |> read,
        β = linear_interpolation(gr_C_D["β"]["β"] |> read, gr_C_D["β"]["data"] |> read, extrapolation_bc = Flat()),
        δe = linear_interpolation(gr_C_D["δe"]["δe"] |> read, gr_C_D["δe"]["data"] |> read, extrapolation_bc = Flat()),
        δf = linear_interpolation(gr_C_D["δf"]["δf"] |> read, gr_C_D["δf"]["data"] |> read, extrapolation_bc = Flat()),
        α_δf = linear_interpolation((gr_C_D["α_δf"]["α"] |> read,  gr_C_D["α_δf"]["δf"] |> read), gr_C_D["α_δf"]["data"] |> read, extrapolation_bc = Flat()),
        ge = linear_interpolation(gr_C_D["ge"]["Δh_nd"] |> read, gr_C_D["ge"]["data"] |> read, extrapolation_bc = Flat())
    )

    C_Y = (
        δr = gr_C_Y["δr"] |> read,
        δa = gr_C_Y["δa"] |> read,
        β_δf = linear_interpolation((gr_C_Y["β_δf"]["β"] |> read,  gr_C_Y["β_δf"]["δf"] |> read), gr_C_Y["β_δf"]["data"] |> read, extrapolation_bc = Flat()),
        p = linear_interpolation((gr_C_Y["p"]["α"] |> read,  gr_C_Y["p"]["δf"] |> read), gr_C_Y["p"]["data"] |> read, extrapolation_bc = Flat()),
        r = linear_interpolation((gr_C_Y["r"]["α"] |> read,  gr_C_Y["r"]["δf"] |> read), gr_C_Y["r"]["data"] |> read, extrapolation_bc = Flat()),
    )

    C_L = (
        δe = gr_C_L["δe"] |> read,
        q = gr_C_L["q"] |> read,
        α_dot = gr_C_L["α_dot"] |> read,
        α = linear_interpolation((gr_C_L["α"]["α"] |> read,  gr_C_L["α"]["stall"] |> read), gr_C_L["α"]["data"] |> read, extrapolation_bc = Flat()),
        δf = linear_interpolation(gr_C_L["δf"]["δf"] |> read, gr_C_L["δf"]["data"] |> read, extrapolation_bc = Flat()),
        ge = linear_interpolation(gr_C_L["ge"]["Δh_nd"] |> read, gr_C_L["ge"]["data"] |> read, extrapolation_bc = Flat())
    )

    C_l = (
        δa = gr_C_l["δa"] |> read,
        δr = gr_C_l["δr"] |> read,
        β = gr_C_l["β"] |> read,
        p = gr_C_l["p"] |> read,
        r = linear_interpolation((gr_C_l["r"]["α"] |> read,  gr_C_l["r"]["δf"] |> read), gr_C_l["r"]["data"] |> read, extrapolation_bc = Flat()),
    )

    C_m = (
        z = gr_C_m["zero"] |> read,
        δe = gr_C_m["δe"] |> read,
        α = gr_C_m["α"] |> read,
        q = gr_C_m["q"] |> read,
        α_dot = gr_C_m["α_dot"] |> read,
        δf = linear_interpolation(gr_C_m["δf"]["δf"] |> read, gr_C_m["δf"]["data"] |> read, extrapolation_bc = Flat()),
    )

    C_n = (
        δr = gr_C_n["δr"] |> read,
        δa = gr_C_n["δa"] |> read,
        β = gr_C_n["β"] |> read,
        p = gr_C_n["p"] |> read,
        r = gr_C_n["r"] |> read,
    )

    close(fid)

    return (C_D = C_D, C_Y = C_Y, C_L = C_L, C_l = C_l, C_m = C_m, C_n = C_n)

end

#the aircraft body reference frame fb is arbitrarily chosen to coincide with
#the aerodynamics frame fa, so the frame transform is trivial
const f_ba = FrameTransform()
const aero_lookup = generate_aero_lookup()

# if this weren't the case, and we cared not only about the rotation but also
#about the velocity lever arm, here's the rigorous way of computing v_wOa_a:
# v_wb_b = v_eb_b - v_ew_b
# v_ea_b = v_eb_b + ω_eb_b × r_ba_b
# v_wa_b = v_ea_b - v_ew_b = v_eb_b + ω_eb_b × r_ba_b - v_ew_b
# v_wa_b = v_wb_b + ω_eb_b × r_ba_b
# v_wa_a = q_ba'(v_wa_b)

@kwdef struct AeroCoeffs
    C_D::Float64 = 0.0
    C_Y::Float64 = 0.0
    C_L::Float64 = 0.0
    C_l::Float64 = 0.0
    C_m::Float64 = 0.0
    C_n::Float64 = 0.0
end

#need to pass aerodynamic data as an input, because since Julia 1.9 accessing
#the fields of a global const allocates
function get_aero_coeffs(lookup = aero_lookup; α, β, p_nd, q_nd, r_nd, δa, δr, δe, δf, α_dot_nd, β_dot_nd, Δh_nd, stall)

    #set sensible bounds
    α = clamp(α, -0.1, 0.36) #0.36 is the highest value (post-stall) tabulated for C_L
    β = clamp(β, -0.2, 0.2)
    α_dot_nd = clamp(α_dot_nd, -0.04, 0.04)
    β_dot_nd = clamp(β_dot_nd, -0.2, 0.2)

    @unpack C_D, C_Y, C_L, C_l, C_m, C_n = lookup

    AeroCoeffs(
        C_D = C_D.z + C_D.ge(Δh_nd) * (C_D.α_δf(α,δf) + C_D.δf(δf)) + C_D.δe(δe) + C_D.β(β),
        C_Y = C_Y.δr * δr + C_Y.δa * δa + C_Y.β_δf(β,δf) + C_Y.p(α,δf) * p_nd + C_Y.r(α,δf) * r_nd,
        C_L = C_L.ge(Δh_nd) * (C_L.α(α,stall) + C_L.δf(δf)) + C_L.δe * δe + C_L.q * q_nd + C_L.α_dot * α_dot_nd,
        C_l = C_l.δa * δa + C_l.δr * δr + C_l.β * β + C_l.p * p_nd + C_l.r(α,δf) * r_nd,
        C_m = C_m.z + C_m.δe * δe + C_m.δf(δf) + C_m.α * α + C_m.q * q_nd + C_m.α_dot * α_dot_nd,
        C_n = C_n.δr * δr + C_n.δa * δa + C_n.β * β + C_n.p * p_nd + C_n.r * r_nd,
    )

end

@kwdef struct Aero <: SystemDefinition
    S::Float64 = 16.165 #wing area
    b::Float64 = 10.912 #wingspan
    c::Float64 = 1.494 #mean aerodynamic chord
    δe_range::NTuple{2,Float64} = deg2rad.((-28, 23)) #elevator deflection range (rad)
    δa_range::NTuple{2,Float64} = deg2rad.((-20, 20)) #aileron deflection range (rad)
    δr_range::NTuple{2,Float64} = deg2rad.((-16, 16)) #rudder deflection range (rad)
    δf_range::NTuple{2,Float64} = deg2rad.((0, 30)) #flap deflection range (rad)
    α_stall::NTuple{2,Float64} = (0.09, 0.36) #α values for stall hysteresis switching
    V_min::Float64 = 1.0 #lower airspeed threshold for non-dimensional angle rates
    τ::Float64 = 0.02 #time constant for filtered airflow angle derivatives
end

#e↑ -> δe↑ -> trailing edge down -> Cm↓ -> pitch down
#a↑ -> δa↑ -> left trailing edge down, right up -> Cl↑ -> roll right
#r↑ -> δr↑ -> rudder trailing edge left -> Cn↓ -> yaw left
#f↑ -> δf↑ -> flap trailing edge down -> CL↑

@kwdef mutable struct AeroU
    e::Ranged{Float64, -1., 1.} = 0.0
    a::Ranged{Float64, -1., 1.} = 0.0
    r::Ranged{Float64, -1., 1.} = 0.0
    f::Ranged{Float64, 0., 1.} = 0.0
end

@kwdef mutable struct AeroS #discrete state
    stall::Bool = false
end

@kwdef struct AeroY
    e::Float64 = 0.0 #normalized elevator control input
    a::Float64 = 0.0 #normalized aileron control input
    r::Float64 = 0.0 #normalized rudder control input
    f::Float64 = 0.0 #normalized flap control input
    α::Float64 = 0.0 #clamped AoA, aerodynamic axes
    β::Float64 = 0.0 #clamped AoS, aerodynamic axes
    α_filt::Float64 = 0.0 #filtered AoA
    β_filt::Float64 = 0.0 #filtered AoS
    α_filt_dot::Float64 = 0.0 #filtered AoA derivative
    β_filt_dot::Float64 = 0.0 #filtered AoS derivative
    stall::Bool = false #stall state
    coeffs::AeroCoeffs = AeroCoeffs() #aerodynamic coefficients
    wr_b::Wrench = Wrench() #aerodynamic Wrench, vehicle frame
end

Systems.X(::Aero) = ComponentVector(α_filt = 0.0, β_filt = 0.0) #filtered airflow angles
Systems.Y(::Aero) = AeroY()
Systems.U(::Aero) = AeroU()
Systems.S(::Aero) = AeroS()

#caution: do not confuse the w-frame in kinematics.ω_wb_b, which refers to the
#wander-azimuth frame (w), with the w in air.v_wb_b, which indicates aerodynamic
#(wind-relative) velocity

function Systems.f_ode!(sys::System{Aero}, ::System{<:PistonThruster},
    air::AirData, kinematics::KinData, terrain::AbstractTerrain)

    @unpack ẋ, x, u, s, constants = sys
    @unpack α_filt, β_filt = x
    @unpack e, a, r, f = u
    @unpack S, b, c, δe_range, δa_range, δr_range, δf_range, α_stall, V_min, τ = constants
    @unpack TAS, q, v_wb_b = air
    @unpack ω_wb_b, n_e, h_o = kinematics
    stall = s.stall

    v_wb_a = f_ba.q'(v_wb_b)

    #for near-zero TAS, the airflow angles are likely to chatter between 0, -π
    #and π. to avoid this, we set a minimum TAS for airflow computation. in this
    #scenario dynamic pressure will be close to zero, so forces and moments will
    #vanish anyway.
    α, β = (TAS > 0.1 ? Air.get_airflow_angles(v_wb_a) : (0.0, 0.0))
    V = max(TAS, V_min) #avoid division by zero

    α_filt_dot = 1/τ * (α - α_filt)
    β_filt_dot = 1/τ * (β - β_filt)

    p_nd = ω_wb_b[1] * b / (2V) #non-dimensional roll rate
    q_nd = ω_wb_b[2] * c / (2V) #non-dimensional pitch rate
    r_nd = ω_wb_b[3] * b / (2V) #non-dimensional yaw rate

    α_dot_nd = α_filt_dot * c / (2V)
    β_dot_nd = β_filt_dot * b / (2V)

    δe = linear_scaling(e, δe_range)
    δa = linear_scaling(a, δa_range)
    δr = linear_scaling(r, δr_range)
    δf = linear_scaling(f, δf_range)

    #non-dimensional height above ground
    loc_Oa = n_e #(2D location of aerodynamics frame, Oa = Ob)
    h_a = h_o #orthometric
    h_trn_a = TerrainData(terrain, loc_Oa).altitude #orthometric
    Δh_nd = (h_a - h_trn_a) / b

    # T = get_wr_b(pwp).F[1]
    # C_T = T / (q * S) #thrust coefficient, not used here

    coeffs = get_aero_coeffs(;
        α, β, p_nd, q_nd, r_nd, δa, δr, δe, δf, α_dot_nd, β_dot_nd, Δh_nd, stall)

    @unpack C_D, C_Y, C_L, C_l, C_m, C_n = coeffs

    q_as = Air.get_stability_axes(α)
    F_aero_s = q * S * SVector{3,Float64}(-C_D, C_Y, -C_L)
    F_aero_a = q_as(F_aero_s)
    τ_aero_a = q * S * SVector{3,Float64}(C_l * b, C_m * c, C_n * b)

    # wr_b = wr_a = Wrench(F_aero_a, τ_aero_a)
    wr_b = Wrench(F_aero_a, τ_aero_a)

    ẋ.α_filt = α_filt_dot
    ẋ.β_filt = β_filt_dot

    sys.y = AeroY(; α, α_filt, α_filt_dot, β, β_filt, β_filt_dot,
        e, a, r, f, stall, coeffs, wr_b)

    return nothing

end


function Systems.f_step!(sys::System{Aero})
    #stall hysteresis
    α = sys.y.α
    α_stall = sys.constants.α_stall
    if α > α_stall[2]
        sys.s.stall = true
    elseif α < α_stall[1]
        sys.s.stall = false
    end
end

Dynamics.get_mp_b(::System{<:Aero}) = MassProperties()
Dynamics.get_hr_b(::System{<:Aero}) = zeros(SVector{3})
Dynamics.get_wr_b(sys::System{Aero}) = sys.y.wr_b


################################# GUI ##########################################


function GUI.draw(sys::System{<:Aero}, p_open::Ref{Bool} = Ref(true),
                window_label::String = "Cessna 172R Aerodynamics")

    @unpack e, a, r, f, α, β, α_filt, β_filt, stall, coeffs, wr_b = sys.y
    @unpack C_D, C_Y, C_L, C_l, C_m, C_n = coeffs

    CImGui.Begin(window_label, p_open)

        CImGui.Text(@sprintf("Elevator Input: %.7f", e))
        CImGui.Text(@sprintf("Aileron Input: %.7f", a))
        CImGui.Text(@sprintf("Rudder Input: %.7f", r))
        CImGui.Text(@sprintf("Flap Setting: %.7f", f))
        CImGui.Text(@sprintf("AoA [Aero]: %.7f deg", rad2deg(α)))
        CImGui.Text(@sprintf("Filtered AoA [Aero]: %.7f deg", rad2deg(α_filt)))
        CImGui.Text(@sprintf("AoS [Aero]: %.7f deg", rad2deg(β)))
        CImGui.Text(@sprintf("Filtered AoS [Aero]: %.7f deg", rad2deg(β_filt)))
        CImGui.Text("Stall Status: $stall")

        if CImGui.TreeNode("Aerodynamic Coefficients")

            CImGui.Text(@sprintf("C_D: %.7f", C_D))
            CImGui.Text(@sprintf("C_Y: %.7f", C_Y))
            CImGui.Text(@sprintf("C_L: %.7f", C_L))
            CImGui.Text(@sprintf("C_l: %.7f", C_l))
            CImGui.Text(@sprintf("C_m: %.7f", C_m))
            CImGui.Text(@sprintf("C_n: %.7f", C_n))

            CImGui.TreePop()
        end

        GUI.draw(wr_b.F, "Aerodynamic Force (O) [Body]", "N")
        GUI.draw(wr_b.τ, "Aerodynamic Torque (O) [Body]", "N*m")

    CImGui.End()

end


###############################################################################
############################# Landing Gear ####################################

struct Ldg <: SystemDefinition
    left::LandingGearUnit{NoSteering, DirectBraking, Strut{SimpleDamper}}
    right::LandingGearUnit{NoSteering, DirectBraking, Strut{SimpleDamper}}
    nose::LandingGearUnit{DirectSteering, NoBraking, Strut{SimpleDamper}}
end

function Ldg()

    mlg_damper = SimpleDamper(k_s = 39404, #2700 lbf/ft
                              k_d_ext = 9340, #640 lbf/(ft/s)
                              k_d_cmp = 9340,
                              F_max = 50000)
    nlg_damper = SimpleDamper(k_s = 26269, #1800 lbf/ft
                              k_d_ext = 3503, #240 lbf/(ft/s)
                              k_d_cmp = 3503,
                              F_max = 50000)

    left = LandingGearUnit(
        strut = Strut(
            t_bs = FrameTransform(r = [-0.381, -1.092, 1.902], q = RQuat() ),
            l_0 = 0.0,
            damper = mlg_damper),
        braking = DirectBraking())

    right = LandingGearUnit(
        strut = Strut(
            t_bs = FrameTransform(r = [-0.381, 1.092, 1.902], q = RQuat() ),
            l_0 = 0.0,
            damper = mlg_damper),
        braking = DirectBraking())

    nose = LandingGearUnit(
        strut = Strut(
            t_bs = FrameTransform(r = [1.27, 0, 1.9] , q = RQuat()),
            l_0 = 0.0,
            damper = nlg_damper),
        steering = DirectSteering())

    Ldg(left, right, nose)

end

Dynamics.get_mp_b(::System{Ldg}) = MassProperties()
Dynamics.get_hr_b(::System{Ldg}) = zeros(SVector{3})

function GUI.draw(sys::System{<:Ldg}, p_open::Ref{Bool} = Ref(true),
                window_label::String = "Cessna 172R Landing Gear")

    @unpack left, right, nose = sys

    CImGui.Begin(window_label, p_open)

        @cstatic(c_left=false, c_right=false, c_nose = false, begin
            @c CImGui.Checkbox("Left Main", &c_left)
            @c CImGui.Checkbox("Right Main", &c_right)
            @c CImGui.Checkbox("Nose", &c_nose)
            c_left && @c GUI.draw(left, &c_left)
            c_right && @c GUI.draw(right, &c_right)
            c_nose && @c GUI.draw(nose, &c_nose)
        end)

    CImGui.End()

end

################################################################################
################################# Payload ######################################

@kwdef struct Payload <: SystemDefinition
    pilot_slot::FrameTransform = FrameTransform(r = SVector{3}(0.183, -0.356, 0.899))
    copilot_slot::FrameTransform = FrameTransform(r = SVector{3}(0.183, 0.356, 0.899))
    lpass_slot::FrameTransform = FrameTransform(r = SVector{3}(-0.681, -0.356, 0.899))
    rpass_slot::FrameTransform = FrameTransform(r = SVector{3}(-0.681, 0.356, 0.899))
    baggage_slot::FrameTransform = FrameTransform(r = SVector{3}(-1.316, 0, 0.899))
end

@kwdef mutable struct PayloadU
    m_pilot::Ranged{Float64, 0., 100.} = 75.0
    m_copilot::Ranged{Float64, 0., 100.} = 75.0
    m_lpass::Ranged{Float64, 0., 100.} = 0.0
    m_rpass::Ranged{Float64, 0., 100.} = 0.0
    m_baggage::Ranged{Float64, 0., 100.} = 50.0
end

Systems.U(::Payload) = PayloadU()

function Dynamics.get_mp_b(sys::System{Payload})
    @unpack m_pilot, m_copilot, m_lpass, m_rpass, m_baggage = sys.u
    @unpack pilot_slot, copilot_slot, lpass_slot, rpass_slot, baggage_slot = sys.constants

    pilot = MassProperties(PointDistribution(m_pilot), pilot_slot)
    copilot = MassProperties(PointDistribution(m_copilot), copilot_slot)
    lpass = MassProperties(PointDistribution(m_lpass), lpass_slot)
    rpass = MassProperties(PointDistribution(m_rpass), rpass_slot)
    baggage = MassProperties(PointDistribution(m_baggage), baggage_slot)

    mp_b = MassProperties() + pilot + copilot + lpass + rpass + baggage
    return mp_b
end

Dynamics.get_wr_b(::System{Payload}) = Wrench()
Dynamics.get_hr_b(::System{Payload}) = zeros(SVector{3})

#################################### GUI #######################################

function GUI.draw!(sys::System{<:Payload},
                    p_open::Ref{Bool} = Ref(true),
                    label::String = "Cessna 172R Payload")

    u = sys.u

    CImGui.Begin(label, p_open)

    CImGui.PushItemWidth(-250)

    u.m_pilot = GUI.safe_slider("Pilot Mass (kg)", u.m_pilot, "%.3f", true)
    u.m_copilot = GUI.safe_slider("Copilot Mass (kg)", u.m_copilot, "%.3f", true)
    u.m_lpass = GUI.safe_slider("Left Passenger Mass (kg)", u.m_lpass, "%.3f", true)
    u.m_rpass = GUI.safe_slider("Right Passenger Mass (kg)", u.m_rpass, "%.3f", true)
    u.m_baggage = GUI.safe_slider("Baggage Mass (kg)", u.m_baggage, "%.3f", true)

    CImGui.PopItemWidth()

    CImGui.End()

end


################################################################################
################################# Fuel #########################################

#assumes fuel is drawn equally from both tanks, no need to model them
#individually for now
@kwdef struct Fuel <: SystemDefinition
    m_full::Float64 = 114.4 #maximum fuel mass (42 gal * 6 lb/gal * 0.454 kg/lb)
    m_res::Float64 = 1.0 #residual fuel mass
end

@kwdef struct FuelY
    x_avail::Float64 = 0.0 #available fuel fraction
    m_total::Float64 = 0.0 #total fuel mass
    m_avail::Float64 = 0.0 #available fuel mass
end

#normalized fuel content (0: residual, 1: full)
Systems.X(::Fuel) = [0.5] #cannot be a scalar, need an AbstractVector{<:Real}
Systems.Y(::Fuel) = FuelY()

function Systems.f_ode!(sys::System{Fuel}, pwp::System{<:PistonThruster})

    @unpack m_full, m_res = sys.constants #no need for subsystems
    x_avail = sys.x[1]
    m_total = m_res + x_avail * (m_full - m_res) #current mass
    m_avail = m_total - m_res
    sys.ẋ .= -pwp.y.engine.ṁ / (m_full - m_res)
    sys.y = FuelY(; x_avail, m_total, m_avail)

end

function Dynamics.get_mp_b(fuel::System{Fuel})

    #in case x becomes negative (fuel consumed beyond x=0 before the engine
    #dies)
    m_fuel = max(0.0, fuel.y.m_total)

    m_left = PointDistribution(0.5m_fuel)
    m_right = PointDistribution(0.5m_fuel)

    #fuel tanks reference frames
    frame_left = FrameTransform(r = SVector{3}(0.325, -2.845, 0))
    frame_right = FrameTransform(r = SVector{3}(0.325, 2.845, 0))

    mp_b = MassProperties()
    mp_b += MassProperties(m_left, frame_left)
    mp_b += MassProperties(m_right, frame_right)

    return mp_b
end

Dynamics.get_wr_b(::System{Fuel}) = Wrench()
Dynamics.get_hr_b(::System{Fuel}) = zeros(SVector{3})

is_fuel_available(sys::System{<:Fuel}) = (sys.y.m_avail > 0)

function GUI.draw(sys::System{Fuel}, p_open::Ref{Bool} = Ref(true),
                window_label::String = "Cessna 172R Fuel System")

    @unpack m_total, m_avail = sys.y

    CImGui.Begin(window_label, p_open)

        CImGui.Text(@sprintf("Total Fuel: %.6f kg", m_total))
        CImGui.Text(@sprintf("Available Fuel: %.6f kg", m_avail))

    CImGui.End()

end


################################################################################
################################# Actuation ####################################

#any Actuation system suitable for the C172 Components
abstract type Actuation <: SystemDefinition end

function assign!(aero::System{<:Aero}, ldg::System{<:Ldg},
                pwp::System{<:PistonThruster}, act::System{<:Actuation})
    throw(MethodError(C172.assign!, (aero, ldg, pwp, act)))
end

Dynamics.get_mp_b(::System{<:Actuation}) = MassProperties()
Dynamics.get_wr_b(::System{<:Actuation}) = Wrench()
Dynamics.get_hr_b(::System{<:Actuation}) = zeros(SVector{3})

################################################################################
################################ Components ######################################

#P is introduced as a type parameter, because PistonThruster is itself a
#parametric type, and therefore not concrete
struct Components{P <: PistonThruster, A <: Actuation} <: AbstractComponents
    afm::Airframe
    aero::Aero
    ldg::Ldg
    fuel::Fuel
    pld::Payload
    pwp::P
    act::A
end

function Components(pwp::PistonThruster, act::Actuation)
    Components(Airframe(), Aero(), Ldg(), Fuel(), Payload(), pwp, act)
end


############################# Update Methods ###################################

function Systems.f_ode!(components::System{<:Components},
                        kin::KinData,
                        air::AirData,
                        trn::AbstractTerrain)

    @unpack act, aero, pwp, ldg, fuel, pld = components

    f_ode!(act) #update actuation system outputs
    assign!(aero, ldg, pwp, act) #assign actuation system outputs to components subsystems
    f_ode!(aero, pwp, air, kin, trn) #update aerodynamics continuous state & outputs
    f_ode!(ldg, kin, trn) #update landing gear continuous state & outputs
    f_ode!(pwp, air, kin) #update powerplant continuous state & outputs
    f_ode!(fuel, pwp) #update fuel system

    Systems.update_y!(components)

end

function Systems.f_step!(components::System{<:Components})
    @unpack aero, ldg, pwp, fuel = components

    f_step!(aero)
    f_step!(ldg)
    f_step!(pwp, is_fuel_available(fuel))

end


#################################### GUI #######################################

function GUI.draw!( components::System{<:Components}, ::System{A},
                    p_open::Ref{Bool} = Ref(true),
                    label::String = "Cessna 172 Components") where {A<:AbstractAvionics}

    @unpack act, pwp, ldg, aero, fuel, pld = components

    CImGui.Begin(label, p_open)

        @cstatic(c_act=false, c_aero=false, c_ldg=false, c_pwp=false, c_fuel=false, c_pld=false,
        begin
            @c CImGui.Checkbox("Actuation", &c_act)
            @c CImGui.Checkbox("Aerodynamics", &c_aero)
            @c CImGui.Checkbox("Landing Gear", &c_ldg)
            @c CImGui.Checkbox("Power Plant", &c_pwp)
            @c CImGui.Checkbox("Fuel", &c_fuel)
            @c CImGui.Checkbox("Payload", &c_pld)
            if c_act
                if A === NoAvionics
                    @c GUI.draw!(act, &c_act)
                else
                    @c GUI.draw(act, &c_act)
                end
            end
            c_aero && @c GUI.draw(aero, &c_aero)
            c_ldg && @c GUI.draw(ldg, &c_ldg)
            c_pwp && @c GUI.draw(pwp, &c_pwp)
            c_fuel && @c GUI.draw(fuel, &c_fuel)
            c_pld && @c GUI.draw!(pld, &c_pld)
        end)

    CImGui.End()

end


################################################################################
################################# Templates ####################################

const Vehicle{F, K, T} = AircraftBase.Vehicle{F, K, T} where {F <: Components, K <: AbstractKinematicDescriptor, T <: AbstractTerrain}

############################### Trimming #######################################
################################################################################

#first 2 are aircraft-agnostic
@kwdef struct TrimState <: AbstractTrimState{7}
    α_a::Float64 = 0.1 #angle of attack, aerodynamic axes
    φ_nb::Float64 = 0.0 #bank angle
    n_eng::Float64 = 0.75 #normalized engine speed (ω/ω_rated)
    throttle::Float64 = 0.47
    aileron::Float64 = 0.014
    elevator::Float64 = -0.0015
    rudder::Float64 = 0.02 #rudder↑ -> aero.u.r↓ -> right yaw
end

@kwdef struct TrimParameters <: AbstractTrimParameters
    Ob::Geographic{NVector, Ellipsoidal} = Geographic(NVector(), HEllip(1050)) #3D location of vehicle frame origin
    ψ_nb::Float64 = 0.0 #geographic heading
    EAS::Float64 = 50.0 #equivalent airspeed
    γ_wb_n::Float64 = 0.0 #wind-relative flight path angle
    ψ_wb_dot::Float64 = 0.0 #WA-relative turn rate
    θ_wb_dot::Float64 = 0.0 #WA-relative pitch rate
    β_a::Float64 = 0.0 #sideslip angle measured in the aerodynamic reference frame
    x_fuel::Ranged{Float64, 0., 1.} = 0.5 #normalized fuel load
    mixture::Ranged{Float64, 0., 1.} = 0.5 #engine mixture control
    flaps::Ranged{Float64, 0., 1.} = 0.0 #flap setting
    payload::C172.PayloadU = C172.PayloadU(m_pilot = 75, m_copilot = 75, m_baggage = 50)
end


function Kinematics.Initializer(trim_state::TrimState,
                                trim_params::TrimParameters,
                                atm_data::AtmData)

    @unpack EAS, β_a, γ_wb_n, ψ_nb, ψ_wb_dot, θ_wb_dot, Ob = trim_params
    @unpack α_a, φ_nb = trim_state

    TAS = Air.EAS2TAS(EAS; ρ = ISAData(Ob, atm_data).ρ)
    v_wb_a = Air.get_velocity_vector(TAS, α_a, β_a)
    v_wb_b = C172.f_ba.q(v_wb_a) #wind-relative aircraft velocity, body frame

    θ_nb = AircraftBase.θ_constraint(; v_wb_b, γ_wb_n, φ_nb)
    e_nb = REuler(ψ_nb, θ_nb, φ_nb)
    q_nb = RQuat(e_nb)

    e_wb = e_nb #initialize WA arbitrarily to NED
    ė_wb = SVector(ψ_wb_dot, θ_wb_dot, 0.0)
    ω_wb_b = Attitude.ω(e_wb, ė_wb)

    loc = NVector(Ob)
    h = HEllip(Ob)

    v_wb_n = q_nb(v_wb_b) #wind-relative aircraft velocity, NED frame
    v_ew_n = atm_data.v_ew_n
    v_eb_n = v_ew_n + v_wb_n

    Kinematics.Initializer(; q_nb, loc, h, ω_wb_b, v_eb_n, Δx = 0.0, Δy = 0.0)

end


function cost(vehicle::System{<:C172.Vehicle})

    @unpack ẋ, y = vehicle

    v_nd_dot = SVector{3}(ẋ.dynamics.v_eb_b) / norm(y.kinematics.data.v_eb_b)
    ω_dot = SVector{3}(ẋ.dynamics.ω_eb_b) #ω should already of order 1
    n_eng_dot = ẋ.components.pwp.engine.ω / vehicle.components.pwp.engine.constants.ω_rated

    sum(v_nd_dot.^2) + sum(ω_dot.^2) + n_eng_dot^2

end

function get_f_target(vehicle::System{<:C172.Vehicle},
                      trim_params::TrimParameters)

    let vehicle = vehicle, trim_params = trim_params
        function (x::TrimState)
            AircraftBase.assign!(vehicle, trim_params, x)
            return cost(vehicle)
        end
    end

end

function Systems.init!(vehicle::System{<:C172.Vehicle}, trim_params::TrimParameters)

    trim_state = TrimState() #could provide initial condition as an optional input

    f_target = get_f_target(vehicle, trim_params)

    #wrapper with the interface expected by NLopt
    f_opt(x::Vector{Float64}, ::Vector{Float64}) = f_target(TrimState(x))

    n = length(trim_state)
    x0 = zeros(n); lower_bounds = similar(x0); upper_bounds = similar(x0); initial_step = similar(x0)

    x0[:] .= trim_state

    lower_bounds[:] .= TrimState(
        α_a = -π/12,
        φ_nb = -π/3,
        n_eng = 0.4,
        throttle = 0,
        aileron = -1,
        elevator = -1,
        rudder = -1)

    upper_bounds[:] .= TrimState(
        α_a = vehicle.components.aero.constants.α_stall[2], #critical AoA is 0.28 < 0.36
        φ_nb = π/3,
        n_eng = 1.1,
        throttle = 1,
        aileron = 1,
        elevator = 1,
        rudder = 1)

    initial_step[:] .= 0.05 #safe value for all optimization variables

    #any of these three algorithms works
    # opt = Opt(:LN_NELDERMEAD, length(x0))
    opt = Opt(:LN_BOBYQA, length(x0))
    # opt = Opt(:GN_CRS2_LM, length(x0))
    opt.min_objective = f_opt
    opt.maxeval = 100000
    opt.stopval = 1e-16
    opt.lower_bounds = lower_bounds
    opt.upper_bounds = upper_bounds
    opt.initial_step = initial_step

    # @btime optimize($opt, $x0)

    (minf, minx, exit_flag) = optimize(opt, x0)

    success = (exit_flag === :STOPVAL_REACHED)
    if !success
        @warn("Trimming optimization failed with exit_flag $exit_flag")
    end
    trim_state_opt = TrimState(minx)
    AircraftBase.assign!(vehicle, trim_params, trim_state_opt)
    return (success = success, trim_state = trim_state_opt)

end

function AircraftBase.trim!(ac::System{<:AircraftBase.Aircraft{<:C172.Vehicle}},
                            trim_params::TrimParameters = TrimParameters())
    Systems.init!(ac, trim_params)
end

function AircraftBase.linearize!(ac::System{<:AircraftBase.Aircraft{<:C172.Vehicle}},
                                trim_params::TrimParameters = TrimParameters())
    linearize!(ac.vehicle, trim_params)
end

################################################################################
############################### C172 Variants ##################################

include(normpath("c172r/c172r.jl")); @reexport using .C172R
include(normpath("c172fbw/c172fbw.jl")); @reexport using .C172FBW
include(normpath("c172rpa/c172rpa.jl")); @reexport using .C172RPA

end