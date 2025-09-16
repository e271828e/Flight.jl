module C172

using LinearAlgebra, StaticArrays, ComponentArrays, UnPack, HDF5, Interpolations
using Logging
using Reexport
using NLopt

using Flight.FlightCore
using Flight.FlightLib

export Cessna172


################################################################################
################################ Airframe #####################################

struct Airframe <: ModelDefinition end

# This component represents the platform's structure, together with any
# components rigidly attached to it, such as the power plant or landing gear,
# but excluding the payload and fuel load. Its mass corresponds roughly to the
# aircraft's Standard Empty Weight

#compute Airframe mass properties in the vehicle reference frame b
const mp_b_afm = let
    #define the airframe as a RigidBodyDistribution; a RigidBodyDistribution is
    #specified in a reference frame c with origin at its center of mass.
    airframe_c = RigidBodyDistribution(767.0, SA[820.0 0 0; 0 1164.0 0; 0 0 1702.0])
    #now, define the transform from vehicle reference frame b to airframe
    #reference frame c (pure translation)
    t_bc = FrameTransform(r = SVector{3}(0.056, 0, 0.582))
    #get the airframe's mass properties to frame b
    MassProperties(airframe_c, t_bc)
end

#the airframe itself receives no external actions nor has any internal angular
#momentum. these are considered separately in the vehicle's aerodynamics, power
#plant and landing gear
Dynamics.get_mp_b(::Model{Airframe}) = mp_b_afm
Dynamics.get_wr_b(::Model{Airframe}) = Wrench()
Dynamics.get_hr_b(::Model{Airframe}) = zeros(SVector{3})


################################################################################
############################ Aerodynamics ######################################


function generate_aero_lookup()

    ################################ C_D data ##################################
    C_D = (
        zero = 0.027,
        δe = (
            x = [-1.0, 0.0, 1.0],
            y = [0.06, 0.0, 0.06]
        ),
        β = (
            x = [-1.0, 0.0, 1.0],
            y = [0.17, 0.0, 0.17]
        ),
        ge = (
            x = [0.0000, 0.1000, 0.1500, 0.2000, 0.3000, 0.4000, 0.5000, 0.6000, 0.7000, 0.8000, 0.9000, 1.0000, 1.1000], #non-dimensional height
            y = [0.4800, 0.5150, 0.6290, 0.7090, 0.8150, 0.8820, 0.9280, 0.9620, 0.9880, 1.0000, 1.0000, 1.0000, 1.0000]
        ),
        δf = (
            x = deg2rad.(range(0, 30, step = 10)),
            y = [0.0000, 0.0070, 0.0120, 0.0180]
        ),
        α_δf = (
            x_α = [ -0.0873, -0.0698, -0.0524, -0.0349, -0.0175, 0.0000, 0.0175, 0.0349, 0.0524, 0.0698, 0.0873, 0.1047, 0.1222, 0.1396, 0.1571, 0.1745, 0.1920, 0.2094, 0.2269, 0.2443, 0.2618, 0.2793, 0.2967, 0.3142, 0.3316, 0.3491],
            x_δf = deg2rad.(range(0, 30, step = 10)),
            y = [
                0.0041 0.0013 0.0001 0.0003 0.0020 0.0052 0.0099 0.0162 0.0240 0.0334 0.0442 0.0566 0.0706 0.0860 0.0962 0.1069 0.1180 0.1298 0.1424 0.1565 0.1727 0.1782 0.1716 0.1618 0.1475 0.1097;
                0.0000 0.0004 0.0023 0.0057 0.0105 0.0168 0.0248 0.0342 0.0452 0.0577 0.0718 0.0874 0.1045 0.1232 0.1353 0.1479 0.1610 0.1746 0.1892 0.2054 0.2240 0.2302 0.2227 0.2115 0.1951 0.1512;
                0.0005 0.0025 0.0059 0.0108 0.0172 0.0251 0.0346 0.0457 0.0583 0.0724 0.0881 0.1053 0.1240 0.1442 0.1573 0.1708 0.1849 0.1995 0.2151 0.2323 0.2521 0.2587 0.2507 0.2388 0.2214 0.1744;
                0.0014 0.0041 0.0084 0.0141 0.0212 0.0299 0.0402 0.0521 0.0655 0.0804 0.0968 0.1148 0.1343 0.1554 0.1690 0.1830 0.1975 0.2126 0.2286 0.2464 0.2667 0.2735 0.2653 0.2531 0.2351 0.1866
            ]'
        )
    )

    ############################# C_Y data #####################################
    C_Y = (
        δr = 0.1870,
        δa = 0.0,
        β_δf = (
            x_β = [-0.3490, 0, 0.3490],
            x_δf = deg2rad.([0, 30]),
            y = [0.1370 0.1060;
                    0.0000 0.0000;
                    -0.1370 -0.1060]
        ),
        p = (
            x_α = [0.0, 0.094],
            x_δf = deg2rad.([0, 30]),
            y = [-0.0750 -0.1610;
                    -0.1450 -0.2310]
        ),
        r = (
            x_α = [0.0, 0.094],
            x_δf = deg2rad.([0, 30]),
            y = [0.2140 0.1620;
                    0.2670 0.2150]
        )
    )

    ############################# C_L data #####################################
    C_L = (
        δe = 0.4300,
        q = 3.900,
        α_dot = 1.700,
        ge = (
            x = [0.0000, 0.1000, 0.1500, 0.2000, 0.3000, 0.4000, 0.5000, 0.6000, 0.7000, 0.8000, 0.9000, 1.0000, 1.1000],
            y = [1.2030, 1.1270, 1.0900, 1.0730, 1.0460, 1.0550, 1.0190, 1.0130, 1.0080, 1.0060, 1.0030, 1.0020, 1.0000]
        ),
        α = (
            x_α = [-0.0900, 0.0000, 0.0900, 0.1000, 0.1200, 0.1400, 0.1600, 0.1700, 0.1900, 0.2100, 0.2400, 0.2600, 0.2800, 0.3000, 0.3200, 0.3400, 0.3600],
            x_stall = [0.0, 1.0],
            y = [
                -0.2200  0.2500  0.7300  0.8300  0.9200  1.0200  1.0800  1.1300  1.1900  1.2500  1.3500  1.4400  1.4700  1.4300  1.3800  1.3000  1.1500;
                -0.2200  0.2500  0.7300  0.7800  0.7900  0.8100  0.8200  0.8300  0.8500  0.8600  0.8800  0.9000  0.9200  0.9500  0.9900  1.0500  1.1500
            ]'
        ),
        δf = (  x = deg2rad.(range(0, 30, step = 10)),
                y = [0.0000, 0.2, 0.3, 0.35])
    )

    ################################## C_l data ################################
    C_l = ( δa = 0.229,
            δr = 0.0147,
            β = -0.09226,
            p = -0.4840,
            r = (   x_α = [0.0, 0.094],
                    x_δf = deg2rad.([0, 30]),
                    y = [0.0798 0.1246;
                        0.1869 0.2317])
            )

    ################################## C_m data ################################
    C_m = ( zero = 0.100,
            δe = -1.1220,
            α = -1.8000,
            q = -12.400,
            α_dot = -7.2700,
            δf = (
                x = deg2rad.([0, 10, 20, 30]),
                y = [0.0000, -0.0654, -0.0981, -0.1140]
            ))

    ################################## C_n data ################################
    C_n = ( δr = -0.0430,
            δa = -0.0053,
            β = 0.05874,
            p = -0.0278,
            r = -0.0937)

    ############################## Interpolations ##############################
    C_D_int = (
        z = C_D.zero,
        β = linear_interpolation(C_D.β.x, C_D.β.y, extrapolation_bc = Flat()),
        δe = linear_interpolation(C_D.δe.x, C_D.δe.y, extrapolation_bc = Flat()),
        δf = linear_interpolation(C_D.δf.x, C_D.δf.y, extrapolation_bc = Flat()),
        α_δf = linear_interpolation((C_D.α_δf.x_α, C_D.α_δf.x_δf), C_D.α_δf.y, extrapolation_bc = Flat()),
        ge = linear_interpolation(C_D.ge.x, C_D.ge.y, extrapolation_bc = Flat())
    )

    C_Y_int = (
        δr = C_Y.δr, δa = C_Y.δa,
        β_δf = linear_interpolation((C_Y.β_δf.x_β, C_Y.β_δf.x_δf), C_Y.β_δf.y, extrapolation_bc = Flat()),
        p = linear_interpolation((C_Y.p.x_α, C_Y.p.x_δf), C_Y.p.y, extrapolation_bc = Flat()),
        r = linear_interpolation((C_Y.r.x_α, C_Y.r.x_δf), C_Y.r.y, extrapolation_bc = Flat())
    )

    C_L_int = (
        δe = C_L.δe, q = C_L.q, α_dot = C_L.α_dot,
        α = linear_interpolation((C_L.α.x_α, C_L.α.x_stall), C_L.α.y, extrapolation_bc = Flat()),
        δf = linear_interpolation(C_L.δf.x, C_L.δf.y, extrapolation_bc = Flat()),
        ge = linear_interpolation(C_L.ge.x, C_L.ge.y, extrapolation_bc = Flat())
    )

    C_l_int = (
        δa = C_l.δa, δr = C_l.δr, β = C_l.β, p = C_l.p,
        r = linear_interpolation((C_l.r.x_α, C_l.r.x_δf), C_l.r.y, extrapolation_bc = Flat())
    )

    C_m_int = (
        z = C_m.zero, δe = C_m.δe, α = C_m.α, q = C_m.q, α_dot = C_m.α_dot,
        δf = linear_interpolation(C_m.δf.x, C_m.δf.y, extrapolation_bc = Flat())
    )

    C_n_int = ( δr = C_n.δr, δa = C_n.δa, β = C_n.β, p = C_n.p, r = C_n.r)

    return (C_D = C_D_int, C_Y = C_Y_int, C_L = C_L_int, C_l = C_l_int, C_m = C_m_int, C_n = C_n_int)

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

@kwdef struct Aero <: ModelDefinition
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

@kwdef mutable struct AeroS
    stall::Bool = false
end

@kwdef struct AeroY
    e::Float64 = 0.0 #normalized elevator deflection
    a::Float64 = 0.0 #normalized aileron deflection
    r::Float64 = 0.0 #normalized rudder deflection
    f::Float64 = 0.0 #normalized flap deflection
    δe::Float64 = 0.0 #elevator deflection (rad)
    δa::Float64 = 0.0 #aileron deflection (rad)
    δr::Float64 = 0.0 #rudder deflection (rad)
    δf::Float64 = 0.0 #flap deflection (rad)
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

Modeling.X(::Aero) = ComponentVector(α_filt = 0.0, β_filt = 0.0) #filtered airflow angles
Modeling.Y(::Aero) = AeroY()
Modeling.U(::Aero) = AeroU()
Modeling.S(::Aero) = AeroS()

#*caution: do not confuse the w-frame in kinematics.ω_wb_b, which refers to the
#wander-azimuth frame (w), with the w in air.v_wb_b, which indicates aerodynamic
#(wind-relative) velocity

function Modeling.f_ode!(mdl::Model{Aero}, terrain::Model{<:AbstractTerrain},
                        air_data::AirData, kin_data::KinData)

    @unpack ẋ, x, u, s, constants = mdl
    @unpack α_filt, β_filt = x
    @unpack e, a, r, f = u
    @unpack S, b, c, δe_range, δa_range, δr_range, δf_range, α_stall, V_min, τ = constants
    @unpack TAS, q, v_wb_b = air_data
    @unpack ω_wb_b, n_e, h_o = kin_data
    stall = s.stall

    v_wb_a = f_ba.q'(v_wb_b)

    #for near-zero TAS, the airflow angles are likely to chatter between 0, -π
    #and π. avoid this by setting set a minimum TAS for airflow computation. in
    #this scenario, dynamic pressure would be close to zero, so forces and
    #moments would vanish anyway.
    α, β = (TAS > 0.1 ? Atmosphere.get_airflow_angles(v_wb_a) : (0.0, 0.0))
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
    h_trn_a = TerrainData(terrain, loc_Oa).elevation #orthometric
    Δh_nd = (h_a - h_trn_a) / b

    # T = get_wr_b(pwp).F[1]
    # C_T = T / (q * S) #thrust coefficient, not used here

    coeffs = get_aero_coeffs(;
        α, β, p_nd, q_nd, r_nd, δa, δr, δe, δf, α_dot_nd, β_dot_nd, Δh_nd, stall)

    @unpack C_D, C_Y, C_L, C_l, C_m, C_n = coeffs

    q_as = Atmosphere.get_stability_axes(α)
    F_aero_s = q * S * SVector{3,Float64}(-C_D, C_Y, -C_L)
    F_aero_a = q_as(F_aero_s)
    τ_aero_a = q * S * SVector{3,Float64}(C_l * b, C_m * c, C_n * b)

    # wr_b = wr_a = Wrench(F_aero_a, τ_aero_a)
    wr_b = Wrench(F_aero_a, τ_aero_a)

    ẋ.α_filt = α_filt_dot
    ẋ.β_filt = β_filt_dot

    mdl.y = AeroY(; e, a, r, f, δe, δa, δr, δf,
                    α, α_filt, α_filt_dot, β, β_filt, β_filt_dot,
                    stall, coeffs, wr_b)

    return nothing

end

function Modeling.f_step!(mdl::Model{Aero})
    #stall hysteresis
    α = mdl.y.α
    α_stall = mdl.α_stall
    if α > α_stall[2]
        mdl.s.stall = true
    elseif α < α_stall[1]
        mdl.s.stall = false
    end
end

Dynamics.get_mp_b(::Model{<:Aero}) = MassProperties()
Dynamics.get_hr_b(::Model{<:Aero}) = zeros(SVector{3})
Dynamics.get_wr_b(mdl::Model{Aero}) = mdl.y.wr_b


################################# GUI ##########################################


function GUI.draw(mdl::Model{<:Aero}, p_open::Ref{Bool} = Ref(true),
                window_label::String = "Cessna 172S Aerodynamics")

    @unpack e, a, r, f, α, β, α_filt, β_filt, stall, coeffs, wr_b = mdl.y
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

struct Ldg <: ModelDefinition
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

Dynamics.get_mp_b(::Model{Ldg}) = MassProperties()
Dynamics.get_hr_b(::Model{Ldg}) = zeros(SVector{3})

#delegate continuous dynamics to submodels
@sm_ode Ldg

#skip output update
Modeling.f_step!(mdl::Model{<:Ldg}) = foreach(f_step!, mdl.submodels)

#approximate height of the aircraft frame origin to the ground
const Δh_to_gnd = 1.81

function GUI.draw(mdl::Model{<:Ldg}, p_open::Ref{Bool} = Ref(true),
                window_label::String = "Cessna 172S Landing Gear")

    @unpack left, right, nose = mdl

    CImGui.Begin(window_label, p_open)

        @cstatic(c_left=false, c_right=false, c_nose = false, begin
            @c CImGui.Checkbox("Left Main", &c_left)
            @c CImGui.Checkbox("Right Main", &c_right)
            @c CImGui.Checkbox("Nose", &c_nose)
            c_left && @c GUI.draw(left, &c_left, "Left Landing Gear Unit")
            c_right && @c GUI.draw(right, &c_right, "Right Landing Gear Unit")
            c_nose && @c GUI.draw(nose, &c_nose, "Nose Landing Gear Unit")
        end)

    CImGui.End()

end

################################################################################
################################# Payload ######################################

@kwdef struct Payload <: ModelDefinition
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

@kwdef struct PayloadY
    m_pilot::Ranged{Float64, 0., 100.} = 75.0
    m_copilot::Ranged{Float64, 0., 100.} = 75.0
    m_lpass::Ranged{Float64, 0., 100.} = 0.0
    m_rpass::Ranged{Float64, 0., 100.} = 0.0
    m_baggage::Ranged{Float64, 0., 100.} = 50.0
end

Modeling.U(::Payload) = PayloadU()
Modeling.Y(::Payload) = PayloadY()

function Dynamics.get_mp_b(mdl::Model{Payload})
    @unpack m_pilot, m_copilot, m_lpass, m_rpass, m_baggage = mdl.u
    @unpack pilot_slot, copilot_slot, lpass_slot, rpass_slot, baggage_slot = mdl.constants

    pilot = MassProperties(PointDistribution(m_pilot), pilot_slot)
    copilot = MassProperties(PointDistribution(m_copilot), copilot_slot)
    lpass = MassProperties(PointDistribution(m_lpass), lpass_slot)
    rpass = MassProperties(PointDistribution(m_rpass), rpass_slot)
    baggage = MassProperties(PointDistribution(m_baggage), baggage_slot)

    mp_b = MassProperties() + pilot + copilot + lpass + rpass + baggage
    return mp_b
end

Dynamics.get_wr_b(::Model{Payload}) = Wrench()
Dynamics.get_hr_b(::Model{Payload}) = zeros(SVector{3})

#################################### GUI #######################################

function GUI.draw!(mdl::Model{<:Payload},
                    p_open::Ref{Bool} = Ref(true),
                    label::String = "Cessna 172S Payload")

    u = mdl.u

    CImGui.Begin(label, p_open)

    CImGui.PushItemWidth(-250)

    u.m_pilot = GUI.safe_slider("Pilot Mass (kg)", u.m_pilot, "%.3f")
    u.m_copilot = GUI.safe_slider("Copilot Mass (kg)", u.m_copilot, "%.3f")
    u.m_lpass = GUI.safe_slider("Left Passenger Mass (kg)", u.m_lpass, "%.3f")
    u.m_rpass = GUI.safe_slider("Right Passenger Mass (kg)", u.m_rpass, "%.3f")
    u.m_baggage = GUI.safe_slider("Baggage Mass (kg)", u.m_baggage, "%.3f")

    CImGui.PopItemWidth()

    CImGui.End()

end


################################################################################
################################# Fuel #########################################

#assumes fuel is drawn equally from both tanks, no need to model them
#individually for now
@kwdef struct Fuel <: ModelDefinition
    m_full::Float64 = 114.4 #maximum fuel mass (42 gal * 6 lb/gal * 0.454 kg/lb)
    m_res::Float64 = 1.0 #residual fuel mass
end

@kwdef struct FuelY
    x_avail::Float64 = 0.0 #available fuel fraction
    m_total::Float64 = 0.0 #total fuel mass
    m_avail::Float64 = 0.0 #available fuel mass
end

#normalized fuel content (0: residual, 1: full)
Modeling.X(::Fuel) = [0.5] #cannot be a scalar, need an AbstractVector{<:Real}
Modeling.Y(::Fuel) = FuelY()

function Modeling.f_ode!(mdl::Model{Fuel}, pwp::Model{<:PistonThruster})

    @unpack m_full, m_res = mdl.constants #no need for submodels
    x_avail = mdl.x[1]
    m_total = m_res + x_avail * (m_full - m_res) #current mass
    m_avail = m_total - m_res
    mdl.ẋ .= -pwp.y.engine.ṁ / (m_full - m_res)
    mdl.y = FuelY(; x_avail, m_total, m_avail)

end

function Dynamics.get_mp_b(fuel::Model{Fuel})

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

Dynamics.get_wr_b(::Model{Fuel}) = Wrench()
Dynamics.get_hr_b(::Model{Fuel}) = zeros(SVector{3})

is_fuel_available(mdl::Model{<:Fuel}) = (mdl.y.m_avail > 0)

function GUI.draw(mdl::Model{Fuel}, p_open::Ref{Bool} = Ref(true),
                window_label::String = "Cessna 172S Fuel System")

    @unpack m_total, m_avail = mdl.y

    CImGui.Begin(window_label, p_open)

        CImGui.Text(@sprintf("Total Fuel: %.6f kg", m_total))
        CImGui.Text(@sprintf("Available Fuel: %.6f kg", m_avail))

    CImGui.End()

end


################################################################################
############################## AbstractActuation ###############################

#any Actuation system suitable for the C172 platform
abstract type AbstractActuation <: ModelDefinition end

function assign!(aero::Model{<:Aero}, ldg::Model{<:Ldg},
                pwp::Model{<:PistonThruster}, act::Model{<:AbstractActuation})
    throw(MethodError(C172.assign!, (aero, ldg, pwp, act)))
end

Dynamics.get_mp_b(::Model{<:AbstractActuation}) = MassProperties()
Dynamics.get_wr_b(::Model{<:AbstractActuation}) = Wrench()
Dynamics.get_hr_b(::Model{<:AbstractActuation}) = zeros(SVector{3})



################################################################################
################################ Systems ######################################

struct Systems{P <: PistonThruster, A <: AbstractActuation} <: AbstractVehicleSystems
    afm::Airframe
    aero::Aero
    ldg::Ldg
    fuel::Fuel
    pld::Payload
    pwp::P
    act::A
end

function Systems(pwp::PistonThruster, act::AbstractActuation)
    Systems(Airframe(), Aero(), Ldg(), Fuel(), Payload(), pwp, act)
end


############################# Update Methods ###################################

function Modeling.f_ode!(systems::Model{<:Systems},
                        terrain::Model{<:AbstractTerrain},
                        kin_data::KinData,
                        air_data::AirData)

    @unpack act, aero, pwp, ldg, fuel, pld = systems

    f_ode!(act) #update actuation system outputs
    assign!(aero, ldg, pwp, act) #assign actuation model outputs to systems submodels
    f_ode!(aero, terrain, air_data, kin_data) #update aerodynamics continuous state & outputs
    f_ode!(ldg, terrain, kin_data) #update landing gear continuous state & outputs
    f_ode!(pwp, air_data, kin_data) #update powerplant continuous state & outputs
    f_ode!(fuel, pwp) #update fuel system

    f_output!(systems)

end

function Modeling.f_step!(systems::Model{<:Systems},
                        ::Model{<:AbstractAtmosphere},
                        ::Model{<:AbstractTerrain})
    @unpack aero, ldg, pwp, fuel = systems

    f_step!(aero)
    f_step!(ldg)
    f_step!(pwp, is_fuel_available(fuel))

end

@no_periodic Systems



#################################### GUI #######################################

function GUI.draw!( systems::Model{<:Systems}, ::Model{<:AbstractAvionics},
                    p_open::Ref{Bool} = Ref(true),
                    label::String = "Systems")

    @unpack act, pwp, ldg, aero, fuel, pld = systems

    CImGui.Begin(label, p_open)

        @cstatic(c_act=false, c_aero=false, c_ldg=false, c_pwp=false, c_fuel=false, c_pld=false,
        begin
            @c CImGui.Checkbox("Actuation", &c_act)
            @c CImGui.Checkbox("Aerodynamics", &c_aero)
            @c CImGui.Checkbox("Landing Gear", &c_ldg)
            @c CImGui.Checkbox("Power Plant", &c_pwp)
            @c CImGui.Checkbox("Fuel", &c_fuel)
            @c CImGui.Checkbox("Payload", &c_pld)
            c_act && @c GUI.draw!(act, &c_act)
            c_aero && @c GUI.draw(aero, &c_aero)
            c_ldg && @c GUI.draw(ldg, &c_ldg)
            c_pwp && @c GUI.draw(pwp, &c_pwp, "Power Plant")
            c_fuel && @c GUI.draw(fuel, &c_fuel)
            c_pld && @c GUI.draw!(pld, &c_pld)
        end)

    CImGui.End()

end


################################################################################
################################# Templates ####################################

const Vehicle = AircraftBase.Vehicle{<:C172.Systems}
const Aircraft = AircraftBase.Aircraft{<:C172.Vehicle}
const Cessna172 = C172.Aircraft

########################## Explicit Initialization #############################
################################################################################

@kwdef struct SystemsInitializer <: AbstractVehicleSystemsInitializer
    engine_state::Piston.EngineStateEnum = Piston.EngineState.off
    mixture_ctl::Piston.MixtureControlEnum = Piston.MixtureControl.auto
    n_eng::Float64 = 0.0 #normalized engine speed
    throttle::Ranged{Float64, 0., 1.} = 0
    mixture::Ranged{Float64, 0., 1.} = 0.5
    elevator::Ranged{Float64, -1., 1.} = 0
    aileron::Ranged{Float64, -1., 1.} = 0
    rudder::Ranged{Float64, -1., 1.} = 0
    flaps::Ranged{Float64, 0., 1.} = 0
    brake_left::Ranged{Float64, 0., 1.} = 0
    brake_right::Ranged{Float64, 0., 1.} = 0
    fuel_load::Ranged{Float64, 0., 1.} = 0.5 #normalized
    payload::C172.PayloadY = C172.PayloadY()
    stall::Bool = false
    α_a_filt::Float64 = 0 #only needed for trim assignments
    β_a_filt::Float64 = 0 #only needed for trim assignments
end

function Init(kin::KinInit = KinInit(), sys::SystemsInitializer = SystemsInitializer())
    AircraftBase.VehicleInitializer(kin, sys)
end


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
    fuel_load::Ranged{Float64, 0., 1.} = 0.5 #normalized fuel load
    mixture::Ranged{Float64, 0., 1.} = 0.5 #engine mixture control
    flaps::Ranged{Float64, 0., 1.} = 0.0 #flap setting
    payload::PayloadY = PayloadY()
end

function Atmosphere.AtmosphericData(mdl::Model{<:AbstractAtmosphere},
                                    trim_params::TrimParameters)
    AtmosphericData(mdl, trim_params.Ob)
end

function Kinematics.Initializer(trim_state::TrimState,
                                trim_params::TrimParameters,
                                atmosphere::Model{<:AbstractAtmosphere})

    @unpack EAS, β_a, γ_wb_n, ψ_nb, ψ_wb_dot, θ_wb_dot, Ob = trim_params
    @unpack α_a, φ_nb = trim_state

    atm_data = AtmosphericData(atmosphere, Ob)
    TAS = Atmosphere.EAS2TAS(EAS; ρ = atm_data.ρ)
    v_wb_a = Atmosphere.get_velocity_vector(TAS, α_a, β_a)
    v_wb_b = C172.f_ba.q(v_wb_a) #wind-relative aircraft velocity, body frame

    θ_nb = AircraftBase.θ_constraint(; v_wb_b, γ_wb_n, φ_nb)
    e_nb = REuler(ψ_nb, θ_nb, φ_nb)
    q_nb = RQuat(e_nb)

    e_wb = e_nb #initialize WA arbitrarily to NED
    ė_wb = SVector(ψ_wb_dot, θ_wb_dot, 0.0)
    ω_wb_b = Attitude.ω(e_wb, ė_wb)

    location = NVector(Ob)
    h = HEllip(Ob)

    v_ew_n = atm_data.v
    v_wb_n = q_nb(v_wb_b) #wind-relative aircraft velocity, NED frame
    v_eb_n = v_ew_n + v_wb_n

    Kinematics.Initializer(; q_nb, location, h, ω_wb_b, v_eb_n)

end


function cost(vehicle::Model{<:C172.Vehicle})

    @unpack ẋ, y = vehicle

    v_nd_dot = SVector{3}(ẋ.dynamics.v_eb_b) / norm(y.kinematics.v_eb_b)
    ω_dot = SVector{3}(ẋ.dynamics.ω_eb_b) #ω should already of order 1
    n_eng_dot = ẋ.systems.pwp.engine.ω / vehicle.systems.pwp.engine.ω_rated

    sum(v_nd_dot.^2) + sum(ω_dot.^2) + n_eng_dot^2

end

function get_f_target(vehicle::Model{<:C172.Vehicle},
                      trim_params::TrimParameters,
                      atmosphere::Model{<:AbstractAtmosphere},
                      terrain::Model{<:AbstractTerrain})

    let vehicle = vehicle, trim_params = trim_params
        function (x::TrimState)
            AircraftBase.assign!(vehicle, trim_params, x, atmosphere, terrain)
            return cost(vehicle)
        end
    end

end

function Modeling.init!(
            vehicle::Model{<:C172.Vehicle},
            trim_params::TrimParameters,
            atmosphere::Model{<:AbstractAtmosphere} = Model(SimpleAtmosphere()),
            terrain::Model{<:AbstractTerrain} = Model(HorizontalTerrain()))

    trim_state = TrimState() #could provide initial condition as an optional input

    f_target = get_f_target(vehicle, trim_params, atmosphere, terrain)

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
        α_a = vehicle.systems.aero.α_stall[2], #critical AoA is 0.28 < 0.36
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
    AircraftBase.assign!(vehicle, trim_params, trim_state_opt, atmosphere, terrain)
    return (success = success, trim_state = trim_state_opt)

end

function AircraftBase.trim!( aircraft::Model{<:Cessna172},
                            params::TrimParameters = TrimParameters(), args...)
    Modeling.init!(aircraft, params, args...)
end

function AircraftBase.linearize!( aircraft::Model{<:Cessna172},
                            params::TrimParameters = TrimParameters(), args...)
    AircraftBase.linearize!(aircraft.vehicle, params, args...)
end

################################################################################
############################### XPlane12Control ###################################

function IODevices.extract_output(aircraft::Model{<:Cessna172}, ::XPlane12ControlMapping)

    t = aircraft.t[]
    @unpack δe, δa, δr, δf = aircraft.y.vehicle.systems.aero
    ψ_sw = aircraft.y.vehicle.systems.ldg.nose.strut.ψ_sw
    ω_prop = aircraft.y.vehicle.systems.pwp.propeller.ω

    ϕ_prop = mod(ω_prop * t, 2π)
    prop_is_disc = (ω_prop > 10 ? true : false)

    drefs = (
        elev_left_pos = "sim/flightmodel2/wing/elevator1_deg[8]",
        elev_right_pos = "sim/flightmodel2/wing/elevator1_deg[9]",
        flap_left_pos = "sim/flightmodel2/wing/flap1_deg[0]",
        flap_right_pos = "sim/flightmodel2/wing/flap1_deg[1]",
        rudder_pos = "sim/flightmodel2/wing/rudder1_deg[10]",
        ail_left_pos = "sim/flightmodel2/wing/aileron1_deg[2]", #not [0]!
        ail_right_pos = "sim/flightmodel2/wing/aileron1_deg[3]", #not [1]!
        prop_is_disc = "sim/flightmodel2/engines/prop_is_disc[0]",
        prop_angle = "sim/flightmodel2/engines/prop_rotation_angle_deg[0]",
        nws_angle = "sim/flightmodel2/gear/tire_steer_actual_deg[0]")

    msgs = (
        Network.xpmsg_set_dref(drefs.elev_left_pos, rad2deg(δe)),
        Network.xpmsg_set_dref(drefs.elev_right_pos, rad2deg(δe)),
        Network.xpmsg_set_dref(drefs.ail_left_pos, rad2deg(δa)),
        Network.xpmsg_set_dref(drefs.ail_right_pos, rad2deg(-δa)),
        Network.xpmsg_set_dref(drefs.flap_left_pos, rad2deg(δf)),
        Network.xpmsg_set_dref(drefs.flap_right_pos, rad2deg(δf)),
        Network.xpmsg_set_dref(drefs.rudder_pos, rad2deg(δr)),
        Network.xpmsg_set_dref(drefs.prop_is_disc, prop_is_disc),
        Network.xpmsg_set_dref(drefs.prop_angle, rad2deg(ϕ_prop)),
        Network.xpmsg_set_dref(drefs.nws_angle, rad2deg(ψ_sw)),
        Network.xpmsg_set_pose(XPlanePose(KinData(aircraft)))
    )

    return msgs

end


################################################################################
############################# Utility Functions ################################

function is_on_gnd(mdl::Model{<:Vehicle}; all::Bool = false)
    wow = SVector{3}(leg.strut.wow for leg in mdl.y.systems.ldg)
    all ? all(wow) : any(wow)
end

################################################################################
############################### C172 Variants ##################################

include(normpath("c172s/c172s.jl")); @reexport using .C172S
include(normpath("c172x/c172x.jl")); @reexport using .C172X

end