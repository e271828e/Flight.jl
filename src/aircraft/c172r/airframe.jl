module C172RAirframe

using StaticArrays
using ComponentArrays
using UnPack
using Printf
using CImGui, CImGui.CSyntax, CImGui.CSyntax.CStatic

using Flight.FlightCore
using Flight.FlightPhysics

using Flight.FlightAircraft.LandingGear
using Flight.FlightAircraft.Propellers
using Flight.FlightAircraft.Piston
using Flight.FlightAircraft.Aircraft

export Airframe

################################################################################
############################ Airframe Subsystems ################################

################################ Structure #####################################

struct Structure <: Component end

# This component represents the airframe structure, together with any components
# rigidly attached to it, such as powerplant or landing gear, but not payload or
# fuel contents. Its mass corresponds roughly to the aircraft's Standard Empty
# Weight

#Structure mass properties computed in the vehicle reference frame b
const mp_Ob_str = let
    #define the structure as a RigidBodyDistribution
    str_G = RigidBodyDistribution(767.0, SA[820.0 0 0; 0 1164.0 0; 0 0 1702.0])
    #define the transform from the origin of the vehicle reference frame (Ob)
    #to the structure's center of mass (G)
    t_Ob_G = FrameTransform(r = SVector{3}(0.056, 0, 0.582))
    #compute the structure's mass properties at Ob
    MassProperties(str_G, t_Ob_G)
end

RigidBody.MassTrait(::System{Structure}) = HasMass()

#the structure itself receives no external actions. these are considered to act
#upon the vehicle's aerodynamics, power plant and landing gear. the same goes
#for rotational angular momentum.
RigidBody.WrenchTrait(::System{Structure}) = GetsNoExternalWrench()
RigidBody.AngMomTrait(::System{Structure}) = HasNoAngularMomentum()

RigidBody.get_mp_Ob(::System{Structure}) = mp_Ob_str


################################################################################
############################ Aerodynamics ######################################

include("data/aero.jl")

#the aircraft body reference frame fb is arbitrarily chosen to coincide with
#the aerodynamics frame fa, so the frame transform is trivial
const aero_lookup = generate_aero_lookup()
const f_ba = FrameTransform()

# if this weren't the case, and we cared not only about the rotation but also
#about the velocity lever arm, here's the rigorous way of computing v_wOa_a:
# v_wOb_b = v_eOb_b - v_ew_b
# v_eOa_b = v_eOb_b + ω_eb_b × r_ObOa_b
# v_wOa_b = v_eOa_b - v_ew_b = v_eOb_b + ω_eb_b × r_ObOa_b - v_ew_b
# v_wOa_b = v_wOb_b + ω_eb_b × r_ObOa_b
# v_wOa_a = q_ba'(v_wOa_b)

Base.@kwdef struct AeroCoeffs
    C_D::Float64 = 0.0
    C_Y::Float64 = 0.0
    C_L::Float64 = 0.0
    C_l::Float64 = 0.0
    C_m::Float64 = 0.0
    C_n::Float64 = 0.0
end

function get_aero_coeffs(; α, β, p_nd, q_nd, r_nd, δa, δr, δe, δf, α_dot_nd, β_dot_nd, Δh_nd, stall)

    #set sensible bounds
    α = clamp(α, -0.1, 0.36) #0.36 is the highest value (post-stall) tabulated for C_L
    β = clamp(β, -0.2, 0.2)
    α_dot_nd = clamp(α_dot_nd, -0.04, 0.04)
    β_dot_nd = clamp(β_dot_nd, -0.2, 0.2)

    @unpack C_D, C_Y, C_L, C_l, C_m, C_n = aero_lookup

    AeroCoeffs(
        C_D = C_D.z + C_D.ge(Δh_nd) * (C_D.α_δf(α,δf) + C_D.δf(δf)) + C_D.δe(δe) + C_D.β(β),
        C_Y = C_Y.δr * δr + C_Y.δa * δa + C_Y.β_δf(β,δf) + C_Y.p(α,δf) * p_nd + C_Y.r(α,δf) * r_nd,
        C_L = C_L.ge(Δh_nd) * (C_L.α(α,stall) + C_L.δf(δf)) + C_L.δe * δe + C_L.q * q_nd + C_L.α_dot * α_dot_nd,
        C_l = C_l.δa * δa + C_l.δr * δr + C_l.β * β + C_l.p * p_nd + C_l.r(α,δf) * r_nd,
        C_m = C_m.z + C_m.δe * δe + C_m.δf(δf) + C_m.α * α + C_m.q * q_nd + C_m.α_dot * α_dot_nd,
        C_n = C_n.δr * δr + C_n.δa * δa + C_n.β * β + C_n.p * p_nd + C_n.r * r_nd,
    )

end

Base.@kwdef struct Aero <: Component
    S::Float64 = 16.165 #wing area
    b::Float64 = 10.912 #wingspan
    c::Float64 = 1.494 #mean aerodynamic chord
    δe_range::NTuple{2,Float64} = deg2rad.((-28, 23)) #elevator deflection range (rad)
    δa_range::NTuple{2,Float64} = deg2rad.((-20, 20)) #aileron deflection range (rad)
    δr_range::NTuple{2,Float64} = deg2rad.((-16, 16)) #rudder deflection range (rad)
    δf_range::NTuple{2,Float64} = deg2rad.((0, 30)) #flap deflection range (rad)
    α_stall::NTuple{2,Float64} = (0.09, 0.36) #α values for stall hysteresis switching
    V_min::Float64 = 1.0 #lower airspeed threshold for non-dimensional angle rates
    τ::Float64 = 0.05 #time constant for filtered airflow angle derivatives
end

#e↑ -> δe↑ -> trailing edge down -> Cm↓ -> pitch down
#a↑ -> δa↑ -> left trailing edge down, right up -> Cl↑ -> roll right
#r↑ -> δr↑ -> rudder trailing edge left -> Cn↓ -> yaw left
#f↑ -> δf↑ -> flap trailing edge down -> CL↑

Base.@kwdef mutable struct AeroU
    e::Ranged{Float64, -1, 1} = 0.0
    a::Ranged{Float64, -1, 1} = 0.0
    r::Ranged{Float64, -1, 1} = 0.0
    f::Ranged{Float64, 0, 1} = 0.0
end

Base.@kwdef mutable struct AeroS #discrete state
    stall::Bool = false
end

Base.@kwdef struct AeroY
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

Systems.init(::SystemX, ::Aero) = ComponentVector(α_filt = 0.0, β_filt = 0.0) #filtered airflow angles
Systems.init(::SystemY, ::Aero) = AeroY()
Systems.init(::SystemU, ::Aero) = AeroU()
Systems.init(::SystemS, ::Aero) = AeroS()

RigidBody.MassTrait(::System{<:Aero}) = HasNoMass()
RigidBody.AngMomTrait(::System{<:Aero}) = HasNoAngularMomentum()
RigidBody.WrenchTrait(::System{<:Aero}) = GetsExternalWrench()

function Systems.f_ode!(sys::System{Aero}, ::System{<:Piston.Thruster},
    air::AirData, kinematics::KinematicData, terrain::System{<:AbstractTerrain})

    @unpack ẋ, x, u, s, params = sys
    @unpack α_filt, β_filt = x
    @unpack e, a, r, f = u
    @unpack S, b, c, δe_range, δa_range, δr_range, δf_range, α_stall, V_min, τ = params
    @unpack TAS, q, v_wOb_b = air
    @unpack ω_lb_b, n_e, h_o = kinematics
    stall = s.stall

    v_wOb_a = f_ba.q'(v_wOb_b)

    #for near-zero TAS, the airflow angles are likely to chatter between 0, -π
    #and π. to avoid this, we set a minimum TAS for airflow computation. in this
    #scenario dynamic pressure will be close to zero, so forces and moments will
    #vanish anyway.
    α, β = (TAS > 0.1 ? get_airflow_angles(v_wOb_a) : (0.0, 0.0))
    V = max(TAS, V_min) #avoid division by zero

    α_filt_dot = 1/τ * (α - α_filt)
    β_filt_dot = 1/τ * (β - β_filt)

    p_nd = ω_lb_b[1] * b / (2V) #non-dimensional roll rate
    q_nd = ω_lb_b[2] * c / (2V) #non-dimensional pitch rate
    r_nd = ω_lb_b[3] * b / (2V) #non-dimensional yaw rate

    α_dot_nd = α_filt_dot * c / (2V)
    β_dot_nd = β_filt_dot * b / (2V)

    δe = linear_scaling(e, δe_range)
    δa = linear_scaling(a, δa_range)
    δr = linear_scaling(r, δr_range)
    δf = linear_scaling(f, δf_range)

    #non-dimensional height above ground
    l2d_Oa = n_e #Oa = Ob
    h_Oa = h_o #orthometric
    h_trn_Oa = TerrainData(terrain, l2d_Oa).altitude #orthometric
    Δh_nd = (h_Oa - h_trn_Oa) / b

    # T = get_wr_b(pwp).F[1]
    # C_T = T / (q * S) #thrust coefficient, not used here

    coeffs = get_aero_coeffs(;
        α, β, p_nd, q_nd, r_nd, δa, δr, δe, δf, α_dot_nd, β_dot_nd, Δh_nd, stall)

    @unpack C_D, C_Y, C_L, C_l, C_m, C_n = coeffs

    q_as = get_stability_axes(α)
    F_aero_s = q * S * SVector{3,Float64}(-C_D, C_Y, -C_L)
    F_aero_a = q_as(F_aero_s)
    M_aero_a = q * S * SVector{3,Float64}(C_l * b, C_m * c, C_n * b)

    # wr_b = wr_a = Wrench(F_aero_a, M_aero_a)
    wr_b = Wrench(F_aero_a, M_aero_a)

    ẋ.α_filt = α_filt_dot
    ẋ.β_filt = β_filt_dot

    sys.y = AeroY(; α, α_filt, α_filt_dot, β, β_filt, β_filt_dot,
        e, a, r, f, stall, coeffs, wr_b)

    return nothing

end

RigidBody.get_wr_b(sys::System{Aero}) = sys.y.wr_b

function Systems.f_step!(sys::System{Aero})
    #stall hysteresis
    α = sys.y.α
    α_stall = sys.params.α_stall
    if α > α_stall[2]
        sys.s.stall = true
    elseif α < α_stall[1]
        sys.s.stall = false
    end
    return false
end

# # splt_α = thplot(t, rad2deg.(α_b);
# #     title = "Angle of Attack", ylabel = L"$\alpha \ (deg)$",
# #     label = "", kwargs...)

# # splt_β = thplot(t, rad2deg.(β_b);
# #     title = "Angle of Sideslip", ylabel = L"$\beta \ (deg)$",
# #     label = "", kwargs...)

# # pd["05_α_β"] = plot(splt_α, splt_β;
# #     plot_title = "Airflow Angles [Airframe]",
# #     layout = (1,2),
# #     kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs


################################# GUI ##########################################


function GUI.draw(sys::System{<:Aero}, window_label::String = "C172R Aerodynamics")

    @unpack e, a, r, f, α, β, α_filt, β_filt, stall, coeffs, wr_b = sys.y
    @unpack C_D, C_Y, C_L, C_l, C_m, C_n = coeffs

    CImGui.Begin(window_label)

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
        GUI.draw(wr_b.M, "Aerodynamic Torque (O) [Body]", "N*m")

    CImGui.End()

end


############################# Landing Gear ####################################

struct Ldg <: Component
    left::LandingGearUnit{NoSteering, DirectBraking, Strut{SimpleDamper}}
    right::LandingGearUnit{NoSteering, DirectBraking, Strut{SimpleDamper}}
    nose::LandingGearUnit{DirectSteering, NoBraking, Strut{SimpleDamper}}
end

RigidBody.MassTrait(::System{Ldg}) = HasNoMass()
RigidBody.WrenchTrait(::System{Ldg}) = GetsExternalWrench()
RigidBody.AngMomTrait(::System{Ldg}) = HasNoAngularMomentum()

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

function GUI.draw(sys::System{<:Ldg}, window_label::String = "Cessna 172R Landing Gear")

    @unpack left, right, nose = sys

    CImGui.Begin(window_label)

        show_left = @cstatic check=false @c CImGui.Checkbox("Left Main", &check)
        show_right = @cstatic check=false @c CImGui.Checkbox("Right Main", &check)
        show_nose = @cstatic check=false @c CImGui.Checkbox("Nose", &check)

    CImGui.End()

    show_left && GUI.draw(left, "Left Main")
    show_right && GUI.draw(right, "Right Main")
    show_nose && GUI.draw(nose, "Nose")

end


################################################################################
################################# Payload ######################################

Base.@kwdef struct Payload <: Component
    pilot_slot::FrameTransform = FrameTransform(r = SVector{3}(0.183, -0.356, 0.899))
    copilot_slot::FrameTransform = FrameTransform(r = SVector{3}(0.183, 0.356, 0.899))
    lpass_slot::FrameTransform = FrameTransform(r = SVector{3}(-0.681, -0.356, 0.899))
    rpass_slot::FrameTransform = FrameTransform(r = SVector{3}(-0.681, 0.356, 0.899))
    baggage_slot::FrameTransform = FrameTransform(r = SVector{3}(-1.316, 0, 0.899))
end

Base.@kwdef mutable struct PayloadU
    m_pilot::Ranged{Float64, 0, 100} = 75.0
    m_copilot::Ranged{Float64, 0, 100} = 75.0
    m_lpass::Ranged{Float64, 0, 100} = 0.0
    m_rpass::Ranged{Float64, 0, 100} = 0.0
    m_baggage::Ranged{Float64, 0, 100} = 50.0
end

Systems.init(::SystemU, ::Payload) = PayloadU()

RigidBody.MassTrait(::System{Payload}) = HasMass()
RigidBody.WrenchTrait(::System{Payload}) = GetsNoExternalWrench()
RigidBody.AngMomTrait(::System{Payload}) = HasNoAngularMomentum()

function RigidBody.get_mp_Ob(sys::System{Payload})
    @unpack m_pilot, m_copilot, m_lpass, m_rpass, m_baggage = sys.u
    @unpack pilot_slot, copilot_slot, lpass_slot, rpass_slot, baggage_slot = sys.params

    pilot = MassProperties(PointDistribution(m_pilot), pilot_slot)
    copilot = MassProperties(PointDistribution(m_copilot), copilot_slot)
    lpass = MassProperties(PointDistribution(m_lpass), lpass_slot)
    rpass = MassProperties(PointDistribution(m_rpass), rpass_slot)
    baggage = MassProperties(PointDistribution(m_baggage), baggage_slot)

    mp_Ob = MassProperties() + pilot + copilot + lpass + rpass + baggage
    return mp_Ob
end

#################################### GUI #######################################

function GUI.draw!(sys::System{<:Payload}, label::String = "Cessna 172R Payload")

    u = sys.u

    CImGui.Begin(label)

    CImGui.PushItemWidth(-60)

    u.m_pilot = GUI.safe_slider("Pilot Mass (kg)", u.m_pilot, 0, 100, "%.3f")
    u.m_copilot = GUI.safe_slider("Copilot Mass (kg)", u.m_copilot, 0, 100, "%.3f")
    u.m_lpass = GUI.safe_slider("Left Passenger Mass (kg)", u.m_lpass, 0, 100, "%.3f")
    u.m_rpass = GUI.safe_slider("Right Passenger Mass (kg)", u.m_rpass, 0, 100, "%.3f")
    u.m_baggage = GUI.safe_slider("Baggage Mass (kg)", u.m_baggage, 0, 100, "%.3f")

    CImGui.PopItemWidth()

    CImGui.End()

end


############################## Fuel #########################################

#assumes fuel is drawn equally from both tanks, no need to model them
#individually for now
Base.@kwdef struct Fuel <: Piston.AbstractFuelSupply
    m_full::Float64 = 114.4 #maximum fuel mass (42 gal * 6 lb/gal * 0.454 kg/lb)
    m_res::Float64 = 1.0 #residual fuel mass
end

Base.@kwdef struct FuelY
    m_total::Float64 = 0.0 #total fuel mass
    m_avail::Float64 = 0.0 #available fuel mass
end

#normalized fuel content (0: residual, 1: full)
Systems.init(::SystemX, ::Fuel) = [0.5] #cannot be a scalar, need an AbstractVector{<:Real}
Systems.init(::SystemY, ::Fuel) = FuelY()

function Systems.f_ode!(sys::System{Fuel}, pwp::System{<:Piston.Thruster})

    @unpack m_full, m_res = sys.params #no need for subsystems
    m_total = m_res + sys.x[1] * (m_full - m_res) #current mass
    m_avail = m_total - m_res
    sys.ẋ .= -pwp.y.engine.ṁ / (m_full - m_res)
    sys.y = FuelY(; m_total, m_avail)

end

Piston.fuel_available(sys::System{<:Fuel}) = (sys.y.m_avail > 0)

function RigidBody.get_mp_Ob(fuel::System{Fuel})

    #in case x becomes negative (fuel consumed beyond x=0 before the engine
    #dies)
    m_fuel = max(0.0, fuel.y.m_total)

    m_left = PointDistribution(0.5m_fuel)
    m_right = PointDistribution(0.5m_fuel)

    #fuel tanks reference frames
    frame_left = FrameTransform(r = SVector{3}(0.325, -2.845, 0))
    frame_right = FrameTransform(r = SVector{3}(0.325, 2.845, 0))

    mp_Ob = MassProperties()
    mp_Ob += MassProperties(m_left, frame_left)
    mp_Ob += MassProperties(m_right, frame_right)

    return mp_Ob
end

function GUI.draw(sys::System{Fuel}, window_label::String = "C172R Fuel System")

    @unpack m_total, m_avail = sys.y

    CImGui.Begin(window_label)

        CImGui.Text(@sprintf("Total Fuel: %.6f kg", m_total))
        CImGui.Text(@sprintf("Available Fuel: %.6f kg", m_avail))

    CImGui.End()

end

################################ Powerplant ####################################

Pwp() = Piston.Thruster(propeller = Propeller(t_bp = FrameTransform(r = [2.055, 0, 0.833])))


################################################################################
################################ Airframe ######################################

#P is introduced as a type parameter, because Piston.Thruster is itself a
#parametric type, and therefore not concrete
Base.@kwdef struct Airframe{P} <: AbstractAirframe
    str::Structure = Structure()
    aero::Aero = Aero()
    ldg::Ldg = Ldg()
    fuel::Fuel = Fuel()
    pld::Payload = Payload()
    pwp::P = Pwp()
end

############################# Update Methods ###################################


function Systems.f_ode!(airframe::System{<:Airframe},
                        kin::KinematicData, air::AirData, trn::System{<:AbstractTerrain})

    @unpack aero, pwp, ldg, fuel, pld = airframe

    f_ode!(aero, pwp, air, kin, trn)
    f_ode!(ldg, kin, trn) #update landing gear continuous state & outputs
    f_ode!(pwp, air, kin) #update powerplant continuous state & outputs
    f_ode!(fuel, pwp) #update fuel system

    update_y!(airframe)

end

function Systems.f_step!(airframe::System{<:Airframe}, ::KinematicSystem)
    @unpack aero, pwp, fuel, ldg, fuel, pld = airframe

    x_mod = false
    x_mod |= f_step!(aero)
    x_mod |= f_step!(ldg)
    x_mod |= f_step!(pwp, fuel)
    return x_mod

end



################################################################################
#################################### GUI #######################################


function GUI.draw(sys::System{<:Airframe}, window_label::String = "Cessna 172R Airframe")

    @unpack pwp, ldg, aero, fuel, pld = sys

    CImGui.Begin(window_label)

        show_aero = @cstatic check=false @c CImGui.Checkbox("Aerodynamics", &check)
        show_ldg = @cstatic check=false @c CImGui.Checkbox("Landing Gear", &check)
        show_pwp = @cstatic check=false @c CImGui.Checkbox("Powerplant", &check)
        show_fuel = @cstatic check=false @c CImGui.Checkbox("Fuel", &check)
        show_pld = @cstatic check=false @c CImGui.Checkbox("Payload", &check)

    CImGui.End()

    show_aero && GUI.draw(aero)
    show_ldg && GUI.draw(ldg)
    show_pwp && GUI.draw(pwp)
    show_fuel && GUI.draw(fuel)
    show_pld && GUI.draw!(pld)

end


end #module