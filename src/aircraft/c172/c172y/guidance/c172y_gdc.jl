


@enumx T=ModeGuidanceLonEnum ModeGuidanceLon begin
    off = 0
    alt = 1
end

@enumx T=ModeGuidanceLatEnum ModeGuidanceLat begin
    off = 0
    seg = 1
end

@enumx T=FlightPhaseEnum FlightPhase begin
    gnd = 0
    air = 1
end

using .FlightPhase: FlightPhaseEnum
using .ModeGuidanceLon: ModeGuidanceLonEnum
using .ModeGuidanceLat: ModeGuidanceLatEnum


#enable JSON parsing of integers as ModeGuidanceLonEnum
StructTypes.StructType(::Type{ModeGuidanceLonEnum}) = StructTypes.CustomStruct()
StructTypes.lowertype(::Type{ModeGuidanceLonEnum}) = Int32 #default enum type
StructTypes.lower(x::ModeGuidanceLonEnum) = Int32(x)

#enable JSON parsing of integers as ModeGuidanceLatEnum
StructTypes.StructType(::Type{ModeGuidanceLatEnum}) = StructTypes.CustomStruct()
StructTypes.lowertype(::Type{ModeGuidanceLatEnum}) = Int32 #default enum type
StructTypes.lower(x::ModeGuidanceLatEnum) = Int32(x)


################################################################################
############################### Guidance #######################################

#a discrete model implementing a specific longitudinal or lateral guidance mode
abstract type AbstractGuidanceChannel <: ModelDefinition end


################################################################################
############################### Vertical Guidance ##############################

@enumx T=AltTrackingStateEnum AltTrackingState begin
    acquire = 0
    hold = 1
end

using .AltTrackingState: AltTrackingStateEnum

@kwdef struct AltitudeTracking <: AbstractGuidanceChannel
    k_h2c::Float64 = 0.2 #h_err to climb rate gain
    h_thr::Float64 = 10.0 #altitude threshold for state switching
end

@kwdef struct AltitudeTrackingY
    state::AltTrackingStateEnum = AltTrackingState.hold
    mode_ctl_lon::ModeControlLonEnum = ModeControlLon.EAS_clm
    h_err::Float64 = 0.0 #current altitude error
    throttle_ref::Float64 = 0.0
    clm_ref::Float64 = 0.0
end

Modeling.U(::AltitudeTracking) = Ref{Union{HEllip, HOrth}}(HEllip(0.0))
Modeling.S(::AltitudeTracking) = Ref{AltTrackingStateEnum}(AltTrackingState.hold)
Modeling.Y(::AltitudeTracking) = AltitudeTrackingY()

altitude_error(h_target::HEllip, kin_data::KinData) = h_target - kin_data.h_e
altitude_error(h_target::HOrth, kin_data::KinData) = h_target - kin_data.h_o

function Modeling.f_periodic!(::NoScheduling, mdl::Model{<:AltitudeTracking},
                        vehicle::Model{<:C172Y.Vehicle})

    state = mdl.s[]
    h_target = mdl.u[]
    @unpack k_h2c, h_thr = mdl.constants

    h_err = altitude_error(h_target, vehicle.y.kinematics)

    if state === AltTrackingState.acquire

        mode_ctl_lon = ModeControlLon.thr_EAS
        throttle_ref = h_err > 0 ? 1.0 : 0.0 #full throttle to climb, idle to descend
        clm_ref = 0.0 #no effect
        (abs(h_err) < h_thr - 1) && (mdl.s[] = AltTrackingState.hold)

    else #AltTrackingState.hold

        mode_ctl_lon = ModeControlLon.EAS_clm
        throttle_ref = 0.0 #no effect
        clm_ref = k_h2c * h_err
        (abs(h_err) > h_thr + 1) && (mdl.s[] = AltTrackingState.acquire)

    end

    mdl.y = AltitudeTrackingY(; state, mode_ctl_lon, h_err, throttle_ref, clm_ref)

end

function Control.reset!(mdl::Model{<:AltitudeTracking})
    mdl.u[] = HEllip(0.0)
    mdl.s[] = AltTrackingState.hold
end

function GUI.draw(mdl::Model{<:AltitudeTracking}, p_open::Ref{Bool} = Ref(true))

    @unpack state, mode_ctl_lon, h_err, throttle_ref, clm_ref = mdl.y
    @unpack h_target[] = mdl.u
    @unpack h_thr = mdl.constants

    Begin("Altitude Guidance", p_open)
        Text("State: $state")
        Text("Control Mode: $mode_ctl_lon")
        Text("Altitude Ref: $(Float64(h_target))")
        Text("Altitude Datum: $(typeof(h_target))")
        Text("Altitude Error: $h_err")
        Text("Altitude Threshold: $h_thr")
        Text("Throttle Setpoint: $throttle_ref")
        Text("Climb Rate Setpoint: $clm_ref")
    End()
end



################################################################################
############################# HorizontalGuidance ###############################

struct Segment
    p1::Geographic{LatLon, Ellipsoidal}
    p2::Geographic{LatLon, Ellipsoidal}
    u_12_h::SVector{3,Float64}
    l_12_h::Float64
    χ_12::Float64
    γ_12::Float64
    q_en::RQuat

    function Segment(ap1::Abstract3DPosition, ap2::Abstract3DPosition)

        p1 = Geographic{LatLon,Ellipsoidal}(ap1)
        p2 = Geographic{LatLon,Ellipsoidal}(ap2)

        r_e1_e = Cartesian(ap1)
        r_e2_e = Cartesian(ap2)
        r_12_e = r_e2_e - r_e1_e

        n_e1 = NVector(p1)
        n_e2 = NVector(p2)
        n_e = NVector(0.5*(n_e1[:] + n_e2[:]))
        q_en = ltf(n_e) #LTF at midpoint

        r_12_n = q_en'(r_12_e)
        r_12_h = SVector{3,Float64}(r_12_n[1], r_12_n[2], 0.0)
        l_12_h = norm(r_12_h)
        @assert l_12_h > 0 "Invalid guidance segment" #must not be vertical or degenerate

        Δh_12 = HEllip(p2) - HEllip(p1)
        u_12_h = r_12_h / l_12_h
        χ_12 = azimuth(u_12_h) #azimuth
        γ_12 = atan(Δh_12, l_12_h) #flight path angle (0 for equal altitude)

        new(p1, p2, u_12_h, l_12_h, χ_12, γ_12, q_en)

    end

end

#construct a Segment from:
#1) initial point p1
#2) an azimuth χ measured on NED(p1)
#3) a distance l along χ measured on the horizontal plane of NED(p1)
#4) an altitude increment Δh

#note: l is the straight-line distance from p1 to p2t, NOT from p1 to p2.
#however, as long as p1 and p2 are relatively close, both their straight-line
#and geodesic distances will be close to l

function Segment(p1::Abstract3DPosition; χ::Real, l::Real, Δh::Real = 0.0)
    geo1 = Geographic{LatLon, Ellipsoidal}(p1)
    q_en1 = ltf(geo1)
    r_12_n1 = SVector{3,Float64}(l*cos(χ), l*sin(χ), 0)
    r_12_e = q_en1(r_12_n1)
    r_e1_e = Cartesian(p1)
    r_e2_e = r_e1_e + r_12_e
    geo2 = Geographic(LatLon(r_e2_e), HEllip(geo1) + Δh)
    Segment(geo1, geo2)
end


Segment() = Segment(Geographic(LatLon(), HEllip()),
                    Geographic(LatLon(0.001, 0), HEllip()))

StructTypes.StructType(::Type{Segment}) = StructTypes.CustomStruct()
StructTypes.lowertype(::Type{Segment}) = NTuple{2, Geographic{LatLon, Ellipsoidal}}
StructTypes.lower(seg::Segment) = (seg.p1, seg.p2)
function StructTypes.construct(::Type{Segment}, ps::NTuple{2, Geographic{LatLon, Ellipsoidal}})
    Segment(ps[1], ps[2])
end


################################################################################


@kwdef struct SegmentGuidance <: AbstractGuidanceChannel
    Δχ_inf::Ranged{Float64, 0., π/2.} = π/2 #intercept angle for x2d_1p → ∞ (rad)
    e_sf::Float64 = 1000 #cross-track distance scale factor (m)
end

@kwdef struct SegmentGuidanceY
    target::Segment = Segment()
    l_1b_h::Float64 = 0.0 #along-track horizontal distance
    e_1b_h::Float64 = 0.0 #cross-track horizontal distance, positive right
    χ_ref::Float64 = 0.0 #course angle reference
end

Modeling.U(::SegmentGuidance) = Ref(Segment())
Modeling.Y(::SegmentGuidance) = SegmentGuidanceY()

function Modeling.f_periodic!(::NoScheduling, mdl::Model{<:SegmentGuidance},
                        vehicle::Model{<:C172Y.Vehicle})

    @unpack target = mdl.u[]
    @unpack Δχ_inf, e_sf = mdl.constants
    @unpack p1, q_en, χ_12, u_12_h = target

    r_eb_e = Cartesian(vehicle.y.kinematics.r_eb_e)
    r_e1_e = Cartesian(p1)
    r_1b_e = r_eb_e - r_e1_e
    r_1b_n = q_en'(r_1b_e)
    r_1b_h = SVector{3,Float64}(r_1b_n[1], r_1b_n[2], 0)
    l_1b_h = u_12_h ⋅ r_1b_h
    e_1b_h = (u_12_h × r_1b_h)[3]
    Δχ = -Float64(Δχ_inf)/(π/2) * atan(e_1b_h / e_sf)
    χ_ref = wrap_to_π(χ_12 + Δχ)

    mdl.y = SegmentGuidanceY(; target, l_1b_h, e_1b_h, χ_ref)

end
    # #these go in the longitudinal guidance law
    # r_eb_e = Cartesian(vehicle.y.kinematics.r_eb_e)
    # r_1b_e = r_eb_e - r_e1_e #position vector from p1 to p
    # r_1b_n1 = q_n1e(r_1b_e)
    # r_1b_h = SVector{3,Float64}(r_1b_n1[1], r_1b_n1[2], 0)
    # l_1b_h = u_12_h ⋅ r_1b_h
    # h_ref = h1 + l_1b_h * tan(γ_12) #nominal altitude at our location
    # clm_ff = V_gnd * sin(γ_12) #nominal climb rate at our current V_gnd
    # clm_ref = clm_ff + mdl.k_h2c * (h_ref - h)

    #we need to compute a normalized parameter p_1b.
    #If p_1b<0, we're behind p1, and we set clm_ff = 0, h_nom =
    #if p_1b > 1, we're beyond p2 and we set clm = 0, h_nom = h2.

################################################################################
################################################################################
################################################################################

function Modeling.f_periodic!(::NoScheduling, mdl::Model{<:ControlLaws},
                        vehicle::Model{<:C172Y.Vehicle})
    @unpack mode_gdc_lon_req, mode_gdc_lat_req, mode_ctl_lon_req, mode_ctl_lat_req,
            eng_start, eng_stop, mixture, flaps, brake_left, brake_right,
            throttle_axis, aileron_axis, elevator_axis, rudder_axis,
            throttle_offset, aileron_offset, elevator_offset, rudder_offset,
            q_ref, EAS_ref, θ_ref, clm_ref, p_ref, φ_ref, χ_ref, β_ref,
            h_target, seg_target = mdl.u

    @unpack ctl_lon, ctl_lat, gdc_lon_alt, gdc_lat_seg = mdl.submodels

    throttle_ref = throttle_axis + throttle_offset
    elevator_ref = elevator_axis + elevator_offset
    aileron_ref = aileron_axis + aileron_offset
    rudder_ref = rudder_axis + rudder_offset

    any_wow = any(SVector{3}(leg.strut.wow for leg in vehicle.y.systems.ldg))
    flight_phase = any_wow ? FlightPhase.gnd : FlightPhase.air

    if flight_phase === FlightPhase.gnd

        mode_gdc_lon = ModeGuidanceLon.off
        mode_gdc_lat = ModeGuidanceLat.off
        mode_ctl_lon = ModeControlLon.direct
        mode_ctl_lat = ModeControlLat.direct

    elseif flight_phase === FlightPhase.air

        mode_gdc_lon = mode_gdc_lon_req
        mode_gdc_lat = mode_gdc_lat_req

        if mode_gdc_lon === ModeGuidanceLon.off

            mode_ctl_lon = mode_ctl_lon_req

        else #mode_gdc_lon === ModeGuidanceLon.alt

            gdc_lon_alt.u[] = h_target
            f_periodic!(gdc_lon_alt, vehicle)

            mode_ctl_lon = gdc_lon_alt.y.mode_ctl_lon
            throttle_ref = gdc_lon_alt.y.throttle_ref
            clm_ref = gdc_lon_alt.y.clm_ref

        end

        if mode_gdc_lat === ModeGuidanceLat.off

            mode_ctl_lat = mode_ctl_lat_req

            #below a v_gnd threshold, override χ mode and revert to φ
            if (mode_ctl_lat === ModeControlLat.χ_β) && (vehicle.y.kinematics.v_gnd < 10.0)
                mode_ctl_lat = ModeControlLat.ModeControlLat.φ_β
            end

        else #mode_gdc_lat === ModeGuidanceLat.seg

            gdc_lat_seg.u[] = seg_target
            f_periodic!(gdc_lat_seg, vehicle)
            mode_ctl_lat = ModeControlLat.χ_β
            χ_ref = gdc_lat_seg.y.χ_ref

        end

    end

    ctl_lon.u.mode = mode_ctl_lon
    @pack! ctl_lon.u = throttle_ref, elevator_ref, q_ref, θ_ref, EAS_ref, clm_ref
    f_periodic!(ctl_lon, vehicle)

    ctl_lat.u.mode = mode_ctl_lat
    @pack! ctl_lat.u = aileron_ref, rudder_ref, p_ref, φ_ref, β_ref, χ_ref
    f_periodic!(ctl_lat, vehicle)

    mdl.y = ControlLawsY(; flight_phase,
        mode_gdc_lon, mode_gdc_lat, mode_ctl_lon, mode_ctl_lat,
        ctl_lon = ctl_lon.y, ctl_lat = ctl_lat.y,
        gdc_lon_alt = gdc_lon_alt.y, gdc_lat_seg = gdc_lat_seg.y)

end


    ############################### Guidance ###################################

    @cstatic h_datum=Int32(0) begin #ellipsoidal: 0, orthometric: 1
    if CImGui.CollapsingHeader("Vertical Guidance")

        if BeginTable("Longitudinal Guidance", 3, CImGui.ImGuiTableFlags_SizingStretchProp )#| CImGui.ImGuiTableFlags_Resizable)# | CImGui.ImGuiTableFlags_BordersInner)
            TableNextRow()
                TableNextColumn();
                    Text("Mode")
                TableNextColumn();
                    mode_button("Off", ModeGuidanceLon.off, u.mode_gdc_lon_req, y.mode_gdc_lon)
                    IsItemActive() && (u.mode_gdc_lon_req = ModeGuidanceLon.off); SameLine()
                    mode_button("Altitude", ModeGuidanceLon.alt, u.mode_gdc_lon_req, y.mode_gdc_lon)
                    if IsItemActive()
                        u.mode_gdc_lon_req = ModeGuidanceLon.alt
                        u.h_target = (h_datum == 0 ? h_e : h_o)
                        u.EAS_ref = EAS
                    end
                TableNextColumn();
            TableNextRow()
                TableNextColumn(); AlignTextToFramePadding(); Text("Altitude (m)")
                TableNextColumn();
                    PushItemWidth(-10)
                    h_target_f64 = safe_input("Altitude Setpoint", Float64(u.h_target), 1, 1.0, "%.3f")
                    u.h_target = (h_datum == 0 ? HEllip(h_target_f64) : HOrth(h_target_f64))
                    PopItemWidth()
                TableNextColumn();
                    h_datum == 0 && Text(@sprintf("%.3f", Float64(h_e)))
                    h_datum == 1 && Text(@sprintf("%.3f", Float64(h_o)))
            TableNextRow()
                TableNextColumn();
                TableNextColumn();
                    @c RadioButton("Ellipsoidal", &h_datum, 0); SameLine()
                    IsItemActive() && (u.h_target = h_e)
                    @c RadioButton("Orthometric", &h_datum, 1)
                    IsItemActive() && (u.h_target = h_o)
            EndTable()
        end #table

        Separator()

    end #header
    end #cstatic


#original C172X1 Avionics, reference for multi-component Avionics

################################################################################
############################### Avionics #######################################


@kwdef struct Avionics{C} <: AbstractAvionics
    ctl::C = Subsampled(ControlLaws(), 2)
end

@sm_updates Avionics

############################### Update methods #################################

function AircraftBase.assign!(systems::Model{<:C172Y.Systems},
                          avionics::Model{<:C172Yv1.Avionics})
    AircraftBase.assign!(systems, avionics.ctl)
end

################################# Trimming #####################################

function Modeling.init!(avionics::Model{<:C172Yv1.Avionics},
                            vehicle::Model{<:C172Y.Vehicle})

    Modeling.init!(avionics.ctl, vehicle)
    update_output!(avionics)

end

################################## GUI #########################################

function GUI.draw!(avionics::Model{<:C172Yv1.Avionics},
                    vehicle::Model{<:C172Y.Vehicle},
                    p_open::Ref{Bool} = Ref(true),
                    label::String = "Cessna172Yv1 Avionics")

    CImGui.Begin(label, p_open)

    @cstatic c_ctl=false begin
        @c CImGui.Checkbox("Flight Control", &c_ctl)
        c_ctl && @c GUI.draw!(avionics.ctl, vehicle, &c_ctl)
    end

    CImGui.End()

end
