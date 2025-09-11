module C172YGuidance

using UnPack, StaticArrays, LinearAlgebra
using StructTypes
using EnumX

using CImGui: Begin, End, PushItemWidth, PopItemWidth, AlignTextToFramePadding,
    Dummy, SameLine, NewLine, IsItemActive, IsItemActivated, Separator, Text,
    Bullet, TextUnformatted, Checkbox, Combo, RadioButton, TableNextColumn,
    TableNextRow, BeginTable, EndTable, TableSetupColumn

using Flight.FlightCore
using Flight.FlightLib

using ...AircraftBase
using ...C172
using ..C172Y: Vehicle, Systems
using ..C172Y.C172YControl: ControlLaws, ModeControlLat, ModeControlLon, is_on_gnd


################################################################################
#################################### Modes #####################################

@enumx T=ModeGuidanceEnum ModeGuidance begin
    direct = 0
    segment = 1
    circular = 2
end

using .ModeGuidance: ModeGuidanceEnum

#enable JSON parsing of integers as ModeGuidanceEnum
StructTypes.StructType(::Type{ModeGuidanceEnum}) = StructTypes.CustomStruct()
StructTypes.lowertype(::Type{ModeGuidanceEnum}) = Int32 #default enum type
StructTypes.lower(x::ModeGuidanceEnum) = Int32(x)


################################################################################
################################## Segment #####################################

struct Segment
    p1::Geographic{LatLon, Ellipsoidal}
    p2::Geographic{LatLon, Ellipsoidal}
    u_12::SVector{2,Float64}
    l_12::Float64
    χ_12::Float64
    γ_12::Float64

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
        r_12_h = SVector{2,Float64}(r_12_n[1], r_12_n[2])
        l_12 = norm(r_12_h)

        l_12_min = 1e-6 #minimum horizontal length
        (l_12 < l_12_min) && throw(ArgumentError("Segment is vertical or degenerate"))

        Δh_12 = HEllip(p2) - HEllip(p1)
        u_12 = r_12_h / l_12
        χ_12 = azimuth(u_12) #azimuth
        γ_12 = atan(Δh_12, l_12) #flight path angle (0 for equal altitude)

        new(p1, p2, u_12, l_12, χ_12, γ_12)

    end

end

Segment() = Segment(Geographic(LatLon()), Geographic(LatLon(1e-3, 0)))

"""
construct a Segment from:
1) initial point p1
2) an azimuth χ on NED(p1)
3) a length l along χ measured on the horizontal plane of NED(p1)
4) a flight path angle γ or altitude increment Δh
note: l is the straight-line distance from p1 to p2t, NOT from p1 to p2.
however, as long as p1 and p2 are relatively close, both their straight-line
and geodesic distances will be close to l
"""
function Segment(p1::Abstract3DPosition; l::Real, χ::Real, kw...)

    if :γ ∈ keys(kw) && :Δh ∈ keys(kw)
        error("Can only provide either γ (flight path angle) or Δh
              (altitude increment) as keyword arguments")
    elseif :γ ∈ keys(kw)
        γ = kw[:γ]
        Δh = l * tan(γ)
    elseif :Δh ∈ keys(kw)
        Δh = kw[:Δh]
    else
        error("Must provide either γ (flight path angle) or Δh
              (altitude increment) as keyword arguments")
    end

    geo1 = Geographic{LatLon, Ellipsoidal}(p1)
    q_en1 = ltf(geo1)
    r_12_n1 = SVector{3,Float64}(l*cos(χ), l*sin(χ), 0)
    r_12_e = q_en1(r_12_n1)
    r_e1_e = Cartesian(p1)
    r_e2_e = r_e1_e + r_12_e
    geo2 = Geographic(LatLon(r_e2_e), HEllip(geo1) + Δh)
    Segment(geo1, geo2)
end

Base.:-(seg::Segment) = Segment(seg.p2, seg.p1)

StructTypes.StructType(::Type{Segment}) = StructTypes.CustomStruct()
StructTypes.lowertype(::Type{Segment}) = NTuple{2, Geographic{LatLon, Ellipsoidal}}
StructTypes.lower(seg::Segment) = (seg.p1, seg.p2)
function StructTypes.construct(::Type{Segment}, ps::NTuple{2, Geographic{LatLon, Ellipsoidal}})
    Segment(ps[1], ps[2])
end



#segment-relative coordinates
@kwdef struct SegmentCoords
    l_1b::Float64 = 0.0 #signed along-track horizontal distance
    e_1b::Float64 = 0.0 #signed cross-track horizontal distance, positive right
    v_1b::Float64 = 0.0 #signed vertical distance, positive up
    h_1b::HEllip = HEllip(0.0) #nominal segment altitude at l_1b
end

function SegmentCoords(s::Segment, p::Abstract3DPosition)

    @unpack p1, χ_12, γ_12, u_12 = s

    q_en = ltf(p)
    r_eb_e = Cartesian(p)
    r_e1_e = Cartesian(p1) #ECEF position of Segment's origin p1
    r_1b_e = r_eb_e - r_e1_e #3D vector from p1 to b, ECEF coordinates
    r_1b_n = q_en'(r_1b_e) #3D vector from p1 to b, Segment's NED coordinates
    r_1b_h = SVector{3,Float64}(r_1b_n[1], r_1b_n[2], 0.0) #horizontal components
    u_12_h = SVector{3, Float64}(u_12[1], u_12[2], 0.0)
    l_1b = u_12_h ⋅ r_1b_h #signed along-track horizontal distance
    e_1b = (u_12_h × r_1b_h)[3] #signed cross-track signed horizontal distance
    h_1b = HEllip(p1) + l_1b * tan(γ_12) #vertical distance l_1b_h
    v_1b = h_1b - HEllip(p)

    return SegmentCoords(; l_1b, e_1b, v_1b, h_1b)

end

function GUI.draw(seg::Segment)

    @unpack p1, p2, u_12, l_12, χ_12, γ_12 = seg

    if BeginTable("Vector", 3, CImGui.ImGuiTableFlags_Resizable | CImGui.ImGuiTableFlags_RowBg)
        TableNextRow()
        TableNextRow()
        TableNextColumn(); Text("Origin Latitude (deg)")
        TableNextColumn(); Text("Origin Longitude (deg)")
        TableNextColumn(); Text("Origin Altitude (m)")
        TableNextRow()
        TableNextColumn(); Text(@sprintf("%.6f deg", rad2deg(p1.loc.ϕ)))
        TableNextColumn(); Text(@sprintf("%.6f deg", rad2deg(p1.loc.λ)))
        TableNextColumn(); Text(@sprintf("%.3f m", Float64(p1.h)))
        TableNextRow()
        TableNextColumn(); Text("End Latitude (deg)")
        TableNextColumn(); Text("End Longitude (deg)")
        TableNextColumn(); Text("End Altitude (m)")
        TableNextRow()
        TableNextColumn(); Text(@sprintf("%.6f deg", rad2deg(p2.loc.ϕ)))
        TableNextColumn(); Text(@sprintf("%.6f deg", rad2deg(p2.loc.λ)))
        TableNextColumn(); Text(@sprintf("%.3f m", Float64(p2.h)))
        TableNextRow()
        TableNextColumn(); Text("Azimuth (deg)")
        TableNextColumn(); Text("Inclination (deg)")
        TableNextColumn(); Text("Horizontal Length (m)")
        TableNextRow()
        TableNextColumn(); Text(@sprintf("%.6f deg", rad2deg(χ_12)))
        TableNextColumn(); Text(@sprintf("%.6f deg", rad2deg(γ_12)))
        TableNextColumn(); Text(@sprintf("%.3f m", l_12))
        EndTable()
    end

end

################################################################################

#e_sf: cross-track error for which intercept angle drops to Δχ_inf / 2
@kwdef struct SegmentGuidance <: ModelDefinition
    Δχ_inf::Ranged{Float64, 0., π/2.} = π/2 #intercept angle for x2d_1p → ∞ (rad)
    e_sf::Float64 = 500 #cross-track error scaling parameter for lateral guidance (m)
    e_thr = 1000 #cross-track error threshold for vertical guidance (m)
end

@kwdef mutable struct SegmentGuidanceU
    target::Segment = Segment()
    horizontal_req::Bool = false #request horizontal guidance
    vertical_req::Bool = false #request vertical guidance
end

@kwdef struct SegmentGuidanceY
    target::Segment = Segment()
    coords::SegmentCoords = SegmentCoords()
    Δχ::Float64 = 0.0 #intercept angle (rad)
    χ_ref::Float64 = 0.0 #course angle reference (rad)
    h_ref::HEllip = 0.0 #altitude reference (m)
    horizontal::Bool = false #horizontal guidance active
    vertical::Bool = false #vertical guidance active
end

Modeling.U(::SegmentGuidance) = SegmentGuidanceU()
Modeling.Y(::SegmentGuidance) = SegmentGuidanceY()

function Modeling.f_periodic!(::NoScheduling, mdl::Model{<:SegmentGuidance},
                                ctl::Model{<:ControlLaws},
                                vehicle::Model{<:Vehicle})

    @unpack target, horizontal_req, vertical_req = mdl.u
    @unpack Δχ_inf, e_sf, e_thr = mdl.constants
    @unpack r_eb_e = vehicle.kinematics.y

    coords = SegmentCoords(target, Cartesian(r_eb_e))
    @unpack l_1b, e_1b, h_1b = coords

    Δχ = -Float64(Δχ_inf)/(π/2) * atan(e_1b / e_sf)
    χ_ref = wrap_to_π(target.χ_12 + Δχ)
    h_ref = h_1b

    horizontal = horizontal_req
    vertical = (abs(e_1b) < e_thr ? vertical_req : false)

    if horizontal
        ctl.lat.u.mode_req = ModeControlLat.χ_β
        ctl.lat.u.χ_ref = χ_ref
        #β_ref unchanged
    end

    if vertical
        ctl.lon.u.mode_req = ModeControlLon.EAS_alt
        ctl.lon.u.h_ref = h_ref
        #EAS_ref unchanged
    end

    mdl.y = SegmentGuidanceY(; target, coords, Δχ, χ_ref, h_ref, horizontal, vertical)

end


################################################################################
########################### CircularGuidance ###################################

@kwdef struct CircularGuidance <: ModelDefinition end
@kwdef struct CircularGuidanceY end
function Modeling.f_periodic!(::NoScheduling, mdl::Model{<:CircularGuidance},
                                ctl::Model{<:ControlLaws},
                                vehicle::Model{<:Vehicle})
end

###############################################################################
############################# GuidanceLaws #####################################

@kwdef struct GuidanceLaws <: ModelDefinition
    seg::SegmentGuidance = SegmentGuidance()
    crc::CircularGuidance = CircularGuidance()
end

@no_ode GuidanceLaws
@no_step GuidanceLaws
@sm_periodic GuidanceLaws

@kwdef mutable struct GuidanceLawsU
    mode_req::ModeGuidanceEnum = ModeGuidance.direct
end

@kwdef struct GuidanceLawsY
    mode::ModeGuidanceEnum = ModeGuidance.direct
    seg::SegmentGuidanceY = SegmentGuidanceY()
    crc::CircularGuidanceY = CircularGuidanceY()
end

Modeling.U(::GuidanceLaws) = GuidanceLawsU()
Modeling.Y(::GuidanceLaws) = GuidanceLawsY()


########################### Update Methods #####################################

function Modeling.f_periodic!(::NoScheduling, mdl::Model{<:GuidanceLaws},
                                ctl::Model{<:ControlLaws},
                                vehicle::Model{<:Vehicle})

    @unpack mode_req = mdl.u
    @unpack seg, crc = mdl.submodels

    mode = (is_on_gnd(vehicle) ? ModeGuidance.direct : mode_req)
    mode === ModeGuidance.segment && f_periodic!(seg, ctl, vehicle)
    mode === ModeGuidance.circular && f_periodic!(crc, ctl, vehicle)

    mdl.y = GuidanceLawsY(; mode, seg = seg.y, crc = crc.y)

end


################################## GUI #########################################

function GUI.draw!(gdc::Model{<:GuidanceLaws}, vehicle::Model{<:Vehicle},
                    p_open::Ref{Bool} = Ref(true))

    @unpack u, y = gdc

    Begin("Cessna172Y Guidance", p_open)

    AlignTextToFramePadding(); Text("Guidance Mode"); SameLine()
    mode_button("Direct", ModeGuidance.direct, u.mode_req, y.mode); SameLine()
    IsItemActive() && (u.mode_req = ModeGuidance.direct)
    mode_button("Segment", ModeGuidance.segment, u.mode_req, y.mode); SameLine()
    IsItemActive() && (u.mode_req = ModeGuidance.segment)
    mode_button("Circular", ModeGuidance.circular, u.mode_req, y.mode)
    IsItemActive() && (u.mode_req = ModeGuidance.circular)

    # CImGui.CollapsingHeader("Segment Guidance") && GUI.draw!(gdc.seg, vehicle)
    # CImGui.CollapsingHeader("Circular Guidance") && GUI.draw!(gdc.crc, vehicle)
    Separator()
    # gdc.y.mode === ModeGuidance.segment && GUI.draw!(gdc.seg, vehicle)
    GUI.draw!(gdc.seg, vehicle)
    gdc.y.mode === ModeGuidance.circular && GUI.draw!(gdc.crc, vehicle)

    End()

end

function GUI.draw!(gdc::Model{<:SegmentGuidance}, vehicle::Model{<:Vehicle})

    @unpack u, constants, y = gdc
    @unpack ϕ_λ, h_e = vehicle.y.kinematics
    @unpack ϕ, λ = ϕ_λ
    @unpack target, coords, Δχ, χ_ref, h_ref, horizontal, vertical = y

    @cstatic(
        spec_from = Cint(1),
        spec_to = Cint(1),
        p1 = Geographic(LatLon(), HEllip()),
        p2 = Geographic(LatLon(1e-3, 0), HEllip()),
        χ_12 = 0.0,
        γ_12 = 0.0,
        l_12 = 1e4,
    begin #cstatic block

    if BeginTable("Segment Builder", 2) #, CImGui.ImGuiTableFlags_SizingStretchProp)# | CImGui.ImGuiTableFlags_Resizable)# | CImGui.ImGuiTableFlags_BordersInner)

        TableNextRow()
            TableNextColumn()
                AlignTextToFramePadding(); TextUnformatted("From"); SameLine()
                @c Combo("##ComboFrom", &spec_from, "Coordinates\0Current Position\0")
            TableNextColumn()
                AlignTextToFramePadding(); TextUnformatted("To"); SameLine()
                @c Combo("##ComboTo", &spec_to, "Coordinates\0Vector (NED)\0")

        TableNextRow()
            TableNextColumn();
            if BeginTable("Build From", 2, CImGui.ImGuiTableFlags_SizingStretchProp)# | CImGui.ImGuiTableFlags_Resizable)# | CImGui.ImGuiTableFlags_BordersInner)
                TableNextRow()
                    TableNextColumn(); AlignTextToFramePadding(); TextUnformatted("Latitude (deg)")
                    TableNextColumn()
                    PushItemWidth(-10)
                    ϕ1 = safe_input("##LatitudeFrom", rad2deg(p1.loc.ϕ), 1e-4, 1e-4, "%.6f") |> deg2rad
                    PopItemWidth()
                TableNextRow()
                    TableNextColumn(); AlignTextToFramePadding(); TextUnformatted("Latitude (deg)")
                    TableNextColumn()
                    PushItemWidth(-10)
                    λ1 = safe_input("##LongitudeFrom", rad2deg(p1.loc.λ), 1e-4, 1e-4, "%.6f") |> deg2rad
                    PopItemWidth()
                TableNextRow()
                    TableNextColumn(); AlignTextToFramePadding(); TextUnformatted("Ellipsoidal Altitude (m)"); SameLine()
                    TableNextColumn()
                    PushItemWidth(-10)
                    h1 = safe_input("##Altitude (m)", Float64(p1.h), 1.0, 1.0, "%.3f")
                    PopItemWidth()
                EndTable()
                spec_from == 0 && (p1 = Geographic(LatLon(ϕ1, λ1), HEllip(h1)))
                spec_from == 1 && (p1 = Geographic(ϕ_λ, h_e))
            end

            TableNextColumn();
            if spec_to == 0 #coordinates
                if BeginTable("To Coordinates", 2, CImGui.ImGuiTableFlags_SizingStretchProp)# | CImGui.ImGuiTableFlags_Resizable)# | CImGui.ImGuiTableFlags_BordersInner)
                    TableNextRow()
                        TableNextColumn(); AlignTextToFramePadding(); TextUnformatted("Latitude (deg)")
                        TableNextColumn()
                        PushItemWidth(-10)
                        ϕ2 = safe_input("##LatitudeFrom", rad2deg(p2.loc.ϕ), 1e-3, 1e-3, "%.6f") |> deg2rad
                        PopItemWidth()
                    TableNextRow()
                        TableNextColumn(); AlignTextToFramePadding(); TextUnformatted("Latitude (deg)")
                        TableNextColumn()
                        PushItemWidth(-10)
                        λ2 = safe_input("##LongitudeFrom", rad2deg(p2.loc.λ), 1e-3, 1e-3, "%.6f") |> deg2rad
                        PopItemWidth()
                    TableNextRow()
                        TableNextColumn(); AlignTextToFramePadding(); TextUnformatted("Ellipsoidal Altitude (m)"); SameLine()
                        TableNextColumn()
                        PushItemWidth(-10)
                        h2 = safe_input("##Altitude (m)", Float64(p2.h), 1.0, 1.0, "%.3f")
                        PopItemWidth()
                    EndTable()
                    p2 = Geographic(LatLon(ϕ2, λ2), HEllip(h2))
                end #table
            elseif spec_to == 1 #vector
                if BeginTable("To NED Coordinates", 2, CImGui.ImGuiTableFlags_SizingStretchProp)# | CImGui.ImGuiTableFlags_Resizable)# | CImGui.ImGuiTableFlags_BordersInner)
                    TableNextRow()
                        TableNextColumn(); AlignTextToFramePadding(); TextUnformatted("Azimuth (deg)")
                        TableNextColumn()
                        PushItemWidth(-10)
                        χ_12 = safe_input("##AzimuthTo", rad2deg(χ_12), 1, 1, "%.3f") |> deg2rad
                        PopItemWidth()
                    TableNextRow()
                        TableNextColumn(); AlignTextToFramePadding(); TextUnformatted("Inclination (deg)")
                        TableNextColumn()
                        PushItemWidth(-10)
                        γ_12 = safe_input("##InclinationTo", rad2deg(γ_12), 1e-1, 1e-1, "%.3f") |> deg2rad
                        PopItemWidth()
                    TableNextRow()
                        TableNextColumn(); AlignTextToFramePadding(); TextUnformatted("Horizontal Length (m)"); SameLine()
                        TableNextColumn()
                        PushItemWidth(-10)
                        l_12 = safe_input("##LengthTo", l_12, 1.0, 1.0, "%.3f")
                        PopItemWidth()
                    EndTable()
                end #table
            end

        EndTable()
    end #table

    CImGui.Button("Set Target"); SameLine()
    if IsItemActivated()
        try
            spec_to == 0 && (u.target = Segment(p1, p2))
            spec_to == 1 && (u.target = Segment(p1; χ = χ_12, l = l_12, γ = γ_12))
        catch
            "Segment construction failed"
        end
    end
    mode_button("Horizontal Guidance", true, u.horizontal_req, y.horizontal); SameLine()
    IsItemActivated() && (u.horizontal_req = !u.horizontal_req)
    mode_button("Vertical Guidance", true, u.vertical_req, y.vertical)
    IsItemActivated() && (u.vertical_req = !u.vertical_req)

    end #cstatic block
    ) #cstatic call

    Separator()

    GUI.draw(y.target)
    if BeginTable("SegmentGuidanceData", 2) #, CImGui.ImGuiTableFlags_SizingStretchProp)# | CImGui.ImGuiTableFlags_Resizable)# | CImGui.ImGuiTableFlags_BordersInner)
        TableNextRow();
            TableNextColumn();
                Text(@sprintf("Along-track distance: %.3f m", coords.l_1b))
                Text(@sprintf("Cross-track distance: %.3f m", coords.e_1b))
                Text(@sprintf("Vertical distance: %.3f m", coords.v_1b))
                Text(@sprintf("Intercept angle: %.3f deg", rad2deg(Δχ)))
            TableNextColumn();
                Text(@sprintf("Course angle reference: %.3f deg", rad2deg(χ_ref)))
                Text(@sprintf("Altitude reference: %.3f m", Float64(h_ref)))
                Text("Horizontal guidance: $horizontal")
                Text("Vertical guidance: $vertical")
        EndTable()
    end

end


end #module