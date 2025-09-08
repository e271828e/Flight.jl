module C172YGuidance

using UnPack, StaticArrays, LinearAlgebra
using StructTypes
using EnumX

# using CImGui: Begin, End, PushItemWidth, PopItemWidth, AlignTextToFramePadding,
#     Dummy, SameLine, NewLine, IsItemActive, IsItemActivated, Separator, Text,
#     Checkbox, RadioButton, TableNextColumn, TableNextRow, BeginTable, EndTable

using Flight.FlightCore
using Flight.FlightLib

using ...AircraftBase
using ...C172
using ..C172Y: Vehicle, Systems
using ..C172Y.C172YControl: ControlLaws, ControlLawsY, ModeControlLat, ModeControlLon, is_on_gnd

@enumx T=ModeGuidanceEnum ModeGuidance begin
    direct = 0
    segment = 1
    circle = 2
end

using .ModeGuidance: ModeGuidanceEnum

#enable JSON parsing of integers as ModeGuidanceEnum
StructTypes.StructType(::Type{ModeGuidanceEnum}) = StructTypes.CustomStruct()
StructTypes.lowertype(::Type{ModeGuidanceEnum}) = Int32 #default enum type
StructTypes.lower(x::ModeGuidanceEnum) = Int32(x)


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

Segment() = Segment(Geographic(LatLon()), Geographic(LatLon(1e-3, 0)))

"""
construct a Segment from:
1) initial point p1
2) an azimuth χ measured on NED(p1)
3) a distance l along χ measured on the horizontal plane of NED(p1)
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
        @assert abs(γ) < π/2
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


################################################################################

#e_sf: cross-track error for which intercept angle is Δχ_inf / 2
@kwdef struct SegmentGuidance <: ModelDefinition
    Δχ_inf::Ranged{Float64, 0., π/2.} = π/2 #intercept angle for x2d_1p → ∞ (rad)
    e_sf::Float64 = 1000 #cross-track error scaling parameter for lateral guidance (m)
    e_thr = 100 #cross-track error threshold for vertical guidance (m)
end

@kwdef mutable struct SegmentGuidanceU
    target::Segment = Segment()
    horizontal_req::Bool = true #request horizontal guidance
    vertical_req::Bool = true #request vertical guidance
end

@kwdef struct SegmentGuidanceY
    target::Segment = Segment()
    l_1b_h::Float64 = 0.0 #along-track horizontal distance
    e_1b_h::Float64 = 0.0 #cross-track horizontal distance, positive right
    χ_ref::Float64 = 0.0 #course angle reference
    h_ref::Float64 = 0.0 #altitude reference
    horizontal::Bool = true #horizontal guidance active
    vertical::Bool = true # guidance active
end

Modeling.U(::SegmentGuidance) = SegmentGuidanceU()
Modeling.Y(::SegmentGuidance) = SegmentGuidanceY()

function Modeling.f_periodic!(::NoScheduling, mdl::Model{<:SegmentGuidance},
                                ctl::Model{<:ControlLaws},
                                vehicle::Model{<:Vehicle})

    @unpack target, horizontal_req, vertical_req = mdl
    @unpack Δχ_inf, e_sf = mdl.constants

    @unpack e_1b_h, h_ref = get_guidance_inputs(target, vehicle.y.kinematics.r_eb_e)
    Δχ = -Float64(Δχ_inf)/(π/2) * atan(e_1b_h / e_sf)
    χ_ref = wrap_to_π(χ_12 + Δχ)

    horizontal = horizontal_req
    vertical = (e_1b_h < e_thr ? vertical_req : false)

    if horizontal
        ctl.lat.u.mode_req = ModeControlLat.χ_β
        ctl.lat.u.χ_ref = χ_ref
        #β_ref unchanged
    end

    if vertical
        ctl.lon.u.mode_req = ModeControlLat.EAS_alt
        ctl.lon.u.h_ref = h_ref
        #EAS_ref unchanged
    end

    mdl.y = SegmentGuidanceY(; target, l_1b_h, e_1b_h, χ_ref, h_ref,
                               horizontal, vertical)

end

function get_guidance_inputs(s::Segment, r_eb_e::Cartesian)

    @unpack p1, q_en, χ_12, γ_12, u_12_h = s

    r_e1_e = Cartesian(p1) #ECEF position of Segment's origin p1
    r_1b_e = r_eb_e - r_e1_e #3D vector from p1 to b, ECEF coordinates
    r_1b_n = q_en'(r_1b_e) #3D vector from p1 to b, Segment's NED coordinates
    r_1b_h = SVector{3,Float64}(r_1b_n[1], r_1b_n[2], 0) #horizontal components
    l_1b_h = u_12_h ⋅ r_1b_h #along-track signed distance
    e_1b_h = (u_12_h × r_1b_h)[3] #cross-track signed distance
    h_ref = p1.h + l_1b_h * tan(γ_12) #nominal altitude at l_1b_h

    return (e_1b_h = e_1b_h, h_ref = h_ref)

end

################################################################################
############################# GuidanceLaws #####################################

@kwdef struct GuidanceLaws{C <: ControlLaws} <: AbstractAvionics
    ctl::C = ControlLaws()
    seg::SegmentGuidance = SegmentGuidance()
end

@no_ode GuidanceLaws
@no_step GuidanceLaws
@sm_periodic GuidanceLaws

@kwdef mutable struct GuidanceLawsU
    mode_req::ModeGuidanceEnum = ModeGuidance.direct
end

@kwdef struct GuidanceLawsY
    mode::ModeGuidanceEnum = ModeGuidance.direct
    ctl::ControlLawsY = ControlLawsY()
    seg::SegmentGuidanceY = SegmentGuidanceY()
end

Modeling.U(::GuidanceLaws) = GuidanceLawsU()
Modeling.Y(::GuidanceLaws) = GuidanceLawsY()


########################### Update Methods #####################################

function Modeling.f_periodic!(::NoScheduling, mdl::Model{<:GuidanceLaws},
                                vehicle::Model{<:Vehicle})

    @unpack mode_req = mdl.u
    @unpack ctl, seg = mdl.submodels

    mode = (is_on_gnd(vehicle) ? ModeGuidance.direct : mode_req)
    mode === ModeGuidance.segment && f_periodic!(seg, ctl, vehicle)

    f_periodic!(ctl, vehicle)
    f_output!(mdl)

end

function Modeling.f_output!(mdl::Model{<:GuidanceLaws})
    @unpack ctl, seg = mdl.submodels
    mdl.y = GuidanceLawsY(mode = mdl.y.mode, ctl = ctl.y, seg = seg.y)
end


function AircraftBase.assign!(systems::Model{<:Systems}, avionics::Model{<:GuidanceLaws})
    AircraftBase.assign!(systems, avionics.ctl)
end


############################# Initialization ###################################

function Modeling.init!(avionics::Model{<:GuidanceLaws}, vehicle::Model{<:Vehicle})
    @unpack ctl, seg = avionics
    Modeling.init!(ctl, vehicle)
    f_output!(avionics)
end

################################## GUI #########################################

# function GUI.draw!(gdc::Model{<:GuidanceLaws},
#                     vehicle::Model{<:Vehicle},
#                     p_open::Ref{Bool} = Ref(true),
#                     label::String = "Cessna172Y Guidance Laws")

#     CImGui.Begin(label, p_open)

#     @cstatic c_ctl=false begin
#         @c CImGui.Checkbox("Flight Control", &c_ctl)
#         c_ctl && @c GUI.draw!(avionics.ctl, vehicle, &c_ctl)
#     end

#     CImGui.End()

# end


end #module