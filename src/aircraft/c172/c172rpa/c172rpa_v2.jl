module C172RPAv2

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays
using StructTypes

using CImGui: Begin, End, PushItemWidth, PopItemWidth, AlignTextToFramePadding,
    Dummy, SameLine, NewLine, IsItemActive, Separator, Text, Checkbox, RadioButton

using Flight.FlightCore
using Flight.FlightLib

using ...AircraftBase
using ...C172
using ..C172RPA
using ..C172RPA.C172RPAControl: Controller, ControllerY
using ..C172RPA.Navigation: Navigator

export Cessna172RPAv2


################################################################################
################################ Sensors #######################################

@kwdef struct Sensors <: SystemDefinition end

function Systems.f_ode!(::System{<:C172RPAv2.Sensors}, ::System{<:C172RPA.Vehicle})
    nothing
end

function AircraftBase.trim!(sensors::System{<:C172RPAv2.Sensors},
                            vehicle::System{<:C172RPA.Vehicle})

    Systems.reset!(sensors)
    @warn "Sensors trim not implemented"

end


################################################################################
############################## Computing #######################################

@kwdef struct Computing{C, N} <: SystemDefinition
    nav::N = Subsampled(Navigator(), 1)
    ctl::C = Subsampled(Controller(), 2)
end

Systems.f_ode!(::System{<:C172RPAv2.Computing}, ::System{<:C172RPA.Vehicle}) = nothing

function AircraftBase.assign!(components::System{<:C172RPA.Components},
                          computing::System{<:C172RPAv2.Computing})
    AircraftBase.assign!(components, computing.ctl)
end

function AircraftBase.trim!(computing::System{<:C172RPAv2.Computing},
                            vehicle::System{<:C172RPA.Vehicle})

    @unpack ctl, nav = computing

    Systems.reset!(computing)
    trim!(nav, vehicle)
    trim!(ctl, vehicle)
    assemble_y!(computing)

end

#################################### GUI #######################################

function GUI.draw!(computing::System{<:C172RPAv2.Computing},
                    vehicle::System{<:C172RPA.Vehicle},
                    p_open::Ref{Bool} = Ref(true),
                    label::String = "Cessna 172 RPAv2 Computing")

    CImGui.Begin(label, p_open)

    @cstatic c_ctl=false c_nav=false begin
        @c CImGui.Checkbox("Navigator", &c_nav)
        @c CImGui.Checkbox("Controller", &c_ctl)
        c_nav && @c GUI.draw!(computing.nav, vehicle, &c_nav)
        c_ctl && @c GUI.draw!(computing.ctl, vehicle, &c_ctl)
    end

    CImGui.End()

end


################################################################################
############################### Avionics #######################################

@kwdef struct Avionics{S, C} <: AbstractAvionics
    sen::S = Sensors()
    cmp::C = Computing()
end

function AircraftBase.assign!(components::System{<:C172RPA.Components},
                          avionics::System{<:C172RPAv2.Avionics})
    AircraftBase.assign!(components, avionics.cmp)
end

function AircraftBase.trim!(avionics::System{<:C172RPAv2.Avionics},
                            vehicle::System{<:C172RPA.Vehicle})

    trim!(avionics.sen, vehicle)
    trim!(avionics.cmp, vehicle)

    assemble_y!(avionics)

end


################################## GUI #########################################

function GUI.draw!(avionics::System{<:C172RPAv2.Avionics},
                    vehicle::System{<:C172RPA.Vehicle},
                    p_open::Ref{Bool} = Ref(true),
                    label::String = "Cessna 172 RPAv2 Avionics")

    CImGui.Begin(label, p_open)

    @cstatic c_sen=false c_cmp=false begin
        @c CImGui.Checkbox("Sensors", &c_sen)
        @c CImGui.Checkbox("Computing", &c_cmp)
        c_sen && @c GUI.draw!(avionics.sen, vehicle, &c_sen)
        c_cmp && @c GUI.draw!(avionics.cmp, vehicle, &c_cmp)
    end

    CImGui.End()

end


################################################################################
############################# Cessna172RPAv2 ###################################

const Cessna172RPAv2{K, T, A} = C172RPA.Aircraft{K, T, A} where {
    K <: AbstractKinematicDescriptor, T <: AbstractTerrain, A <: C172RPAv2.Avionics}

function Cessna172RPAv2(kinematics = LTF(), terrain = HorizontalTerrain())
    C172RPA.Aircraft(kinematics, terrain, C172RPAv2.Avionics())
end


################################################################################
############################ Joystick Mappings #################################

function Systems.assign_input!(sys::System{<:Cessna172RPAv2},
                                data::JoystickData,
                                mapping::IOMapping)
    Systems.assign_input!(sys.avionics.cmp.ctl, data, mapping)
end



end #module