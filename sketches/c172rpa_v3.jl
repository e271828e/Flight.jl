module C172Xv3

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays
using StructTypes

using CImGui: Begin, End, PushItemWidth, PopItemWidth, AlignTextToFramePadding,
    Dummy, SameLine, NewLine, IsItemActive, Separator, Text, Checkbox, RadioButton

using Flight.FlightCore
using Flight.FlightLib
using Flight.FlightAircraft

using Flight.FlightAircraft.C172X.C172XControl: Controller, ControllerY
using Flight.FlightAircraft.C172X.Navigation: Navigator

export Cessna172Xv3


################################################################################
################################ Sensors #######################################

@kwdef struct Sensors <: SystemDefinition end

function Systems.f_ode!(::System{<:C172Xv3.Sensors}, ::System{<:C172X.Vehicle})
    nothing
end

function Systems.init!(sensors::System{<:C172Xv3.Sensors},
                            vehicle::System{<:C172X.Vehicle})

    Systems.reset!(sensors)
    @warn "Sensors init! not implemented"

end


################################################################################
############################## Computing #######################################

@kwdef struct Computing{C, N} <: SystemDefinition
    nav::N = Subsampled(Navigator(), 1)
    ctl::C = Subsampled(Controller(), 2)
end

Systems.f_ode!(::System{<:C172Xv3.Computing}, ::System{<:C172X.Vehicle}) = nothing

function AircraftBase.assign!(components::System{<:C172X.Components},
                          computing::System{<:C172Xv3.Computing})
    AircraftBase.assign!(components, computing.ctl)
end

function Systems.init!(computing::System{<:C172Xv3.Computing},
                            vehicle::System{<:C172X.Vehicle})

    @unpack ctl, nav = computing

    Systems.reset!(computing)
    Systems.init!(nav, vehicle)
    Systems.init!(ctl, vehicle)
    Systems.update_y!(computing)

end

#################################### GUI #######################################

function GUI.draw!(computing::System{<:C172Xv3.Computing},
                    vehicle::System{<:C172X.Vehicle},
                    p_open::Ref{Bool} = Ref(true),
                    label::String = "Cessna 172Xv3 Computing")

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

function AircraftBase.assign!(components::System{<:C172X.Components},
                          avionics::System{<:C172Xv3.Avionics})
    AircraftBase.assign!(components, avionics.cmp)
end

function Systems.init!(avionics::System{<:C172Xv3.Avionics},
                            vehicle::System{<:C172X.Vehicle})

    Systems.init!(avionics.sen, vehicle)
    Systems.init!(avionics.cmp, vehicle)

    Systems.update_y!(avionics)

end


################################## GUI #########################################

function GUI.draw!(avionics::System{<:C172Xv3.Avionics},
                    vehicle::System{<:C172X.Vehicle},
                    p_open::Ref{Bool} = Ref(true),
                    label::String = "Cessna 172Xv3 Avionics")

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
############################# Cessna172Xv3 ###################################

const Cessna172Xv3{K, T, A} = C172X.Aircraft{K, T, A} where {
    K <: AbstractKinematicDescriptor, T <: AbstractTerrain, A <: C172Xv3.Avionics}

function Cessna172Xv3(kinematics = WA(), terrain = HorizontalTerrain())
    C172X.Aircraft(kinematics, terrain, C172Xv3.Avionics())
end


################################################################################
############################ Joystick Mappings #################################

function Systems.assign_input!(sys::System{<:Cessna172Xv3},
                                data::JoystickData,
                                mapping::IOMapping)
    Systems.assign_input!(sys.avionics.cmp.ctl, data, mapping)
end



end #module