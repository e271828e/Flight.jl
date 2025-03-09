module C172Xv2

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays
using StructTypes

using Flight.FlightCore
using Flight.FlightLib

using ...C172
using ..C172X
using ..C172X.C172XControl: Controller, ControllerY
using ..C172X.Navigation: Navigator

export Cessna172Xv2


################################################################################
################################ Sensors #######################################

@kwdef struct Sensors <: SystemDefinition end

function Systems.f_ode!(::System{<:C172Xv2.Sensors}, ::System{<:C172X.Vehicle})
    nothing
end

function Systems.init!(sensors::System{<:C172Xv2.Sensors},
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

Systems.f_ode!(::System{<:C172Xv2.Computing}, ::System{<:C172X.Vehicle}) = nothing

function AircraftBase.assign!(components::System{<:C172X.Components},
                          computing::System{<:C172Xv2.Computing})
    AircraftBase.assign!(components, computing.ctl)
end

function Systems.init!(computing::System{<:C172Xv2.Computing},
                            vehicle::System{<:C172X.Vehicle})

    @unpack ctl, nav = computing

    Systems.reset!(computing)
    Systems.init!(nav, vehicle)
    Systems.init!(ctl, vehicle)
    Systems.update_y!(computing)

end

#################################### GUI #######################################

function GUI.draw!(computing::System{<:C172Xv2.Computing},
                    vehicle::System{<:C172X.Vehicle},
                    p_open::Ref{Bool} = Ref(true),
                    label::String = "Cessna 172Xv2 Computing")

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
                          avionics::System{<:C172Xv2.Avionics})
    AircraftBase.assign!(components, avionics.cmp)
end

function Systems.init!(avionics::System{<:C172Xv2.Avionics},
                            vehicle::System{<:C172X.Vehicle})

    Systems.init!(avionics.sen, vehicle)
    Systems.init!(avionics.cmp, vehicle)

    Systems.update_y!(avionics)

end


################################## GUI #########################################

function GUI.draw!(avionics::System{<:C172Xv2.Avionics},
                    vehicle::System{<:C172X.Vehicle},
                    p_open::Ref{Bool} = Ref(true),
                    label::String = "Cessna 172Xv2 Avionics")

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
############################# Cessna172Xv2 ###################################

const Cessna172Xv2{K, A} = Cessna172X{K, A} where {
    K <: AbstractKinematicDescriptor, A <: C172Xv2.Avionics}

function Cessna172Xv2(kinematics = WA())
    AircraftBase.Aircraft(C172X.Vehicle(kinematics), C172Xv2.Avionics())
end


################################################################################
############################ Joystick Mappings #################################

function Systems.assign_input!(sys::System{<:Cessna172Xv2},
                                data::JoystickData,
                                mapping::IOMapping)
    Systems.assign_input!(sys.avionics.cmp.ctl, data, mapping)
end



end #module