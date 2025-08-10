module C172Xv2

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays
using StructTypes

using Flight.FlightCore
using Flight.FlightLib

using ...C172
using ..C172X
using ..C172X.C172XControl: Controller, ControllerY

export Cessna172Xv2


################################################################################
################################ Sensors #######################################

@kwdef struct Sensors <: ModelDefinition end

@no_ode C172Xv2.Sensors
@no_step C172Xv2.Sensors
@ss_disc C172Xv2.Sensors

function Modeling.init!(sensors::System{<:C172Xv2.Sensors},
                            vehicle::System{<:C172X.Vehicle})

    Modeling.reset!(sensors)
    @warn "Sensors init! not implemented"

end


################################################################################
################################## Navigator ###################################

@kwdef struct Navigator <: ModelDefinition end
@no_ode C172Xv2.Navigator
@no_step C172Xv2.Navigator

function Modeling.f_disc!(::NoScheduling, nav::System{<:Navigator},
                        ::System{<:C172X.Vehicle})

end

function Modeling.init!(nav::System{<:Navigator},
                            vehicle::System{<:C172X.Vehicle})

    @warn "Navigator init! not implemented"

end


################################################################################
############################## Computing #######################################

@kwdef struct Computing{C, N} <: ModelDefinition
    nav::N = Subsampled(Navigator(), 1)
    ctl::C = Subsampled(Controller(), 2)
end

@no_ode C172Xv2.Computing
@no_step C172Xv2.Computing
@ss_disc C172Xv2.Computing

function AircraftBase.assign!(components::System{<:C172X.Components},
                          computing::System{<:C172Xv2.Computing})
    AircraftBase.assign!(components, computing.ctl)
end

function Modeling.init!(computing::System{<:C172Xv2.Computing},
                            vehicle::System{<:C172X.Vehicle})

    @unpack ctl, nav = computing

    Modeling.reset!(computing)
    Modeling.init!(nav, vehicle)
    Modeling.init!(ctl, vehicle)
    Modeling.update_y!(computing)

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

@no_ode C172Xv2.Avionics
@no_step C172Xv2.Avionics
@ss_disc C172Xv2.Avionics

function AircraftBase.assign!(components::System{<:C172X.Components},
                          avionics::System{<:C172Xv2.Avionics})
    AircraftBase.assign!(components, avionics.cmp)
end

function Modeling.init!(avionics::System{<:C172Xv2.Avionics},
                            vehicle::System{<:C172X.Vehicle})

    Modeling.init!(avionics.sen, vehicle)
    Modeling.init!(avionics.cmp, vehicle)

    Modeling.update_y!(avionics)

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

function Modeling.assign_input!(sys::System{<:Cessna172Xv2},
                                data::JoystickData,
                                mapping::IOMapping)
    Modeling.assign_input!(sys.avionics.cmp.ctl, data, mapping)
end



end #module