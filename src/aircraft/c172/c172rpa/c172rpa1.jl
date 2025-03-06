module C172RPA1

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays
using StructTypes

using Flight.FlightCore
using Flight.FlightLib

using ...C172
using ..C172RPA
using ..C172RPA.C172RPAControl: Controller, ControllerY

export Cessna172RPA1


################################################################################
############################### Avionics #######################################

@kwdef struct Avionics{C} <: AbstractAvionics
    ctl::C = Subsampled(Controller(), 2)
end

############################### Update methods #################################

Systems.f_ode!(::System{<:C172RPA1.Avionics}, ::System{<:C172RPA.Vehicle}) = nothing

function AircraftBase.assign!(components::System{<:C172RPA.Components},
                          avionics::System{<:C172RPA1.Avionics})
    AircraftBase.assign!(components, avionics.ctl)
end

################################# Trimming #####################################

function Systems.init!(avionics::System{<:C172RPA1.Avionics},
                            vehicle::System{<:C172RPA.Vehicle})

    Systems.reset!(avionics)
    Systems.init!(avionics.ctl, vehicle)
    Systems.update_y!(avionics)

end

################################## GUI #########################################

function GUI.draw!(avionics::System{<:C172RPA1.Avionics},
                    vehicle::System{<:C172RPA.Vehicle},
                    p_open::Ref{Bool} = Ref(true),
                    label::String = "Cessna172X1 Avionics")

    CImGui.Begin(label, p_open)

    @cstatic c_ctl=false begin
        @c CImGui.Checkbox("Controller", &c_ctl)
        c_ctl && @c GUI.draw!(avionics.ctl, vehicle, &c_ctl)
    end

    CImGui.End()

end


################################################################################
############################# Cessna172RPA1 ###################################

const Cessna172RPA1{K, T, A} = Cessna172RPA{K, T, A} where {
    K <: AbstractKinematicDescriptor, T <: AbstractTerrain, A <: C172RPA1.Avionics}

function Cessna172RPA1(kinematics = WA(), terrain = HorizontalTerrain())
    AircraftBase.Aircraft(C172RPA.Vehicle(kinematics, terrain), C172RPA1.Avionics())
end


################################################################################
############################ Joystick Mappings #################################

function Systems.assign_input!(sys::System{<:Cessna172RPA1},
                                data::JoystickData,
                                mapping::IOMapping)
    Systems.assign_input!(sys.avionics.ctl, data, mapping)
end



end #module