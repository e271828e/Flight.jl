module C172Xv1

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays
using StructTypes

using Flight.FlightCore
using Flight.FlightLib

using ...C172
using ..C172X
using ..C172X.C172XControl: Controller, ControllerY

export Cessna172Xv1


################################################################################
############################### Avionics #######################################

@kwdef struct Avionics{C} <: AbstractAvionics
    ctl::C = Subsampled(Controller(), 2)
end

############################### Update methods #################################

Systems.f_ode!(::System{<:C172Xv1.Avionics}, args...) = nothing
Systems.f_step!(::System{<:C172Xv1.Avionics}, args...) = nothing

@ss_disc C172Xv1.Avionics

function AircraftBase.assign!(components::System{<:C172X.Components},
                          avionics::System{<:C172Xv1.Avionics})
    AircraftBase.assign!(components, avionics.ctl)
end

################################# Trimming #####################################

function Systems.init!(avionics::System{<:C172Xv1.Avionics},
                            vehicle::System{<:C172X.Vehicle})

    Systems.reset!(avionics)
    Systems.init!(avionics.ctl, vehicle)
    Systems.update_y!(avionics)

end

################################## GUI #########################################

function GUI.draw!(avionics::System{<:C172Xv1.Avionics},
                    vehicle::System{<:C172X.Vehicle},
                    p_open::Ref{Bool} = Ref(true),
                    label::String = "Cessna172Xv1 Avionics")

    CImGui.Begin(label, p_open)

    @cstatic c_ctl=false begin
        @c CImGui.Checkbox("Controller", &c_ctl)
        c_ctl && @c GUI.draw!(avionics.ctl, vehicle, &c_ctl)
    end

    CImGui.End()

end


################################################################################
############################# Cessna172Xv1 ###################################

const Cessna172Xv1{K, A} = Cessna172X{K, A} where {
    K <: AbstractKinematicDescriptor, A <: C172Xv1.Avionics}

function Cessna172Xv1(kinematics = WA())
    AircraftBase.Aircraft(C172X.Vehicle(kinematics), C172Xv1.Avionics())
end


################################################################################
############################ Joystick Mappings #################################

function Systems.assign_input!(sys::System{<:Cessna172Xv1},
                                data::JoystickData,
                                mapping::IOMapping)
    Systems.assign_input!(sys.avionics.ctl, data, mapping)
end



end #module