module C172RPAv1

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays
using StructTypes

using Flight.FlightCore
using Flight.FlightCore.Utils

using Flight.FlightPhysics
using Flight.FlightComponents

using ...AircraftBase
using ...C172
using ..C172RPA
using ..C172RPA.FlightControl: FlightController, FlightControllerY

export Cessna172RPAv1


################################################################################
############################### Avionics #######################################

@kwdef struct Avionics{C <: FlightController} <: AbstractAvionics
    N_fcl::Int = 1
    fcl::C = FlightController()
end

@kwdef struct AvionicsY
    fcl::FlightControllerY = FlightControllerY()
end

Systems.Y(::Avionics) = AvionicsY()

################################## Update ######################################

function Systems.f_disc!(avionics::System{<:C172RPAv1.Avionics},
                        vehicle::System{<:C172RPA.Vehicle}, Δt::Real)

    #CAUTION: When scheduling is added, ensure call scheduling is consistent
    #with the Δt we are passing
    f_disc!(avionics.fcl, vehicle, avionics.constants.N_fcl*Δt)
end

function AircraftBase.assign!(components::System{<:C172RPA.Components},
                          avionics::System{<:C172RPAv1.Avionics})
    AircraftBase.assign!(components, avionics.fcl)
end


################################# Trimming #####################################

function AircraftBase.trim!(avionics::System{<:C172RPAv1.Avionics},
                            vehicle::System{<:C172RPA.Vehicle})
    trim!(avionics.fcl, vehicle)
    avionics.y = AvionicsY(avionics.fcl.y)
end


################################## GUI #########################################

function GUI.draw!(avionics::System{<:C172RPAv1.Avionics},
                    vehicle::System{<:C172RPA.Vehicle},
                    p_open::Ref{Bool} = Ref(true),
                    label::String = "Cessna 172 RPAv1 Avionics")

    CImGui.Begin(label, p_open)

    @cstatic c_fcl=false begin
        @c CImGui.Checkbox("Flight Control", &c_fcl)
        c_fcl && @c GUI.draw!(avionics.fcl, vehicle, &c_fcl)
    end

    CImGui.End()

end


################################################################################
############################# Cessna172RPAv1 ###################################

const Cessna172RPAv1{K, T, A} = C172RPA.Aircraft{K, T, A} where {
    K <: AbstractKinematicDescriptor, T <: AbstractTerrain, A <: C172RPAv1.Avionics}

function Cessna172RPAv1(kinematics = LTF(), terrain = HorizontalTerrain())
    C172RPA.Aircraft(kinematics, terrain, C172RPAv1.Avionics())
end


################################################################################
############################ Joystick Mappings #################################

function IODevices.assign_input!(sys::System{<:Cessna172RPAv1}, joystick::Joystick,
                           mapping::InputMapping)
    IODevices.assign_input!(sys.avionics.fcl, joystick, mapping)
end



end #module