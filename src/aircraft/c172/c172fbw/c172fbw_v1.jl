module C172FBWv1

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays
using StructTypes

using Flight.FlightCore
using Flight.FlightCore.Utils

using Flight.FlightPhysics
using Flight.FlightComponents

using ...AircraftBase
using ...C172
using ..C172FBW
using ..C172FBW.FlightControl: Controller, ControllerY

export Cessna172FBWv1

################################################################################
############################### Avionics #######################################

@kwdef struct Avionics{C <: Controller} <: AbstractAvionics
    N_fcl::Int = 1
    fcl::C = Controller()
end

@kwdef mutable struct AvionicsS
    n_disc::Int = 0 #number of f_disc! epochs, for scheduling purposes
end

@kwdef struct AvionicsY{C <: ControllerY}
    n_disc::Int = 0
    fcl::C = ControllerY()
end

Systems.S(::Avionics) = AvionicsS()
Systems.Y(::Avionics) = AvionicsY()


############################### Update methods #################################

function Systems.assemble_y!(sys::System{<:C172FBWv1.Avionics})
    @unpack fcl = sys.subsystems
    sys.y = AvionicsY(; n_disc = sys.s.n_disc, fcl = fcl.y)
end


function Systems.reset!(sys::System{<:C172FBWv1.Avionics})
    sys.s.n_disc = 0 #reset scheduler
    foreach(ss -> Systems.reset!(ss), sys.subsystems) #reset algorithms
end


function Systems.f_disc!(avionics::System{<:C172FBWv1.Avionics},
                        vehicle::System{<:C172FBW.Vehicle}, Δt::Real)

    @unpack N_fcl = avionics.constants
    @unpack fcl = avionics.subsystems
    @unpack n_disc = avionics.s

    (n_disc % N_fcl == 0) && f_disc!(fcl, vehicle, N_fcl*Δt)

    avionics.s.n_disc += 1

    assemble_y!(avionics)

end

function AircraftBase.assign!(components::System{<:C172FBW.Components},
                          avionics::System{<:C172FBWv1.Avionics})
    AircraftBase.assign!(components, avionics.fcl)
end


################################# Trimming #####################################

function AircraftBase.trim!(avionics::System{<:C172FBWv1.Avionics},
                            vehicle::System{<:C172FBW.Vehicle})

    @unpack fcl = avionics

    Systems.reset!(avionics)
    trim!(fcl, vehicle)
    assemble_y!(avionics)

end


################################## GUI #########################################

using CImGui: Begin, End, PushItemWidth, PopItemWidth, AlignTextToFramePadding,
        Dummy, SameLine, NewLine, IsItemActive, Separator, Text, Checkbox, RadioButton

function GUI.draw!(avionics::System{<:C172FBWv1.Avionics},
                    vehicle::System{<:C172FBW.Vehicle},
                    p_open::Ref{Bool} = Ref(true),
                    label::String = "Cessna 172 FBWv1 Avionics")

    CImGui.Begin(label, p_open)

    Text(@sprintf("Epoch: %d", avionics.y.n_disc))
    @cstatic c_fcl=false begin
        @c CImGui.Checkbox("Flight Controller", &c_fcl)
        c_fcl && @c GUI.draw!(avionics.fcl, vehicle, &c_fcl)
    end

    CImGui.End()

end


################################################################################
############################# Cessna172FBWv1 ##################################

const Cessna172FBWv1{K, T, A} = C172FBW.Aircraft{K, T, A} where {
    K <: AbstractKinematicDescriptor, T <: AbstractTerrain, A <: C172FBWv1.Avionics}

function Cessna172FBWv1(kinematics = LTF(), terrain = HorizontalTerrain())
    C172FBW.Aircraft(kinematics, terrain, C172FBWv1.Avionics())
end



################################################################################
############################ Joystick Mappings #################################

function IODevices.assign_input!(sys::System{<:Cessna172FBWv1}, joystick::Joystick,
                           mapping::IOMapping)
    IODevices.assign_input!(sys.avionics.fcl, joystick, mapping)
end


end #module