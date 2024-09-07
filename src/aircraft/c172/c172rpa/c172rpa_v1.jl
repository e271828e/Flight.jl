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
using ..C172RPA.FlightControl: Controller, ControllerY
using ..C172RPA.Navigation: Navigator, NavigatorY

export Cessna172RPAv1


################################################################################
############################### Avionics #######################################

@kwdef struct Avionics{C, N} <: AbstractAvionics
    ctl::C = Subsampled(Controller(), 2)
    nav::N = Subsampled(Navigator(), 1)
end

############################### Update methods #################################

function Systems.f_disc!(::NoScheduling, avionics::System{<:C172RPAv1.Avionics},
                        vehicle::System{<:C172RPA.Vehicle})

    @unpack ctl, nav = avionics.subsystems

    f_disc!(ctl, vehicle)
    f_disc!(nav)

    assemble_y!(avionics)

end

function AircraftBase.assign!(components::System{<:C172RPA.Components},
                          avionics::System{<:C172RPAv1.Avionics})
    AircraftBase.assign!(components, avionics.ctl)
end


################################# Trimming #####################################

function AircraftBase.trim!(avionics::System{<:C172RPAv1.Avionics},
                            vehicle::System{<:C172RPA.Vehicle})

    @unpack ctl, nav = avionics

    Systems.reset!(avionics)
    trim!(ctl, vehicle)
    # trim!(nav, vehicle)
    assemble_y!(avionics)

end


################################## GUI #########################################

using CImGui: Begin, End, PushItemWidth, PopItemWidth, AlignTextToFramePadding,
        Dummy, SameLine, NewLine, IsItemActive, Separator, Text, Checkbox, RadioButton

function GUI.draw!(avionics::System{<:C172RPAv1.Avionics},
                    vehicle::System{<:C172RPA.Vehicle},
                    p_open::Ref{Bool} = Ref(true),
                    label::String = "Cessna 172 RPAv1 Avionics")

    CImGui.Begin(label, p_open)

    Text(@sprintf("Epoch: %d", avionics.y.n_disc))
    @cstatic c_ctl=false begin
        @c CImGui.Checkbox("Flight Controller", &c_ctl)
        c_ctl && @c GUI.draw!(avionics.ctl, vehicle, &c_ctl)
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

function Systems.assign_input!(sys::System{<:Cessna172RPAv1},
                                data::JoystickData,
                                mapping::IOMapping)
    Systems.assign_input!(sys.avionics.ctl, data, mapping)
end



end #module