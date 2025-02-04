module C172FBWv1

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays
using StructTypes

using Flight.FlightCore
using Flight.FlightLib

using ...C172
using ..C172FBW
using ..C172FBW.C172FBWControl: Controller, ControllerY

export Cessna172FBWv1


################################################################################
############################### Avionics #######################################

@kwdef struct Avionics{C} <: AbstractAvionics
    ctl::C = Subsampled(Controller(), 2)
end

############################### Update methods #################################

Systems.f_ode!(::System{<:C172FBWv1.Avionics}, ::System{<:C172FBW.Vehicle}) = nothing

function AircraftBase.assign!(components::System{<:C172FBW.Components},
                          avionics::System{<:C172FBWv1.Avionics})
    AircraftBase.assign!(components, avionics.ctl)
end

################################# Trimming #####################################

function AircraftBase.trim!(avionics::System{<:C172FBWv1.Avionics},
                            vehicle::System{<:C172FBW.Vehicle})

    Systems.reset!(avionics)
    trim!(avionics.ctl, vehicle)
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

    @cstatic c_ctl=false begin
        @c CImGui.Checkbox("Controller", &c_ctl)
        c_ctl && @c GUI.draw!(avionics.ctl, vehicle, &c_ctl)
    end

    CImGui.End()

end


################################################################################
############################# Cessna172FBWv1 ##################################

const Cessna172FBWv1{K, T, A} = C172FBW.Aircraft{K, T, A} where {
    K <: AbstractKinematicDescriptor, T <: AbstractTerrain, A <: C172FBWv1.Avionics}

function Cessna172FBWv1(kinematics = WA(), terrain = HorizontalTerrain())
    C172FBW.Aircraft(kinematics, terrain, C172FBWv1.Avionics())
end



################################################################################
############################ Joystick Mappings #################################

function Systems.assign_input!(sys::System{<:Cessna172FBWv1},
                                data::JoystickData,
                                mapping::IOMapping)

    Systems.assign_input!(sys.avionics.ctl, data, mapping)
end


end #module