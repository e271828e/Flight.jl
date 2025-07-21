module World

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays

using Flight.FlightCore
using Flight.FlightLib
using Flight.FlightLib.AircraftBase: Aircraft

export AbstractWorld, SimpleWorld


################################################################################
################################### World ######################################

abstract type AbstractWorld <: ModelDefinition end

################################################################################
############################## SimpleWorld #####################################

@kwdef struct SimpleWorld{C <: Aircraft, A <: AbstractAtmosphere, T <: AbstractTerrain} <: AbstractWorld
    ac::C = Aircraft()
    atm::A = SimpleAtmosphere()
    trn::T = HorizontalTerrain()
end

function Modeling.f_ode!(world::Model{<:SimpleWorld})
    @unpack ac, atm, trn = world.submodels
    f_ode!(atm)
    f_ode!(trn)
    f_ode!(ac, atm, trn)
    update_output!(world)
end

function Modeling.f_disc!(::NoScheduling, world::Model{<:SimpleWorld})
    @unpack ac, atm, trn = world.submodels
    f_disc!(atm)
    f_disc!(trn)
    f_disc!(ac, atm, trn)
    update_output!(world)
end

function Modeling.f_step!(world::Model{<:SimpleWorld})
    @unpack ac, atm, trn = world.submodels
    f_step!(atm)
    f_step!(trn)
    f_step!(ac, atm, trn)
end

function Modeling.init!( world::Model{<:SimpleWorld},
                        init::Union{<:AircraftBase.VehicleInitializer,
                                    <:AircraftBase.AbstractTrimParameters})
    @unpack ac, atm, trn = world.submodels
    Modeling.init!(atm)
    Modeling.init!(trn)
    Modeling.init!(ac, init, atm, trn)
    update_output!(world) #!
end

################################################################################
############################### XPlane12Control #################################

function IODevices.extract_output(world::Model{<:SimpleWorld},
                                mapping::XPlane12ControlMapping)
    IODevices.extract_output(world.ac, mapping)
end


################################################################################
############################ Joystick Mappings #################################

function IODevices.assign_input!(world::Model{<:SimpleWorld},
                                mapping::IOMapping,
                                data::AbstractJoystickData)
    IODevices.assign_input!(world.ac, mapping, data)
end

################################################################################
#################################### GUI #######################################

function GUI.draw!(world::Model{<:SimpleWorld};
                    p_open::Ref{Bool} = Ref(true), label::String = "Simple World")

    CImGui.Begin(label, p_open)

    @cstatic c_ac=false c_atm=false c_trn=false begin
        @c CImGui.Checkbox("Aircraft", &c_ac)
        c_ac && @c GUI.draw!(world.ac, &c_ac)
        @c CImGui.Checkbox("Atmosphere", &c_atm)
        c_atm && @c GUI.draw!(world.atm, &c_atm)
        @c CImGui.Checkbox("Terrain", &c_trn)
        c_trn && @c GUI.draw!(world.trn, &c_trn)
    end

    CImGui.End()

end

end
