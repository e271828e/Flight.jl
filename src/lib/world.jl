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
    aircraft::C = Aircraft()
    atmosphere::A = SimpleAtmosphere()
    terrain::T = HorizontalTerrain()
end

function Modeling.f_ode!(world::Model{<:SimpleWorld})
    @unpack aircraft, atmosphere, terrain = world
    f_ode!(atmosphere)
    f_ode!(terrain)
    f_ode!(aircraft, atmosphere, terrain)
    update_output!(world)
end

function Modeling.f_disc!(::NoScheduling, world::Model{<:SimpleWorld})
    @unpack aircraft, atmosphere, terrain = world
    f_disc!(atmosphere)
    f_disc!(terrain)
    f_disc!(aircraft, atmosphere, terrain)
    update_output!(world)
end

function Modeling.f_step!(world::Model{<:SimpleWorld})
    @unpack aircraft, atmosphere, terrain = world
    f_step!(atmosphere)
    f_step!(terrain)
    f_step!(aircraft, atmosphere, terrain)
end

function Modeling.init!( world::Model{<:SimpleWorld},
                        init::Union{<:AircraftBase.VehicleInitializer,
                                    <:AircraftBase.AbstractTrimParameters})
    @unpack aircraft, atmosphere, terrain = world
    Modeling.init!(atmosphere)
    Modeling.init!(terrain)
    Modeling.init!(aircraft, init, atmosphere, terrain)
    update_output!(world) #!
end

################################################################################
############################### XPlane12Control #################################

function IODevices.extract_output(world::Model{<:SimpleWorld},
                                mapping::XPlane12ControlMapping)
    IODevices.extract_output(world.aircraft, mapping)
end


################################################################################
############################ Joystick Mappings #################################

function IODevices.assign_input!(world::Model{<:SimpleWorld},
                                mapping::IOMapping,
                                data::AbstractJoystickData)
    IODevices.assign_input!(world.aircraft, mapping, data)
end

################################################################################
#################################### GUI #######################################

function GUI.draw!(world::Model{<:SimpleWorld};
                    p_open::Ref{Bool} = Ref(true), label::String = "Simple World")

    CImGui.Begin(label, p_open)

    @cstatic c_ac=false c_atm=false c_trn=false begin
        @c CImGui.Checkbox("Aircraft", &c_ac)
        c_ac && @c GUI.draw!(world.aircraft, &c_ac)
        @c CImGui.Checkbox("Atmosphere", &c_atm)
        c_atm && @c GUI.draw!(world.atmosphere, &c_atm)
        @c CImGui.Checkbox("Terrain", &c_trn)
        c_trn && @c GUI.draw!(world.terrain, &c_trn)
    end

    CImGui.End()

end

end
