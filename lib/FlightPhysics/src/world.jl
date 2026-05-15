module World

using FlightCore
using FlightCore.GUI
using ..Atmosphere
using ..Terrain
using ..AircraftBase: AircraftBase, Aircraft

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
    (; aircraft, atmosphere, terrain) = world
    f_ode!(atmosphere)
    f_ode!(terrain)
    f_ode!(aircraft, atmosphere, terrain)
    f_output!(world)
end

function Modeling.f_step!(world::Model{<:SimpleWorld})
    (; aircraft, atmosphere, terrain) = world
    f_step!(atmosphere)
    f_step!(terrain)
    f_step!(aircraft, atmosphere, terrain)
end

function Modeling.f_periodic!(::Unconditional, world::Model{<:SimpleWorld})
    (; aircraft, atmosphere, terrain) = world
    f_periodic!(atmosphere)
    f_periodic!(terrain)
    f_periodic!(aircraft, atmosphere, terrain)
    f_output!(world)
end

function Modeling.f_init!( world::Model{<:SimpleWorld},
                        init::Union{<:AircraftBase.VehicleInitializer,
                                    <:AircraftBase.AbstractTrimParameters})
    (; aircraft, atmosphere, terrain) = world
    f_init!(atmosphere)
    f_init!(terrain)
    f_init!(aircraft, init, atmosphere, terrain)
    f_output!(world) #!
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

    BeginWindow(label, p_open)

    @cstatic c_ac=false c_atm=false c_trn=false begin
        @c Checkbox("Aircraft", &c_ac)
        c_ac && @c GUI.draw!(world.aircraft, &c_ac)
        @c Checkbox("Atmosphere", &c_atm)
        c_atm && @c GUI.draw!(world.atmosphere, &c_atm)
        @c Checkbox("Terrain", &c_trn)
        c_trn && @c GUI.draw!(world.terrain, &c_trn)
    end

    EndWindow()

end

end
