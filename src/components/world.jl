module World

using UnPack

using Flight.FlightCore.Systems
using Flight.FlightCore.Plotting
using Flight.FlightCore.GUI
using Flight.FlightCore.IODevices
using Flight.FlightCore.Joysticks
using Flight.FlightCore.XPC

using Flight.FlightPhysics.Environment

using ..Aircraft

export SimpleWorld

abstract type AbstractWorld <: SystemDefinition end


################################################################################
############################### SimpleWorld #####################################

@kwdef struct SimpleWorld{A <: Aircraft.Template, E <: AbstractEnvironment} <: AbstractWorld
    ac::A
    env::E
end

SimpleWorld(ac::Aircraft.Template) = SimpleWorld(ac, SimpleEnvironment())

struct SimpleWorldY{A, E}
    ac::A
    env::E
end

Systems.init(::SystemY, world::SimpleWorld) = SimpleWorldY(init_y(world.ac), init_y(world.env))

get_aircraft(world::System{<:SimpleWorld}) = world.ac
get_environment(world::System{<:SimpleWorld}) = world.env

Aircraft.init_kinematics!(world::System{<:SimpleWorld}, args...) = init_kinematics!(world.ac, args...)

function Systems.f_ode!(sys::System{<:SimpleWorld})
    @unpack ac, env = sys.subsystems
    f_ode!(env)
    f_ode!(ac, env)
    sys.y = SimpleWorldY(ac.y, env.y)
    return nothing
end

function Systems.f_disc!(sys::System{<:SimpleWorld}, Δt::Real)
    @unpack ac, env = sys.subsystems
    x_mod = false
    x_mod |= f_disc!(env, Δt)
    x_mod |= f_disc!(ac, Δt, env)
    sys.y = SimpleWorldY(ac.y, env.y)
    return x_mod
end

#f_step! can rely on its fallback

function Aircraft.trim!(world::System{<:SimpleWorld}, args...)
    trim!(world.ac, args...)
end

############################### Plotting #######################################

function Plotting.make_plots(th::TimeHistory{<:SimpleWorldY}; kwargs...)

    return OrderedDict(
        :aircraft => make_plots(th.ac; kwargs...),
    )

end

############################ Joystick Mappings #################################

function IODevices.assign!(sys::System{<:SimpleWorld}, joystick::Joystick,
                           mapping::InputMapping)
    IODevices.assign!(sys.ac, joystick, mapping)
end

############################# XPlaneConnect ####################################

XPC.set_position!(xp::XPCDevice, y::SimpleWorldY) = XPC.set_position!(xp, y.ac)


################################### GUI ########################################

function GUI.draw!(sys::System{<:SimpleWorld})

    CImGui.Begin("World")

    show_ac = @cstatic check=false @c CImGui.Checkbox("Aircraft", &check)
    show_env = @cstatic check=false @c CImGui.Checkbox("Environment", &check)

    show_ac && GUI.draw!(sys.ac, "Aircraft")
    show_env && GUI.draw!(sys.env, "Environment")

    CImGui.End()

end


end #module