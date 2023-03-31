module World

using UnPack
using CImGui, CImGui.CSyntax, CImGui.CSyntax.CStatic

using Flight.FlightCore
using Flight.FlightPhysics
using ..Aircraft
using ..C172RDirect

export SimpleWorld

abstract type AbstractWorld <: Component end


################################################################################
############################### SimpleWorld #####################################

Base.@kwdef struct SimpleWorld{A <: AircraftTemplate, E <: AbstractEnvironment} <: AbstractWorld
    ac::A = Cessna172R()
    env::E = SimpleEnvironment()
end

struct SimpleWorldU{A, E}
    ac::A
    env::E
end

struct SimpleWorldY{A, E}
    ac::A
    env::E
end

Systems.init(::SystemU, world::SimpleWorld) = SimpleWorldU(init_u(world.ac), init_u(world.env))
Systems.init(::SystemY, world::SimpleWorld) = SimpleWorldY(init_y(world.ac), init_y(world.env))

get_aircraft(world::System{<:SimpleWorld}) = world.ac
get_environment(world::System{<:SimpleWorld}) = world.env

Aircraft.init_kinematics!(world::System{<:SimpleWorld}, args...) = init_kinematics!(world.ac, args...)

function Systems.f_ode!(sys::System{<:SimpleWorld})
    @unpack ac, env = sys
    f_ode!(env)
    f_ode!(ac, env)
    sys.y = SimpleWorldY(ac.y, env.y)
    return nothing
end

function Systems.f_disc!(sys::System{<:SimpleWorld}, Δt::Real)
    @unpack ac, env = sys
    x_mod = false
    x_mod |= f_disc!(env, Δt)
    x_mod |= f_disc!(ac, Δt, env)
    sys.y = SimpleWorldY(ac.y, env.y)
    return x_mod
end

#f_step! can rely on its fallback

############################### Plotting #######################################

function Plotting.make_plots(th::TimeHistory{<:SimpleWorldY}; kwargs...)

    return OrderedDict(
        :aircraft => make_plots(th.ac; kwargs...),
    )

end

############################ Joystick Mappings #################################

function IODevices.assign!(u::SimpleWorldU, joystick::Joystick{XBoxControllerID},
                           mapping::InputMapping)
    IODevices.assign!(u.ac, joystick, mapping)
end

############################# XPlaneConnect ####################################

XPlane.set_position!(xp::XPConnect, y::SimpleWorldY) = XPlane.set_position!(xp, y.ac)


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