module Flight

using Reexport
@reexport using BenchmarkTools

include("core/utils.jl")
include("core/systems.jl")
include("core/sim.jl")
include("core/input.jl")
include("core/output.jl")
include("core/plotting.jl")

include("physics/quaternions.jl")
include("physics/attitude.jl")
include("physics/geodesy.jl")
include("physics/kinematics.jl")
include("physics/rigidbody.jl")

include("components/common/common.jl")

include("components/environment/atmosphere.jl")
include("components/environment/terrain.jl")
include("components/environment/environment.jl")

include("components/aircraft/landinggear.jl")
include("components/aircraft/propellers.jl")
include("components/aircraft/piston.jl")
include("components/aircraft/aircraft.jl")
include("components/aircraft/c172r/c172r.jl")

@reexport using .Utils
@reexport using .Systems
@reexport using .Sim
@reexport using .Input
@reexport using .Output
@reexport using .Plotting

@reexport using .Quaternions
@reexport using .Attitude
@reexport using .Geodesy
@reexport using .Kinematics
@reexport using .RigidBody

@reexport using .Common

@reexport using .Atmosphere
@reexport using .Terrain
@reexport using .Environment

@reexport using .LandingGear
@reexport using .Propellers
@reexport using .Piston
@reexport using .Aircraft
@reexport using .C172R


#force precompilation
# function f_precompile(; save::Bool = true)

#     env = SimpleEnvironment() |> System
#     ac = Cessna172R() |> System
#     kin_init = KinematicInit( h = HOrth(2000))
#     Aircraft.init!(ac, kin_init)

#     sim = Simulation(ac; args_ode = (env, ), t_end = 150, adaptive = true)
#     Sim.run!(sim, verbose = false)
#     plots = make_plots(TimeHistory(sim); Plotting.defaults...)
#     save ? save_plots(plots, save_folder = joinpath("tmp", "nrt_sim_test")) : nothing

# end

#much longer precompilation time, and similar package loading times, but does
#improve first execution times noticeably. not worth it during development.

# f_precompile(save = false)

end
