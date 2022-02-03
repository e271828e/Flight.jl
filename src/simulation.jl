module Simulation

using Flight
using Flight.Modeling
using Flight.Plotting
using Flight.Input
using Flight.Output

using Base.Iterators
using OrdinaryDiffEq
using SciMLBase
using GLFW

export SimulationRun

# I <: Union{Nothing, NTuple{N,AbstractInputInterface} where N},
struct SimulationRun{   M <: Model,
                        I <: Vector{<:AbstractInputInterface},
                        O <: Vector{<:AbstractOutputInterface}}
    model::M
    inputs::I
    outputs::O
    realtime::Bool
    plot_enable::Bool
    plot_path::Union{Nothing, String}
    plot_settings::Union{Nothing, NamedTuple}
end

function SimulationRun(;
    model::Model,
    inputs = Vector{AbstractInputInterface}(),
    outputs = Vector{AbstractOutputInterface}(),
    realtime = false,
    plot_enable = false,
    plot_path = joinpath("tmp", "plots"),
    plot_settings = (linewidth=2, margin = 10mm, guidefontsize = 12))
    SimulationRun(model, inputs, outputs, realtime, plot_enable, plot_path, plot_settings)
end


function run!(sim::SimulationRun)

    mdl = sim.model
    output_div = 1

    # #with this we can close the simulation at any time
    # window = GLFW.CreateWindow(640, 480, "GLFW Callback Test")
    # GLFW.MakeContextCurrent(window)

    for out in sim.outputs
        Output.init!(out)
        Output.update!(out, sim.model.sys)
    end

    t_wall = time()
    t_wall_0 = t_wall

    println("Starting simulation...")
    for i in mdl.integrator

        #the integrator steps automatically at the beginning of each iteration

        #retrieve the dt just taken by the integrator
        dt = mdl.dt

        #compute the wall time epoch corresponding to the simulation time epoch
        #we just reached
        t_wall_next = t_wall + dt

        if sim.realtime

            #only do outputs if running in real time
            if mdl.success_iter % output_div == 0
                for output in sim.outputs
                    Output.update!(output, mdl.sys)
                    # mdl.sys.d.airframe.aero.stall |> println
                    # mdl.sys.y.airframe.aero.α |> println
                end
            end

            #busy wait while wall time catches up
            while (time() < t_wall_next) end
            t_wall = t_wall_next
            # println(time()-t_wall_next)

            #only apply inputs if running in real time
            for input in sim.inputs
                Input.update!(input)
                Input.assign!(mdl.sys, input)
            end

        end

        #move this to an AbstractOutputInterface
        # if mdl.success_iter % 40 == 0
        #     pwf(mdl.sys.u.controls)
        # end

        # Swap front and back buffers
        # GLFW.SwapBuffers(window)

        #when GLFW is implemented, build the model with t_end = Inf and use this to break
        # if !GLFW.WindowShouldClose(window)
        #     break
        # end

    end

    println("Simulation finished in $(time() - t_wall_0) seconds")

    if sim.plot_enable
        plots(mdl; save_path = sim.plot_path, sim.plot_settings...)
    end

end

end