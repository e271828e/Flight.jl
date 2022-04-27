module Simulation

using Flight.Modeling
using Flight.Input
using Flight.Output

using Base.Iterators
using GLFW

export SimulationRun

struct SimulationRun{   M <: Model,
                        I <: Vector{<:AbstractInputInterface},
                        O <: Vector{<:AbstractOutputInterface}}
    model::M
    inputs::I
    outputs::O
    realtime::Bool
end

function SimulationRun(;
    model::Model,
    inputs = Vector{AbstractInputInterface}(),
    outputs = Vector{AbstractOutputInterface}(),
    realtime = false)
    SimulationRun(model, inputs, outputs, realtime)
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
        #     pwf(mdl.sys.u.avionics)
        # end

        # Swap front and back buffers
        # GLFW.SwapBuffers(window)

        #when GLFW is implemented, build the model with t_end = Inf and use this to break
        # if !GLFW.WindowShouldClose(window)
        #     break
        # end

    end

    println("Simulation finished in $(time() - t_wall_0) seconds")

end

end