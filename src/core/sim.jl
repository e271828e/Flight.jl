module Sim

using UnPack
using StructArrays
using SciMLBase: ODEProblem, u_modified!, init as init_problem
using OrdinaryDiffEq: OrdinaryDiffEqAlgorithm, RK4, ODEIntegrator
using DiffEqCallbacks: SavingCallback, DiscreteCallback, CallbackSet, SavedValues

using ..Systems

import SciMLBase: step!, solve!, reinit!, get_proposed_dt
import OrdinaryDiffEq: add_tstop!

import Flight.Utils: TimeHistory

export Simulation, step!


################################################################################
############################# Simulation #######################################

#in this design, the t and x fields of m.sys behave only as temporary
#storage for f_cont! and f_disc! calls, so we have no guarantees about their
#status after a certain step. the only valid sources for t and x at any
#given moment is the integrator's t and u
struct Simulation{S <: System, I <: ODEIntegrator, L <: SavedValues}
    sys::S
    integrator::I
    log::L
    realtime::Bool

    function Simulation(
        sys::System;
        args_c::Tuple = (), #externally supplied arguments to f_cont!
        args_d::Tuple = (), #externally supplied arguments to f_disc!
        sim_callback::Function = no_sim_callback!, #simulation control callback
        realtime::Bool = false,
        solver::OrdinaryDiffEqAlgorithm = RK4(),
        adaptive::Bool = false, #can be overridden in integrator_kwargs
        dt::Real = 0.02, #can be overridden in integrator_kwargs
        t_start::Real = 0.0,
        t_end::Real = 10.0,
        y_save_on::Bool = true,
        y_saveat::Union{Real, AbstractVector{<:Real}} = dt, #save at every dt
        save_on::Bool = false,
        integrator_kwargs...)

        #save_on is set to false because we are not usually interested in the
        #naked System's state vector. everything we need should be available in
        #the System's output struct saved by the SavingCallback
        saveat_arr = (y_saveat isa Real ? (t_start:y_saveat:t_end) : y_saveat)

        params = (sys = sys, args_c = args_c, args_d = args_d, sim_callback = sim_callback)

        log = SavedValues(Float64, typeof(sys.y))

        cb_step = DiscreteCallback((u, t, integrator)->true, f_step!)
        cb_sim = DiscreteCallback((u, t, integrator)->true, f_sim!)
        cb_save = SavingCallback(f_save, log, saveat = saveat_arr)

        if y_save_on
            cb_set = CallbackSet(cb_step, cb_sim, cb_save)
            # cb_set = CallbackSet(cb_step, cb_save)
        else
            cb_set = CallbackSet(cb_step, cb_sim)
            # cb_set = CallbackSet(cb_step)
        end

        #the current System's x value is used as initial condition. a copy is
        #needed, otherwise it's just a reference and will get overwritten during
        #simulation
        x0 = copy(sys.x)
        problem = ODEProblem{true}(f_ode!, x0, (t_start, t_end), params)
        integrator = init_problem(problem, solver; callback = cb_set, save_on,
                                  adaptive, dt, integrator_kwargs...)

        new{typeof(sys), typeof(integrator), typeof(log)}(
            sys, integrator, log, realtime)
    end

end

Base.getproperty(sim::Simulation, s::Symbol) = getproperty(sim, Val(s))

@generated function Base.getproperty(sim::Simulation, ::Val{S}) where {S}
    if S === :t
        return :(getproperty(getfield(sim, :sys), $(QuoteNode(S)))[])
    elseif S ∈ fieldnames(System)
        return :(getproperty(getfield(sim, :sys), $(QuoteNode(S))))
    else
        return :(getfield(sim, $(QuoteNode(S))))
    end
end

#neither f_ode!, f_step! or f_save should allocate

#integrator ODE function. basically, a wrapper around the System's continuous
#dynamics function. as a side effect, modifies the System's ẋ and y.
function f_ode!(u̇, u, p, t)

    @unpack sys, args_c = p

    sys.x .= u #assign current integrator solution to System continuous state
    sys.t[] = t #ditto for time

    f_cont!(sys, args_c...) #call continuous dynamics, updates sys.ẋ and sys.y

    u̇ .= sys.ẋ #update the integrator's derivative

    return nothing

end

#DiscreteCallback function, called on every integration step. brings the
#System's internal x and y up to date with the last integrator's solution step,
#then calls the System's discrete dynamics. as it is, if the System's y
#depends on x or d, and these are modified by f_disc!, the change will not be
#reflected on y until the following integration step
function f_step!(integrator)

    @unpack t, u, p = integrator
    @unpack sys, args_c, args_d = p

    sys.x .= u #assign the updated integrator's state to the System's continuous state
    sys.t[] = t #ditto for time

    #at this point sys.y and sys.ẋ hold the values from the last solver
    #evaluation of f_cont!, not the one corresponding to the updated x. with x
    #up to date, we can now compute the final sys.y and sys.ẋ for this epoch
    f_cont!(sys, args_c...)

    #with the System's ẋ and y up to date, call the discrete dynamics function
    x_modified = f_disc!(sys, args_d...)

    u .= sys.x #assign the (potentially modified) sys.x back to the integrator

    u_modified!(integrator, x_modified)

end

#DiscreteCallback function, calls the user-specified simulation control callback
#after every integration step
function f_sim!(integrator)

    @unpack sys, sim_callback = integrator.p

    sim_callback(sys.u, sys.t[], sys.y, sys.params)
    u_modified!(integrator, false) #avoids needless performance loss

end

#function signature a simulation control callback must adhere to.
#u is the (mutable) System's control input, which the callback may modify
#t is the (dereferenced) System's t field
#y is the (immutable) System's output
#params is the (immutable) System's params field
no_sim_callback!(u, t::Float64, y, params) = nothing

#SavingCallback function, gets called at the end of each step after f_step!
f_save(x, t, integrator) = deepcopy(integrator.p.sys.y)


############################## Method Extensions ###############################

step!(sim::Simulation, args...) = step!(sim.integrator, args...)

get_proposed_dt(sim::Simulation) = get_proposed_dt(sim.integrator)

add_tstop!(sim::Simulation, t) = add_tstop!(sim.integrator, t)

function reinit!(sim::Simulation, args...; kwargs...)

    #for an ODEIntegrator, the optional args... is simply a new initial
    #condition. if not specified, the original initial condition is used.
    reinit!(sim.integrator, args...; kwargs...)

    #apparently, reinit! calls f_ode!, so all System fields are updated in
    #the process, no need to do it manually

    resize!(sim.log.t, 1)
    resize!(sim.log.saveval, 1)
    return nothing
end

function run!(sim::Simulation)

    @unpack sys, integrator, realtime = sim

    # #with this we can close the simulation at any time
    #if realtime
    # window = GLFW.CreateWindow(640, 480, "GLFW Callback Test")
    # GLFW.MakeContextCurrent(window)
    #end

    if isempty(integrator.opts.tstops)
        println("The simulation has hit its end time, reset it using reinit! ",
                " or add further tstops using add_tstop!")
        return
    end

    t_wall = time()
    t_wall_0 = t_wall

    println("Starting simulation...")
    for _ in integrator

        #integrator steps automatically at the beginning of each iteration

        #retrieve the dt just taken by the integrator
        dt = integrator.dt

        #compute the wall time epoch corresponding to the simulation time epoch
        #we just reached
        t_wall_next = t_wall + dt

        if realtime
            #busy wait while wall time catches up
            while (time() < t_wall_next) end
            t_wall = t_wall_next
            # println(time()-t_wall_next)
        end

        # Swap front and back buffers
        # GLFW.SwapBuffers(window)

        #use this to abort a real time simulation
        # if !GLFW.WindowShouldClose(window)
        #     break
        # end

    end

    println("Simulation finished in $(time() - t_wall_0) seconds")

end

TimeHistory(sim::Simulation) = TimeHistory(sim.log.t, sim.log.saveval)

end #module