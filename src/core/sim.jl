module Sim

using UnPack
using StructArrays
using SciMLBase: ODEProblem, u_modified!, init as init_integrator
using OrdinaryDiffEq: OrdinaryDiffEqAlgorithm, ODEIntegrator, RK4
using DiffEqCallbacks: SavingCallback, DiscreteCallback, PeriodicCallback, CallbackSet, SavedValues
using Flight.Utils

using ..Systems

import SciMLBase: step!, reinit!, get_proposed_dt
import OrdinaryDiffEq: add_tstop!

import Flight.Utils: TimeHistory

export Simulation, step!, reinit!, get_proposed_dt, add_tstop!


################################################################################
############################# Simulation #######################################

#in this design, the t and x fields of m.sys behave only as temporary
#storage for f_ode! and f_step! calls, so we have no guarantees about their
#status after a certain step. the only valid sources for t and x at any
#given moment is the integrator's t and u
struct Simulation{S <: System, I <: ODEIntegrator, L <: SavedValues}
    sys::S
    integrator::I
    log::L
    realtime::Bool

    function Simulation(
        sys::System;
        args_ode::Tuple = (), #externally supplied arguments to System's f_ode!
        args_step::Tuple = (), #externally supplied arguments to System's f_step!
        args_disc::Tuple = (), #externally supplied arguments to System's f_disc!
        sys_init!::Function = no_sys_init!, #System initialization function
        sys_io!::Function = no_sys_io!, #System I/O function
        realtime::Bool = false,
        algorithm::OrdinaryDiffEqAlgorithm = RK4(),
        adaptive::Bool = false,
        dt::Real = 0.02, #continuous dynamics integration step
        Δt::Real = Inf, #discrete dynamics execution period
        t_start::Real = 0.0,
        t_end::Real = 10.0,
        save_on::Bool = true,
        saveat::Union{Real, AbstractVector{<:Real}} = Float64[], #defers to save_everystep
        save_everystep::Bool = isempty(saveat),
        sys_init_kwargs = NamedTuple())

        sys_init!(sys; sys_init_kwargs...)

        params = (sys = sys, sys_init! = sys_init!, sys_io! = sys_io!, Δt = Δt,
                  args_ode = args_ode, args_step = args_step, args_disc = args_disc)

        cb_step = DiscreteCallback((u, t, integrator)->true, f_cb_step!)
        cb_disc = PeriodicCallback(f_cb_disc!, Δt)
        cb_io = DiscreteCallback((u, t, integrator)->true, f_cb_io!)

        f_ode!(sys, args_ode...) #update y so that the first log entry is correct
        log = SavedValues(Float64, typeof(sys.y))
        saveat_arr = (saveat isa Real ? (t_start:saveat:t_end) : saveat)
        cb_save = SavingCallback(f_cb_save, log; saveat = saveat_arr, save_everystep)

        if save_on
            cb_set = CallbackSet(cb_step, cb_disc, cb_io, cb_save)
        else
            cb_set = CallbackSet(cb_step, cb_disc, cb_io)
        end

        #the current System's x value is used as initial condition. a copy is
        #needed, otherwise it's just a reference and will get overwritten during
        #simulation
        x0 = copy(sys.x)
        problem = ODEProblem{true}(f_ode_wrapper!, x0, (t_start, t_end), params)
        integrator = init_integrator(
            problem, algorithm; callback = cb_set, save_on = false, adaptive, dt)
        #save_on is set to false because we are not usually interested in the
        #naked System's state vector. everything we need should be made
        #available in the System's output struct saved by the SavingCallback

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

#wrapper around the System's continuous dynamics function f_ode! it modifies the
#System's ẋ and y as a side effect
function f_ode_wrapper!(u̇, u, p, t)

    @unpack sys, args_ode = p

    sys.x .= u #assign current integrator solution to System continuous state
    sys.t[] = t #likewise for time

    f_ode!(sys, args_ode...) #call continuous dynamics, updates sys.ẋ and sys.y

    u̇ .= sys.ẋ #update the integrator's derivative

    return nothing

end

#DiscreteCallback function, called on every integration step. brings the
#System's internal x and y up to date with the last integrator's solution step,
#then calls the System's own f_step!. if the System's y depends on x or s, and
#these are modified by f_step!, the change will not be reflected on y until
#f_ode! is executed again on the following integration step
function f_cb_step!(integrator)

    @unpack t, u, p = integrator
    @unpack sys, args_ode, args_step = p

    sys.x .= u #assign the updated integrator's state to the System's continuous state
    sys.t[] = t #ditto for time

    #at this point sys.y and sys.ẋ hold the values from the last algorithm
    #evaluation of f_ode!, not the one corresponding to the updated x. with x
    #up to date, we can now compute the final sys.y and sys.ẋ for this epoch
    f_ode!(sys, args_ode...)

    #with the System's ẋ and y up to date, call the discrete dynamics function
    x_modified = f_step!(sys, args_step...)

    u .= sys.x #assign the (potentially modified) sys.x back to the integrator

    u_modified!(integrator, x_modified)

end

#PeriodicCallback function, calls the System's discrete dynamics update
#with the period Δt given to the Simulation constructor
function f_cb_disc!(integrator)

    @unpack sys, Δt, args_disc = integrator.p

    x_modified = f_disc!(sys, Δt, args_disc...)

    u_modified!(integrator, x_modified)

end

#DiscreteCallback function, calls the user-specified System I/O function after
#every integration step
function f_cb_io!(integrator)

    @unpack sys, sys_io! = integrator.p

    sys_io!(sys.u, sys.y, sys.t[], sys.params)

    #a System control function will never modify the System's continuous state,
    #so we may as well tell the integrator to avoid any performance hit
    u_modified!(integrator, false)

end

#SavingCallback function, gets called at the end of each step after f_disc!
#and/or f_step!
f_cb_save(x, t, integrator) = deepcopy(integrator.p.sys.y)

#function signature a System initialization function must adhere to.
#x: System's continuous state, which the function may modify
#u: System's control input, which the function may modify
#s: System's discrete state, which the function may modify
#t: (dereferenced) System's t field
#params: (immutable) System's params field
no_sys_init!(sys::System; kwargs...) = nothing

#function signature a System I/O function must adhere to.
#u: (mutable) System's control input, which the function may modify
#y: (immutable) System's output
#t: (dereferenced) System's t field
#params: (immutable) System's params field
no_sys_io!(u, y, t::Float64, params) = nothing


############################## Method Extensions ###############################

step!(sim::Simulation, args...) = step!(sim.integrator, args...)

get_proposed_dt(sim::Simulation) = get_proposed_dt(sim.integrator)

add_tstop!(sim::Simulation, t) = add_tstop!(sim.integrator, t)

#in order to be actually useful, this would need to restore also u and s to
#their original values, which is non-trivial, because it may break the coupling
#between a parent system and its subsystems
function reinit!(sim::Simulation; sys_init_kwargs...)

    @unpack integrator, log = sim
    @unpack p = integrator

    #initialize the System's x, u and s
    println(sys_init_kwargs)
    p.sys_init!(p.sys; sys_init_kwargs...)

    #initialize the ODEIntegrator with the System's initial x. ODEIntegrator's
    #reinit! calls f_ode_wrapper!, so the System's ẋ and y are updated in the
    #process, no need to do it explicitly
    reinit!(integrator, p.sys.x)

    resize!(log.t, 1)
    resize!(log.saveval, 1)
    return nothing
end

function run!(sim::Simulation; verbose = false)

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

    verbose && println("Starting simulation...")
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

    verbose && println("Simulation finished in $(time() - t_wall_0) seconds")
    return nothing

end

TimeHistory(sim::Simulation) = TimeHistory(sim.log.t, sim.log.saveval)

end #module