module Sim

using UnPack, Reexport, StructArrays
using OrdinaryDiffEq: OrdinaryDiffEq, OrdinaryDiffEqAlgorithm, ODEProblem,
                      ODEIntegrator, Heun, RK4, u_modified!, init as init_integrator
using DiffEqCallbacks: SavingCallback, DiscreteCallback, PeriodicCallback,
                       CallbackSet, SavedValues
using RecursiveArrayTools
using Logging

using ..Systems
using ..IODevices
using ..GUI

@reexport using OrdinaryDiffEq: step!, reinit!, add_tstop!, get_proposed_dt
export Simulation, enable_gui!, disable_gui!, attach!
export TimeSeries, get_time, get_data, get_components, get_child_names


################################################################################
################################# SimInfo ########################################

@kwdef struct SimInfo
    algorithm::String = "Undefined"
    t_start::Float64 = 0.0
    t_end::Float64 = 0.0
    dt::Float64 = 0.0 #last time step
    iter::Int64 = 0 #total iterations
    t::Float64 = 0.0 #simulation time
    τ::Float64 = 0.0 #wall-clock time
end


################################################################################
############################### SimData #########################################

@kwdef struct SimData{Y}
    t::Float64
    y::Y #System's y
end

################################################################################
############################# Simulation #######################################

struct Simulation{S <: System, I <: ODEIntegrator, L <: SavedValues, Y}
    sys::S
    integrator::I
    log::L
    gui::Renderer
    info::Ref{SimInfo}
    inputs::Vector{SimInput}
    outputs::Vector{SimOutput}
    started::Base.Event #signals that simulation execution has started
    running::ReentrantLock #while locked it signals that simulation execution is in progress
    stepping::ReentrantLock #must be acquired by simulation inputs to modify the system

    #set sync = 0 so that SwapBuffers returns immediately and doesn't slow down
    #the sim loop; for paced simulations, scheduling is taken care of
    function Simulation(
        sys::System;
        args_ode::Tuple = (), #externally supplied arguments to System's f_ode!
        args_step::Tuple = (), #externally supplied arguments to System's f_step!
        args_disc::Tuple = (), #externally supplied arguments to System's f_disc!
        sys_init!::Function = Systems.init!, #default System initialization function
        user_callback!::Function = user_callback!, #user-specified callback
        algorithm::OrdinaryDiffEqAlgorithm = RK4(),
        adaptive::Bool = false,
        dt::Real = 0.02, #continuous dynamics integration step
        Δt::Real = 0.02, #discrete dynamics execution period (do not set to Inf!)
        t_start::Real = 0.0,
        t_end::Real = 10.0,
        save_on::Bool = true,
        saveat::Union{Real, AbstractVector{<:Real}} = Float64[], #defers to save_everystep
        save_everystep::Bool = isempty(saveat),
        save_start::Bool = false, #initial System's outputs might not be up to date
        save_end::Bool = false,
        )

        @assert (t_end - t_start >= Δt) "Simulation timespan cannot be shorter "* "
                                        than the discrete dynamics update period"

        params = (sys = sys, sys_init! = sys_init!, user_callback! = user_callback!, Δt = Δt,
                  args_ode = args_ode, args_step = args_step, args_disc = args_disc)

        cb_cont = DiscreteCallback((u, t, integrator)->true, f_cb_cont!)
        cb_step = DiscreteCallback((u, t, integrator)->true, f_cb_step!)
        cb_disc = PeriodicCallback(f_cb_disc!, Δt)
        cb_user = DiscreteCallback((u, t, integrator)->true, f_cb_user!)

        log = SavedValues(Float64, typeof(sys.y))
        saveat = (saveat isa Real ? (t_start:saveat:t_end) : saveat)
        cb_save = SavingCallback(f_cb_save, log; saveat, save_everystep, save_start, save_end)

        if save_on
            cb_set = CallbackSet(cb_cont, cb_step, cb_disc, cb_user, cb_save)
        else
            cb_set = CallbackSet(cb_cont, cb_step, cb_disc, cb_user)
        end

        #the current System's x value is used as initial condition. a copy is
        #needed, otherwise it's just a reference and will get overwritten during
        #simulation. if the System has no continuous state, provide a dummy one
        x0 = (has_x(sys) ? copy(sys.x) : [0.0])
        problem = ODEProblem{true}(f_ode_wrapper!, x0, (t_start, t_end), params)
        integrator = init_integrator(problem, algorithm; save_on = false,
                                     callback = cb_set, adaptive, dt)

        #save_on = false because we are not interested in logging the plain
        #System's state vector; everything we need to know about the System
        #should be made available in its (immutable) output y, logged via
        #SavingCallback

        #the GUI Renderer's refresh rate must be uncapped so that calls to
        #update!() return immediately without blocking, so that they do not
        #interfere with simulation scheduling

        info = Ref(SimInfo())

        f_draw = let sys = sys, info = info
            () -> GUI.draw!(sys, info[])
        end

        gui = Renderer(; label = "Simulation", sync = UInt8(0), f_draw)

        new{typeof(sys), typeof(integrator), typeof(log), typeof(sys.y)}(
            sys, integrator, log, gui, info, SimInput[], SimOutput[],
            Base.Event(), ReentrantLock(), ReentrantLock())
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

#function signature for a user callback function, called after every integration
#step. this function MUST NEVER modify the System's ẋ, x or t
user_callback!(::System) = nothing

############################ Stepping functions ################################

@inline has_x(sys::System) = !isnothing(sys.x)

#wrapper around the System's continuous dynamics function f_ode! it modifies the
#System's ẋ and y as a side effect
function f_ode_wrapper!(u̇, u, p, t)

    @unpack sys, args_ode = p

    #assign current integrator solution to System's continuous state
    has_x(sys) && (sys.x .= u)

    #same with time
    sys.t[] = t

    f_ode!(sys, args_ode...) #call continuous dynamics, updates sys.ẋ and sys.y

    has_x(sys) && (u̇ .= sys.ẋ) #update the integrator's derivative

end

#System update function, called on every integration step. brings the
#System's internal x and y up to date with the last integrator's solution step
function f_cb_cont!(integrator)

    @unpack t, u, p = integrator
    @unpack sys, args_ode = p

    #assign final integrator solution for this epoch to System's continuous state
    has_x(sys) && (sys.x .= u)

    #same with time
    sys.t[] = t

    #at this point sys.y and sys.ẋ hold the values from the last algorithm
    #evaluation of f_ode!, not the one corresponding to the updated x. with x
    #up to date, we can now compute the final sys.y and sys.ẋ for this epoch
    f_ode!(sys, args_ode...)

    u_modified!(integrator, false)

end

#DiscreteCallback function, called on every integration step. calls the System's
#own f_step!. modifications to x or s will not propagate to y until the next
#integration step
function f_cb_step!(integrator)

    @unpack u, p = integrator
    @unpack sys, args_step = p

    f_step!(sys, args_step...)

    #assign the (potentially) modified sys.x back to the integrator
    has_x(sys) && (u .= sys.x)

end

#PeriodicCallback function, calls the System's discrete dynamics update with the
#period Δt given to the Simulation constructor. modifications to x or s will not
#propagate to y until the next integration step
function f_cb_disc!(integrator)

    @unpack u, p = integrator
    @unpack sys, Δt, args_disc = p

    f_disc!(sys, Δt, args_disc...)

    #assign the (potentially) modified sys.x back to the integrator
    has_x(sys) && (u .= sys.x)

end

#DiscreteCallback function, calls the user-specified System I/O function after
#every integration step
function f_cb_user!(integrator)

    @unpack u, p = integrator
    @unpack sys, user_callback! = p

    user_callback!(sys)

    #assign the (potentially) modified sys.x back to the integrator
    has_x(sys) && (u .= sys.x)

end

#SavingCallback function, gets called at the end of each step after f_disc!
#and/or f_step!
f_cb_save(x, t, integrator) = deepcopy(integrator.p.sys.y)


####################### OrdinaryDiffEq extensions ##############################

OrdinaryDiffEq.step!(sim::Simulation, args...) = step!(sim.integrator, args...)

OrdinaryDiffEq.get_proposed_dt(sim::Simulation) = get_proposed_dt(sim.integrator)

OrdinaryDiffEq.add_tstop!(sim::Simulation, t) = add_tstop!(sim.integrator, t)

function OrdinaryDiffEq.reinit!(sim::Simulation, sys_init_args...; sys_init_kwargs...)

    @unpack integrator, log = sim
    @unpack p = integrator

    #initialize the System's x, u and s
    p.sys_init!(p.sys, sys_init_args...; sys_init_kwargs...)

    #initialize the ODEIntegrator with the System's initial x. ODEIntegrator's
    #reinit! calls f_ode_wrapper!, so the System's ẋ and y are updated in the
    #process, no need to do it explicitly
    if has_x(p.sys)
        OrdinaryDiffEq.reinit!(integrator, p.sys.x)
    else
        OrdinaryDiffEq.reinit!(integrator)
    end

    resize!(log.t, 1)
    resize!(log.saveval, 1)
    return nothing
end


################################# GUI ##########################################

function GUI.draw!(sys::System, info::SimInfo)
    GUI.draw(info)
    GUI.draw!(sys)
end

function GUI.draw(info::SimInfo)

    @unpack algorithm, t_start, t_end, dt, iter, t, τ = info

    CImGui.Begin("Simulation")
        CImGui.Text("Algorithm: " * algorithm)
        CImGui.Text("Step size: $dt")
        CImGui.Text("Iterations: $iter")
        CImGui.Text(@sprintf("Simulation framerate: %.3f ms/frame (%.1f FPS)",
                            1000 / unsafe_load(CImGui.GetIO().Framerate),
                            unsafe_load(CImGui.GetIO().Framerate)))
        CImGui.Text(@sprintf("Simulation time: %.3f s", t) * " [$t_start, $t_end]")
        CImGui.Text(@sprintf("Wall-clock time: %.3f s", τ))
        CImGui.Text(@sprintf("Simulation pace: x%.3f", (t - t_start) / τ))
    CImGui.End()

end


################################### I/O #######################################

#attaches an input device to the Simulation, enabling it to safely modify the
#System's input during paced or full speed execution
function attach!(sim::Simulation, device::InputDevice;
                    mapping::InputMapping = DefaultMapping())

    input = SimInput(device, sim.sys, mapping, sim.started, sim.running, sim.stepping)
    push!(sim.inputs, input)

end

#attaches an output device to the Simulation, linking it to a dedicated channel
#for simulation output data

# the channel must be buffered. an unbuffered channel returns isready only if
#there is another task waiting on it. however, to prevent the simulation loop
#from blocking, we need the channel to return isready when it is holding data,
#regardless of whether the output interface is currently waiting on it
function attach!(sim::Simulation, device::OutputDevice)
    channel = Channel{SimData{typeof(sim.sys.y)}}(1)
    output = SimOutput(device, channel, sim.started, sim.running)
    push!(sim.outputs, output)
end

function IODevices.put_no_wait!(sim::Simulation)
    data = SimData(sim.t[], sim.y)
    for output in sim.outputs
        put_no_wait!(output, data)
    end
end


################################ Execution #####################################

#compare the following:
# @time sleep(2)

# @time @async sleep(2)
# @time Threads.@spawn sleep(2)

# wait(@time @async sleep(2))
# wait(@time Threads.@spawn sleep(2))

# @time wait(@async sleep(2))
# @time wait(Threads.@spawn @sleep(2))

isdone(sim::Simulation) = isempty(sim.integrator.opts.tstops)

function isdone_err(sim::Simulation)
    isdone(sim) && @error("Simulation has hit its end time, call reinit! ",
                        "or add further tstops using add_tstop!")
end

#for non-paced simulations, we ignore input devices and only start output
#interfaces
function run!(sim::Simulation)
    @sync begin
        for output in sim.outputs
            Threads.@spawn IODevices.start!(output)
        end
        Threads.@spawn sim_loop!(sim)
    end
end

function sim_loop!(sim::Simulation)

    @unpack integrator, started, running, stepping = sim

    t_end = integrator.sol.prob.tspan[2]
    τ = 0

    try

        isdone_err(sim)
        @info("Simulation: Starting on thread $(Threads.threadid())...")

        lock(running)
        notify(started)

        τ = @elapsed begin
            while sim.t[] < t_end
                lock(stepping)
                    step!(sim)
                unlock(stepping)
                put_no_wait!(sim)
            end
        end

        @info("Simulation: Finished in $τ seconds")

    catch ex

        st = stacktrace(catch_backtrace())
        @error("Simulation: Terminated with $ex in $(st[1])")

    finally

        sim_cleanup!(sim)

    end

end

#apparently, if a task is launched from the main thread and it doesn't ever
#block, no other thread will get CPU time until it's done. threfore, we should
#never call _run_paced! directly from the main thread, because it will starve
#all IO threads. it must always be run from a Threads.@spawn'ed thread
function run_paced!(sim::Simulation; kwargs...)
    @sync begin
        for input in sim.inputs
            Threads.@spawn IODevices.start!(input)
        end
        for output in sim.outputs
            Threads.@spawn IODevices.start!(output)
        end
        Threads.@spawn sim_loop_paced!(sim; kwargs...)
    end
end


function sim_loop_paced!(sim::Simulation;
                    pace::Real = 1)

    @unpack sys, integrator, gui, info, started, running, stepping = sim

    t_start = sim.t[]
    t_end = integrator.sol.prob.tspan[2]
    algorithm = sim.integrator.alg |> typeof |> string

    try

        @assert gui.sync == 0
        GUI.init!(gui)

        τ = let wall_time_ref = time()
            ()-> time() - wall_time_ref
        end

        τ_last = τ()

        isdone_err(sim)
        @info("Simulation: Starting on thread $(Threads.threadid())...")

        lock(running)
        notify(started)

        while sim.t[] < t_end

            #simulation time t and wall-clock time τ are related by:
            #t_next = t_start + pace * τ_next

            t_next = sim.t[] + get_proposed_dt(sim)
            τ_next = (t_next - t_start) / pace
            while τ_next > τ() end #busy wait (ugly, should do better but it's not easy)

            lock(stepping)
                step!(sim)
                info[] = SimInfo(; algorithm, t_start, t_end, dt = integrator.dt,
                              iter = integrator.iter, t = sim.t[], τ = τ_last)
                GUI.update!(gui)
            unlock(stepping)

            put_no_wait!(sim)

            if GUI.should_close(gui)
                @info("Simulation: Aborted at t = $(sim.t[])")
                break
            end

        end

        @info("Simulation: Finished in $τ_last seconds")

    catch ex

        st = stacktrace(catch_backtrace())
        @error("Simulation: Terminated with $ex in $(st[1])")

    finally

        sim_cleanup!(sim)

    end

end

function sim_cleanup!(sim::Simulation)

    @unpack gui, started, running = sim

    #signal all interfaces to shut down
    islocked(running) && unlock(running)

    #unblock any interfaces waiting for simulation to start in case we
    #exited with error and didn't get to notify them
    notify(started)

    #reset the event for the next sim execution
    reset(started)

    #unblock any output interfaces waiting for a take!
    put_no_wait!(sim)

    gui._initialized && GUI.shutdown!(gui)
end



################################################################################
############################### TimeSeries ####################################

mutable struct TimeSeries{V, T <: AbstractVector{Float64}, D <: AbstractVector{V}}
    _t::T
    _data::D
    function TimeSeries(t::T, data::D) where {T <: AbstractVector{<:AbstractFloat}, D <: AbstractVector{V}} where {V}
        @assert length(t) == length(data)
        new{V, T, D}(t, data)
    end
end

TimeSeries(t::Real, data) = TimeSeries([Float64(t)], [data])

function TimeSeries(t::AbstractVector, M::Matrix)
    #each Matrix column interpreted as one Vector value
    TimeSeries(t, [M[:, i] for i in 1:size(M,2)])
end

Base.length(ts::TimeSeries) = length(ts._t)

function Base.getproperty(ts::TimeSeries, s::Symbol)
    t = getfield(ts, :_t)
    data = getfield(ts, :_data)
    if s === :_t
        return t
    elseif s === :_data
        return data
    else
        return TimeSeries(t, getproperty(StructArray(data), s))
    end
end

get_time(ts::TimeSeries) = getfield(ts, :_t)
get_data(ts::TimeSeries) = getfield(ts, :_data)

Base.getindex(ts::TimeSeries, i) = TimeSeries(ts._t[i], ts._data[i])

function Base.getindex(ts::TimeSeries{<:AbstractArray{T, N}}, time_ind, comp_ind::Vararg{Any, N}) where {T, N}

    t = ts._t[time_ind]
    data_voa = VectorOfArray(ts._data)[comp_ind..., time_ind]
    data = ndims(data_voa) == 1 ? data_voa : eachslice(data_voa; dims = ndims(data_voa))

    TimeSeries(t, data)

end

function Base.lastindex(ts::TimeSeries)
    lastindex(ts._t)
end
Base.view(ts::TimeSeries, i) = TimeSeries(view(ts._t, i), view(ts._data, i))

#for inspection
get_child_names(::T) where {T <: TimeSeries} = get_child_names(T)
get_child_names(::Type{<:TimeSeries{V}}) where {V} = fieldnames(V)

#could be rewritten as @generated to avoid allocation
function get_components(ts::TimeSeries{<:AbstractVector{T}}) where {T<:Real}
    (TimeSeries(ts._t, y) for y in ts._data |> StructArray |> StructArrays.components)
end

TimeSeries(sim::Simulation) = TimeSeries(sim.log.t, sim.log.saveval)


end #module