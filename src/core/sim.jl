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
export Simulation, enable_gui!, disable_gui!, attach_io!
export TimeSeries, get_time, get_data, get_components, get_child_names


################################################################################
################################# SimInfo ########################################

@kwdef struct SimInfo
    algorithm::String = "No info"
    t_start::Float64 = 0.0
    t_end::Float64 = 0.0
    dt::Float64 = 0.0 #last time step
    iter::Int64 = 0 #total iterations
    t::Float64 = 0.0 #simulation time
    τ::Float64 = 0.0 #wall-clock time
end


################################################################################
############################### Output #########################################

@kwdef struct Output{Y}
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
    channels::Vector{Channel{Output{Y}}}
    started::Base.Event #signals that execution has started
    stepping::ReentrantLock #must be acquired by IO interfaces to modify the system

    #set sync = 0 so that SwapBuffers returns immediately and doesn't slow down
    #the sim loop; for paced simulations, scheduling is taken care of
    function Simulation(
        sys::System;
        args_ode::Tuple = (), #externally supplied arguments to System's f_ode!
        args_step::Tuple = (), #externally supplied arguments to System's f_step!
        args_disc::Tuple = (), #externally supplied arguments to System's f_disc!
        sys_init!::Function = Systems.init!, #default System initialization function
        sys_io!::Function = no_sys_io!, #System I/O function
        gui::Renderer = Renderer(label = "Simulation", sync = 0),
        started::Base.Event = Base.Event(), #to be shared with input interfaces
        stepping::ReentrantLock = ReentrantLock(), #to be shared with input interfaces
        algorithm::OrdinaryDiffEqAlgorithm = RK4(),
        adaptive::Bool = false,
        dt::Real = 0.02, #continuous dynamics integration step
        Δt::Real = 1.0, #discrete dynamics execution period (do not set to Inf!)
        t_start::Real = 0.0,
        t_end::Real = 10.0,
        save_on::Bool = true,
        saveat::Union{Real, AbstractVector{<:Real}} = Float64[], #defers to save_everystep
        save_start::Bool = false, #initial System's outputs might not be up to date
        save_everystep::Bool = isempty(saveat),)

        @assert (t_end - t_start >= Δt) "Simulation timespan cannot be shorter "* "
                                        than the discrete dynamics update period"

        params = (sys = sys, sys_init! = sys_init!, sys_io! = sys_io!, Δt = Δt,
                  args_ode = args_ode, args_step = args_step, args_disc = args_disc)

        cb_cont = DiscreteCallback((u, t, integrator)->true, f_cb_cont!)
        cb_step = DiscreteCallback((u, t, integrator)->true, f_cb_step!)
        cb_disc = PeriodicCallback(f_cb_disc!, Δt)
        cb_io = DiscreteCallback((u, t, integrator)->true, f_cb_io!)

        log = SavedValues(Float64, typeof(sys.y))
        saveat_arr = (saveat isa Real ? (t_start:saveat:t_end) : saveat)
        cb_save = SavingCallback(f_cb_save, log; save_start = save_start,
                                    saveat = saveat_arr, save_everystep)

        if save_on
            cb_set = CallbackSet(cb_cont, cb_step, cb_disc, cb_io, cb_save)
        else
            cb_set = CallbackSet(cb_cont, cb_step, cb_disc, cb_io)
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

        channels = Channel{Output{typeof(sys.y)}}[]

        new{typeof(sys), typeof(integrator), typeof(log), typeof(sys.y)}(
            sys, integrator, log, gui, channels, started, stepping)
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

#get underlying component type
component_type(::Simulation{<:System{C}}) where {C} = C

#get ODE algorithm type
algorithm_type(sim::Simulation) = sim.integrator.alg |> typeof


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

    x_mod = f_step!(sys, args_step...)

    @assert x_mod isa Bool

    u_modified!(integrator, x_mod)

    #assign the modified sys.x back to the integrator
    has_x(sys) && x_mod && (u .= sys.x)

    return nothing

end

#PeriodicCallback function, calls the System's discrete dynamics update with the
#period Δt given to the Simulation constructor. modifications to x or s will not
#propagate to y until the next integration step
function f_cb_disc!(integrator)

    @unpack u, p = integrator
    @unpack sys, Δt, args_disc = p

    x_mod = f_disc!(sys, Δt, args_disc...)

    @assert x_mod isa Bool

    u_modified!(integrator, x_mod)

    #assign the modified sys.x back to the integrator
    has_x(sys) && x_mod && (u .= sys.x)

end

#DiscreteCallback function, calls the user-specified System I/O function after
#every integration step
function f_cb_io!(integrator)

    @unpack sys, sys_io! = integrator.p

    sys_io!(sys)

    #a System I/O function should never modify the System's continuous state, so
    #we may as well tell the integrator to avoid any performance hit
    u_modified!(integrator, false)

end

#SavingCallback function, gets called at the end of each step after f_disc!
#and/or f_step!
f_cb_save(x, t, integrator) = deepcopy(integrator.p.sys.y)

#function signature for a System I/O function. note: an I/O function MUST NOT
#modify a System's ẋ, x or t
no_sys_io!(::System) = nothing


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
    OrdinaryDiffEq.reinit!(integrator, p.sys.x)

    resize!(log.t, 1)
    resize!(log.saveval, 1)
    return nothing
end


################################# GUI ##########################################

enable_gui!(sim::Simulation) = GUI.enable!(sim.gui)
disable_gui!(sim::Simulation) = GUI.disable!(sim.gui)

function update!(sys::System, info::SimInfo, gui::Renderer)
    try
        GUI.render(gui, GUI.draw!, sys, info)
    catch e
        @error "Error while updating window" exception=e
        Base.show_backtrace(stderr, catch_backtrace())
        shutdown!(renderer)
    end
end

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
#System's input during paced or full speed execution, and returns its interface
function attach_io!(sim::Simulation, device::IODevice;
                    mapping::InputMapping = DefaultMapping(),
                    ext_shutdown::Bool = true)

    channel = add_output_channel!(sim)
    IODevices.Interface(device, sim.sys, mapping, channel,
                        sim.started, sim.stepping, ext_shutdown)
end

add_output_channel!(sim::Simulation) = push!(sim.channels, eltype(sim.channels)(1))[end]

pop_output_channel!(sim::Simulation) = pop!(sim.channels)

#for buffered channels, isready returns 1 if there is at least one value in the
#channel. if there is one value and we try to put! another, the simulation loop
#will block. to avoid this, we only put! if !isready. we wait for no one!
@inline function put_no_block!(channel::Channel{T}, data::T) where {T}
    (isopen(channel) && !isready(channel)) && put!(channel, data)
end

@inline function put_no_block!(channels::Vector{Channel{T}}, data::T) where {T}
    for channel in channels
        put_no_block!(channel, data)
    end
end


################################ Execution #####################################

isdone(sim::Simulation) = isempty(sim.integrator.opts.tstops)

function isdone_err(sim::Simulation)
    sim_done = isdone(sim)
    sim_done && @error("Simulation has hit its end time, call reinit! ",
                        "or add further tstops using add_tstop!")
    return sim_done
end

run_thr!(sim::Simulation; kwargs...) = wait(Threads.@spawn(run!(sim; kwargs...)))

function run!(sim::Simulation)

    @unpack integrator, channels, started, stepping = sim
    t_end = integrator.sol.prob.tspan[2]

    isdone_err(sim) && return

    notify(started)

    τ = @elapsed begin
        while sim.t[] < t_end
            lock(stepping)
                step!(sim)
            unlock(stepping)
            output = Output(sim.t[], sim.y)
            put_no_block!(channels, output)
        end
    end

    @info("Simulation: Finished in $τ seconds")
end

#apparently, if a task is launched from the main thread and it doesn't ever
#block, no other thread will get CPU time until it's done. threfore, we should
#never call _run_paced! directly from the main thread, because it will starve
#all IO threads. it must always be run from a Threads.@spawn'ed thread
run_paced_thr!(sim::Simulation; kwargs...) = wait(Threads.@spawn(run_paced!(sim; kwargs...)))

function run_paced!(sim::Simulation;
                    pace::Real = 1)

    isdone_err(sim) && return
    @info("Simulation: Starting on thread $(Threads.threadid())...")

    @unpack sys, integrator, gui, channels, started, stepping = sim

    t_start, t_end = integrator.sol.prob.tspan
    algorithm = algorithm_type(sim) |> string

    gui.sync == 0 || @warn("GUI Renderer sync ",
        "should be set to 0 to avoid interfering with simulation scheduling")

    GUI.init!(gui)

    τ = let wall_time_ref = time()
            ()-> time() - wall_time_ref
        end

    τ_last = τ()

    try

        notify(started)

        while sim.t[] < t_end

            #simulation time t and wall-clock time τ are related by:
            #t_next = t_start + pace * τ_next

            t_next = sim.t[] + get_proposed_dt(sim)
            τ_next = (t_next - t_start) / pace
            while τ_next > τ() end #busy wait (disgusting, needs to do better)

            lock(stepping)
                step!(sim)
                τ_last = τ()
                info = SimInfo(; algorithm, t_start, t_end, dt = integrator.dt,
                              iter = integrator.iter, t = sim.t[], τ = τ_last)
                update!(sys, info, gui)
            unlock(stepping)

            output = Output(sim.t[], sim.y)
            put_no_block!(sim.channels, output)

            if GUI.should_close(gui)
                @info("Simulation: Aborted at t = $(sim.t[])")
                break
            end

        end

    catch ex

        @error("Simulation: Error during execution: $ex")
        Base.show_backtrace(stderr, catch_backtrace())

    finally

        close.(channels)
        GUI.shutdown!(gui)

    end

    @info("Simulation: Finished in $τ_last seconds")

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