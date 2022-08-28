module Sim

using UnPack
using StructArrays
using SciMLBase: ODEProblem, u_modified!, init as init_integrator
using OrdinaryDiffEq: OrdinaryDiffEqAlgorithm, ODEIntegrator, RK4
using DiffEqCallbacks: SavingCallback, DiscreteCallback, PeriodicCallback, CallbackSet, SavedValues

using CImGui
using Printf

using ..Systems
using ..IODevices
using ..GUI

import SciMLBase: step!, reinit!, solve!, get_proposed_dt
import OrdinaryDiffEq: add_tstop!

export Simulation, step!, reinit!, solve!, get_proposed_dt, add_tstop!
export TimeHistory, get_components, get_child_names


################################################################################
################################# Info ########################################

Base.@kwdef struct Info
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

Base.@kwdef struct Output{Y}
    t::Float64
    y::Y #System's y
end

################################################################################
############################# Simulation #######################################

struct Simulation{S <: System, I <: ODEIntegrator, L <: SavedValues, Y}
    sys::S
    integrator::I
    log::L
    output_channels::Vector{Channel{Output{Y}}}
    started::Base.Event #signals that execution has started
    stepping::ReentrantLock #must be acquired by IO interfaces to modify the system

    function Simulation(
        sys::System;
        args_ode::Tuple = (), #externally supplied arguments to System's f_ode!
        args_step::Tuple = (), #externally supplied arguments to System's f_step!
        args_disc::Tuple = (), #externally supplied arguments to System's f_disc!
        sys_init!::Function = no_sys_init!, #System initialization function
        sys_io!::Function = no_sys_io!, #System I/O function
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
        save_everystep::Bool = isempty(saveat),
        sys_init_kwargs = NamedTuple())

        @assert !(sys.x === nothing) "Purely discrete Systems aren't currently supported"

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
        #available in the System's y, saved by the SavingCallback

        output_channels = Channel{Output{typeof(sys.y)}}[]

        new{typeof(sys), typeof(integrator), typeof(log), typeof(sys.y)}(
            sys, integrator, log, output_channels, started, stepping)
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


################# SciMLBase & OrdinaryDiffEq extensions ########################

step!(sim::Simulation, args...) = step!(sim.integrator, args...)

get_proposed_dt(sim::Simulation) = get_proposed_dt(sim.integrator)

add_tstop!(sim::Simulation, t) = add_tstop!(sim.integrator, t)

function reinit!(sim::Simulation; sys_init_kwargs...)

    @unpack integrator, log = sim
    @unpack p = integrator

    if p.sys_init! === no_sys_init!
        println("Warning: Simulation has no sys_init! method; the System's ",
        "states and inputs will remain at their latest values")
    end

    #initialize the System's x, u and s
    p.sys_init!(p.sys; sys_init_kwargs...)

    #initialize the ODEIntegrator with the System's initial x. ODEIntegrator's
    #reinit! calls f_ode_wrapper!, so the System's ẋ and y are updated in the
    #process, no need to do it explicitly
    reinit!(integrator, p.sys.x)

    resize!(log.t, 1)
    resize!(log.saveval, 1)
    return nothing
end


################################### I/O #######################################

#attaches an input device to the Simulation, enabling it to safely modify the
#System's input during paced or full speed execution, and returns its interface
function attach_io!(sim::Simulation, device::IODevice;
                    mapping::InputMapping = DefaultMapping(),
                    ext_shutdown::Bool = true)

    channel = add_output_channel!(sim)
    IODevices.Interface(device, sim.u, mapping, channel,
                        sim.started, sim.stepping, ext_shutdown)
end

add_output_channel!(sim::Simulation) = push!(sim.output_channels,
                                            eltype(sim.output_channels)(1))[end]

pop_output_channel!(sim::Simulation) = pop!(sim.output_channels)

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

function run!(sim::Simulation; rate::Real = Inf, verbose::Bool = false)

    #rate: rate of execution relative to wall time (Inf = unrestricted, 1 ≈ real time )

    if isempty(sim.integrator.opts.tstops)
        println("Simulation has hit t_end, reset it using reinit! ",
                "or add further tstops using add_tstop!")
        return
    end

    rate === Inf ? run_fullspeed!(sim; verbose) : run_paced!(sim; rate, verbose)

end

function run_fullspeed!(sim::Simulation; verbose::Bool)
    τ = @elapsed begin
        for _ in sim.integrator #integrator steps automatically at the beginning of each iteration
            output = Output(sim.t[], sim.y)
            put_no_block!(sim.output_channels, output)
        end
    end
    verbose && println("Simulation: Finished in $τ seconds")
end

#apparently, if a task is launched from the main thread and it doesn't ever
#block, no other thread will get CPU time until it's done. threfore, we should
#never call _run_paced! directly from the main thread, because it will starve
#all IO threads. it must always be run from a Threads.@spawn'ed thread

function run_paced!(sim::Simulation; rate::Real, verbose::Bool)
    wait(Threads.@spawn(_run_paced!(sim; rate, verbose)))
end

function _run_paced!(sim::Simulation; rate::Real, verbose::Bool)

    @unpack sys, integrator, output_channels, started, stepping = sim

    verbose && println("Simulation: Starting on thread $(Threads.threadid())...")

    t_start, t_end = integrator.sol.prob.tspan
    algorithm = algorithm_type(sim) |> string

    #set refresh = 0 so that SwapBuffers returns immediately and doesn't slow
    #down the sim loop
    dashboard = Renderer(refresh = 0, wsize = (320, 240), label = "Simulation")

    GUI.init!(dashboard)

    # GLFW.SetWindowPos(dashboard.renderer._window, 100, 100)

    τ = let wall_time_ref = time()
            ()-> time() - wall_time_ref
        end

    τ_last = τ()

    try

        notify(started)

        while sim.t[] < t_end

            #simulation time t and wall-clock time τ should satisfy:
            #t_next = t_start + rate * τ_next

            t_next = sim.t[] + get_proposed_dt(sim)
            τ_next = (t_next - t_start) / rate
            while τ_next > τ() end #busy wait

            lock(stepping)
                step!(sim)
                τ_last = τ()
                info = Info(; algorithm, t_start, t_end, dt = integrator.dt,
                              iter = integrator.iter, t = sim.t[], τ = τ_last)
                update!(sys, info, dashboard)
            unlock(stepping)

            output = Output(sim.t[], sim.y)
            put_no_block!(sim.output_channels, output)

            if GUI.should_close(dashboard)
                println("Simulation: Aborted at t = $(sim.t[])")
                break
            end

        end

    catch ex

        println("Simulation: Error during execution: $ex")
        Base.show_backtrace(stderr, catch_backtrace())

    finally

        close.(output_channels)
        GUI.shutdown!(dashboard)

    end

    verbose && println("Simulation: Finished in $τ_last seconds")

end

function update!(sys::System, info::Info, dashboard::Renderer)
    GUI.render(dashboard, GUI.draw!, sys.u, sys.y, info)
end

function GUI.draw!(u, y, info::Info)
    GUI.draw!(info)
    GUI.draw!(u, y)
end

function GUI.draw!(info::Info)

    @unpack algorithm, t_start, t_end, dt, iter, t, τ = info

    begin
        CImGui.Begin("Info")
            CImGui.Text("Algorithm: " * algorithm)
            CImGui.Text("Step size: $dt")
            CImGui.Text("Iterations: $iter")
            CImGui.Text(@sprintf("Simulation time: %.3f s", t) * " [$t_start, $t_end]")
            CImGui.Text(@sprintf("Wall-clock time: %.3f s", τ))
            #replace with simulation rate: (t - t_start) / τ
            # CImGui.Text(@sprintf("Simulation: %.3f ms/frame (%.1f FPS)",
            #                     1000 / CImGui.GetIO().Framerate,
            #                     CImGui.GetIO().Framerate))
        CImGui.End()
    end

end


################################################################################
############################### TimeHistory ####################################

mutable struct TimeHistory{V, T <: AbstractVector{Float64}, D <: AbstractVector{V}}
    _t::T
    _data::D
    function TimeHistory(t::T, data::D) where {T, D <: AbstractVector{V}} where {V}
        @assert length(t) == length(data)
        new{V, T, D}(t, data)
    end
end

TimeHistory(t::Real, data) = TimeHistory([Float64(t)], [data])

function TimeHistory(t::AbstractVector, M::Matrix)
    #each Matrix column interpreted as one Vector value
    TimeHistory(t, [M[:, i] for i in 1:size(M,2)])
end

Base.length(th::TimeHistory) = length(th._t)

function Base.getproperty(th::TimeHistory, s::Symbol)
    t = getfield(th, :_t)
    y = getfield(th, :_data)
    if s === :_t
        return t
    elseif s === :_data
        return y
    else
        return TimeHistory(t, getproperty(StructArray(y), s))
    end
end

timestamps(th::TimeHistory) = getfield(th, :_t)

Base.getindex(th::TimeHistory, i) = TimeHistory(th._t[i], th._data[i])
Base.view(th::TimeHistory, i) = TimeHistory(view(th._t, i), view(th._data, i))

#for inspection
get_child_names(::T) where {T <: TimeHistory} = get_child_names(T)
get_child_names(::Type{<:TimeHistory{V}}) where {V} = fieldnames(V)

#could be rewritten as @generated to avoid allocation
function get_components(th::TimeHistory{<:AbstractVector{T}}) where {T<:Real}
    [TimeHistory(th._t, y) for y in th._data |> StructArray |> StructArrays.components]
end

TimeHistory(sim::Simulation) = TimeHistory(sim.log.t, sim.log.saveval)


end #module