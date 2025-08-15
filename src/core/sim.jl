module Sim

using UnPack, Reexport, StructArrays
using OrdinaryDiffEq: OrdinaryDiffEq, OrdinaryDiffEqAlgorithm, ODEProblem,
                      ODEIntegrator, Heun, RK4, u_modified!, init as init_integrator
using DiffEqCallbacks: SavingCallback, DiscreteCallback, PeriodicCallback,
                       CallbackSet, SavedValues
using RecursiveArrayTools
using Logging
using CImGui.lib: ImGuiSliderFlags_Logarithmic

using ..Modeling
using ..IODevices
using ..GUI

@reexport using OrdinaryDiffEq: step!, reinit!, add_tstop!, get_proposed_dt
export Simulation, enable_gui!, disable_gui!, attach!
export TimeSeries, get_time, get_data, get_components
export take_nonblocking!, put_nonblocking!


################################################################################
############################ SimControl ########################################

#this struct is accessed concurrently by the simulation (within the main loop)
#and GUI (within GUI.update!), so it must be guarded by io_lock
@kwdef mutable struct SimControl
    running::Bool = false #to be checked on each loop iteration for termination
    paused::Bool = false #to pause or unpause the simulation
    pace::Float64 = 1.0 #to be set by the SimControl GUI
    algorithm::String = ""
    t_start::Float64 = 0.0
    t_end::Float64 = 0.0
    Δt::Float64 = 0.0 #root system's periodic update interval
    dt::Float64 = 0.0 #last continuous integration step size
    iter::Int64 = 0 #total iterations
    t::Float64 = 0.0 #simulation time
    τ::Float64 = 0.0 #wall-clock time
end

function GUI.draw!(control::SimControl)

    CImGui.Begin("Simulation Control")

        mode_button("Pause", true, false, control.paused; HSV_active = HSV_amber)
        CImGui.IsItemClicked() && (control.paused = !control.paused); CImGui.SameLine()

        control.pace = safe_slider("Pace", control.pace, 0.1, 20.0, "%.3f",
        ImGuiSliderFlags_Logarithmic; show_label = true)

        @unpack algorithm, t_start, t_end, Δt, dt, iter, t, τ = control

        CImGui.Text("Algorithm: " * algorithm)
        CImGui.Text("Continuous step size: $dt")
        CImGui.Text("Periodic update interval: $Δt")
        CImGui.Text("Iterations: $iter")
        CImGui.Text(@sprintf("Simulation time: %.3f s", t) * " [$t_start, $t_end]")
        CImGui.Text(@sprintf("Wall-clock time: %.3f s", τ))
        CImGui.Text(@sprintf("GUI framerate: %.3f ms/frame (%.1f FPS)",
                            1000 / unsafe_load(CImGui.GetIO().Framerate),
                            unsafe_load(CImGui.GetIO().Framerate)))


    CImGui.End()

end

#to control the Simulation, an InputDevice that should define its own IOMapping
#subtype and use it to extend this function
IODevices.assign_input!(::SimControl, ::IOMapping, ::Any) = nothing

################################################################################
############################### SimInterface ###################################

struct SimInterface{D <: IODevice, S <: Model, M <: IOMapping}
    device::D
    mdl::S
    mapping::M
    control::SimControl
    io_start::Base.Event
    io_lock::ReentrantLock
    should_abort::Bool #Simulation should abort when this interface shuts down
end

const SimInput{D} = SimInterface{D} where {D <: InputDevice}
const SimOutput{D} = SimInterface{D} where {D <: OutputDevice}
const SimGUI{D} = SimInterface{D} where {D <: Renderer}

#called from the SimInterface loop thread. this may block, either on get_data
#because the InputDevice is itself blocked, on or lock(io_lock) because the
#Simulation loop is currently stepping
function update!(interface::SimInput)

    @unpack device, mdl, mapping, control, io_lock = interface
    data = IODevices.get_data!(device)

    #the data we got from the InputDevice might be something we're not able to
    #map to the Model (for example, an EOT character)
    try
        lock(io_lock)
        IODevices.assign_input!(mdl, mapping, data)
        IODevices.assign_input!(control, mapping, data)
    catch ex
        @warn("Failed to assign input data $data to Model")
    finally
        unlock(io_lock)
    end

end

#called from the SimInterface loop thread. this may block, either on handle_data
#because the OutputDevice is itself blocked, or on lock(io_lock) because the
#Simulation loop is currently stepping
function update!(interface::SimOutput)

    @unpack device, mdl, mapping, io_lock = interface

    lock(io_lock)
        #this call should never block and always return some usable output
        data = IODevices.extract_output(mdl, mapping)
    unlock(io_lock)

    IODevices.handle_data!(device, data)
end


function update!(gui::SimGUI)

    @unpack device, io_lock = gui

    #the GUI may modify the Model or SimControl, so we need to grab io_lock
    @lock io_lock GUI.render!(device)

    #once io_lock has been released and the simulation thread can proceed, put
    #the GUI thread to sleep for a while to limit framerate and save CPU time
    sleep(0.016)
end


########################## Currently unused ####################################

#unlike lock, trylock does not block if the Channel is already locked, while
#also ensuring it is not modified while we're checking its state

function take_nonblocking!(channel::Channel)
    if trylock(channel)
        data = (isready(channel) ? take!(channel) : nothing)
        unlock(channel)
        return data
    else
        return nothing
    end
end

function put_nonblocking!(channel::Channel{T}, data::T) where {T}
    if trylock(channel)
        (isopen(channel) && !isready(channel)) && put!(channel, data)
        unlock(channel)
    end
end


################################################################################
############################# Simulation #######################################

struct Simulation{D <: ModelDefinition, Y, I <: ODEIntegrator, G <: SimGUI}
    mdl::Model{D, Y}
    integrator::I
    log::SavedValues{Float64, Y}
    gui::G
    control::SimControl
    io_start::Base.Event #to be waited on by IO devices and GUI before entering their update loops
    io_lock::ReentrantLock #to be acquired before modifying the simulated Model or SimControl
    interfaces::Vector{SimInterface}

    function Simulation(
        mdl::Model{D, Y};
        user_callback!::Function = user_callback!, #user-specified callback
        algorithm::OrdinaryDiffEqAlgorithm = RK4(),
        adaptive::Bool = false,
        dt::Real = 0.02, #continuous dynamics integration step size
        Δt::Real = dt, #periodic update interval (do not set to Inf!)
        t_start::Real = 0.0,
        t_end::Real = 10000.0,
        save_on::Bool = true,
        saveat::Union{Real, AbstractVector{<:Real}} = Float64[], #defers to save_everystep
        save_everystep::Bool = isempty(saveat),
        ) where {D, Y}

        t_end - t_start < Δt && @warn(
        "Simulation timespan is shorter than the periodic update interval")

        #need to store a reference to mdl and user_callback! in the integrator
        #params so they can be used within the integrator callbacks
        params = (mdl = mdl, user_callback! = user_callback!)

        cb_step = DiscreteCallback((u, t, integrator)->true, f_cb_step!)
        cb_user = DiscreteCallback((u, t, integrator)->true, f_cb_user!)

        #we don't want to force a call to f_periodic! at t=0, so we set
        #initial_affect=false. if needed, such call can be included in the
        #user-defined Modeling.init!
        cb_periodic = PeriodicCallback(f_cb_periodic!, Δt; initial_affect = false)
        mdl._Δt_root[] = Δt

        #the optional call to f_periodic! at t=0 would correspond to mdl._n[] = 0.
        #therefore, in preparation for the first scheduled periodic update, we
        #should set mdl._n[] = 1
        mdl._n[] = 1

        log = SavedValues(Float64, Y)
        saveat = (saveat isa Real ? (t_start:saveat:t_end) : saveat)
        cb_save = SavingCallback(f_cb_save, log; saveat, save_everystep)

        if save_on
            cb_set = CallbackSet(cb_step, cb_periodic, cb_user, cb_save)
        else
            cb_set = CallbackSet(cb_step, cb_periodic, cb_user)
        end

        #the current Model's x is used as an initial condition. a copy is
        #needed, otherwise x0 will get overwritten during simulation. if the
        #Model has no continuous state, provide a dummy one
        x0 = (has_x(mdl) ? copy(mdl.x) : [0.0])

        #save_on = false because we are not interested in logging the plain
        #Model's state vector; everything we need from the Model should be
        #made available at its output, logged via SavingCallback
        problem = ODEProblem{true}(f_ode_wrapper!, x0, (t_start, t_end), params)
        integrator = init_integrator(problem, algorithm; save_on = false,
                                     callback = cb_set, adaptive, dt)

        control = SimControl()
        io_start = Base.Event()
        io_lock = ReentrantLock()

        f_draw = let control = control, mdl = mdl
            () -> GUI.draw!(control, mdl)
        end

        #sync should be set to 0 to disable VSync. otherwise, the call to
        #GUI.render! will block for a whole display refresh interval while
        #holding the io_lock, leaving the sim loop unable to step in the
        #meantime. the drawback of disabling VSync is that the GUI framerate is
        #uncapped, which unnecessarily taxes the CPU core running the main
        #thread. to prevent this, we can send the GUI to sleep for a while
        #within GUI.update! AFTER io_lock has been released
        renderer = Renderer(; label = "Simulation", sync = UInt8(0), f_draw)
        gui = SimInterface(renderer, mdl, GenericInputMapping(),
                           control, io_start, io_lock, true)

        new{D, Y, typeof(integrator), typeof(gui)}(
            mdl, integrator, log, gui, control, io_start, io_lock, SimInterface[])
    end

end

function Simulation(md::ModelDefinition, args...; kwargs...)
    Simulation(Model(md), args...; kwargs...)
end

function Base.propertynames(::Simulation)
    (fieldnames(Simulation)..., :x, :ẋ, :s, :u, :y, :t, :Δt)
end

Base.getproperty(sim::Simulation, s::Symbol) = getproperty(sim, Val(s))

@generated function Base.getproperty(sim::Simulation, ::Val{S}) where {S}
    if S === :t
        return :(getproperty(getfield(sim, :mdl), $(QuoteNode(S)))[])
    elseif S ∈ (:x, :ẋ, :s, :u, :y, :Δt)
        return :(getproperty(getfield(sim, :mdl), $(QuoteNode(S))))
    else
        return :(getfield(sim, $(QuoteNode(S))))
    end
end

#function signature for a user callback function, called after every integration
#step.
user_callback!(::Model) = nothing


################################ I/O ###########################################

function attach!(sim::Simulation, device::IODevice,
                mapping::IOMapping = get_default_mapping(device);
                should_abort = false)

    interface = SimInterface(device, sim.mdl, mapping,
                            sim.control, sim.io_start, sim.io_lock, should_abort)

    push!(sim.interfaces, interface)

end


############################ Stepping functions ################################

@inline has_x(mdl::Model) = !isnothing(mdl.x)

#wrapper around the root Model's continuous dynamics function
function f_ode_wrapper!(u̇, u, p, t)

    @unpack mdl = p

    #assign current integrator solution to Model's continuous state
    has_x(mdl) && (mdl.x .= u)

    #idem for time
    mdl.t[] = t

    f_ode!(mdl) #call continuous dynamics (updates mdl.ẋ and mdl.y)

    has_x(mdl) && (u̇ .= mdl.ẋ) #update the integrator's derivative

end

#DiscreteCallback wrapper around the root Model's post-integration step function
function f_cb_step!(integrator)

    @unpack u, p = integrator
    @unpack mdl = p

    f_step!(mdl) #potentially modifies x, u, s or y

    #assign the (potentially) modified mdl.x back to the integrator
    has_x(mdl) && (u .= mdl.x)

end

#DiscreteCallback wrapper around the user-defined post-integration step callback
function f_cb_user!(integrator)

    @unpack u, p = integrator
    @unpack mdl, user_callback! = p

    user_callback!(mdl) #potentially modifies x, u, s or y

    #assign the (potentially) modified mdl.x back to the integrator
    has_x(mdl) && (u .= mdl.x)

end

#PeriodicCallback wrapper around the root Model's periodic update function
function f_cb_periodic!(integrator)

    @unpack u, p = integrator
    @unpack mdl = p

    f_periodic!(mdl) #call scheduled periodic update, potentially updates mdl.s and mdl.y

    #increment the periodic update counter
    mdl._n[] += 1

end

#SavingCallback function, gets called at the end of each step after f_periodic!
#and/or f_step!
function f_cb_save(x, t, integrator)
    deepcopy(integrator.p.mdl.y)
end


####################### OrdinaryDiffEq extensions ##############################

OrdinaryDiffEq.step!(sim::Simulation, args...) = step!(sim.integrator, args...)

OrdinaryDiffEq.get_proposed_dt(sim::Simulation) = get_proposed_dt(sim.integrator)

function OrdinaryDiffEq.reinit!(sim::Simulation, init_args...; init_kwargs...)

    @unpack mdl, integrator, log = sim

    #drop the log entries from the last run
    resize!(log.t, 0)
    resize!(log.saveval, 0)

    #reset scheduling counter, so f_periodic! is guaranteed to execute if called by
    #Modeling.init!
    mdl._n[] = 0

    #initialize the Model's x, u and s
    Modeling.init!(mdl, init_args...; init_kwargs...)

    #initialize the integrator with the Model's initial x. within the
    #integrator's reinit! f_cb_save and f_ode_wrapper! are called, in that order
    if has_x(mdl)
        OrdinaryDiffEq.reinit!(integrator, mdl.x)
    else
        OrdinaryDiffEq.reinit!(integrator)
    end

    #prepare scheduling counter for the first integration step
    mdl._n[] = 1

    return nothing
end

function init!(sim::Simulation, args...; kwargs...)
    OrdinaryDiffEq.reinit!(sim, args...; kwargs...)
end


################################# GUI ##########################################

function GUI.draw!(control::SimControl, mdl::Model)
    GUI.draw!(control)
    GUI.draw!(mdl)
end


################################################################################
############################# Threaded loops ###################################

################################ SimInterface ##################################

function start!(interface::SimInterface{D}) where {D <: IODevice}

    @unpack device, control, io_start, io_lock, should_abort = interface

    @info("$D: Starting on thread $(Threads.threadid())...")

    try

        IODevices.init!(device)

        @info("$D: Waiting for Simulation...")
        wait(io_start)
        @info("$D: Running...")

        while true

            !(@lock io_lock control.running) && break

            if IODevices.should_close(device)
                @info("$D: Shutdown request received...")
                break
            end

            #we need to yield somewhere because this loop is not guaranteed to
            #block at any point, so it could slow down the simulation thread
            #significantly. in the case of SimGUI, this is by design: the
            #framerate is uncapped precisely to prevent the renderer from
            #blocking for VSync while update! is holding the io_lock
            yield()
            update!(interface)

        end

    catch ex

        @error("$D: Error during execution: $ex")

    finally
        if should_abort
            @lock io_lock begin
                control.paused = false
                control.running = false
            end
        end
        @info("$D: Shutting down...")
        IODevices.shutdown!(device)
        @info("$D: Closed")
    end

end


################################ Simulation ####################################

function start!(sim::Simulation)

    @unpack mdl, integrator, gui, control, io_start, io_lock, interfaces = sim

    try

        if isempty(sim.integrator.opts.tstops)
            @error("Simulation has hit its end time, call reinit! to reset it")
        end

        τ = let wall_time_ref = time()
            ()-> time() - wall_time_ref
        end

        t_start = integrator.sol.prob.tspan[1]
        t_end = integrator.sol.prob.tspan[2]

        lock(io_lock)
            control.running = true
            control.paused = false
            control.algorithm = sim.integrator.alg |> typeof |> string
            control.t_start = t_start
            control.t_end = t_end
            control.Δt = mdl._Δt_root[]
        unlock(io_lock)

        notify(io_start)

        τ_last = τ()

        @info("Simulation: Starting on thread $(Threads.threadid())...")

        Δτ = @elapsed begin

            while sim.t[] < t_end

                lock(io_lock)
                    @unpack running, paused, pace = control
                    control.dt = integrator.dt
                    control.iter = integrator.iter
                    control.t = sim.t[]
                    control.τ = τ()
                unlock(io_lock)

                if !running
                    @info("Simulation: Aborted at t = $(sim.t[])")
                    break
                end

                if control.paused
                    τ_last = τ()
                    continue #skips next Simulation step
                end

                #compute wall-clock time at the end of next simulation step
                τ_next = τ_last + get_proposed_dt(sim) / pace

                #update simulation
                @lock io_lock step!(sim)

                #wait for wall-clock time to catch up. busy wait is not ideal,
                #but it is more responsive than sleep(τ_next - τ()). sleep()
                #only guarantees a minimum sleep duration, and its granularity
                #is not great
                while τ_next > τ() end

                τ_last = τ_next

            end

        end

        @info("Simulation: Finished in $Δτ seconds")

    catch ex

        st = stacktrace(catch_backtrace())
        @error("Simulation: Terminated with $ex in $(st[1])")

    finally

        sim_cleanup!(sim)

    end

end

function sim_cleanup!(sim::Simulation)

    @unpack control, io_start, io_lock = sim

    #if the simulation ran to conclusion, signal IO threads to shut down
    @lock io_lock begin
        control.paused = false
        control.running = false
    end

    #unblock any IO threads still waiting for Simulation start
    notify(io_start)

    #prepare for another run
    reset(io_start)

end


################################################################################
############################### Execution ######################################

#* caution: if a task is launched from the main thread and it doesn't
#* ever block or yield, no other thread will get CPU time until it's done.

#* caution: on non-Windows platforms, CImGui needs to run on the main thread!

function run!(sim::Simulation; pace = Inf)

    req_threads = 1 + length(sim.interfaces) #sim loop + interfaces
    threads = Threads.nthreads()
    req_threads <= threads || error(
        "Running this simulation requires at least $req_threads available
        threads, the current Julia session only has $threads")

    sim.control.pace = pace

    @sync begin
        for interface in sim.interfaces
            Threads.@spawn start!(interface)
        end
        Threads.@spawn start!(sim)
    end

end

function run_interactive!(sim::Simulation; pace = 1.0)

    req_threads = 2 + length(sim.interfaces) #GUI + sim loop + interfaces
    threads = Threads.nthreads()
    req_threads <= threads || error(
        "Running this simulation requires at least $req_threads available
        threads, the current Julia session only has $threads")

    sim.control.pace = pace

    @sync begin
        for interface in sim.interfaces
            Threads.@spawn start!(interface)
        end
        Threads.@spawn start!(sim)
        start!(sim.gui) #CImGui must run on the main thread
    end

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

function get_components(ts::TimeSeries{<:AbstractArray{T}}) where {T<:Real}
    (TimeSeries(ts._t, y) for y in ts._data |> StructArray |> StructArrays.components)
end

Base.propertynames(ts::TimeSeries) = propertynames(getfield(ts, :_data)[1])

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

Base.eltype(::TimeSeries{V}) where {V} = V

Base.view(ts::TimeSeries, i) = TimeSeries(view(ts._t, i), view(ts._data, i))

TimeSeries(sim::Simulation) = TimeSeries(sim.log.t, sim.log.saveval)


################################################################################
############################### Inspection #####################################

# #prevent monstrous type signatures from flooding the REPL

function Base.show(io::IO, ::Simulation{D}) where {D <: ModelDefinition}
    maxlen = 100 #maximum length for ModelDefinition type parameter
    md_str = sprint(show, D)
    md_str = (length(md_str) < maxlen ? md_str : first(md_str, maxlen) * "...")
    write(io, "Simulation{" * md_str * "}(...)")
end

function Base.show(io::IO, ::TimeSeries{D}) where {D}
    maxlen = 200 #maximum length for ModelDefinition type parameter
    dtype_str = sprint(show, D)
    dtype_str = (length(dtype_str) < maxlen ? dtype_str : first(dtype_str, maxlen) * "...")
    write(io, "TimeSeries{" * dtype_str * "}(...)")
end

end #module
