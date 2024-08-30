module Sim

using UnPack, Reexport, StructArrays
using OrdinaryDiffEq: OrdinaryDiffEq, OrdinaryDiffEqAlgorithm, ODEProblem,
                      ODEIntegrator, Heun, RK4, u_modified!, init as init_integrator
using DiffEqCallbacks: SavingCallback, DiscreteCallback, PeriodicCallback,
                       CallbackSet, SavedValues
using RecursiveArrayTools
using Logging
using CImGui.LibCImGui

using ..Systems
using ..IODevices
using ..GUI

@reexport using OrdinaryDiffEq: step!, reinit!, add_tstop!, get_proposed_dt
export Simulation, enable_gui!, disable_gui!, attach!
export TimeSeries, get_time, get_data, get_components, get_child_names
export take_nonblocking!, put_nonblocking!


################################################################################
############################ SimControl ########################################

#this struct will be accessed concurrently by the simulation (within the main
#loop) and GUI (within the call to GUI.update!) threads. therefore, access to it
#must be guarded by acquisition of io_lock
@kwdef mutable struct SimControl
    running::Bool = false #to be checked on each loop iteration for termination
    paused::Bool = false #to pause or unpause the simulation
    pace::Float64 = 1.0 #to be set by the SimControl GUI
    algorithm::String = ""
    t_start::Float64 = 0.0
    t_end::Float64 = 0.0
    dt::Float64 = 0.0 #last time step
    iter::Int64 = 0 #total iterations
    t::Float64 = 0.0 #simulation time
    τ::Float64 = 0.0 #wall-clock time
end

function GUI.draw!(control::SimControl)

    CImGui.Begin("Simulation Control")

        mode_button("Pause", true, false, control.paused; HSV_active = HSV_amber)
        CImGui.IsItemClicked() && (control.paused = !control.paused); CImGui.SameLine()

        dynamic_button("Abort", HSV_red)
        if CImGui.IsItemClicked()
            control.paused = false
            control.running = false
        end
        CImGui.SameLine()

        control.pace = safe_slider("Pace", control.pace, 0.1, 20.0, "%.3f",
        ImGuiSliderFlags_Logarithmic; show_label = true)

        @unpack algorithm, t_start, t_end, dt, iter, t, τ = control

        CImGui.Text("Algorithm: " * algorithm)
        CImGui.Text("Step size: $dt")
        CImGui.Text("Iterations: $iter")
        CImGui.Text(@sprintf("Simulation time: %.3f s", t) * " [$t_start, $t_end]")
        CImGui.Text(@sprintf("Wall-clock time: %.3f s", τ))
        CImGui.Text(@sprintf("GUI framerate: %.3f ms/frame (%.1f FPS)",
                            1000 / unsafe_load(CImGui.GetIO().Framerate),
                            unsafe_load(CImGui.GetIO().Framerate)))


    CImGui.End()

end

################################################################################
############################### SimInterface ###################################

struct SimInterface{D <: IODevice, T, M <: IOMapping}
    device::D
    channel::Channel{T}
    mapping::M
    control::SimControl
    io_start::Base.Event
    io_lock::ReentrantLock
end

const SimInput = SimInterface{<:InputDevice}
const SimOutput = SimInterface{<:OutputDevice}

#called from the SimInterface loop thread. this may block, either on get_data
#because the InputDevice is itself blocked, or on put! because the Simulation
#loop has yet to take!
function update!(input::SimInput)
    put!(input.channel, IODevices.get_data(input.device))
end

#called from the SimInterface loop thread. this may block, either on handle_data
#because the OutputDevice is itself blocked, or on take! because the Simulation
#loop has yet to put!
function update!(output::SimOutput)
    IODevices.handle_data(output.device, take!(output.channel))
end

#called from the Simulation loop thread, will never block
function update!(input::SimInput, sys::System)
    @unpack device, channel, mapping = input
    Systems.assign_data!(sys, take_nonblocking!(channel), device, mapping)
end

#called from the Simulation loop thread, will never block
function update!(output::SimOutput, sys::System)
    @unpack device, channel, mapping = output
    put_nonblocking!(channel, Systems.extract_data(sys, device, mapping))
end

function take_nonblocking!(channel::Channel)
    #unlike lock, trylock avoids blocking if the Channel is already locked,
    #while also ensuring it is not modified while we're checking its state
    if trylock(channel)
        data = (isready(channel) ? take!(channel) : nothing)
        unlock(channel)
        return data
    else
        return nothing
    end
end

function put_nonblocking!(channel::Channel{T}, data::T) where {T}
    #unlike lock, trylock avoids blocking if the Channel is already locked,
    #while also ensuring it is not modified while we're checking its state
    if trylock(channel) #ensure Channel is not modified while we're checking its state
        (isopen(channel) && !isready(channel)) && put!(channel, data)
        unlock(channel)
    end
end


################################################################################
################################# SimGUI #######################################

struct SimGUI
    renderer::Renderer
    control::SimControl
    io_start::Base.Event
    io_lock::ReentrantLock
end


################################################################################
############################# Simulation #######################################

struct Simulation{D <: SystemDefinition, Y, I <: ODEIntegrator}
    sys::System{D, Y}
    integrator::I
    log::SavedValues{Float64, Y}
    gui::SimGUI
    control::SimControl
    io_start::Base.Event #to be waited on by IO devices and GUI before entering their update loops
    io_lock::ReentrantLock #to be acquired before modifying the simulated System or SimControl
    interfaces::Vector{SimInterface}

    function Simulation(
        sys::System{D, Y};
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
        ) where {D, Y}

        @assert (t_end - t_start >= Δt) "Simulation timespan cannot be shorter "* "
                                        than the discrete dynamics update period"

        params = (sys = sys, sys_init! = sys_init!, user_callback! = user_callback!, Δt = Δt)

        cb_step = DiscreteCallback((u, t, integrator)->true, f_cb_step!)
        cb_disc = PeriodicCallback(f_cb_disc!, Δt)
        cb_user = DiscreteCallback((u, t, integrator)->true, f_cb_user!)

        #before initializing the log, compute and assign y, so its initial value
        #is consistent with x, u and s
        f_ode!(sys)
        log = SavedValues(Float64, Y)
        saveat = (saveat isa Real ? (t_start:saveat:t_end) : saveat)
        cb_save = SavingCallback(f_cb_save, log; saveat, save_everystep)

        if save_on
            cb_set = CallbackSet(cb_step, cb_disc, cb_user, cb_save)
        else
            cb_set = CallbackSet(cb_step, cb_disc, cb_user)
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
        control = SimControl()
        io_start = Base.Event()
        io_lock = ReentrantLock()

        f_draw = let control = control, sys = sys
            () -> GUI.draw!(control, sys)
        end

        #the GUI Renderer's refresh rate must be uncapped (no VSync), so that
        #calls to GUI.update!() return immediately without blocking and
        #therefore do not interfere with simulation stepping
        gui = SimGUI(Renderer(; label = "Simulation", sync = UInt8(0), f_draw),
                    control, io_start, io_lock)

        new{D, Y, typeof(integrator)}(
            sys, integrator, log, gui, control, io_start, io_lock, SimInterface[])
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
#step. this function MUST NOT modify the System's ẋ, x or t
user_callback!(::System) = nothing


################################ I/O ###########################################

#the Channel for an OutputInterface must be buffered (size > 0), so that the
#Simulation thread can write at least one output data item without blocking.
#this is because there is no straightforward way to check whether the output
#thread is blocked on a take! call to the output Channel. in contrast, the
#Channel for an InputInterface could be unbuffered (although it does not have to
#be), because by calling isready the simulation loop can tell whether the input
#thread is already blocked on a put! call on that Channel, and may therefore
#take! from it without blocking

function attach!(sim::Simulation, device::IODevice,
                 mapping::IOMapping = DefaultMapping())

    channel = Channel{IODevices.data_type(device)}(1)
    interface = SimInterface(device, channel, mapping,
                            sim.control, sim.io_start, sim.io_lock)

    push!(sim.interfaces, interface)

end


############################ Stepping functions ################################

@inline has_x(sys::System) = !isnothing(sys.x)

#wrapper around the System's continuous dynamics function f_ode! it modifies the
#System's ẋ and y as a side effect
function f_ode_wrapper!(u̇, u, p, t)

    @unpack sys = p

    #assign current integrator solution to System's continuous state
    has_x(sys) && (sys.x .= u)

    #same with time
    sys.t[] = t

    f_ode!(sys) #call continuous dynamics, updates sys.ẋ and sys.y

    has_x(sys) && (u̇ .= sys.ẋ) #update the integrator's derivative

end

#DiscreteCallback function, called on every integration step. calls the System's
#own f_step!. modifications to x or s will not propagate to y until the next
#integration step
function f_cb_step!(integrator)

    @unpack u, p = integrator
    @unpack sys = p

    f_step!(sys)

    #assign the (potentially) modified sys.x back to the integrator
    has_x(sys) && (u .= sys.x)

end

#PeriodicCallback function, calls the System's discrete dynamics update with the
#period Δt given to the Simulation constructor. modifications to x or s will not
#propagate to y until the next integration step
function f_cb_disc!(integrator)

    @unpack u, p = integrator
    @unpack sys, Δt = p

    f_disc!(sys, Δt)

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

function OrdinaryDiffEq.reinit!(sim::Simulation, sys_init! = Systems.init!)

    @unpack sys, integrator, log = sim
    @unpack p = integrator

    #initialize the System's x, u and s
    sys_init!(sys)

    #let it propagate to y. within its reinit! method, the integrator will call
    #the SavingCallback to set the first entry in the log, which is sys.y
    f_ode!(sys)

    #initialize the ODEIntegrator with the System's initial x. ODEIntegrator's
    #reinit! calls f_ode_wrapper!, so ẋ and y are updated in the process
    if has_x(p.sys)
        OrdinaryDiffEq.reinit!(integrator, p.sys.x)
    else
        OrdinaryDiffEq.reinit!(integrator)
    end

    #drop the log entries from the last run and keep the newly initialized one
    resize!(log.t, 1)
    resize!(log.saveval, 1)

    return nothing
end


################################# GUI ##########################################

function GUI.draw!(control::SimControl, sys::System)
    GUI.draw!(control)
    GUI.draw!(sys)
end


################################################################################
############################# Threaded loops ###################################

################################ SimGUI ########################################

function start!(gui::SimGUI)

    @unpack renderer, control, io_start, io_lock = gui

    @info("SimGUI: Starting on thread $(Threads.threadid())...")

    try

        #ensure VSync is disabled. otherwise, the call to GUI.update! will block
        #while holding the system lock, and therefore the sim loop will
        #not be able to step in the meantime
        @assert renderer.sync == 0
        GUI.init!(renderer)

        @info("SimGUI: Waiting for Simulation...")
        wait(io_start)
        @info("SimGUI: Running...")

        while true

            if !(@lock io_lock control.running) || GUI.should_close(renderer)
                @info("SimGUI: Shutting down...")
                break
            end

            #we need a yield here, because the GUI thread never sleeps or blocks
            #by itself, and therefore tends to slow down the simulation thread.
            #this is because the framerate is uncapped, precisely to prevent the
            #GUI thread from blocking for VSync while holding the sys_lock
            yield()
            @lock io_lock GUI.update!(renderer)
        end

    catch ex

        @error("SimGUI: Error during execution: $ex")

    finally
        #abort simulation upon closing the GUI
        @lock io_lock begin
            control.paused = false
            control.running = false
        end
        renderer._initialized && GUI.shutdown!(renderer)
        @info("SimGUI: Closed")
    end

end


################################ SimInterface ##################################

function start!(interface::SimInterface{D}) where {D <: IODevice}

    @unpack device, control, io_start, io_lock = interface

    @info("$D: Starting on thread $(Threads.threadid())...")

    try

        IODevices.init!(device)

        @info("$D: Waiting for Simulation...")
        wait(io_start)
        @info("$D: Running...")

        while true

            if !(@lock io_lock control.running) || IODevices.should_close(device)
                @info("$D: Shutting down...")
                break
            end

            #we need a yield, because this loop is not guaranteed to block at
            #any point, and therefore could slow down the simulation thread
            yield()
            update!(interface)

        end

    catch ex

        @error("$D: Error during execution: $ex")

    finally
        IODevices.shutdown!(device)
        @info("$D: Closed")
    end

end


################################ Simulation ####################################

#compare the following:
# @time sleep(2)

# @time @async sleep(2)
# @time Threads.@spawn sleep(2)

# wait(@time @async sleep(2))
# wait(@time Threads.@spawn sleep(2))

# @time wait(@async sleep(2))
# @time wait(Threads.@spawn @sleep(2))

function start!(sim::Simulation)

    @unpack sys, integrator, gui, control, io_start, io_lock, interfaces = sim

    try

        if isempty(sim.integrator.opts.tstops)
            @error("Simulation has hit its end time, call reinit! ",
                   "or add further tstops using add_tstop!")
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
                    continue #skip next Simulation step, this is spinning
                end

                τ_next = τ_last + get_proposed_dt(sim) / pace
                while τ_next > τ() end #busy wait (should do better, but it's not that easy)

                @lock io_lock step!(sim)
                for interface in interfaces
                    update!(interface, sys)
                end

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

    @unpack sys, control, io_start, io_lock, interfaces = sim

    #if the simulation ran to conclusion, signal IO threads to shut down
    @lock io_lock begin
        control.paused = false
        control.running = false
    end

    # reset(io_start)
    #maybe we should't call reset(io_start) immediately afterwards, because it
    #might be executed before all waiting threads have had time to unblock

    #make sure all IO Channels are emptied so any IO threads do not block on
    #them or they unblock if they were blocked
    for interface in interfaces
        update!(interface, sys)
    end

    #unblock any IO threads still waiting for Simulation start
    notify(io_start)

end


################################################################################
############################### Execution ######################################

#apparently, if a task is launched from the main thread and it doesn't ever
#block, no other thread will get CPU time until it's done. threfore, we should
#never call start! directly from the main thread, because it will starve all IO
#threads. it must always be run from a Threads.@spawn'ed thread

function run!(sim::Simulation)

    sim.control.pace = Inf

    reset(sim.io_start)

    @sync begin
        for interface in sim.interfaces
            Threads.@spawn start!(interface)
        end
        Threads.@spawn start!(sim)
    end

end

function run_interactive!(sim::Simulation; pace = 1.0)

    sim.control.pace = pace

    reset(sim.io_start)

    @sync begin
        for interface in sim.interfaces
            Threads.@spawn start!(interface)
        end
        Threads.@spawn start!(sim.gui)
        Threads.@spawn start!(sim)
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
