module Sim

using UnPack, Reexport, StructArrays
using OrdinaryDiffEq: OrdinaryDiffEq, OrdinaryDiffEqAlgorithm, ODEProblem,
                      ODEIntegrator, Heun, RK4, u_modified!, init as init_integrator
using DiffEqCallbacks: SavingCallback, DiscreteCallback, PeriodicCallback,
                       CallbackSet, SavedValues
using RecursiveArrayTools
using Logging
using CImGui.lib: ImGuiSliderFlags_Logarithmic

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
#loop) and GUI (within the call to GUI.update!), so access to it must be guarded
#by io_lock
@kwdef mutable struct SimControl
    running::Bool = false #to be checked on each loop iteration for termination
    paused::Bool = false #to pause or unpause the simulation
    pace::Float64 = 1.0 #to be set by the SimControl GUI
    algorithm::String = ""
    t_start::Float64 = 0.0
    t_end::Float64 = 0.0
    Δt::Float64 = 0.0 #discrete step size
    dt::Float64 = 0.0 #last continuous integration step size
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

        @unpack algorithm, t_start, t_end, Δt, dt, iter, t, τ = control

        CImGui.Text("Algorithm: " * algorithm)
        CImGui.Text("Discrete step size: $Δt")
        CImGui.Text("Continuous step size: $dt")
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

struct SimInterface{D <: IODevice, S <: System, M <: IOMapping}
    device::D
    sys::S
    mapping::M
    control::SimControl
    io_start::Base.Event
    io_lock::ReentrantLock
    should_abort::Bool #aborts Simulation when closed
end

const SimInput{D} = SimInterface{D} where {D <: InputDevice}
const SimOutput{D} = SimInterface{D} where {D <: OutputDevice}
const SimGUI{D} = SimInterface{D} where {D <: Renderer}


#called from the SimInterface loop thread. this may block, either on get_data
#because the InputDevice is itself blocked, on or lock(io_lock) because the
#Simulation loop is currently stepping
function update!(interface::SimInput)

    @unpack device, sys, mapping, io_lock = interface
    data = IODevices.get_data!(device)

    #the data we got from the InputDevice might be something we're not able to
    #assign to the System (for example, an EOT character). we need to handle
    #that scenario
    try
        lock(io_lock)
        Systems.assign_input!(sys, data, mapping)
    catch ex
        @warn("Failed to assign input data $data to System")
        # println(ex)
    finally
        unlock(io_lock)
    end

end


#called from the SimInterface loop thread. this may block, either on handle_data
#because the OutputDevice is itself blocked, or on lock(io_lock) because the
#Simulation loop is currently stepping
function update!(interface::SimOutput{<:OutputDevice{T}}) where {T}

    @unpack device, sys, mapping, io_lock = interface

    lock(io_lock)
        #this call should never block and always return some usable output
        data = Systems.extract_output(sys, T, mapping)
    unlock(io_lock)

    IODevices.handle_data!(device, data)
end


function update!(gui::SimGUI)
    @unpack device, io_lock = gui
    #this may modify the System or SimControl, so we need to grab the IO_lock
    @lock io_lock GUI.render!(device)
end



################################################################################
############################# Simulation #######################################

struct Simulation{D <: SystemDefinition, Y, I <: ODEIntegrator, G <: SimGUI}
    sys::System{D, Y}
    integrator::I
    log::SavedValues{Float64, Y}
    gui::G
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

        t_end - t_start < Δt && @warn(
        "Simulation timespan is shorter than the discrete dynamics update period")

        params = (sys = sys, sys_init! = sys_init!, user_callback! = user_callback!)

        #initial_affect=true causes f_disc! to be called before the initial step
        cb_step = DiscreteCallback((u, t, integrator)->true, f_cb_step!)
        cb_disc = PeriodicCallback(f_cb_disc!, Δt; initial_affect = false)
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


        #when the integrator is instantiated, if initial_affect==true for
        #cb_disc, the integrator will call f_cb_disc! and sys.n will be
        #incremented. if we then mak reset System's root Δt and discrete iteration counter. this
        #must be done after the integrator is instantiated, because
        sys.Δt_root[] = Δt
        sys.n[] = 0

        #save_on = false because we are not interested in logging the plain
        #System's state vector; everything we need to know about the System
        #should be made available in its (immutable) output y, logged via
        #SavingCallback
        control = SimControl()
        io_start = Base.Event()
        io_lock = ReentrantLock()

        #the GUI Renderer's refresh rate must be uncapped (no VSync), so that
        #calls to GUI.update!() return immediately without blocking and
        #therefore do not interfere with simulation stepping
        f_draw = let control = control, sys = sys
            () -> GUI.draw!(control, sys)
        end
        #sync must be set to 0 to disable VSync. otherwise, the call to update!
        #will block for a whole display refresh interval while holding the
        #io_lock, leaving the sim loop unable to step in the meantime
        renderer = Renderer(; label = "Simulation", sync = UInt8(0), f_draw)
        gui = SimInterface(renderer, sys, DefaultMapping(),
                           control, io_start, io_lock, true)

        new{D, Y, typeof(integrator), typeof(gui)}(
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
                mapping::IOMapping = DefaultMapping(); should_abort = false)

    interface = SimInterface(device, sim.sys, mapping,
                            sim.control, sim.io_start, sim.io_lock, should_abort)

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
    @unpack sys = p

    f_disc!(sys)

    #increment the discrete iteration counter
    sys.n[] += 1

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

function OrdinaryDiffEq.reinit!(sim::Simulation, init_args...; init_kwargs...)

    @unpack sys, integrator, log = sim
    @unpack p = integrator

    #reset scheduling counter (must be done before Systems.init!, in case the
    #init methods need to call f_disc!)
    sys.n[] = 0

    #initialize the System's x, u and s
    Systems.init!(sys, init_args...; init_kwargs...)

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

    #reset scheduling counter (also AFTER integrator reinit!, which may have
    #called f_disc_cb! if initial_affect==true, incrementing the counter)
    sys.n[] = 0

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
            #block at any point, and therefore could slow down the simulation
            #thread significantly. in the case of SimGUI, this is by design: the
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
            control.Δt = sim.Δt_root[]
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

end


################################################################################
############################### Execution ######################################

#apparently, if a task is launched from the main thread and it doesn't ever
#block, no other thread will get CPU time until it's done. threfore, we should
#never call start! directly from the main thread unless we know the task will
#block or yield at some point.

#important: on platforms other than Windows, CImGui needs to run on the main
#thread!

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
