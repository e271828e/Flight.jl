module IODevices

using StaticArrays
using UnPack
using Logging

export InputDevice, InputMapping, DefaultMapping, InputInterface
export OutputDevice, OutputInterface
export put_no_wait!


################################################################################
################################# IODevice #####################################

abstract type IODevice end

#each device should control its own update rate from within the update! method,
#preferably via calls to blocking functions (such as GLFW.SwapBuffers with
#GLFW.SwapInterval > 0)
init!(device::D) where {D<:IODevice} = MethodError(init!, (device, )) |> throw
shutdown!(device::D) where {D<:IODevice} = MethodError(shutdown!, (device, )) |> throw
should_close(device::D) where {D<:IODevice} = MethodError(should_close, (device, )) |> throw

################################################################################
############################ InputInterface ####################################

abstract type InputDevice <: IODevice end
abstract type InputMapping end
struct DefaultMapping <: InputMapping end

#for every input device it wants to support, the target type should at least
#extend the three argument assign! method below with a ::DefaultMapping
#argument. for alternative mappings, it can define additional subtypes of
#InputMapping and their corresponding assign! methods
function assign!(target::Any, device::InputDevice, mapping::InputMapping)
    @warn("Assigment method for target $(typeof(target)) from device"*
    "$(typeof(device)) with mapping $(typeof(mapping)) not implemented")
end

function assign!(target, device::InputDevice)
    assign!(target, device, DefaultMapping())
end

#an input device must control its own update rate through calls to blocking
#functions from within the update! method (such as GLFW.SwapBuffers with
#GLFW.SwapInterval > 0 or a socket receive)
update!(device::D) where {D<:InputDevice} = MethodError(update!, (device, )) |> throw

struct DummyInputDevice <: InputDevice end

init!(::DummyInputDevice) = println("DummyInputDevice initialized")
assign!(::Any, ::DummyInputDevice, ::DefaultMapping) = nothing
update!(::DummyInputDevice) = (sleep(1); println("DummyInputDevice updated"))
shutdown!(::DummyInputDevice) = println("DummyInputDevice shutting down")
should_close(::DummyInputDevice) = false

################################################################################

struct InputInterface{D <: InputDevice, T,  M <: InputMapping}
    device::D
    target::T #target for input assignment, typically the simulated System
    mapping::M #selected device-to-target mapping, used for dispatch on assign!
    start::Base.Event #to be waited on before entering the update loop
    running::ReentrantLock #to be checked on each loop iteration for termination
    stepping::ReentrantLock #to acquire before modifying the target
end

function start!(interface::InputInterface{D}) where {D}

    @unpack device, target, mapping, start, running, stepping = interface

    @info("$D Interface: Starting on thread $(Threads.threadid())...")

    try

        init!(device)
        wait(start)

        while true

            if !islocked(running) || should_close(device)
                @info("$D Interface: Shutdown requested")
                break
            end

            update!(device)
            lock(stepping) #ensure the Simulation is not currently stepping
                assign!(target, device, mapping)
            unlock(stepping)

        end

    catch ex

        @error("$D Interface: Error during execution: $ex")

    finally
        shutdown!(device)
        @info("$D Interface: Closed")
    end

end

################################################################################
############################ OutputInterface ###################################

abstract type OutputDevice <: IODevice end

update!(device::D, data) where {D<:OutputDevice} = MethodError(update!, (device, data)) |> throw

struct DummyOutputDevice <: OutputDevice end

init!(::DummyOutputDevice) = println("DummyOutputDevice initialized")
update!(::DummyOutputDevice, data) = (sleep(1); println("DummyOutputDevice updated"))
shutdown!(::DummyOutputDevice) = println("DummyOutputDevice shutting down")
should_close(::DummyOutputDevice) = false

################################################################################

struct OutputInterface{D <: OutputDevice, C <: Channel}
    device::D
    channel::C #channel to which Simulation's output will be put!
    start::Base.Event #to be waited on before entering the update loop
    running::ReentrantLock #to be checked on each loop iteration for termination
end

#for buffered channels, isready returns 1 if there is at least one value in the
#channel. if there is one value and we try to put! another, the simulation loop
#will block. to avoid this, we should only put! if our channel !isready
@inline function put_no_wait!(channel::Channel{T}, data::T) where {T}
    (isopen(channel) && !isready(channel)) && put!(channel, data)
end

@inline function put_no_wait!(interface::OutputInterface, data)
    put_no_wait!(interface.channel, data)
end

#the update rate for an output device is implicitly controlled by the simulation
#loop. the call to take! below will block until the simulation writes a new
#output to the channel
function start!(interface::OutputInterface{D}) where {D}

    @unpack device, channel, start, running = interface

    @info("$D Interface: Starting on thread $(Threads.threadid())...")

    try

        init!(device)
        wait(start)

        while true

            #if sim loop is not running, break immediately to avoid getting
            #stuck on take!
            if !islocked(running) || should_close(device)
                @info("$D Interface: Shutdown requested")
                break
            end

            data = take!(channel)
            update!(device, data)

        end

    catch ex

        if ex isa InvalidStateException
            @info("$D Interface: Channel closed")
        else
            @error("$D Interface: Error during execution: $ex")
        end

    finally
        # close(channel)
        shutdown!(device)
        @info("$D Interface: Closed")
    end

end



end