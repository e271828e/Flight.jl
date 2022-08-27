module IODevices

using StaticArrays
using UnPack

export IODevice, InputMapping, DefaultMapping


################################################################################
################################# IODevice #####################################

abstract type IODevice end
abstract type InputMapping end
struct DefaultMapping <: InputMapping end

#each device should control its own update rate from within the update! method,
#preferably via calls to blocking functions (such as GLFW.SwapBuffers with
#GLFW.SwapInterval > 0)
init!(device::D) where {D<:IODevice} = MethodError(init!, (device, ))
shutdown!(device::D) where {D<:IODevice} = MethodError(shutdown!, (device, ))
should_close(device::D) where {D<:IODevice} = MethodError(should_close, (device, ))
update!(device::D, data) where {D<:IODevice} = MethodError(update!, (device, data))

#for every input device it wants to support, the target type should at least
#extend the three argument assign! method below with a ::DefaultMapping
#argument. for alternative mappings, it can define additional subtypes of
#InputMapping and their corresponding assign! methods
function assign!(target, device::IODevice, mapping::InputMapping)
    MethodError(update!, (target, device, mapping))
end

function assign!(target, device::IODevice)
    assign!(target, device, DefaultMapping())
end


################################################################################
################################ Interface #####################################

#add mode: IOMode, Input <: IOMode, Output <: IOMode, InputOutput <: IOMode
struct Interface{D <: IODevice, T,  M <: InputMapping, C <: Channel}
    device::D
    target::T #input target, typically the simulation's u
    mapping::M #selected device to target mapping, used for dispatch on assign!
    channel::C
    start::Base.Event #to be waited on before entering the update loop
    target_lock::ReentrantLock #to acquire before modifying the target
    ext_shutdown::Bool #whether to observe shutdown requests received by the IO device
end

start!(io::Interface; verbose = true) = (Threads.@spawn _start!(io; verbose))

function _start!(io::Interface{D}; verbose = true) where {D}

    @unpack device, target, mapping, channel, start, sim_stepping, ext_shutdown = io

    verbose && println("$D: Starting on thread $(Threads.threadid())...")

    init!(device)
    wait(start)

    try

        while true

            data = take!(channel)
            update!(device, data)

            lock(target_lock) #ensure the Simulation is not currently stepping
                assign!(target, device, mapping)
            unlock(target_lock)

            if ext_shutdown && should_close(device)
                println("$D: Shutdown requested")
                break
            end

        end

    catch ex

        if ex isa InvalidStateException
            println("$D: Channel closed")
        else
            println("$D: Error during execution: $ex")
        end

    finally
        println("$D: Exiting...")
        shutdown!(device)
    end

end


end