module IODevices

using StaticArrays
using UnPack
using Logging

export InputDevice, InputMapping, DefaultMapping
export OutputDevice
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
############################### InputDevice ####################################

abstract type InputDevice <: IODevice end
abstract type InputMapping end
struct DefaultMapping <: InputMapping end

#an input device must control its own update rate through calls to blocking
#functions from within the update_input! method (such as GLFW.SwapBuffers with
#GLFW.SwapInterval > 0 or a socket receive)
function update_input!(device::D) where {D<:InputDevice}
    MethodError(update_input!, (device, )) |> throw
end

#for every input device it wants to support, the target type should at least
#extend the three argument assign! method below with a ::DefaultMapping
#argument. for alternative mappings, it can define additional subtypes of
#InputMapping and their corresponding assign! methods
function assign_input!(target::Any, device::InputDevice, mapping::InputMapping)
    @warn("Assigment method for target $(typeof(target)) from device"*
    "$(typeof(device)) with mapping $(typeof(mapping)) not implemented")
end

function assign_input!(target, device::InputDevice)
    assign_input!(target, device, DefaultMapping())
end

struct DummyInputDevice <: InputDevice end

init!(::DummyInputDevice) = println("DummyInputDevice initialized")
update_input!(::DummyInputDevice) = (sleep(1); println("DummyInputDevice updated"))
assign_input!(::Any, ::DummyInputDevice, ::DefaultMapping) = nothing
shutdown!(::DummyInputDevice) = println("DummyInputDevice shutting down")
should_close(::DummyInputDevice) = false


################################################################################
############################### OutputDevice ###################################

abstract type OutputDevice <: IODevice end

function process_output!(device::D, data) where {D<:OutputDevice}
    MethodError(process_output!, (device, data)) |> throw
end

struct DummyOutputDevice <: OutputDevice end

init!(::DummyOutputDevice) = println("DummyOutputDevice initialized")
process_output!(::DummyOutputDevice, data) = (sleep(1); println("DummyOutputDevice updated"))
shutdown!(::DummyOutputDevice) = println("DummyOutputDevice shutting down")
should_close(::DummyOutputDevice) = false



end