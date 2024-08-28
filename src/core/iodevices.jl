module IODevices

using StaticArrays
using UnPack
using Logging

export InputDevice, OutputDevice, IOMapping, DefaultMapping


################################################################################
################################# IODevice #####################################

abstract type IODevice end

abstract type IOMapping end
struct DefaultMapping <: IOMapping end

#each device should control its own update rate from within the update! method,
#preferably via calls to blocking functions (such as GLFW.SwapBuffers with
#GLFW.SwapInterval > 0)
init!(device::D) where {D<:IODevice} = MethodError(init!, (device, )) |> throw
shutdown!(device::D) where {D<:IODevice} = MethodError(shutdown!, (device, )) |> throw
should_close(device::D) where {D<:IODevice} = MethodError(should_close, (device, )) |> throw

#returns a Channel suitable for the type of data produced or expected by the
#IODevice. this should be called by the SimInput or SimOutput constructors.
function Base.Channel(device::D, size::Int) where {D <: IODevice}
    MethodError(Base.Channel, (device, size)) |> throw
end


################################################################################
############################### InputDevice ####################################

abstract type InputDevice <: IODevice end

#returns a new instance of input data. may block
function get_data(device::D) where {D<:InputDevice}
    MethodError(get_data, (device, )) |> throw
end


################################################################################
############################## OutputDevice ####################################

abstract type OutputDevice <: IODevice end

#processes an instance of output data. may block
function handle_data(device::D, data::Any) where {D<:OutputDevice}
    MethodError(handle_data, (device, data)) |> throw
end


end