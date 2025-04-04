module IODevices

using StaticArrays
using UnPack
using Logging

export IODevice, InputDevice, OutputDevice, IOMapping, GenericMapping
export get_default_mapping


################################################################################
################################# IOMapping ####################################

abstract type IOMapping end
struct GenericMapping <: IOMapping end

################################################################################
################################# IODevice #####################################

abstract type IODevice end

function get_default_mapping(device::IODevice)
    MethodError(get_default_mapping, (device, )) |> throw
end

init!(::D) where {D<:IODevice} = nothing
shutdown!(::D) where {D<:IODevice} = nothing
should_close(::D) where {D<:IODevice} = false

################################################################################
############################### InputDevice ####################################

abstract type InputDevice <: IODevice end

#return a new instance of input data to be assigned. may block
function get_data!(device::D) where {D<:InputDevice}
    MethodError(get_data!, (device, )) |> throw
end

#customize input data assignment through IOMapping subtype dispatch
function assign_input!(target::Any, mapping::IOMapping, data::Any)
    MethodError(assign_input!, (target, mapping, data)) |> throw
end

assign_input!(::Any, ::GenericMapping, ::Any) = nothing

################################################################################
############################## OutputDevice ####################################

abstract type OutputDevice <: IODevice end

#customize output data generation through IOMapping subtype dispatch
function extract_output(source::Any, mapping::IOMapping)
    MethodError(extract_output, (source, mapping)) |> throw
end

extract_output(::Any, ::GenericMapping) = nothing

#process an instance of the extracted output data. may block
function handle_data!(device::D, data::Any) where {D<:OutputDevice}
    MethodError(handle_data!, (device, data)) |> throw
end



end