module IODevices

using Logging

export IODevice, InputDevice, OutputDevice
export IOMapping, GenericInputMapping
export get_default_mapping
export InputMappingError


################################################################################
################################# IOMapping ####################################

abstract type IOMapping end

struct GenericInputMapping <: IOMapping end

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

"""
    InputMappingError(cause)

Signals that an incoming datum could not be mapped onto its target for a reason
that is *not* a programming error (e.g. a malformed or out-of-range network
packet). An `assign_input!` implementation that wants such a datum tolerated
—logged and skipped, leaving the input interface alive— should catch the
underlying parsing/conversion exception and rethrow it as an `InputMappingError`.
Any other exception propagates, terminating the input interface (and the
Simulation, if its `should_abort` flag is set).
"""
struct InputMappingError <: Exception
    cause::Exception
end

Base.showerror(io::IO, e::InputMappingError) =
    (print(io, "InputMappingError caused by: "); showerror(io, e.cause))

################################################################################
############################## OutputDevice ####################################

abstract type OutputDevice <: IODevice end

#customize output data generation through IOMapping subtype dispatch
function extract_output(source::Any, mapping::IOMapping)
    MethodError(extract_output, (source, mapping)) |> throw
end

#process an instance of the extracted output data. may block
function handle_data!(device::D, data::Any) where {D<:OutputDevice}
    MethodError(handle_data!, (device, data)) |> throw
end



end
