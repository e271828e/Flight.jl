module Plotting

export make_plots, save_plots

#default plotting options
const defaults = Dict{Symbol, Any}()

#stubs
function make_plots(args...; kwargs...)
    error("Plotting functionality is not loaded. Please run `using Plots` to enable `make_plots`.")
end

function save_plots(args...; kwargs...)
    error("Plotting functionality is not loaded. Please run `using Plots` to enable `save_plots`.")
end

end #module
