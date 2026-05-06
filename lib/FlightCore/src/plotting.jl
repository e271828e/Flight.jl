module Plotting

export make_plots, save_plots

#default plotting options
const defaults = Dict{Symbol, Any}()

#stubs
function make_plots(args...; kwargs...)
    MethodError(make_plots, args) |> throw
end

function save_plots(args...; kwargs...)
    MethodError(save_plots, args) |> throw
end

end #module
