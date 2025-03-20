module Plotting

using Reexport
using Dates
using UnPack
using StructArrays
using Logging
using RecursiveArrayTools: VectorOfArray
@reexport using Plots
@reexport using LaTeXStrings
@reexport using DataStructures: OrderedDict

@reexport using ..Sim: Simulation, TimeSeries, get_components, get_child_names
using ..Types

export make_plots, save_plots

const defaults = Dict(
    :size => (1920, 1080),
    :dpi => 100,
    :margin => 10Plots.mm,
    :topmargin => 5Plots.mm,
    :linewidth => 2,
    :plot_titlefontfamily => "Computer Modern",
    :titlefontfamily => "Computer Modern",
    :fontfamily_subplot => "Computer Modern", #controls guide, tick and legend font
    :plot_titlefontsize => 20,
    :titlefontsize => 16,
    :guidefontsize => 12,
    :tickfontsize => 10,
    :legendfontsize => 12,
    )

################################################################################
################################# Recipes ######################################

@recipe function f(ts::TimeSeries{<:Real})

    xguide --> L"$t \: (s)$"
    return ts._t, ts._data

end

@recipe function f(ts::TimeSeries{<:Ranged{T}}) where {T}

    return ts._t, T.(ts._data)

end

@recipe function f(ts::TimeSeries{<:AbstractVector{<:Real}}; ts_split = :none)

    #ts._data is a Vector{AbstractVector{<:Real}}; convert it to a matrix
    data = Array(VectorOfArray(ts._data))'

    #number of matrix columns corresponds to the AbstractVector's length
    vlength = size(data)[2]

    xguide --> L"$t \ (s)$"

    label --> (vlength <= 3 ?  ["x" "y" "z"][:, 1:vlength] : (1:vlength)')
    if ts_split === :h
        layout --> (1, vlength)
        link --> :y #alternative: :none
    elseif ts_split === :v
        layout --> (vlength, 1)
    else
        layout --> 1
    end

    return ts._t, data

end

################################################################################
######################### Multi-Plot Specifications ############################

make_plots(::T; kwargs...) where {T<:TimeSeries} = @warn("Method make_plots not extended for $T")

#these yield a single figure so they can be directly handled by the Plots
#pipeline as recipes
make_plots(ts::TimeSeries{<:Real}; kwargs...) = plot(ts; kwargs...)

make_plots(ts::TimeSeries{<:AbstractVector{<:Real}}; kwargs...) = plot(ts; kwargs...)

#complex Systems whose outputs are NamedTuples will typically require multiple
#figures, so we cannot use a @recipe for them. we need to handle the TimeSeries
#recursively
function make_plots(ts::TimeSeries{<:NamedTuple}; kwargs...)

    pd = OrderedDict{Symbol, Any}()
    for name in get_child_names(ts)
        child_plots = make_plots(getproperty(ts, name); kwargs...)::Union{Nothing, OrderedDict, Plots.Plot}
        !isnothing(child_plots) && (pd[name] = child_plots)
    end

    return (!isempty(pd) ? pd : nothing)

end

function make_plots(sim::Simulation;
                    plot_level = :full, #:simplified
                    plot_settings...)
    make_plots(TimeSeries(sim); plot_level, plot_settings...)
end


################################################################################
############################# Plot Saving ######################################

function save_plots(dict::OrderedDict{Symbol, T} where {T};
                    save_folder::Union{String, Nothing} = nothing, format = :png)

    save_folder = mkpath(save_folder === nothing ?
        joinpath("tmp", Dates.format(now(), "yyyy_mm_dd_HHMMSS")) : save_folder)

    n = 0
    for (label, child) in zip(keys(dict), values(dict))

        if isa(child, OrderedDict)
            save_subfolder = mkpath(joinpath(save_folder, String(label)))
            save_plots(dict[label]; save_folder = save_subfolder, format)

        elseif isa(child, Plots.Plot)
            n += 1
            plot_filename = joinpath(save_folder, string(n, pad = 2)*"_"*String(label)*"."*String(format))
            savefig(child, plot_filename)
            @info("Saved figure $plot_filename")

        elseif !isnothing(child)
            @error("Invalid entry type ($(typeof(child))")

        end
    end

end


end #module