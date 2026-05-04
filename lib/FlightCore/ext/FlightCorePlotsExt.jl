module FlightCorePlotsExt

using Plots, LaTeXStrings
using Dates, DataStructures
using RecursiveArrayTools: VectorOfArray

using FlightCore.Sim: Simulation, TimeSeries
using FlightCore.Plotting

function __init__()
    # Populate the defaults dict when Plots is loaded
    merge!(Plotting.defaults, Dict(
        :size => (1920, 1080),
        :dpi => 100,
        :margin => 10Plots.mm,
        :topmargin => 5Plots.mm,
        :linewidth => 2,
        :markersize => 8,
        :plot_titlefontfamily => "Computer Modern",
        :titlefontfamily => "Computer Modern",
        :fontfamily_subplot => "Computer Modern",
        :plot_titlefontsize => 20,
        :titlefontsize => 16,
        :guidefontsize => 12,
        :tickfontsize => 10,
        :legendfontsize => 12,
    ))
end

################################################################################
################################# Recipes ######################################

@recipe function f(ts::TimeSeries{<:Real})

    xguide --> L"$t \: (s)$"
    return ts._t, ts._data

end

@recipe function f(ts::TimeSeries{<:AbstractVector{<:Real}}; ts_split=:none)

    #ts._data is a Vector{AbstractVector{<:Real}}; convert it to a matrix
    data = Array(VectorOfArray(ts._data))'

    #number of matrix columns corresponds to the AbstractVector's length
    vlength = size(data)[2]

    xguide --> L"$t \ (s)$"

    label --> (vlength <= 3 ? ["x" "y" "z"][:, 1:vlength] : (1:vlength)')
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
################################# Extensions ###################################

Plotting.make_plots(sim::Simulation; kwargs...) = make_plots(TimeSeries(sim); kwargs...)
Plotting.make_plots(::T; kwargs...) where {T<:TimeSeries} = @info("make_plots not implemented for $T")

#these yield a single figure so they can be directly handled by the Plots
#pipeline as recipes
Plotting.make_plots(ts::TimeSeries{<:Real}; kwargs...) = plot(ts; kwargs...)
Plotting.make_plots(ts::TimeSeries{<:AbstractVector{<:Real}}; kwargs...) = plot(ts; kwargs...)

#complex Models with NamedTuple outputs will generally require multiple
#figures, so we cannot use a @recipe for them. we need to handle the TimeSeries
#recursively
function Plotting.make_plots(ts::TimeSeries{<:NamedTuple}; kwargs...)

    pd = OrderedDict{Symbol,Any}()
    for name in propertynames(ts)
        child_plots = make_plots(getproperty(ts, name); kwargs...)::Union{Nothing,OrderedDict,Plots.Plot}
        !isnothing(child_plots) && (pd[name] = child_plots)
    end

    return (!isempty(pd) ? pd : nothing)

end

####################################################################################

Plotting.save_plots(sim::Simulation, args...; kwargs...) = save_plots(TimeSeries(sim), args...; kwargs...)

function Plotting.save_plots(ts::TimeSeries, args...; kwargs...)
    pd = make_plots(ts; kwargs...)
    save_plots(pd, args...; kwargs...)
end


function Plotting.save_plots(dict::OrderedDict{Symbol,T} where {T},
    folder::Union{String,Nothing}=nothing,
    format::Symbol=:png;
    kwargs...)

    folder = mkpath(!isnothing(folder) ? folder : joinpath("tmp", Dates.format(now(), "yyyy_mm_dd_HHMMSS")))

    for (index, (label, child)) in enumerate(zip(keys(dict), values(dict)))

        if isa(child, OrderedDict)
            subfolder = mkpath(joinpath(folder, String(label)))
            save_plots(dict[label], subfolder, format; kwargs...)

        elseif isa(child, Plots.Plot)
            plot_filename = joinpath(folder, string(index, pad=2) * "_" * String(label) * "." * String(format))
            savefig(child, plot_filename)
            @info("Saved figure $plot_filename")

        elseif !isnothing(child)
            @error("Invalid entry type: ($(typeof(child))")

        end
    end

end

end #module
