module Plotting

using Reexport
using Dates
@reexport using Plots
@reexport using Measures
@reexport using LaTeXStrings
@reexport using StructArrays
@reexport using RecursiveArrayTools
@reexport using DataStructures: OrderedDict

#this is how relative imports would work
using ..Misc
using ..Modeling

export make_plots, save_plots


################################################################################

make_plots(::T; kwargs...) where {T<:TimeHistory} = println("Method make_plots not extended for $T")

#these yield a single figure so they can be handled directly by the Plots
#pipeline directly as recipes
make_plots(th::TimeHistory{<:Real}; kwargs...) = plot(th; kwargs...)

make_plots(th::TimeHistory{<:AbstractVector{<:Real}}; kwargs...) = plot(th; kwargs...)

#complex Systems whose outputs are NamedTuples will typically require multiple
#figures, so we cannot use a @recipe for them. we need to handle the TimeHistory
#recursively
function make_plots(th::TimeHistory{<:NamedTuple}; kwargs...)

    pd = OrderedDict{Symbol, Any}()
    for name in Modeling.get_child_names(th)
        child_plots = make_plots(getproperty(th, name); kwargs...)::Union{Nothing, OrderedDict, Plots.Plot}
        !isnothing(child_plots) ? pd[name] = child_plots : nothing
    end

    return (!isempty(pd) ? pd : nothing)

end

function make_plots(mdl::Model;
                    plot_level = :full, #:simplified
                    linewidth=2, margin = 10mm, guidefontsize = 1,
                    kwargs...)
    make_plots(TimeHistory(mdl); plot_level, linewidth, margin, guidefontsize, kwargs...)
end


@recipe function plot(th::TimeHistory{<:Real})

    xguide --> L"$t \: (s)$"
    return th._t, th._y

end

@recipe function plot(th::TimeHistory{<:AbstractVector{<:Real}}; th_split = :none)

    #th._y is a Vector{AbstractVector{<:Real}}; convert it to a matrix
    y_matrix = Array(VectorOfArray(th._y))'

    #number of matrix columns corresponds to the AbstractVector's length
    vlength = size(y_matrix)[2]

    xguide --> L"$t \ (s)$"

    label --> (vlength <= 3 ?  ["x" "y" "z"][:, 1:vlength] : (1:vlength)')
    if th_split === :h
        layout --> (1, vlength)
        link --> :y #alternative: :none
    elseif th_split === :v
        layout --> (vlength, 1)
    else
        layout --> 1
    end

    return th._t, y_matrix

end


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
            println("Saved figure $plot_filename")

        elseif !isnothing(child)
            error("Invalid entry type ($(typeof(child))")

        end
    end

end


end #module