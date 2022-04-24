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

make_plots(::T; kwargs...) where {T<:THNew} = println("Method make_plots not extended for $T")

#these produce a single figure so they can be handled by the Plots pipeline
#directly as recipes
make_plots(th::THNew{<:Real}; kwargs...) = plot(th; kwargs...)

make_plots(th::THNew{<:AbstractVector{<:Real}}; kwargs...) = plot(th; kwargs...)

#complex Systems whose outputs are NamedTuples will typically require multiple
#figures, so we cannot use a @recipe for them. we need to handle the THNew
#recursively
function make_plots(th::THNew{<:NamedTuple}; mode, kwargs...)

    pd = Dict{Symbol, Any}()
    for name in Modeling.get_child_names(th)
        child_plots = make_plots(getproperty(th, name); mode, kwargs...)::Union{Nothing, NamedTuple, Plots.Plot}
        !isnothing(child_plots) ? pd[name] = child_plots : nothing
    end

    return NamedTuple(pd)

end

# function make_plots(mdl::Model; mode::Symbol = :basic, kwargs...)
#     make_plots(THNew(mdl); mode, save_path, kwargs...)
# end



@recipe function plot(th::THNew{<:Real})

    xguide --> L"$t \: (s)$"
    return th._t, th._y

end

@recipe function plot(th::THNew{<:AbstractVector{<:Real}}; th_split = :none)

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


function save_plots(nt::NamedTuple{K, T} where {K, T <: NTuple{N, Union{NamedTuple, Plots.Plot}} where {N}};
                    save_folder::Union{String, Nothing} = nothing, format = :png)

    save_folder = mkpath(save_folder === nothing ?
        joinpath("tmp", Dates.format(now(), "yyyy_mm_dd_HHMMSS")) : save_folder)

    n = 0
    for (label, child) in zip(keys(nt), values(nt))

        if isa(child, NamedTuple)
            save_subfolder = mkpath(joinpath(save_folder, String(label)))
            save_plots(nt[label]; save_folder = save_subfolder, format)

        elseif isa(child, Plots.Plot)
            n += 1
            plot_filename = joinpath(save_folder, string(n, pad = 2)*"_"*String(label)*"."*String(format))
            savefig(child, plot_filename)
            println("Saved figure $plot_filename")

        else
            error("Invalid entry type")
        end
    end

end


end #module