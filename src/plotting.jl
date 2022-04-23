module Plotting

using Reexport
using Dates
@reexport using Plots
@reexport using Measures
@reexport using StructArrays
@reexport using RecursiveArrayTools
@reexport using LaTeXStrings

using Flight.Modeling

export TimeHistory
export plots, save_plots, thplot, thplot!
export make_plots

# ########################## plots method ###############################

#our entry plotting entry point cannot be a recipe. a recipe is called within
#the plot() pipeline, which creates a single figure. however, a Vector of System
#outputs will typically need to generate multiple plots from its values. for
#this we define a new function plots().

plots(args...; kwargs...) = println("Not implemented")

#for Systems without outputs
plots(::AbstractVector{<:Real}, ::AbstractVector{Nothing}; kwargs...) = nothing

#takes care of SystemGroups
function plots(t, data::AbstractVector{<:NamedTuple}; mode, save_path, kwargs...)

    c = data |> StructArray |> StructArrays.components
    for (c_label, c_data) in zip(keys(c), values(c))
        println("Generating plots for $c_label")
        save_path_c = mkpath(joinpath(save_path, String(c_label)))
        plots(t, c_data; mode, save_path = save_path_c, kwargs...)
    end

end

function save_plots(d::Dict{String,Plots.Plot}; save_path, format = :png)
    for (id, p) in zip(keys(d), values(d))
        savefig(p, joinpath(save_path, id*"."*String(format)))
    end
end

# ############################ TimeHistory #################################

mutable struct TimeHistory{D}
    t::AbstractVector{<:Real}
    data::D
end

#for all the following TimeHistory subtypes, a single figure is enough, so they
#can be handled by recipes
@recipe function f(th::TimeHistory{<:AbstractVector{<:Real}})

    xguide --> L"$t \: (s)$"
    return th.t, th.data

end

@recipe function f(th::TimeHistory{<:AbstractMatrix{<:Real}}; th_split = :none)

    xguide --> L"$t \ (s)$"

    vlength = size(th.data)[2]
    label --> (vlength <= 3 ?  ["x" "y" "z"][:, 1:vlength] : (1:vlength)')
    if th_split === :h
        layout --> (1, vlength)
        link --> :y #alternative: :none
    elseif th_split === :v
        layout --> (vlength, 1)
    else
        layout --> 1
    end

    return th.t, th.data

end

#convert to TimeHistory{Matrix} and return it to the pipeline for dispatching
@recipe function f(th::TimeHistory{<:AbstractVector{<:AbstractVector{<:Real}}})
    return TimeHistory(th.t, Array(VectorOfArray(th.data))')
end

thplot(t, data; kwargs...) = plot(TimeHistory(t, data); kwargs...)
thplot!(t, data; kwargs...) = plot!(TimeHistory(t, data); kwargs...)


################################################################################

# function make_plots(mdl::Model; mode::Symbol = :basic, kwargs...)
#     make_plots(THNew(mdl); mode, save_path, kwargs...)
# end

# #try to dispatch to the Plots method; this is useful when assembling
# make_plots(args...; kwargs...) = plot(args...; kwargs...)

make_plots(::THNew; kwargs...) = println("Not implemented")

#complex Systems whose outputs are NamedTuples will typically require multiple
#figures, so we cannot use a @recipe for them. we need to handle the THNew
#recursively
function make_plots(th::THNew{<:NamedTuple}; mode, kwargs...)

    pd = Dict{Symbol, Any}()
    for name in Modeling.get_child_names(th)
        println("Generating plots for $name")
        #this will return either nothing, a Dict or a Plots.Plot instance
        child_plots = make_plots(getproperty(th, name); mode, kwargs...)
        !isnothing(child_plots) ? pd[name] = child_plots : nothing
    end

    return pd

end

# #remove this?
# make_plots(::THNew{Nothing}; kwargs...) = nothing

#these produce a single figure so they can be passed to the Plots pipeline
#directly as recipes
make_plots(th::THNew{<:Real}; kwargs...) = plot(th; kwargs...)
make_plots(th::THNew{<:AbstractVector{<:Real}}; kwargs...) = plot(th; kwargs...)

@recipe function f(th::THNew{<:Real})

    xguide --> L"$t \: (s)$"
    return th._t, th._y

end

@recipe function f(th::THNew{<:AbstractVector{<:Real}}; th_split = :none)

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


function save_plots(d::Dict{Symbol, T} where {T <: Any};
                    save_path::Union{String, Nothing} = nothing, format = :png)

    save_path = (save_path === nothing ?
        joinpath("tmp", Dates.format(now(), "yyyy_mm_dd_HHMMSS")) : save_path)

    for (label, child) in zip(keys(d), values(d))
        child_save_path = joinpath(save_path, String(label))
        if isa(child, Dict)
            mkpath(child_save_path)
            save_plots(d[label]; save_path = child_save_path, format)
        elseif isa(child, Plots.Plot)
            savefig(child, child_save_path*"."*String(format))
        else
            error("Invalid Dict entry type")
        end
    end

end

function save_plots(plt::Plots.Plot;
                    save_path::Union{String,Nothing} = nothing, format = :png)

    save_path = (save_path === nothing ?
        joinpath("tmp", Dates.format(now(), "yyyy_mm_dd_HHMMSS")) : save_path)

    savefig(plt, save_path*"."*String(format))
end


end