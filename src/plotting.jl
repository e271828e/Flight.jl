module Plotting

using Reexport
@reexport using Plots
@reexport using StructArrays
@reexport using RecursiveArrayTools
@reexport using LaTeXStrings

using Flight.Utils
using Flight.ModelingTools

export TimeHistory
export plots, save_plots, thplot, thplot!

########################## plots method ###############################

#our entry plotting entry point cannot be a recipe. a recipe is called within
#the plot() pipeline, which creates a single figure. however, a Vector of System
#outputs will typically need to generate multiple plots from its values. for
#this we define a new method plots(). maybe define different detail levels /
#modes for

plots(args...; kwargs...) = println("Not implemented")

function plots(mdl::Model; mode::Symbol = :basic,
    save_path::Union{String,Nothing} = nothing, kwargs...)
    #generate default path tmp/plots/current_date
    save_path = (save_path === nothing ?
        joinpath("tmp", Dates.format(now(), "yyyy_mm_dd_HHMMSS")) : save_path)
    mkpath(save_path)
    plots(mdl.log.t, mdl.log.saveval; mode, save_path, kwargs...)
end

#for Systems without outputs
plots(::AbstractVector{<:Real}, ::AbstractVector{Nothing}; kwargs...) = nothing

#for AbstractAirframeNode and AirframeGroup Systems
function plots(t, data::AbstractVector{<:NamedTuple}; mode, save_path, kwargs...)

    c = data |> StructArray |> StructArrays.components
    for (c_label, c_data) in zip(keys(c), values(c))
        save_path_c = mkpath(joinpath(save_path, String(c_label)))
        plots(t, c_data; mode, save_path = save_path_c, kwargs...)
    end

end

function save_plots(d::Dict{String,Plots.Plot}; save_path, format = :png)
    for (id, p) in zip(keys(d), values(d))
        savefig(p, joinpath(save_path, id*"."*String(format)))
    end
end

############################ TimeHistory #################################

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

end