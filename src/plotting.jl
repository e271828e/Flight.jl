module Plotting

using StructArrays, RecursiveArrayTools
using LaTeXStrings
using Plots
using Plots.PlotMeasures

export TimeHistory
export thplot

mutable struct TimeHistory{D}
    t::AbstractVector{<:Real}
    data::D
end


@recipe function f(v::TimeHistory{<:AbstractVector{<:Real}})

    xguide --> L"$t \: (s)$"

    return v.t, v.data

end

@recipe function f(v::TimeHistory{<:AbstractVector{<:AbstractVector{<:Real}}}; split = :none)

    t = v.t
    data = convert(Array, VectorOfArray(v.data))'
    vlength = size(data)[2]

    xguide --> L"$t \: (s)$"

    label --> (vlength <= 3 ?  ["x" "y" "z"][:, 1:vlength] : (1:vlength)')

    if split === :h
        layout --> (1, vlength)
        link --> :y #alternative: :none
    elseif split === :v
        layout --> (vlength, 1)
    else
        layout --> 1
    end
    delete!(plotattributes, :split)

    return t, data

end

thplot(t, data; kwargs...) = plot(TimeHistory(t, data); kwargs...)

end