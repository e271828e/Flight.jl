using StaticArrays, BlockArrays, LabelledArrays, BenchmarkTools
includet("statevectorfast.jl")
using .StateVectorFast


function test_fast()
    println("Testing StateVectorFast")

    @getlength
end