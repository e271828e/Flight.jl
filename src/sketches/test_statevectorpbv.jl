using StaticArrays, BlockArrays, LabelledArrays, BenchmarkTools
# includet("mylabelledarrays.jl")
includet("statevectorpbv.jl")
# using .MyLabelledArrays
using .StateVectorPBV

function test_pbv()
    println("Hello")


    # @btime (y .= x .+ 2 .* x) setup=(x=rand(10); y = zeros(10))
    # @btime (y .= x .+ 2 .* x) setup=(x=PseudoBlockArray(zeros(10), [4,3,3]); y=similar(x))
    # @btime (x .+ 2 .* x) setup=(x=PseudoBlockArray(zeros(10), [4,3,3]); y=similar(x))

    x_plain = rand(10)
    x_rbd = Node{:rbd}()
    y_rbd = similar(Node{:rbd})
    # @benchmark (y_rbd .= x_rbd .+ 2 .* x_rbd) setup=(x_rbd=Node{:rbd}(); y_rbd=similar(x))
    # @benchmark ($y_rbd .= $x_rbd .+ $x_rbd)
    # @btime ($y_rbd .= sin.($x_rbd) .+ cos.($x_rbd))
    # @btime (y = $y_rbd[:]; x = $x_rbd; y.= sin.(x) .+ cos.(x))
    # @btime ($y_rbd .= sin.($x_rbd) .+ cos.($x_rbd).+ exp.($x_rbd))
    # @btime (sin.(cos.($x_plain)))
    # @btime broadcast!(+, $y_rbd, $x_rbd, $x_rbd)
    @btime $x_rbd + 2.0 * $x_rbd .* $x_rbd
end