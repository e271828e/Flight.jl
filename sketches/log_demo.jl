using ComponentArrays, RecursiveArrayTools

function plotdemo()

    x = ComponentVector(a = ones(2), b = ComponentVector(b1 = fill(2, 5), b2 = zeros(3)))
    v = Vector{typeof(x)}()

    #make a Vector of ComponentVectors to simulate a simulation log
    for i in 1:10
        push!(v, copy(x))
    end

    #this returns a 10-element Vector of ComponentVector which can be indexed as
    #a matrix, while preserving the axis structure of its ComponentVector
    #elements
    v = VectorOfArray(v)

    #this produces an actual ComponentMatrix of size [length(x), 10]. the
    #awesome part is that its row dimension preserves the Axes metadata in the
    #original x. its column dimension has a Flat Axis.
    m = convert(Array, v)

    #we can now easily extract and plot the log of any specific block in the
    #original x
    log_b = m[:b, :]
    @show log_b2 = log_b[:b2, :]

end