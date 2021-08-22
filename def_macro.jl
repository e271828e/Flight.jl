macro hello(L)
    ex = Expr(:tuple) #this is the expression head, all the pushes afterwards add an argument to the tuple
    println(L)
    for label in L
        label = QuoteNode(label)
        ex_ss = quote
            f_cont!(sys.subsystems[$label], args...)
        end
        push!(ex.args, ex_ss)
    end
    # push!(ex.args, :(NamedTuple{$L}()))
    return ex
end