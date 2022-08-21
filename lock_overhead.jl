function put_loop(l::ReentrantLock, r::Ref{Float64})
    a = π/4
    println("Put loop running on thread $(Threads.threadid())")
    while true
        lock(l)
        # try
            r[] += a
        # finally
        unlock(l)
        # end
    end
end

function take_loop(ch::Channel)
    println("Take loop running on thread $(Threads.threadid())")
    for _ in 1:1_000_000
        b = take!(ch)
    end
end

function test_multi()

    mylock = ReentrantLock()
    r = Ref(0.0)

    #single thread
    @async put_loop(mylock, r)
    t = @elapsed wait(@async take_loop(mylock, r))
    close(channel)

    println("Approximate context switching time (μs): $(t/1e6)")

end

function test_single()

    mylock = ReentrantLock()
    r = Ref(0.0)

    #single thread
    @async put_loop(mylock, r)
    t = @elapsed wait(@async take_loop(mylock, r))
    close(channel)

    println("Approximate context switching time (μs): $(t/1e6)")

end