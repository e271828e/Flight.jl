function put_loop(ch::Channel)
    a = π/4
    println("Put loop running on thread $(Threads.threadid())")
    while true
        try
            put!(ch, a)
        catch ex
            break
        end
    end
end

function take_loop(ch::Channel)
    println("Take loop running on thread $(Threads.threadid())")
    for _ in 1:1_000_000
        b = take!(ch)
    end
end

function test_multi()

    channel = Channel{Float64}()

    Threads.@spawn put_loop(channel)
    t = @elapsed wait(Threads.@spawn take_loop(channel))
    close(channel)

    println("Approximate context switching time (μs): $(t/1e6)")

end

function test_single()

    channel = Channel{Float64}()

    #single thread
    @async put_loop(channel)
    t = @elapsed wait(@async take_loop(channel))
    close(channel)

    println("Approximate context switching time (μs): $(t/1e6)")

end


function test_lock_overhead(l::ReentrantLock)
    a = π/4
    # println("Put loop running on thread $(Threads.threadid())")
    for _ in 1:10
        lock(l)
        # try
            a += a
        # finally
        unlock(l)
        # end
    end
    return a
end