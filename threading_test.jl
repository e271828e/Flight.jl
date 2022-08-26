function f_busywait(t)
    time0 = time()
    while time() - time0 < t end
end

function f_background()
    t = @elapsed begin
        println("Starting background task on thread $(Threads.threadid())")
    end
    println("Background task done in $t seconds")
end

function f_nested()
    Threads.@spawn f_background()
    f_busywait(5)
end

function f_parallel()
    @sync begin
        Threads.@spawn f_busywait(5)
        Threads.@spawn f_background()
    end
end

#conclusion: if a background task is started from the main task, but the main
#task doesn't ever block, the background task never gets any CPU time