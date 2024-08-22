function threading_sketch()

    c = Channel{Int}(1)
    @sync begin
        Threads.@spawn begin
            while isopen(c)
                Core.println("Taken $(take!(c))")
            end
        end
        Threads.@spawn begin
            for i in 1:5
                sleep(1)
                @lock c begin
                    if !isready(c)
                        Core.println("Putting $i")
                        put!(c, i)
                    end
                end
            end
            close(c)
            Core.println("Bye")
        end

    end
end