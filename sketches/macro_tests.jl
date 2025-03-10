macro echo(input)
    println(input.args)
    expr = quote
       $input
    end
    return expr

end

macro no_cont(input)
    esc(:(Systems.f_ode!(::System{<:($input)}, args...) = nothing))
end

macro no_disc(input)
    :(Systems.f_disc!(::NoScheduling, ::System{<:($input)}, args...) = nothing)
end

macro no_step(input)
    :(Systems.f_step!(::System{<:($input)}, args...) = nothing)
end

macro no_dynamics(input)
    esc(quote
       @no_cont $input
       @no_step $input
       @no_disc $input
    end)
end

macro ss_cont(input)
    quote
        @inline function Systems.f_ode!(sys::System{<:($input)}, args...)
            for ss in sys.subsystems
                f_ode!(ss, args...)
            end
            update_y!(sys)
            return nothing
        end
    end
end

macro ss_disc(input)
    quote
        @inline function Systems.f_disc!(::NoScheduling, sys::System{<:($input)}, args...)
            for ss in sys.subsystems
                f_disc!(ss, args...)
            end
            update_y!(sys)
            return nothing
        end
    end
end

macro ss_step(input)
    quote
        @inline function Systems.f_step!(sys::System{<:($input)}, args...)
            for ss in sys.subsystems
                f_step!(ss, args...)
            end
            update_y!(sys)
            return nothing
        end
    end
end