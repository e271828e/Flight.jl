using ComponentArrays
using Flight
using UnPack
using OrdinaryDiffEq

Base.@kwdef struct MassSpringDamper <: SystemDescriptor
    ω_n::Float64 = 1.0 #undamped natural frequency
    ζ::Float64 = 1.0 #damping ratio
    Δt_d::Float64 = 2 #spring switching time interval
    v_ϵ::Float64 = 1e-6 #velocity threshold
end

Base.@kwdef struct MassSpringDamperY
    p::Float64 = 0.0
    v::Float64 = 0.0
    a::Float64 = 0.0
end

Base.@kwdef mutable struct MassSpringDamperS
    damper_engaged::Bool = false
    t_last_moving::Float64 = 0.0
    t_last_stopped::Float64 = 0.0
end

Systems.init(::SystemX, ::MassSpringDamper) = ComponentVector(p = 0.0, v = 0.0)
Systems.init(::SystemU, ::MassSpringDamper) = Ref(0.0)
Systems.init(::SystemY, ::MassSpringDamper) = MassSpringDamperY()
Systems.init(::SystemS, ::MassSpringDamper) = MassSpringDamperS()

function Systems.f_ode!(sys::System{<:MassSpringDamper})

    @unpack x, u, s, params = sys
    @unpack p, v = x
    @unpack ω_n, ζ = params

    f = u[] #de-reference u to get the underlying Float64
    ζ *= s.damper_engaged
    a = f - ω_n^2 * p - 2ζ * ω_n * v

    #update sys.ẋ
    sys.ẋ.v = a
    sys.ẋ.p = v

    #update sys.y (cannot be mutated, we need to assign a new instance to it)
    sys.y = MassSpringDamperY(; p, v, a)
end

function Systems.f_step!(sys::System{<:MassSpringDamper})

    @unpack x, s, t, params = sys
    @unpack Δt_d, v_ϵ = params

    #t also needs to be de-referenced
    abs(x.v) > v_ϵ ? s.t_last_moving = t[] : s.t_last_stopped = t[]

    Δt_stopped = t[] - s.t_last_moving
    Δt_moving = t[] - s.t_last_stopped

    if s.damper_engaged
        Δt_stopped > Δt_d ? s.damper_engaged = false : nothing
    else #!s.damper_engaged
        Δt_moving > Δt_d ? s.damper_engaged = true : nothing
    end

    return false
end



function get_sim()

    sys_desc = MassSpringDamper(; ω_n = 20.0, ζ = 0.4, Δt_d = 1)
    sys = System(sys_desc)
    # @show typeof(sys)
    # Utils.showfields(sys) #show the type's fields with their current values

    function sys_init!(sys; p::Real = 0.0, v::Real = 0.0, f::Real = 0.0, )
        sys.x.p = p
        sys.x.v = v
        sys.u[] = f
        sys.s.damper_engaged = false
        sys.s.t_last_moving = 0
        sys.s.t_last_stopped = 0
        @show sys.x
        @show sys.u
    end

    function sys_io!(u, t, y, params)
        # println("System IO called")
        u[] += 0
    end

    sim = Simulation(sys; algorithm = Heun(), sys_init!, sys_io!,
                        init_kwargs = (f = 0.0, p = 0.0, v = 0.0))

    return sim



end