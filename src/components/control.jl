module Control

using ComponentArrays, StaticArrays, UnPack, LinearAlgebra
using ControlSystems: ControlSystemsBase, ControlSystems, ss

using Flight.FlightCore.Systems
using Flight.FlightCore.Plotting
using Flight.FlightCore.GUI

export LinearStateSpace, PIContinuous, PIDDiscrete

################################################################################
########################### LinearStateSpace ###################################

const tV = AbstractVector{<:Float64}
const tM = AbstractMatrix{<:Float64}

struct LinearStateSpace{ LX, LU, LY, #state, input and output vector lengths
                        tX <: tV, tU <: tV, tY <: tV,
                        tA <: tM, tB <: tM, tC <: tM, tD <: tM} <: SystemDefinition

    ẋ0::tX; x0::tX; u0::tU; y0::tY; #reference values (for linearized systems)
    A::tA; B::tB; C::tC; D::tD; #state-space matrices
    x_cache::tX; y_cache::tY; y_cache_out::tY;
    Δx_cache::tX; Δu_cache::tU

    function LinearStateSpace(ẋ0, x0, u0, y0, A, B, C, D)

        lengths = map(length, (x0, u0, y0))
        types = map(typeof, (x0, u0, y0, A, B, C, D))

        vectors = map(copy, (ẋ0, x0, u0, y0))
        matrices = map(copy, (A, B, C, D))
        caches = map(copy, (x0, y0, y0, x0, u0))

        new{lengths..., types...}(vectors..., matrices..., caches...)

    end

end

LinearStateSpace(; ẋ0, x0, u0, y0, A, B, C, D) = LinearStateSpace(ẋ0, x0, u0, y0, A, B, C, D)

ControlSystems.ss(cmp::LinearStateSpace) = ControlSystems.ss(cmp.A, cmp.B, cmp.C, cmp.D)

function LinearStateSpace(sys::ControlSystemsBase.StateSpace{ControlSystemsBase.Continuous, <:AbstractFloat})
    @unpack A, B, C, D, nx, nu, ny = sys
    ẋ0 = zeros(nx); x0 = zeros(nx); u0 = zeros(nu); y0 = zeros(ny)
    LinearStateSpace(; ẋ0, x0, u0, y0, A, B, C, D)
end

Systems.init(::SystemX, cmp::LinearStateSpace) = copy(cmp.x0)
Systems.init(::SystemU, cmp::LinearStateSpace) = copy(cmp.u0)
Systems.init(::SystemY, cmp::LinearStateSpace) = SVector{length(cmp.y0)}(cmp.y0)

function Systems.f_ode!(sys::System{<:LinearStateSpace{LX, LU, LY}}) where {LX, LU, LY}

    @unpack ẋ, x, u, y, params = sys
    @unpack ẋ0, x0, u0, y0, A, B, C, D, x_cache, y_cache, y_cache_out, Δx_cache, Δu_cache = params

    #the following is equivalent to:
    #ẋ = ẋ0 + A * (x - x0) + B * (u - u0)
    #y = y0 + C * (x - x0) + D * (u - u0)
    #... but without allocations

    @. Δx_cache = x - x0
    @. Δu_cache = u - u0

    ẋ .= ẋ0
    mul!(x_cache, A, Δx_cache)
    ẋ .+= x_cache
    mul!(x_cache, B, Δu_cache)
    ẋ .+= x_cache

    y_cache_out .= y0
    mul!(y_cache, C, Δx_cache)
    y_cache_out .+= y_cache
    mul!(y_cache, D, Δu_cache)
    y_cache_out .+= y_cache

    sys.y = SVector{LY}(y_cache_out)

    return nothing

end

function Base.filter(cmp::LinearStateSpace; x = keys(cmp.x0), u = keys(cmp.u0), y = keys(cmp.y0))

    #to do: make it work for scalars

    x_ind = x
    u_ind = u
    y_ind = y

    ẋ0 = cmp.ẋ0[x_ind]
    x0 = cmp.x0[x_ind]
    u0 = cmp.u0[u_ind]
    y0 = cmp.y0[y_ind]
    A = cmp.A[x_ind, x_ind]
    B = cmp.B[x_ind, u_ind]
    C = cmp.C[y_ind, x_ind]
    D = cmp.D[y_ind, u_ind]

    return LinearStateSpace(; ẋ0, x0, u0, y0, A, B, C, D)

end

# #convert a mutable ComponentVector into a (non-allocating) immutable
# #ComponentVector with an underlying StaticVector, preserving the original
# #ComponentVector's axes
# function static_cv(L::Integer, x::ComponentVector)
#     ComponentVector(SVector{L}(getdata(x)), getaxes(x))
# end


####################### Proportional-Integral Compensator ######################
################################################################################

struct PIContinuous{N} <: SystemDefinition #Parallel form
    k_p::SVector{N,Float64} #proportional gain
    k_i::SVector{N,Float64} #integral gain
    k_l::SVector{N,Float64} #integrator leak factor
    β_p::SVector{N,Float64} #proportional path input weighting factor
end

function PIContinuous{N}(; k_p::Real = 1.0, k_i::Real = 0.1, k_l::Real = 0.0,
                           β_p::Real = 1.0) where {N}
    s2v = (x)->fill(x,N)
    PIContinuous{N}(s2v(k_p), s2v(k_i), s2v(k_l), s2v(β_p))
end

Base.@kwdef struct PIContinuousU{N}
    setpoint::MVector{N,Float64} = zeros(N) #commanded setpoint
    feedback::MVector{N,Float64} = zeros(N) #plant feedback (non-inverted)
    bound_lo::MVector{N,Float64} = fill(-Inf, N) #lower output bounds
    bound_hi::MVector{N,Float64} = fill(Inf, N) #higher output bounds
    sat_ext::MVector{N,Int64} = zeros(Int64, N) #external (signed) saturation signal
    anti_windup::MVector{N,Bool} = ones(Bool, N) #enable anti wind-up
    reset::MVector{N,Bool} = zeros(Bool, N) #reset PID states and null outputs
end

Base.@kwdef struct PIContinuousY{N}
    setpoint::SVector{N,Float64} = zeros(SVector{N}) #commanded setpoint
    feedback::SVector{N,Float64} = zeros(SVector{N}) #plant feedback (non-inverted)
    bound_lo::SVector{N,Float64} = fill(-Inf, N) #lower output bounds
    bound_hi::SVector{N,Float64} = fill(Inf, N) #higher output bounds
    sat_ext::SVector{N,Int64} = zeros(Int64, N) #external (signed) saturation signal
    anti_windup::SVector{N,Bool} = zeros(SVector{N, Bool}) #anti wind-up enabled
    reset::SVector{N,Bool} = zeros(SVector{N, Bool}) #reset PID states and null outputs
    u_p::SVector{N,Float64} = zeros(SVector{N}) #proportional path input
    y_p::SVector{N,Float64} = zeros(SVector{N}) #proportional path output
    u_i::SVector{N,Float64} = zeros(SVector{N}) #integral path input
    y_i::SVector{N,Float64} = zeros(SVector{N}) #integral path output
    int_halt::SVector{N,Bool} = zeros(SVector{N, Bool}) #integration halted
    out_free::SVector{N,Float64} = zeros(SVector{N}) #total output, free
    out::SVector{N,Float64} = zeros(SVector{N}) #actual output
    sat_out::SVector{N,Int64} = zeros(SVector{N, Int64}) #current output saturation status
end

Systems.init(::SystemX, ::PIContinuous{N}) where {N} = zeros(N)
Systems.init(::SystemY, ::PIContinuous{N}) where {N} = PIContinuousY{N}()
Systems.init(::SystemU, ::PIContinuous{N}) where {N} = PIContinuousU{N}()

function Systems.f_ode!(sys::System{<:PIContinuous{N}}) where {N}

    @unpack ẋ, x, u, params = sys
    @unpack k_p, k_i, k_l, β_p = params

    x_i = SVector{N, Float64}(x)
    setpoint = SVector(u.setpoint)
    feedback = SVector(u.feedback)
    bound_lo = SVector(u.bound_lo)
    bound_hi = SVector(u.bound_hi)
    anti_windup = SVector(u.anti_windup)
    sat_ext = SVector(u.sat_ext)
    reset = SVector(u.reset)

    u_p = β_p .* setpoint - feedback
    u_i = setpoint - feedback

    y_p = k_p .* u_p .* .!reset
    y_i = k_i .* x_i .* .!reset
    out_free = y_p + y_i #raw output
    out = clamp.(out_free, bound_lo, bound_hi) #clamped output

    sat_hi = out_free .>= bound_hi
    sat_lo = out_free .<= bound_lo
    sat_out = sat_hi - sat_lo
    int_halt = ((sign.(u_i .* sat_out) .> 0) .|| (sign.(u_i .* sat_ext) .> 0)) .&& anti_windup

    ẋ .= (u_i .* .!int_halt - k_l .* x_i) .* .!reset

    sys.y = PIContinuousY(; setpoint, feedback, bound_lo, bound_hi,
                            sat_ext, anti_windup, reset,
                            u_p, y_p, u_i, y_i, int_halt,
                            out_free, out, sat_out)

end

function Systems.f_step!(sys::System{<:PIContinuous{N}}) where {N}

    x = SVector{N,Float64}(sys.x)
    x_new = x .* .!sys.u.reset
    x_mod = any(x .!= x_new)

    sys.x .= x_new
    return x_mod

end


# ############################## Plotting ########################################

function Plotting.make_plots(th::TimeHistory{<:PIContinuousY}; kwargs...)

    pd = OrderedDict{Symbol, Plots.Plot}()

    setpoint = plot(th.setpoint; title = "Setpoint", ylabel = L"$r$", kwargs...)
    feedback = plot(th.feedback; title = "Feedback", ylabel = L"$f$", kwargs...)

    pd[:sf] = plot(setpoint, feedback;
        plot_title = "Setpoint & Feedback",
        layout = (1,2),
        link = :y,
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    sat_ext = plot(th.sat_ext; title = "External Saturation", ylabel = "", kwargs...)
    anti_windup = plot(th.anti_windup; title = "Anti-Windup Enable", ylabel = "", kwargs...)
    reset = plot(th.reset; title = "Reset", ylabel = "", kwargs...)

    pd[:ctl] = plot(anti_windup, sat_ext, reset;
        plot_title = "External Control Signals",
        layout = (1,3),
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    u_p = plot(th.u_p; title = "Input", ylabel = L"$u_p$", kwargs...)
    y_p = plot(th.y_p; title = "Output", ylabel = L"$y_p$", kwargs...)

    pd[:prop] = plot(u_p, y_p;
        plot_title = "Proportional Path",
        layout = (1,2),
        link = :y,
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    u_i = plot(th.u_i; title = "Input", ylabel = L"$u_i$", kwargs...)
    y_i = plot(th.y_i; title = "Output", ylabel = L"$y_i$", kwargs...)
    int_halt = plot(th.int_halt; title = "Integrator Halted", ylabel = "", kwargs...)

    pd[:int] = plot(u_i, y_i, int_halt;
        plot_title = "Integral Path",
        layout = (1,3),
        link = :y,
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    out_free = plot(th.out_free; title = "Free", ylabel = L"$y_{free}$", kwargs...)
    out = plot(th.out; title = "Actual", ylabel = L"$y$", kwargs...)
    sat_out = plot(th.sat_out; title = "Saturation", ylabel = L"$S$", kwargs...)

    pd[:out] = plot(out_free, out, sat_out;
        plot_title = "PID Output",
        layout = (1,3),
        link = :y,
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    return pd

end

#################################### GUI #######################################


function GUI.draw(sys::System{<:PIContinuous{N}}, label::String = "PIContinuous{$N}") where {N}

    @unpack y, params = sys

    CImGui.Begin(label)

    if CImGui.TreeNode("Params")
        @unpack k_p, k_i, k_l, β_p = params
        CImGui.Text("Proportional Gain: $k_p")
        CImGui.Text("Integral Gain: $k_i")
        CImGui.Text("Leak Factor: = $k_l")
        CImGui.Text("Proportional Path Weighting Factor: = $β_p")
        CImGui.TreePop()
    end


    if CImGui.TreeNode("Outputs")
        @unpack setpoint, feedback, bound_lo, bound_hi, sat_ext, anti_windup, reset,
                u_p, y_p, u_i, y_i, int_halt, out_free, out, sat_out = y

        CImGui.Text("Setpoint = $setpoint")
        CImGui.Text("Feedback = $feedback")
        CImGui.Text("Lower Bound = $bound_lo")
        CImGui.Text("Upper Bound = $bound_hi")
        CImGui.Text("External Saturation = $sat_ext")
        CImGui.Text("Anti Wind-up = $anti_windup")
        CImGui.Text("Reset = $reset")
        CImGui.Text("Proportional Path Input = $u_p")
        CImGui.Text("Proportional Path Output = $y_p")
        CImGui.Text("Integral Path Input = $u_i")
        CImGui.Text("Integral Path Output = $y_i")
        CImGui.Text("Integrator Halted = $int_halt")
        CImGui.Text("Free Output = $out_free")
        CImGui.Text("Actual Output = $out")
        CImGui.Text("Output Saturation = $sat_out")
        CImGui.TreePop()
    end

    CImGui.End()

end #function

function GUI.draw!(sys::System{<:PIContinuous{N}}, label::String = "PIContinuous{$N}") where {N}

    GUI.draw(sys, label)

    CImGui.Begin(label)

    #this shows how to get input from widgets without the @cstatic and @c macros
    if CImGui.TreeNode("Inputs")
        for i in 1:N
            if CImGui.TreeNode("[$i]")
                @unpack reset = sys.u
                reset_ref = Ref(reset[i])
                CImGui.Checkbox("Reset", reset_ref) #; CImGui.SameLine()
                reset[i] = reset_ref[]
                CImGui.TreePop()
            end
        end
        CImGui.TreePop()
    end

    CImGui.End()

end

################# Proportional-Integral-Derivative Compensator #################
################################################################################


struct PIDDiscrete{N} <: SystemDefinition #Parallel form
    k_p::SVector{N,Float64} #proportional gain
    k_i::SVector{N,Float64} #integral gain
    k_d::SVector{N,Float64} #derivative gain
    τ_d::SVector{N,Float64} #derivative filter time constant
    β_p::SVector{N,Float64} #proportional path setpoint weighting factor
    β_d::SVector{N,Float64} #derivative path setpoint weighting factor
end

function PIDDiscrete{N}(; k_p::Real = 1.0, k_i::Real = 0.1,
                        k_d::Real = 0, τ_d::Real = 0.05,
                        β_p::Real = 1.0, β_d::Real = 1.0) where {N}
    s2v = (x)->fill(x,N)
    PIDDiscrete{N}( s2v(k_p), s2v(k_i),
                   s2v(k_d), s2v(τ_d),
                   s2v(β_p), s2v(β_d),)
end

Base.@kwdef struct PIDDiscreteU{N}
    setpoint::MVector{N,Float64} = zeros(N) #commanded setpoint
    feedback::MVector{N,Float64} = zeros(N) #plant feedback (non-inverted)
    bound_lo::MVector{N,Float64} = fill(-Inf, N) #lower output bounds
    bound_hi::MVector{N,Float64} = fill(Inf, N) #higher output bounds
    sat_ext::MVector{N,Int64} = zeros(Int64, N) #external (signed) saturation signal
    anti_windup::MVector{N,Bool} = ones(Bool, N) #enable anti wind-up
    reset::MVector{N,Bool} = zeros(Bool, N) #reset PID states and null outputs
end

Base.@kwdef struct PIDDiscreteS{N}
    x_i0::MVector{N,Float64} = zeros(N) #previous integrator path state
    x_d0::MVector{N,Float64} = zeros(N) #previous derivative path state
    sat_out_0::MVector{N,Int64} = zeros(N) #previous output saturation status
end

Base.@kwdef struct PIDDiscreteY{N}
    setpoint::SVector{N,Float64} = zeros(SVector{N}) #commanded setpoint
    feedback::SVector{N,Float64} = zeros(SVector{N}) #plant feedback (non-inverted)
    bound_lo::SVector{N,Float64} = fill(-Inf, N) #lower output bounds
    bound_hi::SVector{N,Float64} = fill(Inf, N) #higher output bounds
    sat_ext::SVector{N,Int64} = zeros(Int64, N) #external (signed) saturation signal
    anti_windup::SVector{N,Bool} = zeros(SVector{N, Bool}) #anti wind-up enabled
    reset::SVector{N,Bool} = zeros(SVector{N, Bool}) #reset PID states and null outputs
    u_p::SVector{N,Float64} = zeros(SVector{N}) #proportional path input
    u_i::SVector{N,Float64} = zeros(SVector{N}) #integral path input
    u_d::SVector{N,Float64} = zeros(SVector{N}) #derivative path input
    int_halt::SVector{N,Bool} = zeros(SVector{N, Bool}) #integration halted
    y_p::SVector{N,Float64} = zeros(SVector{N}) #proportional path output
    y_i::SVector{N,Float64} = zeros(SVector{N}) #integral path output
    y_d::SVector{N,Float64} = zeros(SVector{N}) #derivative path output
    out_free::SVector{N,Float64} = zeros(SVector{N}) #non-clamped PID output
    out::SVector{N,Float64} = zeros(SVector{N}) #actual PID output
    sat_out::SVector{N,Int64} = zeros(SVector{N, Int64}) #current output saturation status
end

Systems.init(::SystemY, ::PIDDiscrete{N}) where {N} = PIDDiscreteY{N}()
Systems.init(::SystemU, ::PIDDiscrete{N}) where {N} = PIDDiscreteU{N}()
Systems.init(::SystemS, ::PIDDiscrete{N}) where {N} = PIDDiscreteS{N}()

function Systems.f_disc!(sys::System{<:PIDDiscrete{N}}, Δt::Real) where {N}

    @unpack s, u, params = sys
    @unpack k_p, k_i, k_d, τ_d, β_p, β_d = params

    setpoint = SVector(u.setpoint)
    feedback = SVector(u.feedback)
    bound_lo = SVector(u.bound_lo)
    bound_hi = SVector(u.bound_hi)
    anti_windup = SVector(u.anti_windup)
    sat_ext = SVector(u.sat_ext)
    reset = SVector(u.reset)

    x_i0 = SVector(s.x_i0)
    x_d0 = SVector(s.x_d0)
    sat_out_0 = SVector(s.sat_out_0)

    u_p = β_p .* setpoint - feedback
    u_d = β_d .* setpoint - feedback
    u_i = setpoint - feedback

    #integration is halted in those components for which:
    #((output saturation has the same sign as integrator path input) OR
    #(external saturation signal has the same sign as integrator path input))
    #AND (anti_windup is enabled)
    int_halt = ((sign.(u_i .* sat_out_0) .> 0) .|| (sign.(u_i .* sat_ext) .> 0)) .&& anti_windup
    α = 1 ./ (τ_d .+ Δt)

    x_i = x_i0 + Δt * k_i .* u_i .* .!int_halt
    x_d = α .* τ_d .* x_d0 + Δt * α .* k_d .* u_d

    y_p = k_p .* u_p
    y_i = x_i
    y_d = α .* (-x_d0 + k_d .* u_d)

    x_i = x_i .* .!reset
    x_d = x_d .* .!reset

    y_p = y_p .* .!reset
    y_i = y_i .* .!reset
    y_d = y_d .* .!reset

    out_free = y_p + y_i + y_d
    out = clamp.(out_free, bound_lo, bound_hi)

    sat_hi = out_free .>= bound_hi
    sat_lo = out_free .<= bound_lo
    sat_out = sat_hi - sat_lo

    s.x_i0 .= x_i
    s.x_d0 .= x_d
    s.sat_out_0 .= sat_out

    sys.y = PIDDiscreteY(; setpoint, feedback, bound_lo, bound_hi,
                           sat_ext, anti_windup, reset,
                           u_p, u_i, u_d, int_halt, y_p, y_i, y_d,
                           out_free, out, sat_out)

    return false

end

function Plotting.make_plots(th::TimeHistory{<:PIDDiscreteY}; kwargs...)

    pd = OrderedDict{Symbol, Plots.Plot}()

    setpoint = plot(th.setpoint; title = "Setpoint", ylabel = L"$r$", kwargs...)
    feedback = plot(th.feedback; title = "Feedback", ylabel = L"$f$", kwargs...)

    pd[:sf] = plot(setpoint, feedback;
        plot_title = "Setpoint & Feedback",
        layout = (1,2),
        link = :y,
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    sat_ext = plot(th.sat_ext; title = "External Saturation", ylabel = "", kwargs...)
    anti_windup = plot(th.anti_windup; title = "Anti-Windup Enable", ylabel = "", kwargs...)
    reset = plot(th.reset; title = "Reset", ylabel = "", kwargs...)

    pd[:ctl] = plot(anti_windup, sat_ext, reset;
        plot_title = "External Control Signals",
        layout = (1,3),
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    u_p = plot(th.u_p; title = "Input", ylabel = L"$u_p$", kwargs...)
    y_p = plot(th.y_p; title = "Output", ylabel = L"$y_p$", kwargs...)

    pd[:prop] = plot(u_p, y_p;
        plot_title = "Proportional Path",
        layout = (1,2),
        link = :y,
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    u_i = plot(th.u_i; title = "Input", ylabel = L"$u_i$", kwargs...)
    y_i = plot(th.y_i; title = "Output", ylabel = L"$y_i$", kwargs...)
    int_halt = plot(th.int_halt; title = "Integrator Halted", ylabel = "", kwargs...)

    pd[:int] = plot(u_i, y_i, int_halt;
        plot_title = "Integral Path",
        layout = (1,3),
        link = :y,
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    u_d = plot(th.u_d; title = "Input", ylabel = L"$u_d$", kwargs...)
    y_d = plot(th.y_d; title = "Output", ylabel = L"$y_d$", kwargs...)

    pd[:der] = plot(u_d, y_d;
        plot_title = "Derivative Path",
        layout = (1,2),
        link = :y,
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    out_free = plot(th.out_free; title = "Free", ylabel = L"$y_{free}$", kwargs...)
    out = plot(th.out; title = "Actual", ylabel = L"$y$", kwargs...)
    sat_out = plot(th.sat_out; title = "Saturation", ylabel = L"$S$", kwargs...)

    pd[:out] = plot(out_free, out, sat_out;
        plot_title = "PID Output",
        layout = (1,3),
        link = :y,
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    return pd

end


end #module