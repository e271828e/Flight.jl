module Control

using ComponentArrays, StaticArrays, UnPack, LinearAlgebra
using ControlSystems: ControlSystemsBase, ControlSystems, ss
using RobustAndOptimalControl

using Flight.FlightCore.Systems
using Flight.FlightCore.Plotting
using Flight.FlightCore.GUI

export LinearStateSpace, submodel
export PIContinuous, PIDDiscreteVector
export PIDParams, PIDDiscrete, IntegratorDiscrete, LeadLagDiscrete

export PIDOpt



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

function RobustAndOptimalControl.named_ss(lss::Control.LinearStateSpace)
    x_labels, u_labels, y_labels = map(collect ∘ propertynames, (lss.x0, lss.u0, lss.y0))
    named_ss(ss(lss), x = x_labels, u = u_labels, y = y_labels)
end

function LinearStateSpace(sys::ControlSystemsBase.StateSpace{ControlSystemsBase.Continuous, <:AbstractFloat})
    @unpack A, B, C, D, nx, nu, ny = sys
    ẋ0 = zeros(nx); x0 = zeros(nx); u0 = zeros(nu); y0 = zeros(ny)
    LinearStateSpace(; ẋ0, x0, u0, y0, A, B, C, D)
end

Systems.init(::SystemX, cmp::LinearStateSpace) = copy(cmp.x0)
Systems.init(::SystemU, cmp::LinearStateSpace) = copy(cmp.u0)
Systems.init(::SystemY, cmp::LinearStateSpace) = SVector{length(cmp.y0)}(cmp.y0)

function Systems.f_ode!(sys::System{<:LinearStateSpace{LX, LU, LY}}) where {LX, LU, LY}

    @unpack ẋ, x, u, y, constants = sys
    @unpack ẋ0, x0, u0, y0, A, B, C, D, x_cache, y_cache, y_cache_out, Δx_cache, Δu_cache = constants

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

function submodel(cmp::LinearStateSpace; x = keys(cmp.x0), u = keys(cmp.u0), y = keys(cmp.y0))

    #to do: generalize for scalars

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

@kwdef struct PIContinuousU{N}
    setpoint::MVector{N,Float64} = zeros(N) #commanded setpoint
    feedback::MVector{N,Float64} = zeros(N) #plant feedback (non-inverted)
    bound_lo::MVector{N,Float64} = fill(-Inf, N) #lower output bounds
    bound_hi::MVector{N,Float64} = fill(Inf, N) #higher output bounds
    sat_ext::MVector{N,Int64} = zeros(Int64, N) #external (signed) saturation signal
    anti_windup::MVector{N,Bool} = ones(Bool, N) #enable anti wind-up
    reset::MVector{N,Bool} = zeros(Bool, N) #reset PID states and null outputs
end

@kwdef struct PIContinuousY{N}
    setpoint::SVector{N,Float64} = zeros(SVector{N}) #commanded setpoint
    feedback::SVector{N,Float64} = zeros(SVector{N}) #plant feedback (non-inverted)
    bound_lo::SVector{N,Float64} = fill(-Inf, SVector{N}) #lower output bounds
    bound_hi::SVector{N,Float64} = fill(Inf, SVector{N}) #higher output bounds
    anti_windup::SVector{N,Bool} = zeros(SVector{N, Bool}) #anti wind-up enabled
    reset::SVector{N,Bool} = zeros(SVector{N, Bool}) #reset PID states and null outputs
    u_p::SVector{N,Float64} = zeros(SVector{N}) #proportional path input
    u_i::SVector{N,Float64} = zeros(SVector{N}) #integral path input
    y_p::SVector{N,Float64} = zeros(SVector{N}) #proportional path output
    y_i::SVector{N,Float64} = zeros(SVector{N}) #integral path output
    out_free::SVector{N,Float64} = zeros(SVector{N}) #total output, free
    out::SVector{N,Float64} = zeros(SVector{N}) #actual output
    sat_ext::SVector{N,Int64} = zeros(SVector{N, Int64}) #external (signed) saturation signal
    sat_out::SVector{N,Int64} = zeros(SVector{N, Int64}) #current output saturation status
    int_halt::SVector{N,Bool} = zeros(SVector{N, Bool}) #integration halted
end

Systems.init(::SystemX, ::PIContinuous{N}) where {N} = zeros(N)
Systems.init(::SystemY, ::PIContinuous{N}) where {N} = PIContinuousY{N}()
Systems.init(::SystemU, ::PIContinuous{N}) where {N} = PIContinuousU{N}()

function Systems.f_ode!(sys::System{<:PIContinuous{N}}) where {N}

    @unpack ẋ, x, u, constants = sys
    @unpack k_p, k_i, k_l, β_p = constants

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

    @unpack setpoint, feedback, bound_lo, bound_hi, sat_ext, anti_windup, reset,
            u_p, y_p, u_i, y_i, int_halt, out_free, out, sat_out = sys.y

    # CImGui.Begin(label)

        CImGui.Text("Setpoint = $setpoint")
        CImGui.Text("Feedback = $feedback")
        CImGui.Text("Lower Bound = $bound_lo")
        CImGui.Text("Upper Bound = $bound_hi")
        CImGui.Text("Anti Wind-up = $anti_windup")
        CImGui.Text("Reset = $reset")
        CImGui.Text("Proportional Path Input = $u_p")
        CImGui.Text("Proportional Path Output = $y_p")
        CImGui.Text("Integral Path Input = $u_i")
        CImGui.Text("Integral Path Output = $y_i")
        CImGui.Text("Free Output = $out_free")
        CImGui.Text("Actual Output = $out")
        CImGui.Text("External Saturation = $sat_ext")
        CImGui.Text("Output Saturation = $sat_out")
        CImGui.Text("Integrator Halted = $int_halt")

    # CImGui.End()

end #function

function GUI.draw!(sys::System{<:PIContinuous{N}}, label::String = "PIContinuous{$N}") where {N}

    CImGui.Begin(label)

    if CImGui.TreeNode("Outputs")
        GUI.draw(sys, label)
        CImGui.TreePop()
    end
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


struct PIDDiscreteVector{N} <: SystemDefinition #Parallel form
    k_p::SVector{N,Float64} #proportional gain
    k_i::SVector{N,Float64} #integral gain
    k_d::SVector{N,Float64} #derivative gain
    τ_d::SVector{N,Float64} #derivative filter time constant
    β_p::SVector{N,Float64} #proportional path setpoint weighting factor
    β_d::SVector{N,Float64} #derivative path setpoint weighting factor
end

function PIDDiscreteVector{N}(; k_p::Real = 1.0, k_i::Real = 0.1,
                        k_d::Real = 0, τ_d::Real = 0.05,
                        β_p::Real = 1.0, β_d::Real = 1.0) where {N}
    s2v = (x)->fill(x,N)
    PIDDiscreteVector{N}( s2v(k_p), s2v(k_i),
                   s2v(k_d), s2v(τ_d),
                   s2v(β_p), s2v(β_d),)
end

@kwdef struct PIDDiscreteVectorU{N}
    setpoint::MVector{N,Float64} = zeros(N) #commanded setpoint
    feedback::MVector{N,Float64} = zeros(N) #plant feedback (non-inverted)
    bound_lo::MVector{N,Float64} = fill(-Inf, N) #lower output bounds
    bound_hi::MVector{N,Float64} = fill(Inf, N) #higher output bounds
    sat_ext::MVector{N,Int64} = zeros(Int64, N) #external (signed) saturation signal
    anti_windup::MVector{N,Bool} = ones(Bool, N) #enable anti wind-up
    reset::MVector{N,Bool} = zeros(Bool, N) #reset PID states and null outputs
end

@kwdef struct PIDDiscreteVectorS{N}
    x_i0::MVector{N,Float64} = zeros(N) #previous integrator path state
    x_d0::MVector{N,Float64} = zeros(N) #previous derivative path state
    sat_out_0::MVector{N,Int64} = zeros(N) #previous output saturation status
end

@kwdef struct PIDDiscreteVectorY{N}
    setpoint::SVector{N,Float64} = zeros(SVector{N}) #commanded setpoint
    feedback::SVector{N,Float64} = zeros(SVector{N}) #plant feedback (non-inverted)
    bound_lo::SVector{N,Float64} = fill(-Inf, SVector{N}) #lower output bounds
    bound_hi::SVector{N,Float64} = fill(Inf, SVector{N}) #higher output bounds
    anti_windup::SVector{N,Bool} = zeros(SVector{N, Bool}) #anti wind-up enabled
    reset::SVector{N,Bool} = zeros(SVector{N, Bool}) #reset PID states and null outputs
    u_p::SVector{N,Float64} = zeros(SVector{N}) #proportional path input
    u_i::SVector{N,Float64} = zeros(SVector{N}) #integral path input
    u_d::SVector{N,Float64} = zeros(SVector{N}) #derivative path input
    y_p::SVector{N,Float64} = zeros(SVector{N}) #proportional path output
    y_i::SVector{N,Float64} = zeros(SVector{N}) #integral path output
    y_d::SVector{N,Float64} = zeros(SVector{N}) #derivative path output
    out_free::SVector{N,Float64} = zeros(SVector{N}) #non-clamped PID output
    out::SVector{N,Float64} = zeros(SVector{N}) #actual PID output
    sat_ext::SVector{N,Int64} = zeros(SVector{N, Int64}) #external (signed) saturation signal
    sat_out::SVector{N,Int64} = zeros(SVector{N, Int64}) #output saturation status
    int_halt::SVector{N,Bool} = zeros(SVector{N, Bool}) #integration halted
end

Systems.init(::SystemY, ::PIDDiscreteVector{N}) where {N} = PIDDiscreteVectorY{N}()
Systems.init(::SystemU, ::PIDDiscreteVector{N}) where {N} = PIDDiscreteVectorU{N}()
Systems.init(::SystemS, ::PIDDiscreteVector{N}) where {N} = PIDDiscreteVectorS{N}()

function reset!(sys::System{<:PIDDiscreteVector{N}}) where {N}
    sys.u.setpoint .= 0
    sys.u.feedback .= 0
    sys.u.sat_ext .= 0
    reset_prev = SVector{N}(sys.u.reset)
    sys.u.reset .= true
    f_disc!(sys, 1.0)
    sys.u.reset .= reset_prev
end

function Systems.f_disc!(sys::System{<:PIDDiscreteVector{N}}, Δt::Real) where {N}

    @unpack s, u, constants = sys
    @unpack k_p, k_i, k_d, τ_d, β_p, β_d = constants

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

    α = 1 ./ (τ_d .+ Δt)

    u_p = β_p .* setpoint - feedback
    u_d = β_d .* setpoint - feedback
    u_i = setpoint - feedback

    #integration is halted in those components for which:
    #((output saturation has the same sign as integrator path output increment) OR
    #(external saturation signal has the same sign as integrator path output increment))
    #AND (anti_windup is enabled)
    Δx_i = Δt * k_i .* u_i
    int_halt = ((sign.(Δx_i .* sat_out_0) .> 0) .|| (sign.(Δx_i .* sat_ext) .> 0)) .&& anti_windup
    x_i = x_i0 + Δx_i .* .!int_halt
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

    sys.y = PIDDiscreteVectorY(; setpoint, feedback, bound_lo, bound_hi, anti_windup,
                           reset, u_p, u_i, u_d, y_p, y_i, y_d,
                           out_free, out, sat_ext, sat_out, int_halt)

    return false

end

function Plotting.make_plots(th::TimeHistory{<:PIDDiscreteVectorY}; kwargs...)

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

    pd[:out] = plot(out_free, output, sat_out;
        plot_title = "PID Output",
        layout = (1,3),
        link = :y,
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    return pd

end

#################################### GUI #######################################


function GUI.draw(sys::System{<:PIDDiscreteVector{N}}, label::String = "PIDDiscreteVector{$N}") where {N}

    @unpack setpoint, feedback, bound_lo, bound_hi, anti_windup, reset,
            u_p, u_i, u_d, y_p, y_i, y_d, out_free, output, sat_ext, sat_out, int_halt = sys.y

    # CImGui.Begin(label)

        CImGui.Text("Setpoint = $setpoint")
        CImGui.Text("Feedback = $feedback")
        CImGui.Text("Lower Bound = $bound_lo")
        CImGui.Text("Upper Bound = $bound_hi")
        CImGui.Text("Anti Wind-up = $anti_windup")
        CImGui.Text("Reset = $reset")
        CImGui.Text("Proportional Path Input = $u_p")
        CImGui.Text("Proportional Path Output = $y_p")
        CImGui.Text("Integral Path Input = $u_i")
        CImGui.Text("Integral Path Output = $y_i")
        CImGui.Text("Derivative Path Input = $u_d")
        CImGui.Text("Derivative Path Output = $y_d")
        CImGui.Text("Free Output = $out_free")
        CImGui.Text("Actual Output = $output")
        CImGui.Text("External Saturation = $sat_ext")
        CImGui.Text("Output Saturation = $sat_out")
        CImGui.Text("Integrator Halted = $int_halt")

    # CImGui.End()

end #function


################# IntegratorDiscrete #################
################################################################################

struct IntegratorDiscrete <: SystemDefinition end

@kwdef mutable struct IntegratorDiscreteU
    input::Float64 = 0 #current input
    sat_ext::Int64 = 0 #external (signed) saturation signal
    bound_lo::Float64 = -Inf #lower output bound
    bound_hi::Float64 = Inf #higher output bound
end

@kwdef mutable struct IntegratorDiscreteS
    x0::Float64 = 0 #previous integrator state
    sat_out_0::Int64 = 0 #previous output saturation state
end

@kwdef struct IntegratorDiscreteY
    input::Float64 = 0
    sat_ext::Float64 = 0
    bound_lo::Float64 = -Inf
    bound_hi::Float64 = Inf
    x1::Float64 = 0 #current state
    output::Float64 = 0 #current output
    sat_out::Int64 = 0 #current output saturation status
    halted::Bool = false #integration halted
end

Systems.init(::SystemY, ::IntegratorDiscrete) = IntegratorDiscreteY()
Systems.init(::SystemU, ::IntegratorDiscrete) = IntegratorDiscreteU()
Systems.init(::SystemS, ::IntegratorDiscrete) = IntegratorDiscreteS()

function reset!(sys::System{<:IntegratorDiscrete})
    sys.u.input = 0
    sys.u.sat_ext = 0
    sys.s.x0 = 0
    sys.s.sat_out_0 = 0
    f_disc!(sys, 1.0)
end

function Systems.f_disc!(sys::System{<:IntegratorDiscrete}, Δt::Real)

    @unpack s, u = sys
    @unpack input, sat_ext, bound_lo, bound_hi = u
    @unpack x0, sat_out_0 = s

    halted = ((sign(input * sat_out_0) > 0) || (sign(input * sat_ext) > 0))
    x1 = x0 + Δt * input * !halted
    output = clamp(x1, bound_lo, bound_hi)

    sat_hi = x1 >= bound_hi
    sat_lo = x1 <= bound_lo
    sat_out = sat_hi - sat_lo

    s.x0 = x1
    s.sat_out_0 = sat_out

    sys.y = IntegratorDiscreteY(; input, sat_ext, bound_lo, bound_hi, x1, output, sat_out, halted)

    return false

end


#################################### GUI #######################################


function GUI.draw(sys::System{<:IntegratorDiscrete}, label::String = "IntegratorDiscrete")

    @unpack x0, sat_out_0 = sys.s
    @unpack input, sat_ext, bound_lo, bound_hi, x1, output, sat_out, halted = sys.y

    # CImGui.Begin(label)

        CImGui.Text("input = $input")
        CImGui.Text("sat_ext = $sat_ext")
        CImGui.Text("bound_lo = $bound_lo")
        CImGui.Text("bound_hi = $bound_hi")
        CImGui.Text("x0 = $x0")
        CImGui.Text("sat_out_0 = $sat_out_0")
        CImGui.Text("x1 = $x1")
        CImGui.Text("output = $output")
        CImGui.Text("sat_out = $sat_out")
        CImGui.Text("halted = $halted")

    # CImGui.End()

end #function


######################### LeadLagDiscreteCompensator ##############################
################################################################################

#Lead / lag compensator with pole p < 0, zero z < 0, gain k. Discretized by
#Tustin transform
#|p| > |z|: lead
#|p| < |z|: lag

struct LeadLagDiscrete <: SystemDefinition end

@kwdef mutable struct LeadLagDiscreteU
    z::Float64 = -1.0 #zero location (z < 0)
    p::Float64 = -10.0 #pole location (p < 0)
    k::Float64 = 1.0 #gain
    u1::Float64 = 0.0 #input
end

@kwdef mutable struct LeadLagDiscreteS
    u0::Float64 = 0.0 #previous input
    x0::Float64 = 0.0 #previous output
end

@kwdef struct LeadLagDiscreteY
    z::Float64 = -1.0
    p::Float64 = -10.0
    k::Float64 = 1.0
    u1::Float64 = 0.0
    y1::Float64 = 0.0 #current output
end

Systems.init(::SystemY, ::LeadLagDiscrete) = LeadLagDiscreteY()
Systems.init(::SystemU, ::LeadLagDiscrete) = LeadLagDiscreteU()
Systems.init(::SystemS, ::LeadLagDiscrete) = LeadLagDiscreteS()

function reset!(sys::System{<:LeadLagDiscrete})
    sys.u.u1 = 0
    sys.s.u0 = 0
    sys.s.x0 = 0
    f_disc!(sys, 1.0) #update outputs (keeps the actual p and z)
end

function Systems.f_disc!(sys::System{<:LeadLagDiscrete}, Δt::Real)

    @unpack s, u = sys
    @unpack z, p, k, u1 = u
    @unpack u0, x0 = s

    a0 = (2+p*Δt)/(2-p*Δt)
    b1 = (2-z*Δt)/(2-p*Δt)
    b0 = (-2-z*Δt)/(2-p*Δt)

    x1 = a0 * x0 + b1 * u1 + b0 * u0
    y1 = k * x1

    sys.y = LeadLagDiscreteY(; z, p, k, u1, y1)

    s.x0 = x1
    s.u0 = u1

    return false

end


#################################### GUI #######################################

function GUI.draw(sys::System{<:LeadLagDiscrete}, label::String = "Discrete Lead Compensator")

    @unpack z, p, k, u1, y1 = sys.y

    # CImGui.Begin(label)

        CImGui.Text("z = $z")
        CImGui.Text("p = $p")
        CImGui.Text("k = $k")
        CImGui.Text("u1 = $u1")
        CImGui.Text("y1 = $y1")

    # CImGui.End()

end #function


####################### Gain-Schedulable PID Compensator #######################
################################################################################

@kwdef struct PIDParams{T} <: FieldVector{4, T} #parallel form
    k_p::T = 1.0
    k_i::T = 0.0
    k_d::T = 0.1
    τ_f::T = 0.01
end

function Base.NamedTuple(params::PIDParams)
    names = fieldnames(typeof(params))
    fields = map(n -> getfield(params, n), names)
    NamedTuple{names}(fields)
end

Base.getproperty(params::PIDParams, name::Symbol) = getproperty(params, Val(name))

@generated function Base.getproperty(params::PIDParams, ::Val{S}) where {S}
    if S ∈ fieldnames(PIDParams)
        return :(getfield(params, $(QuoteNode(S))))
    elseif S === :T_i
        return :(params.k_p / params.k_i)
    elseif S === :T_d
        return :(params.k_d / params.k_p)
    end
end

@kwdef struct PIDDiscrete <: SystemDefinition end

@kwdef mutable struct PIDDiscreteU
    k_p::Float64 = 1.0 #proportional gain
    k_i::Float64 = 0.1 #integral gain
    k_d::Float64 = 0.1 #derivative gain
    τ_f::Float64 = 0.01 #derivative filter time constant
    β_p::Float64 = 1.0 #proportional path setpoint weighting factor
    β_d::Float64 = 1.0 #derivative path setpoint weighting factor
    bound_lo::Float64 = -Inf #lower output bound
    bound_hi::Float64 = Inf #higher output bound
    input::Float64 = 0.0 #current input signal
    sat_ext::Int64 = 0 #external (signed) saturation input signal
end

@kwdef mutable struct PIDDiscreteS
    x_i0::Float64 = 0.0 #previous integrator path state
    x_d0::Float64 = 0.0 #previous derivative path state
    sat_out_0::Int64 = 0 #previous output saturation status
end

@kwdef struct PIDDiscreteY
    k_p::Float64 = 1.0 #proportional gain
    k_i::Float64 = 0.1 #integral gain
    k_d::Float64 = 0.1 #derivative gain
    τ_f::Float64 = 0.01 #derivative filter time constant
    β_p::Float64 = 1.0 #proportional path setpoint weighting factor
    β_d::Float64 = 1.0 #derivative path setpoint weighting factor
    bound_lo::Float64 = -Inf
    bound_hi::Float64 = Inf
    input::Float64 = 0.0
    sat_ext::Float64 = 0.0
    u_p::Float64 = 0.0 #proportional path input
    u_i::Float64 = 0.0 #integral path input
    u_d::Float64 = 0.0 #derivative path input
    y_p::Float64 = 0.0 #proportional path output
    y_i::Float64 = 0.0 #integral path output
    y_d::Float64 = 0.0 #derivative path output
    out_free::Float64 = 0 #unbounded output
    sat_out::Int64 = 0 #output saturation status
    output::Float64 = 0 #actual PID output
    int_halted::Bool = false #integration halted
end

Systems.init(::SystemY, ::PIDDiscrete) = PIDDiscreteY()
Systems.init(::SystemU, ::PIDDiscrete) = PIDDiscreteU()
Systems.init(::SystemS, ::PIDDiscrete) = PIDDiscreteS()

function reset!(sys::System{<:PIDDiscrete})
    sys.u.input = 0
    sys.u.sat_ext = 0
    sys.s.x_i0 = 0
    sys.s.x_d0 = 0
    sys.s.sat_out_0 = 0
    f_disc!(sys, 1.0)
end

function assign!(sys::System{<:PIDDiscrete}, params::PIDParams{<:Real})
    @unpack k_p, k_i, k_d, τ_f = params
    @pack! sys.u = k_p, k_i, k_d, τ_f
end

function Systems.f_disc!(sys::System{<:PIDDiscrete}, Δt::Real)

    @unpack k_p, k_i, k_d, τ_f, β_p, β_d, bound_lo, bound_hi, input, sat_ext = sys.u
    @unpack x_i0, x_d0, sat_out_0 = sys.s

    α = 1 / (τ_f + Δt)

    u_p = β_p * input
    u_d = β_d * input
    u_i = input

    int_halted = ((sign(u_i * sat_out_0) > 0) || (sign(u_i * sat_ext) > 0))

    x_i = x_i0 + Δt * u_i * !int_halted
    x_d = α * τ_f * x_d0 + Δt * α * k_d * u_d

    y_p = k_p * u_p
    y_i = k_i * x_i
    y_d = α * (-x_d0 + k_d * u_d)
    out_free = y_p + y_i + y_d

    sat_hi = out_free >= bound_hi
    sat_lo = out_free <= bound_lo
    sat_out = sat_hi - sat_lo

    output = clamp(out_free, bound_lo, bound_hi)

    sys.y = PIDDiscreteY(; k_p, k_i, k_d, τ_f, β_p, β_d, bound_lo, bound_hi, input, sat_ext,
                u_p, u_i, u_d, y_p, y_i, y_d, out_free, sat_out, output, int_halted)

    sys.s.x_i0 = x_i
    sys.s.x_d0 = x_d
    sys.s.sat_out_0 = sat_out

    return false

end

function Plotting.make_plots(th::TimeHistory{<:PIDDiscreteY}; kwargs...)

    pd = OrderedDict{Symbol, Plots.Plot}()

    input = plot(th.input; title = "Input", ylabel = L"$e$", kwargs...)
    output = plot(th.output; title = "Output", ylabel = L"$y$", kwargs...)

    pd[:sf] = plot(input, output;
        plot_title = "Input & Output",
        layout = (1,2),
        link = :y,
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    k_p = plot(th.β_p; title = "Proportional Gain", ylabel = L"$k_p$", kwargs...)
    k_i = plot(th.β_d; title = "Integral Gain", ylabel = L"$k_i$", kwargs...)
    k_d = plot(th.bound_lo; title = "Derivative Gain", ylabel = L"$k_d$", kwargs...)
    τ_f = plot(th.bound_hi; title = "Derivative Filter Time Constant", ylabel = L"$\tau_f$", kwargs...)

    pd[:p1] = plot(k_p, k_i, k_d, τ_f;
        plot_title = "Parameters",
        layout = (4,1),
        link = :none,
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    β_p = plot(th.β_p; title = "Proportional Input Weighting", ylabel = L"$\beta_p$", kwargs...)
    β_d = plot(th.β_d; title = "Derivative Input Weighting", ylabel = L"$\beta_d$", kwargs...)
    bound_lo = plot(th.bound_lo; title = "Lower Output Bound", ylabel = L"$y_{min}$", kwargs...)
    bound_hi = plot(th.bound_hi; title = "Upper Output Bound", ylabel = L"$y_{max}$", kwargs...)

    pd[:p2] = plot(β_p, β_d, bound_lo, bound_hi;
        plot_title = "Parameters",
        layout = (4,1),
        link = :none,
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    sat_ext = plot(th.sat_ext; title = "External Saturation Input", ylabel = "", kwargs...)
    sat_out = plot(th.sat_out; title = "Output Saturation", ylabel = "", kwargs...)
    int_halted = plot(th.int_halted; title = "Integrator Halted", ylabel = "", kwargs...)

    pd[:awu] = plot(sat_ext, sat_out, int_halted;
        plot_title = "Anti-Windup",
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

    pd[:int] = plot(u_i, y_i, int_halted;
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
    output = plot(th.output; title = "Actual", ylabel = L"$y$", kwargs...)

    pd[:output] = plot(out_free, output, sat_out;
        plot_title = "PID Output",
        layout = (1,3),
        link = :y,
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    return pd

end

#################################### GUI #######################################


function GUI.draw(sys::System{<:PIDDiscrete}, label::String = "PIDDiscrete")

    @unpack k_p, k_i, k_d, τ_f, β_p, β_d, bound_lo, bound_hi, input, sat_ext, u_p, u_i, u_d,
            y_p, y_i, y_d, out_free, sat_out, output, int_halted = sys.y

    # CImGui.Begin(label)

        CImGui.Text("Proportional Gain = $k_p")
        CImGui.Text("Integral Gain = $k_i")
        CImGui.Text("Derivative Gain = $k_d")
        CImGui.Text("Derivative Filter Time Constant = $τ_f")
        CImGui.Text("Proportional Input Weighting = $β_p")
        CImGui.Text("Derivative Input Weighting = $β_d")
        CImGui.Text("Lower Output Bound = $bound_lo")
        CImGui.Text("Upper Output Bound = $bound_hi")
        CImGui.Text("Input = $input")
        CImGui.Text("External Saturation Input = $sat_ext")
        CImGui.Text("Proportional Path Input = $u_p")
        CImGui.Text("Proportional Path Output = $y_p")
        CImGui.Text("Integral Path Input = $u_i")
        CImGui.Text("Integral Path Output = $y_i")
        CImGui.Text("Derivative Path Input = $u_d")
        CImGui.Text("Derivative Path Output = $y_d")
        CImGui.Text("Free Output = $out_free")
        CImGui.Text("Actual Output = $output")
        CImGui.Text("Output Saturation = $sat_out")
        CImGui.Text("Integrator Halted = $int_halted")

    # CImGui.End()

end #function


################################################################################
############################ PID Optimization ##################################
################################################################################

module PIDOpt

using StaticArrays
using UnPack
using NLopt
using ControlSystems
using Trapz: trapz

using ..Control

@kwdef struct Metrics{T} <: FieldVector{5, T}
    Ms::T #maximum sensitivity
    ∫e::T #integrated absolute error
    ef::T #absolute steady-state error
    ∫u::T #integrated absolute control effort
    up::T #absolute peak control effort
end

@kwdef struct Settings
    t_sim::Float64 = 5.0
    maxeval::Int64 = 5000
    lower_bounds::PIDParams = PIDParams(; k_p = 0.0, k_i = 0.0, k_d = 0.0, τ_f = 0.01)
    upper_bounds::PIDParams = PIDParams(; k_p = 50.0, k_i = 50.0, k_d = 10.0, τ_f = 0.01)
    initial_step::PIDParams = PIDParams(; k_p = 0.01, k_i = 0.01, k_d = 0.01, τ_f = 0.01)
end

@kwdef struct Results
    exit_flag::Symbol
    cost::Float64
    metrics::Metrics{Float64}
    params::PIDParams{Float64}
end

function build_PID(params::PIDParams{<:Real})
    @unpack k_p, k_i, k_d, τ_f = params
    (k_p + k_i * tf(1, [1,0]) + k_d * tf([1, 0], [τ_f, 1])) |> ss
end

function Metrics(plant::AbstractStateSpace, pid::AbstractStateSpace,
                       settings::Settings)

    #hinfnorm appears to be quite brittle, so instead we brute force the
    #computation of maximum sensitivity transfer function magnitude
    S = sensitivity(plant, pid) #sensitivity function
    S_tf = tf(S)
    iω_range = ((10^x)*im for x in range(-3, 3, length=1000))
    S_range = [abs(S_tf.(iω)[1]) for iω in iω_range]
    Ms = maximum(S_range)

    T = output_comp_sensitivity(plant, pid) #complementary sensitivity function (AKA closed loop)
    T_step = step(T, settings.t_sim)
    t = T_step.t
    y = T_step.y |> vec
    abs_e = abs.(y .- 1.0)
    #integrated error normalized with respect to the length of the time window
    ∫e = trapz(t, abs_e)/t[end]
    ef = abs_e[end]

    CS = G_CS(plant, pid) #reference to control input
    CS_step = step(CS, settings.t_sim)
    t = CS_step.t
    y = CS_step.y |> vec
    abs_u = abs.(y .- 1.0)
    ∫u = trapz(t, abs_u)/t[end]
    up = maximum(abs_u)

    Metrics(; Ms, ∫e, ef, ∫u, up)

end

function cost(plant::AbstractStateSpace, pid::AbstractStateSpace,
              settings::Settings, weights::Metrics{<:Real})
    costs = Metrics(plant, pid, settings)
    return sum(costs .* weights) / sum(weights)
end


function optimize_PID(  plant::LTISystem;
                    params_0::PIDParams = PIDParams(), #initial condition
                    settings::Settings = Settings(),
                    weights::Metrics{<:Real} = Metrics(ones(5)),
                    global_search::Bool = true)

    x0 = params_0 |> Vector
    lower_bounds = settings.lower_bounds |> Vector
    upper_bounds = settings.upper_bounds |> Vector
    initial_step = settings.initial_step |> Vector
    maxeval = settings.maxeval

    x0 = clamp.(x0, lower_bounds, upper_bounds)

    plant = ss(plant)
    f_opt = let plant = plant, settings = settings, weights = weights
        function (x::Vector{Float64}, ::Vector{Float64})
            pid = build_PID(PIDParams(x...))
            cost(plant, pid, settings, weights)
        end
    end

    minx = x0

    if global_search
        opt_glb = Opt(:GN_DIRECT_L, length(x0))
        opt_glb.min_objective = f_opt
        opt_glb.maxeval = maxeval
        opt_glb.stopval = 1e-5
        opt_glb.lower_bounds = lower_bounds
        opt_glb.upper_bounds = upper_bounds
        opt_glb.initial_step = initial_step

        (minf, minx, exit_flag) = optimize(opt_glb, x0)
    end

    #second pass with local optimizer
    opt_loc = Opt(:LN_BOBYQA, length(x0))
    # opt_loc = Opt(:LN_COBYLA, length(x0))
    opt_loc.min_objective = f_opt
    opt_loc.maxeval = 5000
    opt_loc.stopval = 1e-5
    opt_loc.lower_bounds = lower_bounds
    opt_loc.upper_bounds = upper_bounds
    opt_loc.initial_step = initial_step

    (minf, minx, exit_flag) = optimize(opt_loc, minx)

    params_opt = PIDParams(minx...)
    pid_opt = build_PID(params_opt)
    metrics_opt = Metrics(plant, pid_opt, settings)

    return Results(exit_flag, minf, metrics_opt, params_opt)


end

function check_results(results::Results, thresholds::Metrics{Float64})

    @unpack exit_flag, metrics = results

    success = true
    success && (exit_flag === :ROUNDOFF_LIMITED) | (exit_flag === :STOPVAL_REACHED)
    success && (metrics.Ms < thresholds.Ms) #sensitivity function maximum magnitude
    success && (metrics.∫e < thresholds.∫e) #normalized absolute integrated error
    success && (metrics.ef < thresholds.ef) #remaining error after t_sim

end

end ######################### submodule ########################################

end #module