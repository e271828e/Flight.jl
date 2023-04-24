module Control

using ComponentArrays, StaticArrays, UnPack, LinearAlgebra
using CImGui, CImGui.CSyntax
using ControlSystems: ControlSystems, ss

using Flight.FlightCore

export LinearStateSpace, PIContinuous, PIDDiscrete

################################################################################
########################### LinearStateSpace ###################################

const tV = AbstractVector{<:Float64}
const tM = AbstractMatrix{<:Float64}

struct LinearStateSpace{ LX, LU, LY, #state, input and output vector lengths
                        tX <: tV, tU <: tV, tY <: tV,
                        tA <: tM, tB <: tM, tC <: tM, tD <: tM} <: Component

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

struct PIContinuous{N} <: Component #Parallel form
    k_p::SVector{N,Float64} #proportional gain
    k_i::SVector{N,Float64} #integral gain
    k_l::SVector{N,Float64} #integrator leak factor
    bounds::NTuple{2,MVector{N,Float64}} #output bounds
end

function PIContinuous{N}(; k_p::Real = 1.0, k_i::Real = 0.1, k_l::Real = 0.0,
                            bounds::NTuple{2,Real} = (-1.0, 1.0)) where {N}
    s2v = (x)->fill(x,N)
    PIContinuous{N}(s2v(k_p), s2v(k_i), s2v(k_l), (s2v(bounds[1]), s2v(bounds[2])))
end

Base.@kwdef struct PIContinuousU{N}
    setpoint::MVector{N,Float64} = zeros(N) #commanded setpoint
    feedback::MVector{N,Float64} = zeros(N) #plant feedback (non-inverted)
    hold::MVector{N,Bool} = zeros(Bool, N) #hold integrator state
    reset::MVector{N,Bool} = zeros(Bool, N) #reset integrator state
    sat_enable::MVector{N,Bool} = ones(Bool, N)
end

Base.@kwdef struct PIContinuousY{N}
    hold::SVector{N,Bool} = zeros(SVector{N, Bool}) #hold integrator state
    reset::SVector{N,Bool} = zeros(SVector{N, Bool}) #reset integrator state
    error::SVector{N,Float64} = zeros(SVector{N}) #error = setpoint-feedback
    state::SVector{N,Float64} = zeros(SVector{N}) #integrator state
    out_p::SVector{N,Float64} = zeros(SVector{N}) #proportional term
    out_i::SVector{N,Float64} = zeros(SVector{N}) #integral term
    out_free::SVector{N,Float64} = zeros(SVector{N}) #total output, free
    out::SVector{N,Float64} = zeros(SVector{N}) #actual output
    sat_enable::SVector{N,Bool} = zeros(SVector{N, Bool}) #saturation enabled
    sat_status::SVector{N,Int64} = zeros(SVector{N, Int64}) #saturation status
    int_status::SVector{N,Bool} = zeros(SVector{N, Bool}) #integrator active
end

Systems.init(::SystemX, ::PIContinuous{N}) where {N} = zeros(N)
Systems.init(::SystemY, ::PIContinuous{N}) where {N} = PIContinuousY{N}()
Systems.init(::SystemU, ::PIContinuous{N}) where {N} = PIContinuousU{N}()

function Systems.f_ode!(sys::System{<:PIContinuous{N}}) where {N}

    @unpack ẋ, x, u, params = sys
    @unpack k_p, k_i, k_l, bounds = params

    state = SVector{N, Float64}(x)
    setpoint = SVector(u.setpoint)
    feedback = SVector(u.feedback)
    hold = SVector(u.hold)
    reset = SVector(u.reset)
    sat_enable = SVector(u.sat_enable)

    error = setpoint - feedback
    out_p = k_p .* error
    out_i = k_i .* state
    out_free = out_p + out_i #raw output
    out_clamped = clamp.(out_free, bounds[1], bounds[2]) #clamped output
    out = (out_free .* .!sat_enable) .+ (out_clamped .* sat_enable)

    sat_upper = (out_free .>= bounds[2]) .* sat_enable
    sat_lower = (out_free .<= bounds[1]) .* sat_enable
    sat_status = sat_upper - sat_lower
    int_status = sign.(error .* sat_status) .<= 0 #enable integrator?

    ẋ .= (error .* int_status - k_l .* state) .* .!reset .* .!hold

    sys.y = PIContinuousY(; hold, reset, error, state, out_p, out_i,
                            out_free, out, sat_enable, sat_status, int_status)

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

    splt_e = plot(th.error; title = "Error",
        ylabel = L"$e$", kwargs...)

    splt_s = plot(th.state; title = "Integrator State",
        ylabel = L"$s$", kwargs...)

    splt_y_p = plot(th.out_p; title = "Proportional Term",
        ylabel = L"$y_p$", kwargs...)

    splt_y_i = plot(th.out_i; title = "Integral Term",
        ylabel = L"$y_i$", kwargs...)

    splt_y_free = plot(th.out_free; title = "Unbounded Output",
        ylabel = L"$y_{free}$", kwargs...)

    splt_y = plot(th.out; title = "Actual Output",
        ylabel = L"$y$", kwargs...)

    splt_sat = plot(th.sat_status; title = "Saturation Status",
        ylabel = L"$S$", kwargs...)

    pd[:es] = plot(splt_e, splt_s, splt_sat;
        plot_title = "Error & Integrator State",
        layout = (1,3),
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    pd[:pi] = plot(splt_y_p, splt_y_i, splt_sat;
        plot_title = "Proportional & Integral Terms",
        layout = (1,3),
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    pd[:out] = plot(splt_y_free, splt_y, splt_sat;
        plot_title = "Total Output",
        layout = (1,3),
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    return pd

end

#################################### GUI #######################################


function GUI.draw(sys::System{<:PIContinuous{N}}, label::String = "PIContinuous{$N}") where {N}

    @unpack y, params = sys

    CImGui.Begin(label)

    if CImGui.TreeNode("Params")
        @unpack k_p, k_i, k_l, bounds = params
        CImGui.Text("Proportional Gain: $k_p")
        CImGui.Text("Integral Gain: $k_i")
        CImGui.Text("Leak Factor: = $k_l")
        CImGui.Text("Output Bounds = $bounds")
        CImGui.TreePop()
    end


    if CImGui.TreeNode("Outputs")
        @unpack hold, reset, error, sat_enable, sat_status, out_free, out = y
        CImGui.Text("Hold = $hold")
        CImGui.Text("Reset = $reset")
        CImGui.Text("Error = $error")
        CImGui.Text("Unbounded Output = $out_free")
        CImGui.Text("Output = $out")
        CImGui.Text("Saturation Enabled = $sat_enable")
        CImGui.Text("Saturation Status = $sat_status")
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
                @unpack hold, reset, sat_enable = sys.u
                hold_ref = Ref(hold[i])
                reset_ref = Ref(reset[i])
                sat_ref = Ref(sat_enable[i])
                CImGui.Checkbox("Hold", hold_ref)
                CImGui.SameLine()
                CImGui.Checkbox("Reset", reset_ref)
                CImGui.SameLine()
                CImGui.Checkbox("Enable Saturation", sat_ref)
                hold[i] = hold_ref[]
                reset[i] = reset_ref[]
                sat_enable[i] = sat_ref[]
                CImGui.TreePop()
            end
        end
        CImGui.TreePop()
    end

    CImGui.End()

end

################# Proportional-Integral-Derivative Compensator #################
################################################################################


struct PIDDiscrete{N} <: Component #Parallel form
    k_p::SVector{N,Float64} #proportional gain
    k_i::SVector{N,Float64} #integral gain
    k_d::SVector{N,Float64} #derivative gain
    τ_d::SVector{N,Float64} #derivative filter time constant
    β_p::SVector{N,Float64} #proportional path setpoint weighting factor
    β_d::SVector{N,Float64} #derivative path setpoint weighting factor
    bounds::NTuple{2,MVector{N,Float64}} #output bounds
end

function PIDDiscrete{N}(; k_p::Real = 1.0, k_i::Real = 0.1,
                        k_d::Real = 0, τ_d::Real = 0.05,
                        β_p::Real = 1.0, β_d::Real = 0,
                        bounds::NTuple{2,Real} = (-1.0, 1.0)) where {N}
    s2v = (x)->fill(x,N)
    PIDDiscrete{N}( s2v(k_p), s2v(k_i),
                   s2v(k_d), s2v(τ_d),
                   s2v(β_p), s2v(β_d),
                   (s2v(bounds[1]), s2v(bounds[2])))
end

Base.@kwdef struct PIDDiscreteU{N}
    setpoint::MVector{N,Float64} = zeros(N) #commanded setpoint
    feedback::MVector{N,Float64} = zeros(N) #plant feedback (non-inverted)
    sat_enable::MVector{N,Bool} = ones(Bool, N) #enable PID output saturation
    int_hold::MVector{N,Bool} = zeros(Bool, N) #hold integrator state
    reset::MVector{N,Bool} = zeros(Bool, N) #reset PID
end

Base.@kwdef struct PIDDiscreteS{N}
    x_i0::MVector{N,Float64} = zeros(N) #previous integrator path state
    x_d0::MVector{N,Float64} = zeros(N) #previous derivative path state
    sat_0::MVector{N,Int64} = zeros(N) #previous PID saturation status
end

Base.@kwdef struct PIDDiscreteY{N}
    setpoint::SVector{N,Float64} = zeros(SVector{N}) #commanded setpoint
    feedback::SVector{N,Float64} = zeros(SVector{N}) #plant feedback (non-inverted)
    sat_enable::SVector{N,Bool} = zeros(SVector{N, Bool}) #saturation enabled
    int_hold::SVector{N,Bool} = zeros(SVector{N, Bool}) #integrator state hold
    reset::SVector{N,Bool} = zeros(SVector{N, Bool}) #reset PID
    u_p::SVector{N,Float64} = zeros(SVector{N}) #proportional path input
    u_i::SVector{N,Float64} = zeros(SVector{N}) #integral path input
    u_d::SVector{N,Float64} = zeros(SVector{N}) #derivative path input
    int_active::SVector{N,Bool} = zeros(SVector{N, Bool}) #integrator active
    y_p::SVector{N,Float64} = zeros(SVector{N}) #proportional path output
    y_i::SVector{N,Float64} = zeros(SVector{N}) #integral path output
    y_d::SVector{N,Float64} = zeros(SVector{N}) #derivative path output
    out_free::SVector{N,Float64} = zeros(SVector{N}) #non-clamped PID output
    out::SVector{N,Float64} = zeros(SVector{N}) #actual PID output
    sat::SVector{N,Int64} = zeros(SVector{N, Int64}) #current saturation status
end

Systems.init(::SystemY, ::PIDDiscrete{N}) where {N} = PIDDiscreteY{N}()
Systems.init(::SystemU, ::PIDDiscrete{N}) where {N} = PIDDiscreteU{N}()
Systems.init(::SystemS, ::PIDDiscrete{N}) where {N} = PIDDiscreteS{N}()

function Systems.f_disc!(sys::System{<:PIDDiscrete{N}}, Δt::Real) where {N}

    @unpack s, u, params = sys
    @unpack k_p, k_i, k_d, τ_d, β_p, β_d, bounds = params

    setpoint = SVector(u.setpoint)
    feedback = SVector(u.feedback)
    sat_enable = SVector(u.sat_enable)
    int_hold = SVector(u.int_hold)
    reset = SVector(u.reset)

    x_i0 = SVector(s.x_i0)
    x_d0 = SVector(s.x_d0)
    sat_0 = SVector(s.sat_0)

    u_p = β_p .* setpoint - feedback
    u_d = β_d .* setpoint - feedback
    u_i = setpoint - feedback

    #integration must only operate in those components for which:
    # -integrator path input does not have the same sign as PID output saturation
    # -integrator state is not being held
    int_active = (sign.(u_i .* sat_0) .<= 0) .* .!int_hold
    α = 1 ./ (τ_d .+ Δt)

    x_i = x_i0 + Δt * k_i .* u_i .* int_active
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
    out_clamped = clamp.(out_free, bounds[1], bounds[2])
    out = (out_free .* .!sat_enable) .+ (out_clamped .* sat_enable)

    sat_upper = (out_free .>= bounds[2]) .* sat_enable
    sat_lower = (out_free .<= bounds[1]) .* sat_enable
    sat = sat_upper - sat_lower

    s.x_i0 .= x_i
    s.x_d0 .= x_d
    s.sat_0 .= sat

    sys.y = PIDDiscreteY(; setpoint, feedback, sat_enable, int_hold, reset,
                           u_p, u_i, u_d, int_active, y_p, y_i, y_d,
                           out_free, out, sat)

    return false

end

function Plotting.make_plots(th::TimeHistory{<:PIDDiscreteY}; kwargs...)

    pd = OrderedDict{Symbol, Plots.Plot}()

    setpoint = plot(th.setpoint; title = "Setpoint", ylabel = L"$r$", kwargs...)
    feedback = plot(th.feedback; title = "Feedback", ylabel = L"$f$", kwargs...)

    pd[:sf] = plot(setpoint, feedback;
        plot_title = "Setpoint & Feedback",
        layout = (1,2),
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    sat_enable = plot(th.sat_enable; title = "Saturation Enabled", ylabel = "", kwargs...)
    int_hold = plot(th.int_hold; title = "Integrator Hold", ylabel = "", kwargs...)
    reset = plot(th.reset; title = "Reset", ylabel = "", kwargs...)

    pd[:ctl] = plot(sat_enable, int_hold, reset;
        plot_title = "Control Flags",
        layout = (1,3),
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    u_p = plot(th.u_p; title = "Input", ylabel = L"$u_p$", kwargs...)
    y_p = plot(th.y_p; title = "Output", ylabel = L"$y_p$", kwargs...)

    pd[:prop] = plot(u_p, y_p;
        plot_title = "Proportional Path",
        layout = (1,2),
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    u_i = plot(th.u_i; title = "Input", ylabel = L"$u_i$", kwargs...)
    y_i = plot(th.y_i; title = "Output", ylabel = L"$y_i$", kwargs...)

    pd[:int] = plot(u_i, y_i;
        plot_title = "Integral Path",
        layout = (1,2),
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    u_d = plot(th.u_d; title = "Input", ylabel = L"$u_d$", kwargs...)
    y_d = plot(th.y_d; title = "Output", ylabel = L"$y_d$", kwargs...)

    pd[:der] = plot(u_d, y_d;
        plot_title = "Derivative Path",
        layout = (1,2),
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    y_free = plot(th.out_free; title = "Unclamped", ylabel = L"$y_{free}$", kwargs...)
    y = plot(th.out; title = "Actual", ylabel = L"$y$", kwargs...)
    sat = plot(th.sat; title = "Saturation Status", ylabel = L"$S$", kwargs...)

    pd[:out] = plot(y_free, y, sat;
        plot_title = "PID Output",
        layout = (1,3),
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    return pd

end


end #module