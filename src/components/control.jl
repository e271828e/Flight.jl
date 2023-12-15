module Control


################################################################################
########################## Continuous Systems ##################################
################################################################################

module Continuous

using ComponentArrays, StaticArrays, UnPack, LinearAlgebra
using ControlSystems: ControlSystemsBase, ControlSystems, ss
using RobustAndOptimalControl

using Flight.FlightCore

################################################################################
########################### LinearizedSS ###################################

const tV = AbstractVector{<:Float64}
const tM = AbstractMatrix{<:Float64}

struct LinearizedSS{ LX, LU, LY, #state, input and output vector lengths
                        tX <: tV, tU <: tV, tY <: tV,
                        tA <: tM, tB <: tM, tC <: tM, tD <: tM} <: SystemDefinition

    ẋ0::tX; x0::tX; u0::tU; y0::tY; #reference values (for linearized systems)
    A::tA; B::tB; C::tC; D::tD; #state-space matrices
    x_cache::tX; y_cache::tY; y_cache_out::tY;
    Δx_cache::tX; Δu_cache::tU

    function LinearizedSS(ẋ0, x0, u0, y0, A, B, C, D)

        lengths = map(length, (x0, u0, y0))
        types = map(typeof, (x0, u0, y0, A, B, C, D))

        vectors = map(copy, (ẋ0, x0, u0, y0))
        matrices = map(copy, (A, B, C, D))
        caches = map(copy, (x0, y0, y0, x0, u0))

        new{lengths..., types...}(vectors..., matrices..., caches...)

    end

end

LinearizedSS(; ẋ0, x0, u0, y0, A, B, C, D) = LinearizedSS(ẋ0, x0, u0, y0, A, B, C, D)

ControlSystems.ss(cmp::LinearizedSS) = ControlSystems.ss(cmp.A, cmp.B, cmp.C, cmp.D)

function RobustAndOptimalControl.named_ss(lss::LinearizedSS)
    x_labels, u_labels, y_labels = map(collect ∘ propertynames, (lss.x0, lss.u0, lss.y0))
    named_ss(ss(lss), x = x_labels, u = u_labels, y = y_labels)
end

function LinearizedSS(sys::ControlSystemsBase.StateSpace{ControlSystemsBase.Continuous, <:AbstractFloat})
    @unpack A, B, C, D, nx, nu, ny = sys
    ẋ0 = zeros(nx); x0 = zeros(nx); u0 = zeros(nu); y0 = zeros(ny)
    LinearizedSS(; ẋ0, x0, u0, y0, A, B, C, D)
end

Systems.X(cmp::LinearizedSS) = copy(cmp.x0)
Systems.U(cmp::LinearizedSS) = copy(cmp.u0)
Systems.Y(cmp::LinearizedSS) = SVector{length(cmp.y0)}(cmp.y0)

function Systems.f_ode!(sys::System{<:LinearizedSS{LX, LU, LY}}) where {LX, LU, LY}

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

function submodel(cmp::LinearizedSS; x = keys(cmp.x0), u = keys(cmp.u0), y = keys(cmp.y0))

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

    return LinearizedSS(; ẋ0, x0, u0, y0, A, B, C, D)

end

# #convert a mutable ComponentVector into a (non-allocating) immutable
# #ComponentVector with an underlying StaticVector, preserving the original
# #ComponentVector's axes
# function static_cv(L::Integer, x::ComponentVector)
#     ComponentVector(SVector{L}(getdata(x)), getaxes(x))
# end


####################### Proportional-Integral Compensator ######################
################################################################################

struct PIVector{N} <: SystemDefinition end

@kwdef struct PIVectorInput{N}
    k_p::MVector{N,Float64} = ones(N) #proportional gain
    k_i::MVector{N,Float64} = zeros(N) #integral gain
    k_l::MVector{N,Float64} = zeros(N) #integrator leak factor
    β_p::MVector{N,Float64} = ones(N) #proportional path setpoint weighting factor
    bound_lo::MVector{N,Float64} = fill(-Inf, N) #lower output bounds
    bound_hi::MVector{N,Float64} = fill(Inf, N) #higher output bounds
    input::MVector{N,Float64} = zeros(N) #input (setpoint - feedback)
    sat_ext::MVector{N,Int64} = zeros(Int64, N) #external (signed) saturation signal
    reset::MVector{N,Bool} = zeros(Bool, N) #reset PID continuous state
end

@kwdef struct PIVectorOutput{N}
    k_p::SVector{N,Float64} = ones(SVector{N})
    k_i::SVector{N,Float64} = zeros(SVector{N})
    k_l::SVector{N,Float64} = zeros(SVector{N})
    β_p::SVector{N,Float64} = ones(SVector{N})
    bound_lo::SVector{N,Float64} = fill(-Inf, SVector{N}) #lower output bounds
    bound_hi::SVector{N,Float64} = fill(Inf, SVector{N}) #higher output bounds
    input::SVector{N,Float64} = zeros(SVector{N})
    sat_ext::SVector{N,Int64} = zeros(SVector{N, Int64}) #external (signed) saturation signal
    reset::SVector{N,Bool} = zeros(SVector{N, Bool}) #reset PID states and null outputs
    u_p::SVector{N,Float64} = zeros(SVector{N}) #proportional path input
    u_i::SVector{N,Float64} = zeros(SVector{N}) #integral path input
    y_p::SVector{N,Float64} = zeros(SVector{N}) #proportional path output
    y_i::SVector{N,Float64} = zeros(SVector{N}) #integral path output
    out_free::SVector{N,Float64} = zeros(SVector{N}) #total output, free
    sat_out::SVector{N,Int64} = zeros(SVector{N, Int64}) #current output saturation status
    output::SVector{N,Float64} = zeros(SVector{N}) #actual output
    int_halted::SVector{N,Bool} = zeros(SVector{N, Bool}) #integration halted
end

Systems.X(::PIVector{N}) where {N} = zeros(N)
Systems.Y(::PIVector{N}) where {N} = PIVectorOutput{N}()
Systems.U(::PIVector{N}) where {N} = PIVectorInput{N}()

function Systems.f_ode!(sys::System{<:PIVector{N}}) where {N}

    @unpack ẋ, x, u = sys

    x_i = SVector{N, Float64}(x)
    k_p, k_i, k_l, β_p = map(SVector, (u.k_p, u.k_i, u.k_l, u.β_p))
    input, bound_lo, bound_hi, sat_ext, reset = map(SVector, (
        u.input, u.bound_lo, u.bound_hi, u.sat_ext, u.reset))

    u_p = β_p .* input
    u_i = input

    y_p = k_p .* u_p .* .!reset
    y_i = x_i .* .!reset
    out_free = y_p + y_i #raw output
    output = clamp.(out_free, bound_lo, bound_hi) #clamped output

    sat_hi = out_free .>= bound_hi
    sat_lo = out_free .<= bound_lo
    sat_out = sat_hi - sat_lo
    int_halted = ((sign.(u_i .* sat_out) .> 0) .|| (sign.(u_i .* sat_ext) .> 0))

    ẋ .= (k_i .* u_i .* .!int_halted - k_l .* x_i) .* .!reset

    sys.y = PIVectorOutput(; k_p, k_i, k_l, β_p, bound_lo, bound_hi, input, sat_ext, reset,
                     u_p, u_i, y_p, y_i, out_free, sat_out, output, int_halted)

end

function Systems.f_step!(sys::System{<:PIVector{N}}) where {N}

    x = SVector{N,Float64}(sys.x)
    x_new = x .* .!sys.u.reset
    x_mod = any(x .!= x_new)

    sys.x .= x_new
    return x_mod

end


# ############################## Plotting ########################################

function Plotting.make_plots(ts::TimeSeries{<:PIVectorOutput}; kwargs...)

    pd = OrderedDict{Symbol, Plots.Plot}()

    input = plot(ts.input; title = "Input", ylabel = L"$e$", kwargs...)
    output = plot(ts.output; title = "Output", ylabel = L"$y$", kwargs...)
    reset = plot(ts.reset; title = "Reset", ylabel = L"$r$", kwargs...)

    pd[:sf] = plot(input, output, reset;
        plot_title = "Input, Output & Reset",
        layout = (1,3),
        link = :y,
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    k_p = plot(ts.k_p; title = "Proportional Gain", ylabel = L"$k_p$", kwargs...)
    k_i = plot(ts.k_i; title = "Integral Gain", ylabel = L"$k_i$", kwargs...)
    k_l = plot(ts.k_l; title = "Leak Factor", ylabel = L"$k_d$", kwargs...)

    pd[:p1] = plot(k_p, k_i, k_l;
        plot_title = "Parameters",
        layout = (3,1),
        link = :none,
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    β_p = plot(ts.β_p; title = "Proportional Input Weighting", ylabel = L"$\beta_p$", kwargs...)
    bound_lo = plot(ts.bound_lo; title = "Lower Output Bound", ylabel = L"$y_{min}$", kwargs...)
    bound_hi = plot(ts.bound_hi; title = "Upper Output Bound", ylabel = L"$y_{max}$", kwargs...)

    pd[:p2] = plot(β_p, bound_lo, bound_hi;
        plot_title = "Parameters",
        layout = (3,1),
        link = :none,
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    sat_ext = plot(ts.sat_ext; title = "External Saturation Input", ylabel = "", kwargs...)
    sat_out = plot(ts.sat_out; title = "Output Saturation", ylabel = "", kwargs...)
    int_halted = plot(ts.int_halted; title = "Integrator Halted", ylabel = "", kwargs...)

    pd[:awu] = plot(sat_ext, sat_out, int_halted;
        plot_title = "Anti-Windup",
        layout = (1,3),
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    u_p = plot(ts.u_p; title = "Input", ylabel = L"$u_p$", kwargs...)
    y_p = plot(ts.y_p; title = "Output", ylabel = L"$y_p$", kwargs...)

    pd[:prop] = plot(u_p, y_p;
        plot_title = "Proportional Path",
        layout = (1,2),
        link = :y,
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    u_i = plot(ts.u_i; title = "Input", ylabel = L"$u_i$", kwargs...)
    y_i = plot(ts.y_i; title = "Output", ylabel = L"$y_i$", kwargs...)

    pd[:int] = plot(u_i, y_i, int_halted;
        plot_title = "Integral Path",
        layout = (1,3),
        link = :y,
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    out_free = plot(ts.out_free; title = "Free", ylabel = L"$y_{free}$", kwargs...)
    output = plot(ts.output; title = "Actual", ylabel = L"$y$", kwargs...)

    pd[:output] = plot(out_free, output, sat_out;
        plot_title = "PID Output",
        layout = (1,3),
        link = :y,
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    return pd

end

#################################### GUI #######################################


function GUI.draw(sys::System{<:PIVector{N}}, label::String = "PIVector{$N}") where {N}

    @unpack k_p, k_i, k_l, β_p, bound_lo, bound_hi, input, sat_ext, reset,
            u_p, u_i, y_p, y_i, out_free, sat_out, output, int_halted = sys.y

    # CImGui.Begin(label)

        CImGui.Text("Proportional Gain = $k_p")
        CImGui.Text("Integral Gain = $k_i")
        CImGui.Text("Integrator Leak Factor = $k_l")
        CImGui.Text("Proportional Input Weighting = $β_p")
        CImGui.Text("Lower Output Bound = $bound_lo")
        CImGui.Text("Upper Output Bound = $bound_hi")
        CImGui.Text("Input = $input")
        CImGui.Text("External Saturation Input = $sat_ext")
        CImGui.Text("Reset Input = $reset")
        CImGui.Text("Proportional Path Input = $u_p")
        CImGui.Text("Proportional Path Output = $y_p")
        CImGui.Text("Integral Path Input = $u_i")
        CImGui.Text("Integral Path Output = $y_i")
        CImGui.Text("Free Output = $out_free")
        CImGui.Text("Output Saturation = $sat_out")
        CImGui.Text("Actual Output = $output")
        CImGui.Text("Integrator Halted = $int_halted")

    # CImGui.End()

end #function

function GUI.draw!(sys::System{<:PIVector{N}}, label::String = "PIVector{$N}") where {N}

    CImGui.Begin(label)

    if CImGui.TreeNode("Outputs")
        GUI.draw(sys, label)
        CImGui.TreePop()
    end
    #shows how to get input from widgets without the @cstatic and @c macros
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

end #submodule


################################################################################
############################# Discrete Systems #################################
################################################################################

module Discrete

using StaticArrays, UnPack, LinearAlgebra
using RobustAndOptimalControl

using Flight.FlightCore

using ..Control


############################# Integrator ###############################
################################################################################

struct Integrator <: SystemDefinition end

@kwdef mutable struct IntegratorInput
    input::Float64 = 0 #current input
    sat_ext::Int64 = 0 #external (signed) saturation signal
    bound_lo::Float64 = -Inf #lower output bound
    bound_hi::Float64 = Inf #higher output bound
end

@kwdef mutable struct IntegratorState
    x0::Float64 = 0 #previous integrator state
    sat_out_0::Int64 = 0 #previous output saturation state
end

@kwdef struct IntegratorOutput
    input::Float64 = 0
    sat_ext::Float64 = 0
    bound_lo::Float64 = -Inf
    bound_hi::Float64 = Inf
    x1::Float64 = 0 #current state
    output::Float64 = 0 #current output
    sat_out::Int64 = 0 #current output saturation status
    halted::Bool = false #integration halted
end

Systems.Y(::Integrator) = IntegratorOutput()
Systems.U(::Integrator) = IntegratorInput()
Systems.S(::Integrator) = IntegratorState()

function Systems.reset!(sys::System{<:Integrator}, x0::Real = 0.0)
    sys.u.input = 0
    sys.u.sat_ext = 0
    sys.s.x0 = x0
    sys.s.sat_out_0 = 0
end

function Systems.f_disc!(sys::System{<:Integrator}, Δt::Real)

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

    sys.y = IntegratorOutput(; input, sat_ext, bound_lo, bound_hi, x1, output, sat_out, halted)

    return false

end

################################################################################

struct IntegratorVector{N} <: SystemDefinition end

@kwdef mutable struct IntegratorVectorInput{N}
    input::MVector{N,Float64} = zeros(Float64, N)
    sat_ext::MVector{N,Int64} = zeros(Int64, N) #external (signed) saturation signal
    bound_lo::MVector{N,Float64} = fill(-Inf, N) #lower output bounds
    bound_hi::MVector{N,Float64} = fill(Inf, N) #higher output bounds
end

@kwdef mutable struct IntegratorVectorState{N}
    x0::MVector{N,Float64} = zeros(N) #previous integrator state
    sat_out_0::MVector{N,Int64} = zeros(N) #previous output saturation status
end

@kwdef struct IntegratorVectorOutput{N}
    input::SVector{N,Float64} = zeros(SVector{N})
    sat_ext::SVector{N,Int64} = zeros(SVector{N, Int64}) #external (signed) saturation signal
    bound_lo::SVector{N,Float64} = fill(-Inf, SVector{N}) #lower output bounds
    bound_hi::SVector{N,Float64} = fill(Inf, SVector{N}) #higher output bounds
    x1::SVector{N,Float64} = zeros(N) #current integrator state
    output::SVector{N,Float64} = zeros(SVector{N}) #actual PID output
    sat_out::SVector{N,Int64} = zeros(SVector{N, Int64}) #output saturation status
    halted::SVector{N,Bool} = zeros(SVector{N, Bool}) #integration halted
end

Systems.Y(::IntegratorVector{N}) where {N} = IntegratorVectorOutput{N}()
Systems.U(::IntegratorVector{N}) where {N} = IntegratorVectorInput{N}()
Systems.S(::IntegratorVector{N}) where {N} = IntegratorVectorState{N}()

function Systems.reset!(sys::System{<:IntegratorVector{N}}, x0::AbstractVector{<:Real} = zeros(SVector{N})) where {N}
    sys.u.input .= 0
    sys.u.sat_ext .= 0
    sys.s.x0 .= x0
    sys.s.sat_out_0 .= 0
end

function Systems.f_disc!(sys::System{<:IntegratorVector}, Δt::Real)

    @unpack s, u = sys

    input, bound_lo, bound_hi, sat_ext = map(SVector, (
        u.input, u.bound_lo, u.bound_hi, u.sat_ext))

    x0 = SVector(s.x0)
    sat_out_0 = SVector(s.sat_out_0)

    halted = ((sign.(input .* sat_out_0) .> 0) .|| (sign.(input .* sat_ext) .> 0))
    x1 = x0 + Δt .* input .* .!halted
    output = clamp.(x1, bound_lo, bound_hi)

    sat_hi = x1 .>= bound_hi
    sat_lo = x1 .<= bound_lo
    sat_out = sat_hi - sat_lo

    s.x0 .= x1
    s.sat_out_0 .= sat_out

    sys.y = IntegratorVectorOutput(; input, sat_ext, bound_lo, bound_hi, x1, output, sat_out, halted)

    return false

end


#################################### GUI #######################################


function GUI.draw(sys::Union{System{<:Integrator}, System{<:IntegratorVector}},
                    label::String = "Integrator")

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

struct LeadLag <: SystemDefinition end

@kwdef mutable struct LeadLagInput
    z::Float64 = -1.0 #zero location (z < 0)
    p::Float64 = -10.0 #pole location (p < 0)
    k::Float64 = 1.0 #gain
    u1::Float64 = 0.0 #input
end

@kwdef mutable struct LeadLagState
    u0::Float64 = 0.0 #previous input
    x0::Float64 = 0.0 #previous output
end

@kwdef struct LeadLagOutput
    z::Float64 = -1.0
    p::Float64 = -10.0
    k::Float64 = 1.0
    u1::Float64 = 0.0
    y1::Float64 = 0.0 #current output
end

Systems.Y(::LeadLag) = LeadLagOutput()
Systems.U(::LeadLag) = LeadLagInput()
Systems.S(::LeadLag) = LeadLagState()

function Systems.reset!(sys::System{<:LeadLag})
    sys.u.u1 = 0
    sys.s.u0 = 0
    sys.s.x0 = 0
end

function Systems.f_disc!(sys::System{<:LeadLag}, Δt::Real)

    @unpack s, u = sys
    @unpack z, p, k, u1 = u
    @unpack u0, x0 = s

    a0 = (2+p*Δt)/(2-p*Δt)
    b1 = (2-z*Δt)/(2-p*Δt)
    b0 = (-2-z*Δt)/(2-p*Δt)

    x1 = a0 * x0 + b1 * u1 + b0 * u0
    y1 = k * x1

    sys.y = LeadLagOutput(; z, p, k, u1, y1)

    s.x0 = x1
    s.u0 = u1

    return false

end


#################################### GUI #######################################

function GUI.draw(sys::System{<:LeadLag}, label::String = "Discrete Lead Compensator")

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
    else
        return :(error("PIDParams has no property $S"))
    end
end

################################# Scalar Version ###############################

@kwdef struct PID <: SystemDefinition end

@kwdef mutable struct PIDInput
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

@kwdef mutable struct PIDState
    x_i0::Float64 = 0.0 #previous integrator path state
    x_d0::Float64 = 0.0 #previous derivative path state
    sat_out_0::Int64 = 0 #previous output saturation status
end

@kwdef struct PIDOutput
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

Systems.Y(::PID) = PIDOutput()
Systems.U(::PID) = PIDInput()
Systems.S(::PID) = PIDState()

function Systems.reset!(sys::System{<:PID})
    sys.u.input = 0
    sys.u.sat_ext = 0
    sys.s.x_i0 = 0
    sys.s.x_d0 = 0
    sys.s.sat_out_0 = 0
end

function assign!(sys::System{<:PID}, params::PIDParams{<:Real})
    @unpack k_p, k_i, k_d, τ_f = params
    @pack! sys.u = k_p, k_i, k_d, τ_f
end

function Systems.f_disc!(sys::System{<:PID}, Δt::Real)

    @unpack k_p, k_i, k_d, τ_f, β_p, β_d, bound_lo, bound_hi, input, sat_ext = sys.u
    @unpack x_i0, x_d0, sat_out_0 = sys.s

    α = 1 / (τ_f + Δt)

    u_p = β_p * input
    u_d = β_d * input
    u_i = input

    int_halted = ((sign(u_i * sat_out_0) > 0) || (sign(u_i * sat_ext) > 0))

    x_i = x_i0 + Δt * k_i * u_i * !int_halted
    x_d = α * τ_f * x_d0 + Δt * α * k_d * u_d

    y_p = k_p * u_p
    y_i = x_i
    y_d = α * (-x_d0 + k_d * u_d)
    out_free = y_p + y_i + y_d

    sat_hi = out_free >= bound_hi
    sat_lo = out_free <= bound_lo
    sat_out = sat_hi - sat_lo

    output = clamp(out_free, bound_lo, bound_hi)

    sys.y = PIDOutput(; k_p, k_i, k_d, τ_f, β_p, β_d, bound_lo, bound_hi, input, sat_ext,
                u_p, u_i, u_d, y_p, y_i, y_d, out_free, sat_out, output, int_halted)

    sys.s.x_i0 = x_i
    sys.s.x_d0 = x_d
    sys.s.sat_out_0 = sat_out

    return false

end

############################## Vector Version ##################################

struct PIDVector{N} <: SystemDefinition end

@kwdef struct PIDVectorInput{N}
    k_p::MVector{N,Float64} = ones(N) #proportional gain
    k_i::MVector{N,Float64} = zeros(N) #integral gain
    k_d::MVector{N,Float64} = zeros(N) #derivative gain
    τ_f::MVector{N,Float64} = 0.01 * ones(N) #derivative filter time constant
    β_p::MVector{N,Float64} = ones(N) #proportional path setpoint weighting factor
    β_d::MVector{N,Float64} = ones(N) #derivative path setpoint weighting factor
    bound_lo::MVector{N,Float64} = fill(-Inf, N) #lower output bounds
    bound_hi::MVector{N,Float64} = fill(Inf, N) #higher output bounds
    input::MVector{N,Float64} = zeros(Float64, N) #input
    sat_ext::MVector{N,Int64} = zeros(Int64, N) #external (signed) saturation signal
end

@kwdef struct PIDVectorState{N}
    x_i0::MVector{N,Float64} = zeros(N) #previous integrator path state
    x_d0::MVector{N,Float64} = zeros(N) #previous derivative path state
    sat_out_0::MVector{N,Int64} = zeros(N) #previous output saturation status
end

@kwdef struct PIDVectorOutput{N}
    k_p::SVector{N,Float64} = ones(SVector{N})
    k_i::SVector{N,Float64} = zeros(SVector{N})
    k_d::SVector{N,Float64} = zeros(SVector{N})
    τ_f::SVector{N,Float64} = 0.01ones(SVector{N})
    β_p::SVector{N,Float64} = ones(SVector{N})
    β_d::SVector{N,Float64} = ones(SVector{N})
    bound_lo::SVector{N,Float64} = fill(-Inf, SVector{N}) #lower output bounds
    bound_hi::SVector{N,Float64} = fill(Inf, SVector{N}) #higher output bounds
    input::SVector{N,Float64} = zeros(SVector{N})
    sat_ext::SVector{N,Int64} = zeros(SVector{N, Int64}) #external (signed) saturation signal
    u_p::SVector{N,Float64} = zeros(SVector{N}) #proportional path input
    u_i::SVector{N,Float64} = zeros(SVector{N}) #integral path input
    u_d::SVector{N,Float64} = zeros(SVector{N}) #derivative path input
    y_p::SVector{N,Float64} = zeros(SVector{N}) #proportional path output
    y_i::SVector{N,Float64} = zeros(SVector{N}) #integral path output
    y_d::SVector{N,Float64} = zeros(SVector{N}) #derivative path output
    out_free::SVector{N,Float64} = zeros(SVector{N}) #non-clamped PID output
    sat_out::SVector{N,Int64} = zeros(SVector{N, Int64}) #output saturation status
    output::SVector{N,Float64} = zeros(SVector{N}) #actual PID output
    int_halted::SVector{N,Bool} = zeros(SVector{N, Bool}) #integration halted
end

Systems.Y(::PIDVector{N}) where {N} = PIDVectorOutput{N}()
Systems.U(::PIDVector{N}) where {N} = PIDVectorInput{N}()
Systems.S(::PIDVector{N}) where {N} = PIDVectorState{N}()

function Systems.reset!(sys::System{<:PIDVector{N}}) where {N}
    sys.u.input .= 0
    sys.u.sat_ext .= 0
    sys.s.x_i0 .= 0
    sys.s.x_d0 .= 0
    sys.s.sat_out_0 .= 0
end

function Systems.f_disc!(sys::System{<:PIDVector{N}}, Δt::Real) where {N}

    @unpack s, u = sys

    k_p, k_i, k_d, τ_f, β_p, β_d = map(SVector, (
        u.k_p, u.k_i, u.k_d, u.τ_f, u.β_p, u.β_d))

    input, bound_lo, bound_hi, sat_ext = map(SVector, (
        u.input, u.bound_lo, u.bound_hi, u.sat_ext))

    x_i0, x_d0, sat_out_0 = map(SVector, (
        s.x_i0, s.x_d0, s.sat_out_0))

    α = 1 ./ (τ_f .+ Δt)

    u_p = β_p .* input
    u_d = β_d .* input
    u_i = input

    int_halted = ((sign.(u_i .* sat_out_0) .> 0) .|| (sign.(u_i .* sat_ext) .> 0))

    x_i = x_i0 + Δt * k_i .* u_i .* .!int_halted
    x_d = α .* τ_f .* x_d0 + Δt * α .* k_d .* u_d

    y_p = k_p .* u_p
    y_i = x_i
    y_d = α .* (-x_d0 + k_d .* u_d)
    out_free = y_p + y_i + y_d

    sat_hi = out_free .>= bound_hi
    sat_lo = out_free .<= bound_lo
    sat_out = sat_hi - sat_lo

    output = clamp.(out_free, bound_lo, bound_hi)

    s.x_i0 .= x_i
    s.x_d0 .= x_d
    s.sat_out_0 .= sat_out

    sys.y = PIDVectorOutput(; k_p, k_i, k_d, τ_f, β_p, β_d, bound_lo, bound_hi, input, sat_ext,
                u_p, u_i, u_d, y_p, y_i, y_d, out_free, sat_out, output, int_halted)

    return false

end


################################### Plots ######################################

function Plotting.make_plots(ts::Union{TimeSeries{<:PIDOutput},
                                       TimeSeries{<:PIDVectorOutput}}; kwargs...)

    pd = OrderedDict{Symbol, Plots.Plot}()

    input = plot(ts.input; title = "Input", ylabel = L"$e$", kwargs...)
    output = plot(ts.output; title = "Output", ylabel = L"$y$", kwargs...)

    pd[:sf] = plot(input, output;
        plot_title = "Input & Output",
        layout = (1,2),
        link = :y,
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    k_p = plot(ts.k_p; title = "Proportional Gain", ylabel = L"$k_p$", kwargs...)
    k_i = plot(ts.k_i; title = "Integral Gain", ylabel = L"$k_i$", kwargs...)
    k_d = plot(ts.k_d; title = "Derivative Gain", ylabel = L"$k_d$", kwargs...)
    τ_f = plot(ts.τ_f; title = "Derivative Filter Time Constant", ylabel = L"$\tau_f$", kwargs...)

    pd[:p1] = plot(k_p, k_i, k_d, τ_f;
        plot_title = "Parameters",
        layout = (4,1),
        link = :none,
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    β_p = plot(ts.β_p; title = "Proportional Input Weighting", ylabel = L"$\beta_p$", kwargs...)
    β_d = plot(ts.β_d; title = "Derivative Input Weighting", ylabel = L"$\beta_d$", kwargs...)
    bound_lo = plot(ts.bound_lo; title = "Lower Output Bound", ylabel = L"$y_{min}$", kwargs...)
    bound_hi = plot(ts.bound_hi; title = "Upper Output Bound", ylabel = L"$y_{max}$", kwargs...)

    pd[:p2] = plot(β_p, β_d, bound_lo, bound_hi;
        plot_title = "Parameters",
        layout = (4,1),
        link = :none,
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    sat_ext = plot(ts.sat_ext; title = "External Saturation Input", ylabel = "", kwargs...)
    sat_out = plot(ts.sat_out; title = "Output Saturation", ylabel = "", kwargs...)
    int_halted = plot(ts.int_halted; title = "Integrator Halted", ylabel = "", kwargs...)

    pd[:awu] = plot(sat_ext, sat_out, int_halted;
        plot_title = "Anti-Windup",
        layout = (1,3),
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    u_p = plot(ts.u_p; title = "Input", ylabel = L"$u_p$", kwargs...)
    y_p = plot(ts.y_p; title = "Output", ylabel = L"$y_p$", kwargs...)

    pd[:prop] = plot(u_p, y_p;
        plot_title = "Proportional Path",
        layout = (1,2),
        link = :y,
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    u_i = plot(ts.u_i; title = "Input", ylabel = L"$u_i$", kwargs...)
    y_i = plot(ts.y_i; title = "Output", ylabel = L"$y_i$", kwargs...)

    pd[:int] = plot(u_i, y_i, int_halted;
        plot_title = "Integral Path",
        layout = (1,3),
        link = :y,
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    u_d = plot(ts.u_d; title = "Input", ylabel = L"$u_d$", kwargs...)
    y_d = plot(ts.y_d; title = "Output", ylabel = L"$y_d$", kwargs...)

    pd[:der] = plot(u_d, y_d;
        plot_title = "Derivative Path",
        layout = (1,2),
        link = :y,
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    out_free = plot(ts.out_free; title = "Free", ylabel = L"$y_{free}$", kwargs...)
    output = plot(ts.output; title = "Actual", ylabel = L"$y$", kwargs...)

    pd[:output] = plot(out_free, output, sat_out;
        plot_title = "PID Output",
        layout = (1,3),
        link = :y,
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    return pd

end

#################################### GUI #######################################


function GUI.draw(sys::Union{System{<:PID}, System{<:PIDVector}},
                            label::String = "Discrete PID")

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


############################## LQRTracker ######################################
################################################################################

struct LQRTracker{NX, NU, NZ, NUX, NUZ} <: SystemDefinition end

function LQRTracker{NX, NU, NZ}() where {NX, NU, NZ}
    @assert NZ <= NU "Can't have more command variables than control inputs"
    NUX = NU * NX
    NUZ = NU * NZ
    LQRTracker{NX, NU, NZ, NUX, NUZ}()
end

@kwdef struct LQRTrackerInput{NX, NU, NZ, NUX, NUZ}
    C_fbk::MMatrix{NU, NX, Float64, NUX} = zeros(NU, NX) #state feedback matrix
    C_fwd::MMatrix{NU, NZ, Float64, NUZ} = zeros(NU, NZ) #feedforward matrix
    C_int::MMatrix{NU, NZ, Float64, NUZ} = zeros(NU, NZ) #integrator gain matrix
    x_trim::MVector{NX, Float64} = zeros(NX) #trim point state
    u_trim::MVector{NU, Float64} = zeros(NU) #trim point control input
    z_trim::MVector{NZ, Float64} = zeros(NZ) #trim point command vector
    bound_lo::MVector{NU,Float64} = fill(-Inf, NU) #lower output bounds
    bound_hi::MVector{NU,Float64} = fill(Inf, NU) #upper output bounds
    sat_ext::MVector{NU,Int64} = zeros(Int64, NU) #saturation input signal
    z_sp::MVector{NZ, Float64} = zeros(NZ) #command vector set point
    z::MVector{NZ, Float64} = zeros(NZ) #current command vector value
    x::MVector{NX, Float64} = zeros(NX) #current state vector value
end

function LQRTrackerInput{NX, NU, NZ}(args...; kwargs...) where {NX, NU, NZ}
    NUX = NU * NX
    NUZ = NU * NZ
    LQRTrackerInput{NX, NU, NZ, NUX, NUZ}(args...; kwargs...)
end

@kwdef struct LQRTrackerOutput{NX, NU, NZ, NUX, NUZ}
    C_fbk::SMatrix{NU, NX, Float64, NUX} = zeros(SMatrix{NU, NX}) #state feedback matrix
    C_fwd::SMatrix{NU, NZ, Float64, NUZ} = zeros(SMatrix{NU, NZ}) #feedforward matrix
    C_int::SMatrix{NU, NZ, Float64, NUZ} = zeros(SMatrix{NU, NZ}) #integrator gain matrix
    x_trim::SVector{NX, Float64} = zeros(SVector{NX}) #trim point state
    u_trim::SVector{NU, Float64} = zeros(SVector{NU}) #trim point control input
    z_trim::SVector{NZ, Float64} = zeros(SVector{NZ}) #trim point command variable
    bound_lo::SVector{NU,Float64} = fill(-Inf, SVector{NU}) #lower output bounds
    bound_hi::SVector{NU,Float64} = fill(Inf, SVector{NU}) #upper output bounds
    sat_ext::SVector{NU,Int64} = zeros(SVector{NU, Int64}) #saturation input signal
    z_sp::SVector{NZ, Float64} = zeros(SVector{NZ}) #command variable set point
    z::SVector{NZ, Float64} = zeros(SVector{NZ}) #current command vector value
    x::SVector{NX, Float64} = zeros(SVector{NX}) #current state vector value
    int_in::SVector{NU,Float64} = zeros(SVector{NU}) #integrator input
    int_halted::SVector{NU,Bool} = zeros(SVector{NU, Bool}) #integration halted
    int_out::SVector{NU,Float64} = zeros(SVector{NU}) #integrator output
    out_free::SVector{NU,Float64} = zeros(SVector{NU}) #total output, free
    out_sat::SVector{NU,Int64} = zeros(SVector{NU, Int64}) #current output saturation status
    output::SVector{NU,Float64} = zeros(SVector{NU}) #actual output
end

function LQRTrackerOutput{NX, NU, NZ}(args...; kwargs...) where {NX, NU, NZ}
    NUX = NU * NX
    NUZ = NU * NZ
    LQRTrackerOutput{NX, NU, NZ, NUX, NUZ}(args...; kwargs...)
end

@kwdef struct LQRTrackerState{NX, NU}
    int_out_0::MVector{NU,Float64} = zeros(NU) #previous integrator path state
    out_sat_0::MVector{NU,Int64} = zeros(NU) #previous output saturation status
end

function Systems.Y(::LQRTracker{NX, NU, NZ, NUX, NUZ}) where {NX, NU, NZ, NUX, NUZ}
    LQRTrackerOutput{NX, NU, NZ, NUX, NUZ}()
end

function Systems.U(::LQRTracker{NX, NU, NZ, NUX, NUZ}) where {NX, NU, NZ, NUX, NUZ}
    LQRTrackerInput{NX, NU, NZ, NUX, NUZ}()
end

function Systems.S(::LQRTracker{NX, NU, NZ, NUX, NUZ}) where {NX, NU, NZ, NUX, NUZ}
    LQRTrackerState{NX, NU}()
end

function Systems.reset!(sys::System{<:LQRTracker})
    sys.u.z_sp .= 0
    sys.u.z .= 0
    sys.u.x .= 0
    sys.u.sat_ext .= 0
    sys.s.int_out_0 .= 0
    sys.s.out_sat_0 .= 0
end

function Systems.f_disc!(sys::System{<:LQRTracker}, Δt::Real)

    @unpack s, u = sys

    C_fbk, C_fwd, C_int = map(SMatrix, (u.C_fbk, u.C_fwd, u.C_int))
    x_trim, u_trim, z_trim = map(SVector, (u.x_trim, u.u_trim, u.z_trim))
    bound_lo, bound_hi, sat_ext = map(SVector, (u.bound_lo, u.bound_hi, u.sat_ext))
    z_sp, z, x = map(SVector, (u.z_sp, u.z, u.x))

    int_out_0 = SVector(s.int_out_0)
    out_sat_0 = SVector(s.out_sat_0)

    int_in = C_int * (z_sp - z)
    int_halted = ((sign.(int_in .* out_sat_0) .> 0) .|| (sign.(int_in .* sat_ext) .> 0))
    int_out = int_out_0 + Δt * int_in .* .!int_halted

    out_free = u_trim + int_out + C_fwd * (z_sp - z_trim) - C_fbk * (x - x_trim)

    out_sat = (out_free .>= bound_hi) - (out_free .<= bound_lo)
    output = clamp.(out_free, bound_lo, bound_hi)

    s.int_out_0 .= int_out
    s.out_sat_0 .= out_sat

    sys.y = LQRTrackerOutput(; C_fbk, C_fwd, C_int, x_trim, u_trim, z_trim,
        bound_lo, bound_hi, sat_ext, z_sp, z, x,
        int_in, int_out, int_halted, out_free, out_sat, output)

    return false

end

function GUI.draw(sys::System{<:LQRTracker})

    @unpack C_fbk, C_fwd, C_int, x_trim, u_trim, z_trim, bound_lo, bound_hi,
            sat_ext, z_sp, z, x, int_in, int_halted, int_out,
            out_free, out_sat, output = sys.y

        CImGui.Text("Feedback Gain = $C_fbk")
        CImGui.Text("Forward Gain = $C_fwd")
        CImGui.Text("Integral Gain = $C_int")
        CImGui.Text("Trim State Vector = $x_trim")
        CImGui.Text("Trim Control Vector = $u_trim")
        CImGui.Text("Trim Command Vector = $z_trim")
        CImGui.Text("Lower Output Bound = $bound_lo")
        CImGui.Text("Upper Output Bound = $bound_hi")
        CImGui.Text("External Saturation Input = $sat_ext")
        CImGui.Text("Set Point Command Vector = $z_sp")
        CImGui.Text("Current Command Vector = $z")
        CImGui.Text("Current State Vector = $x")
        CImGui.Text("Integrator Input = $int_in")
        CImGui.Text("Integrator Halted = $int_halted")
        CImGui.Text("Integrator Output = $int_out")
        CImGui.Text("Free Output = $out_free")
        CImGui.Text("Output Saturation = $out_sat")
        CImGui.Text("Actual Output = $output")

end #function

@kwdef struct LQRTrackerParams{CB, CF, CI, X, U, Z}
    C_fbk::CB
    C_fwd::CF
    C_int::CI
    x_trim::X
    u_trim::U
    z_trim::Z
end

const LQRTrackerPoint = LQRTrackerParams{CB, CF, CI, X, U, Z} where {
    CB <: AbstractMatrix,
    CF <: AbstractMatrix,
    CI <: AbstractMatrix,
    X <: AbstractVector,
    U <: AbstractVector,
    Z <: AbstractVector}


function assign!(lqr::System{<:LQRTracker}, params::LQRTrackerPoint)
    @unpack C_fbk, C_fwd, C_int, x_trim, u_trim, z_trim = params
    lqr.u.C_fbk .= C_fbk
    lqr.u.C_fwd .= C_fwd
    lqr.u.C_int .= C_int
    lqr.u.x_trim .= x_trim
    lqr.u.u_trim .= u_trim
    lqr.u.z_trim .= z_trim
end

# const LQRTrackerParamsStatic{NX, NU, NZ, NUX, NUZ} = LQRTrackerParams{
#     SMatrix{NU, NX, Float64, NUX},
#     SMatrix{NU, NZ, Float64, NUZ},
#     SMatrix{NU, NZ, Float64, NUZ},
#     SVector{NX, Float64},
#     SVector{NU, Float64},
#     SVector{NZ, Float64}}

end #submodule


################################################################################
############################ PID Optimization ##################################
################################################################################

module PIDOpt

using StaticArrays, UnPack, NLopt, ControlSystems
using Trapz: trapz
using ..Control.Discrete: PIDParams

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

end #submodule

################################################################################
################################################################################
################################################################################


using Reexport

using .Continuous
using .Discrete
using .PIDOpt

end #module