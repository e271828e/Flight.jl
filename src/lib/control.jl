module Control

############################## Continuous Models ###############################
################################################################################
module Continuous ##############################################################

using ComponentArrays, StaticArrays, LinearAlgebra
using ControlSystems: ControlSystemsBase, ControlSystems, ss
using Plots, LaTeXStrings, DataStructures

using Flight.FlightCore

####################### Proportional-Integral Compensator ######################
################################################################################

@kwdef struct PIVector{N} <: ModelDefinition
    k_p::MVector{N,Float64} = ones(N) #proportional gain
    k_i::MVector{N,Float64} = zeros(N) #integral gain
    k_l::MVector{N,Float64} = zeros(N) #integrator leak factor
    β_p::MVector{N,Float64} = ones(N) #proportional path reference weighting factor
    bound_lo::MVector{N,Float64} = fill(-Inf, N) #lower output bounds
    bound_hi::MVector{N,Float64} = fill(Inf, N) #higher output bounds
end

@kwdef struct PIVectorU{N}
    input::MVector{N,Float64} = zeros(N) #input (reference - feedback)
    sat_ext::MVector{N,Int64} = zeros(Int64, N) #external (signed) saturation signal
end

@kwdef struct PIVectorY{N}
    k_p::SVector{N,Float64} = ones(SVector{N})
    k_i::SVector{N,Float64} = zeros(SVector{N})
    k_l::SVector{N,Float64} = zeros(SVector{N})
    β_p::SVector{N,Float64} = ones(SVector{N})
    bound_lo::SVector{N,Float64} = fill(-Inf, SVector{N}) #lower output bounds
    bound_hi::SVector{N,Float64} = fill(Inf, SVector{N}) #higher output bounds
    input::SVector{N,Float64} = zeros(SVector{N})
    sat_ext::SVector{N,Int64} = zeros(SVector{N, Int64}) #external (signed) saturation signal
    u_p::SVector{N,Float64} = zeros(SVector{N}) #proportional path input
    u_i::SVector{N,Float64} = zeros(SVector{N}) #integral path input
    x_i::SVector{N,Float64} = zeros(SVector{N}) #integrator state
    y_p::SVector{N,Float64} = zeros(SVector{N}) #proportional path output
    y_i::SVector{N,Float64} = zeros(SVector{N}) #integral path output
    out_free::SVector{N,Float64} = zeros(SVector{N}) #total output, free
    sat_out::SVector{N,Int64} = zeros(SVector{N, Int64}) #current output saturation status
    output::SVector{N,Float64} = zeros(SVector{N}) #actual output
    int_halted::SVector{N,Bool} = zeros(SVector{N, Bool}) #integration halted
end

Modeling.X(::PIVector{N}) where {N} = zeros(N)
Modeling.Y(::PIVector{N}) where {N} = PIVectorY{N}()
Modeling.U(::PIVector{N}) where {N} = PIVectorU{N}()

@no_periodic PIVector
@no_step PIVector

function Modeling.f_ode!(mdl::Model{<:PIVector{N}}) where {N}

    (; ẋ, x, u, parameters) = mdl
    (; input, sat_ext) = u
    (; k_p, k_i, k_l, β_p, bound_lo, bound_hi) = parameters

    input, sat_ext, k_p, k_i, k_l, β_p, bound_lo, bound_hi = map(SVector, (
    input, sat_ext, k_p, k_i, k_l, β_p, bound_lo, bound_hi))

    x_i = SVector{N, Float64}(x)

    u_p = β_p .* input
    u_i = input

    y_p = k_p .* u_p
    y_i = x_i
    out_free = y_p + y_i #raw output
    output = clamp.(out_free, bound_lo, bound_hi) #clamped output

    sat_hi = out_free .>= bound_hi
    sat_lo = out_free .<= bound_lo
    sat_out = sat_hi - sat_lo
    int_halted = ((sign.(u_i .* sat_out) .> 0) .|| (sign.(u_i .* sat_ext) .> 0))

    ẋ .= k_i .* u_i .* .!int_halted - k_l .* x_i

    mdl.y = PIVectorY(; k_p, k_i, k_l, β_p, bound_lo, bound_hi, input, sat_ext,
                     u_p, u_i, x_i, y_p, y_i, out_free, sat_out, output, int_halted)

end

function Modeling.init!(mdl::Model{<:PIVector})
    mdl.u.input .= 0
    mdl.u.sat_ext .= 0
    mdl.x .= 0
end



# ############################## Plotting ########################################

function Plotting.make_plots(ts::TimeSeries{<:PIVectorY}; kwargs...)

    pd = OrderedDict{Symbol, Plots.Plot}()

    input = plot(ts.input; title = "Input", ylabel = L"$e$", kwargs...)
    output = plot(ts.output; title = "Output", ylabel = L"$y$", kwargs...)

    pd[:sf] = plot(input, output;
        plot_title = "Input & Output",
        layout = (2,1),
        link = :y,
        kwargs...)

    k_p = plot(ts.k_p; title = "Proportional Gain", ylabel = L"$k_p$", kwargs...)
    k_i = plot(ts.k_i; title = "Integral Gain", ylabel = L"$k_i$", kwargs...)
    k_l = plot(ts.k_l; title = "Leak Factor", ylabel = L"$k_d$", kwargs...)

    pd[:p1] = plot(k_p, k_i, k_l;
        plot_title = "Parameters",
        layout = (3,1),
        link = :none,
        kwargs...)

    β_p = plot(ts.β_p; title = "Proportional Input Weighting", ylabel = L"$\beta_p$", kwargs...)
    bound_lo = plot(ts.bound_lo; title = "Lower Output Bound", ylabel = L"$y_{min}$", kwargs...)
    bound_hi = plot(ts.bound_hi; title = "Upper Output Bound", ylabel = L"$y_{max}$", kwargs...)

    pd[:p2] = plot(β_p, bound_lo, bound_hi;
        plot_title = "Parameters",
        layout = (3,1),
        link = :none,
        kwargs...)

    sat_ext = plot(ts.sat_ext; title = "External Saturation Input", ylabel = "", kwargs...)
    sat_out = plot(ts.sat_out; title = "Output Saturation", ylabel = "", kwargs...)
    int_halted = plot(ts.int_halted; title = "Integrator Halted", ylabel = "", kwargs...)

    pd[:awu] = plot(sat_ext, sat_out, int_halted;
        plot_title = "Anti-Windup",
        layout = (3,1),
        kwargs...)

    u_p = plot(ts.u_p; title = "Input", ylabel = L"$u_p$", kwargs...)
    y_p = plot(ts.y_p; title = "Output", ylabel = L"$y_p$", kwargs...)

    pd[:prop] = plot(u_p, y_p;
        plot_title = "Proportional Path",
        layout = (2,1),
        link = :y,
        kwargs...)

    u_i = plot(ts.u_i; title = "Input", ylabel = L"$u_i$", kwargs...)
    y_i = plot(ts.y_i; title = "Output", ylabel = L"$y_i$", kwargs...)

    pd[:int] = plot(u_i, y_i, int_halted;
        plot_title = "Integral Path",
        layout = (3,1),
        link = :y,
        kwargs...)

    out_free = plot(ts.out_free; title = "Free", ylabel = L"$y_{free}$", kwargs...)
    output = plot(ts.output; title = "Actual", ylabel = L"$y$", kwargs...)

    pd[:output] = plot(out_free, output, sat_out;
        plot_title = "PID Output",
        layout = (3,1),
        link = :y,
        kwargs...)

    return pd

end

#################################### GUI #######################################


function GUI.draw(mdl::Model{<:PIVector{N}}, label::String = "PIVector{$N}") where {N}

    (; k_p, k_i, k_l, β_p, bound_lo, bound_hi, input, sat_ext,
        u_p, u_i, x_i, y_p, y_i, out_free, sat_out, output, int_halted) = mdl.y

    # CImGui.Begin(label)

        CImGui.Text("Proportional Gain = $k_p")
        CImGui.Text("Integral Gain = $k_i")
        CImGui.Text("Integrator Leak Factor = $k_l")
        CImGui.Text("Proportional Input Weighting = $β_p")
        CImGui.Text("Lower Output Bound = $bound_lo")
        CImGui.Text("Upper Output Bound = $bound_hi")
        CImGui.Text("Input = $input")
        CImGui.Text("External Saturation Input = $sat_ext")
        CImGui.Text("Proportional Path Input = $u_p")
        CImGui.Text("Proportional Path Output = $y_p")
        CImGui.Text("Integrator State = $x_i")
        CImGui.Text("Integral Path Input = $u_i")
        CImGui.Text("Integral Path Output = $y_i")
        CImGui.Text("Free Output = $out_free")
        CImGui.Text("Output Saturation = $sat_out")
        CImGui.Text("Actual Output = $output")
        CImGui.Text("Integrator Halted = $int_halted")

    # CImGui.End()

end #function

end #submodule


################################ Discrete Models ###############################
################################################################################
module Discrete ###############################################################

using StaticArrays, LinearAlgebra
using StructArrays, Interpolations, HDF5 #for lookups
using RobustAndOptimalControl
using Plots, LaTeXStrings, DataStructures

using Flight.FlightCore

############################# Integrator ###############################
################################################################################

@kwdef struct Integrator <: ModelDefinition
    bound_lo::Ref{Float64} = Ref(-Inf) #lower output bound
    bound_hi::Ref{Float64} = Ref(Inf) #higher output bound
end

@kwdef mutable struct IntegratorInput
    input::Float64 = 0 #current input
    sat_ext::Int64 = 0 #external (signed) saturation signal
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

Modeling.Y(::Integrator) = IntegratorOutput()
Modeling.U(::Integrator) = IntegratorInput()
Modeling.S(::Integrator) = IntegratorState()

function Modeling.init!(mdl::Model{<:Integrator}, x0::Real = 0.0)
    mdl.u.input = 0
    mdl.u.sat_ext = 0
    mdl.s.x0 = x0
    mdl.s.sat_out_0 = 0
end

@no_step Integrator
@no_ode Integrator

function Modeling.f_periodic!(::NoScheduling, mdl::Model{<:Integrator})

    (; s, u, Δt, parameters) = mdl
    (; input, sat_ext) = u
    (; x0, sat_out_0) = s
    bound_lo = parameters.bound_lo[]
    bound_hi = parameters.bound_hi[]

    halted = ((sign(input * sat_out_0) > 0) || (sign(input * sat_ext) > 0))
    x1 = x0 + Δt * input * !halted
    output = clamp(x1, bound_lo, bound_hi)

    sat_hi = x1 >= bound_hi
    sat_lo = x1 <= bound_lo
    sat_out = sat_hi - sat_lo

    s.x0 = x1
    s.sat_out_0 = sat_out

    mdl.y = IntegratorOutput(; input, sat_ext, bound_lo, bound_hi, x1, output, sat_out, halted)

end

################################################################################

@kwdef struct IntegratorVector{N} <: ModelDefinition
    bound_lo::MVector{N,Float64} = fill(-Inf, N) #lower output bounds
    bound_hi::MVector{N,Float64} = fill(Inf, N) #higher output bounds
end

@kwdef mutable struct IntegratorVectorInput{N}
    input::MVector{N,Float64} = zeros(Float64, N)
    sat_ext::MVector{N,Int64} = zeros(Int64, N) #external (signed) saturation signal
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

Modeling.Y(::IntegratorVector{N}) where {N} = IntegratorVectorOutput{N}()
Modeling.U(::IntegratorVector{N}) where {N} = IntegratorVectorInput{N}()
Modeling.S(::IntegratorVector{N}) where {N} = IntegratorVectorState{N}()

function Modeling.init!(mdl::Model{<:IntegratorVector{N}}, x0::AbstractVector{<:Real} = zeros(SVector{N})) where {N}
    mdl.u.input .= 0
    mdl.u.sat_ext .= 0
    mdl.s.x0 .= x0
    mdl.s.sat_out_0 .= 0
end

@no_ode IntegratorVector
@no_step IntegratorVector

function Modeling.f_periodic!(::NoScheduling, mdl::Model{<:IntegratorVector})

    (; s, u, Δt, parameters) = mdl
    (; input, sat_ext) = u
    (; bound_lo, bound_hi) = parameters

    input, sat_ext, bound_lo, bound_hi = map(SVector, (
    input, sat_ext, bound_lo, bound_hi))

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

    mdl.y = IntegratorVectorOutput(; input, sat_ext, bound_lo, bound_hi, x1, output, sat_out, halted)

end


#################################### GUI #######################################


function GUI.draw(mdl::Union{Model{<:Integrator}, Model{<:IntegratorVector}},
                    label::String = "Integrator")

    (; x0, sat_out_0) = mdl.s
    (; input, sat_ext, bound_lo, bound_hi, x1, output, sat_out, halted) = mdl.y

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

@kwdef struct LeadLag <: ModelDefinition
    z::Ref{Float64} = Ref(-1.0) #zero location (z < 0)
    p::Ref{Float64} = Ref(-10.0) #pole location (p < 0)
    k::Ref{Float64} = Ref(1.0) #gain
end

@kwdef mutable struct LeadLagInput
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
    u1::Float64 = 0.0 #current input
    y1::Float64 = 0.0 #current output
end

Modeling.U(::LeadLag) = LeadLagInput()
Modeling.S(::LeadLag) = LeadLagState()
Modeling.Y(::LeadLag) = LeadLagOutput()

function Modeling.init!(mdl::Model{<:LeadLag})
    mdl.u.u1 = 0
    mdl.s.u0 = 0
    mdl.s.x0 = 0
end

@no_ode LeadLag
@no_step LeadLag

function Modeling.f_periodic!(::NoScheduling, mdl::Model{<:LeadLag})

    (; parameters, s, u, Δt) = mdl
    (; u1) = u
    (; u0, x0) = s
    z = parameters.z[]
    p = parameters.p[]
    k = parameters.k[]

    a0 = (2+p*Δt)/(2-p*Δt)
    b1 = (2-z*Δt)/(2-p*Δt)
    b0 = (-2-z*Δt)/(2-p*Δt)

    x1 = a0 * x0 + b1 * u1 + b0 * u0
    y1 = k * x1

    mdl.y = LeadLagOutput(; z, p, k, u1, y1)

    s.x0 = x1
    s.u0 = u1

end


#################################### GUI #######################################

function GUI.draw(mdl::Model{<:LeadLag}, label::String = "Discrete Lead Compensator")

    (; z, p, k, u1, y1) = mdl.y

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


################################# Scalar Version ###############################

@kwdef struct PID <: ModelDefinition
    k_p::Ref{Float64} = Ref(1.0) #proportional gain
    k_i::Ref{Float64} = Ref(0.0) #integral gain
    k_d::Ref{Float64} = Ref(0.0) #derivative gain
    τ_f::Ref{Float64} = Ref(0.01) #derivative filter time constant
    β_p::Ref{Float64} = Ref(1.0) #proportional path reference weighting factor
    β_d::Ref{Float64} = Ref(1.0) #derivative path reference weighting factor
    bound_lo::Ref{Float64} = Ref(-Inf) #lower output bound
    bound_hi::Ref{Float64} = Ref(Inf) #higher output bound
end

@kwdef mutable struct PIDInput
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
    k_i::Float64 = 0.0 #integral gain
    k_d::Float64 = 0.0 #derivative gain
    τ_f::Float64 = 0.01 #derivative filter time constant
    β_p::Float64 = 1.0 #proportional path reference weighting factor
    β_d::Float64 = 1.0 #derivative path reference weighting factor
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

Modeling.Y(::PID) = PIDOutput()
Modeling.U(::PID) = PIDInput()
Modeling.S(::PID) = PIDState()

function Modeling.init!(mdl::Model{<:PID})
    mdl.u.input = 0
    mdl.u.sat_ext = 0
    mdl.s.x_i0 = 0
    mdl.s.x_d0 = 0
    mdl.s.sat_out_0 = 0
end

@no_ode PID
@no_step PID

function Modeling.f_periodic!(::NoScheduling, mdl::Model{<:PID})

    (; parameters, u, s, Δt) = mdl
    (; k_p, k_i, k_d, τ_f, β_p, β_d, bound_lo, bound_hi) = parameters
    (; input, sat_ext) = u
    (; x_i0, x_d0, sat_out_0) = s

    #extract parameter values from Refs
    k_p, k_i, k_d, τ_f, β_p, β_d, bound_lo, bound_hi = map(r->getproperty(r, :x), (
    k_p, k_i, k_d, τ_f, β_p, β_d, bound_lo, bound_hi))

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

    mdl.y = PIDOutput(; k_p, k_i, k_d, τ_f, β_p, β_d, bound_lo, bound_hi, input, sat_ext,
                u_p, u_i, u_d, y_p, y_i, y_d, out_free, sat_out, output, int_halted)

    mdl.s.x_i0 = x_i
    mdl.s.x_d0 = x_d
    mdl.s.sat_out_0 = sat_out

end

############################## Vector Version ##################################

@kwdef struct PIDVector{N} <: ModelDefinition
    k_p::MVector{N,Float64} = ones(N) #proportional gain
    k_i::MVector{N,Float64} = zeros(N) #integral gain
    k_d::MVector{N,Float64} = zeros(N) #derivative gain
    τ_f::MVector{N,Float64} = 0.01 * ones(N) #derivative filter time constant
    β_p::MVector{N,Float64} = ones(N) #proportional path reference weighting factor
    β_d::MVector{N,Float64} = ones(N) #derivative path reference weighting factor
    bound_lo::MVector{N,Float64} = fill(-Inf, N) #lower output bounds
    bound_hi::MVector{N,Float64} = fill(Inf, N) #higher output bounds
end

@kwdef struct PIDVectorInput{N}
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

Modeling.Y(::PIDVector{N}) where {N} = PIDVectorOutput{N}()
Modeling.U(::PIDVector{N}) where {N} = PIDVectorInput{N}()
Modeling.S(::PIDVector{N}) where {N} = PIDVectorState{N}()

function Modeling.init!(mdl::Model{<:PIDVector{N}}) where {N}
    mdl.u.input .= 0
    mdl.u.sat_ext .= 0
    mdl.s.x_i0 .= 0
    mdl.s.x_d0 .= 0
    mdl.s.sat_out_0 .= 0
end

@no_ode PIDVector
@no_step PIDVector


function Modeling.f_periodic!(::NoScheduling, mdl::Model{<:PIDVector{N}}) where {N}

    (; parameters, s, u, Δt) = mdl
    (; k_p, k_i, k_d, τ_f, β_p, β_d, bound_lo, bound_hi) = parameters
    (; x_i0, x_d0, sat_out_0) = s
    (; input, sat_ext) = u

    k_p, k_i, k_d, τ_f, β_p, β_d, bound_lo, bound_hi = map(SVector, (
    k_p, k_i, k_d, τ_f, β_p, β_d, bound_lo, bound_hi))

    input, sat_ext, x_i0, x_d0, sat_out_0 = map(SVector, (
    input, sat_ext, x_i0, x_d0, sat_out_0))

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

    mdl.y = PIDVectorOutput(; k_p, k_i, k_d, τ_f, β_p, β_d, bound_lo, bound_hi, input, sat_ext,
                u_p, u_i, u_d, y_p, y_i, y_d, out_free, sat_out, output, int_halted)

end


################################### Plots ######################################

function Plotting.make_plots(ts::Union{TimeSeries{<:PIDOutput},
                                       TimeSeries{<:PIDVectorOutput}}; kwargs...)

    pd = OrderedDict{Symbol, Plots.Plot}()

    input = plot(ts.input; title = "Input", ylabel = L"$e$", kwargs...)
    output = plot(ts.output; title = "Output", ylabel = L"$y$", kwargs...)

    pd[:sf] = plot(input, output;
        plot_title = "Input & Output",
        layout = (2,1),
        link = :y,
        kwargs...)

    k_p = plot(ts.k_p; title = "Proportional Gain", ylabel = L"$k_p$", kwargs...)
    k_i = plot(ts.k_i; title = "Integral Gain", ylabel = L"$k_i$", kwargs...)
    k_d = plot(ts.k_d; title = "Derivative Gain", ylabel = L"$k_d$", kwargs...)
    τ_f = plot(ts.τ_f; title = "Derivative Filter Time Constant", ylabel = L"$\tau_f$", kwargs...)

    pd[:p1] = plot(k_p, k_i, k_d, τ_f;
        plot_title = "Parameters",
        layout = (4,1),
        link = :none,
        kwargs...)

    β_p = plot(ts.β_p; title = "Proportional Input Weighting", ylabel = L"$\beta_p$", kwargs...)
    β_d = plot(ts.β_d; title = "Derivative Input Weighting", ylabel = L"$\beta_d$", kwargs...)
    bound_lo = plot(ts.bound_lo; title = "Lower Output Bound", ylabel = L"$y_{min}$", kwargs...)
    bound_hi = plot(ts.bound_hi; title = "Upper Output Bound", ylabel = L"$y_{max}$", kwargs...)

    pd[:p2] = plot(β_p, β_d, bound_lo, bound_hi;
        plot_title = "Parameters",
        layout = (4,1),
        link = :none,
        kwargs...)

    sat_ext = plot(ts.sat_ext; title = "External Saturation Input", ylabel = "", kwargs...)
    sat_out = plot(ts.sat_out; title = "Output Saturation", ylabel = "", kwargs...)
    int_halted = plot(ts.int_halted; title = "Integrator Halted", ylabel = "", kwargs...)

    pd[:awu] = plot(sat_ext, sat_out, int_halted;
        plot_title = "Anti-Windup",
        layout = (3,1),
        kwargs...)

    u_p = plot(ts.u_p; title = "Input", ylabel = L"$u_p$", kwargs...)
    y_p = plot(ts.y_p; title = "Output", ylabel = L"$y_p$", kwargs...)

    pd[:prop] = plot(u_p, y_p;
        plot_title = "Proportional Path",
        layout = (2,1),
        link = :y,
        kwargs...)

    u_i = plot(ts.u_i; title = "Input", ylabel = L"$u_i$", kwargs...)
    y_i = plot(ts.y_i; title = "Output", ylabel = L"$y_i$", kwargs...)

    pd[:int] = plot(u_i, y_i, int_halted;
        plot_title = "Integral Path",
        layout = (3,1),
        link = :y,
        kwargs...)

    u_d = plot(ts.u_d; title = "Input", ylabel = L"$u_d$", kwargs...)
    y_d = plot(ts.y_d; title = "Output", ylabel = L"$y_d$", kwargs...)

    pd[:der] = plot(u_d, y_d;
        plot_title = "Derivative Path",
        layout = (3,1),
        link = :y,
        kwargs...)

    out_free = plot(ts.out_free; title = "Free", ylabel = L"$y_{free}$", kwargs...)
    output = plot(ts.output; title = "Actual", ylabel = L"$y$", kwargs...)

    pd[:output] = plot(out_free, output, sat_out;
        plot_title = "PID Output",
        layout = (3,1),
        link = :y,
        kwargs...)

    return pd

end

#################################### GUI #######################################


function GUI.draw(mdl::Union{Model{<:PID}, Model{<:PIDVector}},
                            label::String = "Discrete PID")

    (; k_p, k_i, k_d, τ_f, β_p, β_d, bound_lo, bound_hi, input, sat_ext,
        u_p, u_i, u_d, y_p, y_i, y_d, out_free, sat_out, output, int_halted) = mdl.y

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


GUI.draw!(mdl::Model{<:PID}, label::String = "Discrete PID") = GUI.draw(mdl, label)


############################## LQR ######################################
################################################################################

@kwdef struct LQR{NX, NU, NZ, NUX, NUZ} <: ModelDefinition
    K_fbk::MMatrix{NU, NX, Float64, NUX} = zeros(NU, NX) #state feedback matrix
    K_fwd::MMatrix{NU, NZ, Float64, NUZ} = zeros(NU, NZ) #feedforward matrix
    K_int::MMatrix{NU, NZ, Float64, NUZ} = zeros(NU, NZ) #integrator gain matrix
    x_trim::MVector{NX, Float64} = zeros(NX) #trim point state
    u_trim::MVector{NU, Float64} = zeros(NU) #trim point control input
    z_trim::MVector{NZ, Float64} = zeros(NZ) #trim point command vector
    bound_lo::MVector{NU,Float64} = fill(-Inf, NU) #lower output bounds
    bound_hi::MVector{NU,Float64} = fill(Inf, NU) #upper output bounds
end

function LQR{NX, NU, NZ}() where {NX, NU, NZ}
    @assert NZ <= NU "Can't have more command variables than control inputs"
    NUX = NU * NX
    NUZ = NU * NZ
    LQR{NX, NU, NZ, NUX, NUZ}()
end

@kwdef struct LQRInput{NX, NU, NZ, NUX, NUZ}
    sat_ext::MVector{NU,Int64} = zeros(Int64, NU) #saturation input signal
    z_ref::MVector{NZ, Float64} = zeros(NZ) #command vector reference
    z::MVector{NZ, Float64} = zeros(NZ) #current command vector value
    x::MVector{NX, Float64} = zeros(NX) #current state vector value
end

function LQRInput{NX, NU, NZ}(args...; kwargs...) where {NX, NU, NZ}
    NUX = NU * NX
    NUZ = NU * NZ
    LQRInput{NX, NU, NZ, NUX, NUZ}(args...; kwargs...)
end

@kwdef struct LQROutput{NX, NU, NZ, NUX, NUZ}
    K_fbk::SMatrix{NU, NX, Float64, NUX} = zeros(SMatrix{NU, NX}) #state feedback matrix
    K_fwd::SMatrix{NU, NZ, Float64, NUZ} = zeros(SMatrix{NU, NZ}) #feedforward matrix
    K_int::SMatrix{NU, NZ, Float64, NUZ} = zeros(SMatrix{NU, NZ}) #integrator gain matrix
    x_trim::SVector{NX, Float64} = zeros(SVector{NX}) #trim point state
    u_trim::SVector{NU, Float64} = zeros(SVector{NU}) #trim point control input
    z_trim::SVector{NZ, Float64} = zeros(SVector{NZ}) #trim point command variable
    bound_lo::SVector{NU,Float64} = fill(-Inf, SVector{NU}) #lower output bounds
    bound_hi::SVector{NU,Float64} = fill(Inf, SVector{NU}) #upper output bounds
    sat_ext::SVector{NU,Int64} = zeros(SVector{NU, Int64}) #saturation input signal
    z_ref::SVector{NZ, Float64} = zeros(SVector{NZ}) #command variable reference
    z::SVector{NZ, Float64} = zeros(SVector{NZ}) #current command vector value
    x::SVector{NX, Float64} = zeros(SVector{NX}) #current state vector value
    int_in::SVector{NU,Float64} = zeros(SVector{NU}) #integrator input
    int_halted::SVector{NU,Bool} = zeros(SVector{NU, Bool}) #integration halted
    int_out::SVector{NU,Float64} = zeros(SVector{NU}) #integrator output
    out_free::SVector{NU,Float64} = zeros(SVector{NU}) #total output, free
    out_sat::SVector{NU,Int64} = zeros(SVector{NU, Int64}) #current output saturation status
    output::SVector{NU,Float64} = zeros(SVector{NU}) #actual output
end

function LQROutput{NX, NU, NZ}(args...; kwargs...) where {NX, NU, NZ}
    NUX = NU * NX
    NUZ = NU * NZ
    LQROutput{NX, NU, NZ, NUX, NUZ}(args...; kwargs...)
end

@kwdef struct LQRState{NX, NU}
    int_out_0::MVector{NU,Float64} = zeros(NU) #previous integrator path state
    out_sat_0::MVector{NU,Int64} = zeros(NU) #previous output saturation status
end

function Modeling.Y(::LQR{NX, NU, NZ, NUX, NUZ}) where {NX, NU, NZ, NUX, NUZ}
    LQROutput{NX, NU, NZ, NUX, NUZ}()
end

function Modeling.U(::LQR{NX, NU, NZ, NUX, NUZ}) where {NX, NU, NZ, NUX, NUZ}
    LQRInput{NX, NU, NZ, NUX, NUZ}()
end

function Modeling.S(::LQR{NX, NU, NZ, NUX, NUZ}) where {NX, NU, NZ, NUX, NUZ}
    LQRState{NX, NU}()
end

function Modeling.init!(mdl::Model{<:LQR})
    mdl.u.z_ref .= 0
    mdl.u.z .= 0
    mdl.u.x .= 0
    mdl.u.sat_ext .= 0
    mdl.s.int_out_0 .= 0
    mdl.s.out_sat_0 .= 0
end

@no_ode LQR
@no_step LQR

function Modeling.f_periodic!(::NoScheduling, mdl::Model{<:LQR})

    (; parameters, s, u, Δt) = mdl
    (; K_fbk, K_fwd, K_int, x_trim, u_trim, z_trim, bound_lo, bound_hi) = parameters
    (; int_out_0, out_sat_0) = s
    (; sat_ext, z_ref, z, x) = u

    K_fbk, K_fwd, K_int = map(SMatrix, (
    K_fbk, K_fwd, K_int))

    x_trim, u_trim, z_trim, bound_lo, bound_hi = map(SVector, (
    x_trim, u_trim, z_trim, bound_lo, bound_hi))

    sat_ext, z_ref, z, x = map(SVector, (
    sat_ext, z_ref, z, x))

    int_out_0, out_sat_0 = map(SVector, (
    int_out_0, out_sat_0))

    int_in = K_int * (z_ref - z)
    int_halted = ((sign.(int_in .* out_sat_0) .> 0) .|| (sign.(int_in .* sat_ext) .> 0))
    int_out = int_out_0 + Δt * int_in .* .!int_halted

    out_free = u_trim + int_out + K_fwd * (z_ref - z_trim) - K_fbk * (x - x_trim)

    out_sat = (out_free .>= bound_hi) - (out_free .<= bound_lo)
    output = clamp.(out_free, bound_lo, bound_hi)

    s.int_out_0 .= int_out
    s.out_sat_0 .= out_sat

    mdl.y = LQROutput(; K_fbk, K_fwd, K_int, x_trim, u_trim, z_trim,
        bound_lo, bound_hi, sat_ext, z_ref, z, x,
        int_in, int_out, int_halted, out_free, out_sat, output)

end

function GUI.draw(mdl::Model{<:LQR})

    (; K_fbk, K_fwd, K_int, x_trim, u_trim, z_trim, bound_lo, bound_hi,
        sat_ext, z_ref, z, x, int_in, int_halted, int_out,
        out_free, out_sat, output) = mdl.y

    CImGui.Text("Feedback Gain = $K_fbk")
    CImGui.Text("Forward Gain = $K_fwd")
    CImGui.Text("Integral Gain = $K_int")
    CImGui.Text("Trim State Vector = $x_trim")
    CImGui.Text("Trim Control Vector = $u_trim")
    CImGui.Text("Trim Command Vector = $z_trim")
    CImGui.Text("Lower Output Bound = $bound_lo")
    CImGui.Text("Upper Output Bound = $bound_hi")
    CImGui.Text("External Saturation Input = $sat_ext")
    CImGui.Text("Reference Command Vector = $z_ref")
    CImGui.Text("Current Command Vector = $z")
    CImGui.Text("Current State Vector = $x")
    CImGui.Text("Integrator Input = $int_in")
    CImGui.Text("Integrator Halted = $int_halted")
    CImGui.Text("Integrator Output = $int_out")
    CImGui.Text("Free Output = $out_free")
    CImGui.Text("Output Saturation = $out_sat")
    CImGui.Text("Actual Output = $output")

end #function


########################## Controller Data Handling ############################
################################################################################

@kwdef struct PIDData{T} <: FieldVector{4, T}
    k_p::T = 1.0
    k_i::T = 0.0
    k_d::T = 0.1
    τ_f::T = 0.01
end

@kwdef struct LQRData{CB, CF, CI, X, U, Z}
    K_fbk::CB
    K_fwd::CF
    K_int::CI
    x_trim::X
    u_trim::U
    z_trim::Z
end

############################## Data Points #####################################

const PIDDataPoint = PIDData{Float64}

Base.getproperty(data::PIDDataPoint, name::Symbol) = getproperty(data, Val(name))

@generated function Base.getproperty(data::PIDDataPoint, ::Val{S}) where {S}
    if S ∈ fieldnames(PIDData)
        return :(getfield(data, $(QuoteNode(S))))
    elseif S === :T_i
        return :(data.k_p / data.k_i)
    elseif S === :T_d
        return :(data.k_d / data.k_p)
    else
        return :(error("PIDData has no property $S"))
    end
end


const LQRDataPoint{NX, NU, NZ, NUX, NUZ} = LQRData{
    SMatrix{NU, NX, Float64, NUX},
    SMatrix{NU, NZ, Float64, NUZ},
    SMatrix{NU, NZ, Float64, NUZ},
    SVector{NX, Float64},
    SVector{NU, Float64},
    SVector{NZ, Float64}}

function LQRDataPoint(; K_fbk::AbstractMatrix{<:Real},
    K_fwd::AbstractMatrix{<:Real}, K_int::AbstractMatrix{<:Real},
    x_trim::AbstractVector{<:Real}, u_trim::AbstractVector{<:Real},
    z_trim::AbstractVector{<:Real})

    NX, NU, NZ = map(length, (x_trim, u_trim, z_trim))

    #extract data to avoid issues with custom Array types
    LQRDataPoint{NX, NU, NZ, NU * NX, NU * NZ}(
        K_fbk[:, :], K_fwd[:, :], K_int[:, :], x_trim[:], u_trim[:], z_trim[:])

end

function assign!(mdl::Model{<:PID}, point::PIDDataPoint)
    (; k_p, k_i, k_d, τ_f) = point
    mdl.parameters.k_p[] = k_p
    mdl.parameters.k_i[] = k_i
    mdl.parameters.k_d[] = k_d
    mdl.parameters.τ_f[] = τ_f
end

function assign!(mdl::Model{<:LQR}, point::LQRDataPoint)
    (; K_fbk, K_fwd, K_int, x_trim, u_trim, z_trim) = point
    mdl.parameters.K_fbk .= K_fbk
    mdl.parameters.K_fwd .= K_fwd
    mdl.parameters.K_int .= K_int
    mdl.parameters.x_trim .= x_trim
    mdl.parameters.u_trim .= u_trim
    mdl.parameters.z_trim .= z_trim
end


############################### Lookup Data ####################################

const LookupBounds = Vector{NTuple{2, Float64}}

function save_lookup_data(points::Union{Array{<:PIDDataPoint}, Array{<:LQRDataPoint}},
                        bounds::LookupBounds, fname::String = joinpath(@__DIR__, "test.h5"))

    @assert ndims(points) == length(bounds)

    #convert to NamedTuple of Arrays
    nt = StructArrays.components(StructArray(points))

    fid = h5open(fname, "w")

    create_group(fid, "data")
    foreach(propertynames(nt)) do name
        array = getproperty(nt, name)
        fid["data"][string(name)] = stack(array)
    end

    fid["bounds"] = stack(bounds) #2xN matrix

    close(fid)

end


#returns (Array{PIDDataPoint}, LookupBounds)
function load_lookup_data_pid(fname::String)

    fid = h5open(fname, "r")
    #read entries as ordered in PIDData
    points_tuple = map(name -> read(fid["data"][string(name)]), fieldnames(PIDData))
    bounds_stacked = read(fid["bounds"])
    close(fid)

    #convert to Array{PIDDataPoint}
    points =  StructArray{PIDDataPoint}(points_tuple) |> Array

    #rearrange stacked bounds into LookupBounds
    bounds = mapslices(x->tuple(x...), bounds_stacked, dims = 1) |> vec

    #number of interpolation dimensions must match bounds vector length
    @assert ndims(points) == length(bounds)

    return (points = points, bounds = bounds)

end


#returns (Array{<:LQRDataPoint}, LookupBounds)
function load_lookup_data_lqr(fname::String)

    fid = h5open(fname, "r")
    #read entries as ordered in LQRData
    points_stacked = LQRData(map(name -> read(fid["data"][string(name)]), fieldnames(LQRData))...)
    bounds_stacked = read(fid["bounds"])
    close(fid)

    #determine number of interpolation dimensions
    D = ndims(points_stacked.x_trim) - 1

    #generate tuple of N-dimensional arrays from stacked LQRData fields
    points_tuple = map(fieldnames(LQRData)) do name
        stacked_field = getproperty(points_stacked, name)
        if name ∈ (:x_trim, :u_trim, :z_trim) #vector parameter
            return map(SVector{size(stacked_field)[1]}, eachslice(stacked_field; dims = Tuple(2:D+1)))
        elseif name ∈ (:K_fbk, :K_fwd, :K_int) #matrix parameter
            return map(SMatrix{size(stacked_field)[1:2]...}, eachslice(stacked_field; dims = Tuple(3:D+2)))
        else
            error("Unknown field")
        end
    end

    #convert tuple of Arrays to Array{<:LQRDataPoint}
    NX, NU, NZ = map(name -> size(getproperty(points_stacked, name), 1), (:x_trim, :u_trim, :z_trim))
    points = StructArray{LQRDataPoint{NX, NU, NZ, NU*NX, NU*NZ}}(points_tuple) |> Array

    #rearrange stacked bounds into LookupBounds
    bounds = mapslices(x->tuple(x...), bounds_stacked, dims = 1) |> vec
    @assert length(bounds) == D

    return (points = points, bounds = bounds)

end

############################### Lookups ########################################

const PIDDataLookup = PIDData{T} where {T <: AbstractInterpolation}

const LQRDataLookup = LQRData{CB, CF, CI, X, U, Z} where {
    CB <: AbstractInterpolation,
    CF <: AbstractInterpolation,
    CI <: AbstractInterpolation,
    X <: AbstractInterpolation,
    U <: AbstractInterpolation,
    Z <: AbstractInterpolation}


function build_interps(points::Union{Array{<:PIDDataPoint}, Array{<:LQRDataPoint}}, bounds::LookupBounds)

    #convert array of points to a NamedTuple of arrays
    @assert ndims(points) == length(bounds)
    itp_lengths = size(points)
    nt_points = StructArrays.components(StructArray(points))

    #define interpolation mode and scaling, handling singleton dimensions
    itp_args = map(bounds, itp_lengths) do b, l
        r = range(b..., length = l)
        (mode, scaling) = length(r) > 1 ? (BSpline(Linear()), r) : (NoInterp(), 1:1)
        return (mode = mode, scaling = scaling)
    end |> collect |> StructArray

    (; mode, scaling) = itp_args
    interps = [extrapolate(scale(interpolate(getproperty(nt_points, p), tuple(mode...)), scaling...), Flat())
                for p in propertynames(nt_points)]

    return interps

end

build_lookup_pid(fname::String) = PIDData(build_interps(load_lookup_data_pid(fname)...)...)
build_lookup_lqr(fname::String) = LQRData(build_interps(load_lookup_data_lqr(fname)...)...)

function (lookup::PIDDataLookup)(args::Vararg{Real, N}) where {N}
    (; k_p, k_i, k_d, τ_f) = lookup
    PIDData(; k_p = k_p(args...),
                k_i = k_i(args...),
                k_d = k_d(args...),
                τ_f = τ_f(args...)
                )
end

function (lookup::LQRDataLookup)(args::Vararg{Real, N}) where {N}
    (; K_fbk, K_fwd, K_int, x_trim, u_trim, z_trim) = lookup
    LQRData(;
        K_fbk = K_fbk(args...),
        K_fwd = K_fwd(args...),
        K_int = K_int(args...),
        x_trim = x_trim(args...),
        u_trim = u_trim(args...),
        z_trim = z_trim(args...),
        )
end


end #submodule

############################ PID Optimization ##################################
################################################################################

module PIDOpt

using StaticArrays, NLopt, ControlSystems
using RobustAndOptimalControl: hinfnorm2
using Trapz: trapz
using ..Control.Discrete: PIDData, PIDDataPoint

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
    lower_bounds::PIDDataPoint = PIDDataPoint(; k_p = 0.0, k_i = 0.0, k_d = 0.0, τ_f = 0.01)
    upper_bounds::PIDDataPoint = PIDDataPoint(; k_p = 50.0, k_i = 50.0, k_d = 10.0, τ_f = 0.01)
    initial_step::PIDDataPoint = PIDDataPoint(; k_p = 0.01, k_i = 0.01, k_d = 0.01, τ_f = 0.01)
end

@kwdef struct Results
    exit_flag::Symbol
    cost::Float64
    metrics::Metrics{Float64}
    data::PIDDataPoint
end

function build_PID(data::PIDDataPoint)
    (; k_p, k_i, k_d, τ_f) = data
    (k_p + k_i * tf(1, [1,0]) + k_d * tf([1, 0], [τ_f, 1])) |> ss
end

function Metrics(plant::AbstractStateSpace, pid::AbstractStateSpace,
                       settings::Settings)

    S = sensitivity(plant, pid) #sensitivity function

    #robust computation of H-∞ norm
    Ms, _ = hinfnorm2(S)
    Ms = min(Ms, 1e3) #allow optimizing unstable systems

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
                    data_0::PIDDataPoint = PIDDataPoint(), #initial condition
                    settings::Settings = Settings(),
                    weights::Metrics{<:Real} = Metrics(ones(5)),
                    global_search::Bool = true)

    x0 = data_0 |> Vector
    lower_bounds = settings.lower_bounds |> Vector
    upper_bounds = settings.upper_bounds |> Vector
    initial_step = settings.initial_step |> Vector
    maxeval = settings.maxeval

    x0 = clamp.(x0, lower_bounds, upper_bounds)

    plant = ss(plant)
    f_opt = let plant = plant, settings = settings, weights = weights
        function (x::Vector{Float64}, ::Vector{Float64})
            pid = build_PID(PIDDataPoint(x...))
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

    data_opt = PIDDataPoint(minx...)
    pid_opt = build_PID(data_opt)
    metrics_opt = Metrics(plant, pid_opt, settings)

    return Results(exit_flag, minf, metrics_opt, data_opt)


end

function check_results(results::Results, thresholds::Metrics{Float64})

    (; exit_flag, metrics) = results

    success = true
    if !((exit_flag === :ROUNDOFF_LIMITED) | (exit_flag === :STOPVAL_REACHED))
        @warn "Unexpected exit flag: $exit_flag"
        success = false
    end
    if !(metrics.Ms < thresholds.Ms) #sensitivity function maximum magnitude
        @warn "Sensitivity function maximum magnitude exceeded: $(metrics.Ms) > $(thresholds.Ms)"
        success = false
    end
    if !(metrics.∫e < thresholds.∫e) #normalized absolute integrated error
        @warn "Normalized absolute integrated error exceeded: $(metrics.∫e) > $(thresholds.∫e)"
        success = false
    end
    if !(metrics.ef < thresholds.ef) #remaining error after t_sim
        @warn "Absolute final error exceeded: $(metrics.ef) > $(thresholds.ef)"
        success = false
    end
    return success

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