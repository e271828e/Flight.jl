module Control

using Flight.FlightCore

############################## Continuous Models ###############################
################################################################################
module Continuous ##############################################################

using ComponentArrays, StaticArrays, UnPack, LinearAlgebra
using ControlSystems: ControlSystemsBase, ControlSystems, ss
using RobustAndOptimalControl
using Plots, LaTeXStrings, DataStructures

using Flight.FlightCore

using ..Control

################################################################################
########################### LinearizedSS ###################################

const tV = AbstractVector{<:Float64}
const tM = AbstractMatrix{<:Float64}

struct LinearizedSS{ LX, LU, LY, #state, input and output vector lengths
                        tX <: tV, tU <: tV, tY <: tV,
                        tA <: tM, tB <: tM, tC <: tM, tD <: tM} <: ModelDefinition

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

function LinearizedSS(mdl::ControlSystemsBase.StateSpace{ControlSystemsBase.Continuous, <:AbstractFloat})
    @unpack A, B, C, D, nx, nu, ny = mdl
    ẋ0 = zeros(nx); x0 = zeros(nx); u0 = zeros(nu); y0 = zeros(ny)
    LinearizedSS(; ẋ0, x0, u0, y0, A, B, C, D)
end

Modeling.X(cmp::LinearizedSS) = copy(cmp.x0)
Modeling.U(cmp::LinearizedSS) = copy(cmp.u0)
Modeling.Y(cmp::LinearizedSS) = SVector{length(cmp.y0)}(cmp.y0)

@no_periodic LinearizedSS
@no_step LinearizedSS

function Modeling.f_ode!(mdl::Model{<:LinearizedSS{LX, LU, LY}}) where {LX, LU, LY}

    @unpack ẋ, x, u, y, parameters = mdl
    @unpack ẋ0, x0, u0, y0, A, B, C, D, x_cache, y_cache, y_cache_out, Δx_cache, Δu_cache = parameters

    #non-allocating equivalent of:
    #ẋ = ẋ0 + A * (x - x0) + B * (u - u0)
    #y = y0 + C * (x - x0) + D * (u - u0)

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

    mdl.y = SVector{LY}(y_cache_out)

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

function submodel(nss::RobustAndOptimalControl.NamedStateSpace;
                  x = nss.x, u = nss.u, y = nss.y)

    #to do: generalize for scalars

    x_axis = Axis(nss.x)
    u_axis = Axis(nss.u)
    y_axis = Axis(nss.y)

    A_nss = ComponentMatrix(nss.A, x_axis, x_axis)
    B_nss = ComponentMatrix(nss.B, x_axis, u_axis)
    C_nss = ComponentMatrix(nss.C, y_axis, x_axis)
    D_nss = ComponentMatrix(nss.D, y_axis, u_axis)

    A_sub = A_nss[x, x]
    B_sub = B_nss[x, u]
    C_sub = C_nss[y, x]
    D_sub = D_nss[y, u]

    ss_sub = ss(A_sub, B_sub, C_sub, D_sub)

    named_ss(ss_sub; x, u, y)

end


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

    @unpack ẋ, x, u, parameters = mdl
    @unpack input, sat_ext = u
    @unpack k_p, k_i, k_l, β_p, bound_lo, bound_hi = parameters

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

    @unpack k_p, k_i, k_l, β_p, bound_lo, bound_hi, input, sat_ext,
            u_p, u_i, x_i, y_p, y_i, out_free, sat_out, output, int_halted = mdl.y

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

using StaticArrays, UnPack, LinearAlgebra
using StructArrays, Interpolations, HDF5 #for lookups
using RobustAndOptimalControl
using Plots, LaTeXStrings, DataStructures

using Flight.FlightCore

using ..Control


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

    @unpack s, u, Δt, parameters = mdl
    @unpack input, sat_ext = u
    @unpack x0, sat_out_0 = s
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

    @unpack s, u, Δt, parameters = mdl
    @unpack input, sat_ext = u
    @unpack bound_lo, bound_hi = parameters

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

    @unpack x0, sat_out_0 = mdl.s
    @unpack input, sat_ext, bound_lo, bound_hi, x1, output, sat_out, halted = mdl.y

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

    @unpack parameters, s, u, Δt = mdl
    @unpack u1 = u
    @unpack u0, x0 = s
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

    @unpack z, p, k, u1, y1 = mdl.y

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

@kwdef struct PID <: ModelDefinition
    k_p::Ref{Float64} = Ref(1.0) #proportional gain
    k_i::Ref{Float64} = Ref(0.1) #integral gain
    k_d::Ref{Float64} = Ref(0.1) #derivative gain
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
    k_i::Float64 = 0.1 #integral gain
    k_d::Float64 = 0.1 #derivative gain
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

    @unpack parameters, u, s, Δt = mdl
    @unpack k_p, k_i, k_d, τ_f, β_p, β_d, bound_lo, bound_hi = parameters
    @unpack input, sat_ext = u
    @unpack x_i0, x_d0, sat_out_0 = s

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

    @unpack parameters, s, u, Δt = mdl
    @unpack k_p, k_i, k_d, τ_f, β_p, β_d, bound_lo, bound_hi = parameters
    @unpack x_i0, x_d0, sat_out_0 = s
    @unpack input, sat_ext = u

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

    @unpack k_p, k_i, k_d, τ_f, β_p, β_d, bound_lo, bound_hi, input, sat_ext, u_p, u_i, u_d,
            y_p, y_i, y_d, out_free, sat_out, output, int_halted = mdl.y

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
    C_fbk::MMatrix{NU, NX, Float64, NUX} = zeros(NU, NX) #state feedback matrix
    C_fwd::MMatrix{NU, NZ, Float64, NUZ} = zeros(NU, NZ) #feedforward matrix
    C_int::MMatrix{NU, NZ, Float64, NUZ} = zeros(NU, NZ) #integrator gain matrix
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
    C_fbk::SMatrix{NU, NX, Float64, NUX} = zeros(SMatrix{NU, NX}) #state feedback matrix
    C_fwd::SMatrix{NU, NZ, Float64, NUZ} = zeros(SMatrix{NU, NZ}) #feedforward matrix
    C_int::SMatrix{NU, NZ, Float64, NUZ} = zeros(SMatrix{NU, NZ}) #integrator gain matrix
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

    @unpack parameters, s, u, Δt = mdl
    @unpack C_fbk, C_fwd, C_int, x_trim, u_trim, z_trim, bound_lo, bound_hi = parameters
    @unpack int_out_0, out_sat_0 = s
    @unpack sat_ext, z_ref, z, x = u

    C_fbk, C_fwd, C_int = map(SMatrix, (
    C_fbk, C_fwd, C_int))

    x_trim, u_trim, z_trim, bound_lo, bound_hi = map(SVector, (
    x_trim, u_trim, z_trim, bound_lo, bound_hi))

    sat_ext, z_ref, z, x = map(SVector, (
    sat_ext, z_ref, z, x))

    int_out_0, out_sat_0 = map(SVector, (
    int_out_0, out_sat_0))

    int_in = C_int * (z_ref - z)
    int_halted = ((sign.(int_in .* out_sat_0) .> 0) .|| (sign.(int_in .* sat_ext) .> 0))
    int_out = int_out_0 + Δt * int_in .* .!int_halted

    out_free = u_trim + int_out + C_fwd * (z_ref - z_trim) - C_fbk * (x - x_trim)

    out_sat = (out_free .>= bound_hi) - (out_free .<= bound_lo)
    output = clamp.(out_free, bound_lo, bound_hi)

    s.int_out_0 .= int_out
    s.out_sat_0 .= out_sat

    mdl.y = LQROutput(; C_fbk, C_fwd, C_int, x_trim, u_trim, z_trim,
        bound_lo, bound_hi, sat_ext, z_ref, z, x,
        int_in, int_out, int_halted, out_free, out_sat, output)

end

function GUI.draw(mdl::Model{<:LQR})

    @unpack C_fbk, C_fwd, C_int, x_trim, u_trim, z_trim, bound_lo, bound_hi,
            sat_ext, z_ref, z, x, int_in, int_halted, int_out,
            out_free, out_sat, output = mdl.y

        CImGui.Text("Feedback Gain = $C_fbk")
        CImGui.Text("Forward Gain = $C_fwd")
        CImGui.Text("Integral Gain = $C_int")
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

@kwdef struct LQRParams{CB, CF, CI, X, U, Z}
    C_fbk::CB
    C_fwd::CF
    C_int::CI
    x_trim::X
    u_trim::U
    z_trim::Z
end


############################### Gain Scheduling ################################
################################################################################

const LookupBounds{N} = NTuple{N, NTuple{2, Real}}

############################## Data Points #####################################

const PIDPoint = PIDParams{<:Real}

const LQRPoint = LQRParams{CB, CF, CI, X, U, Z} where {
    CB <: AbstractMatrix{<:Real},
    CF <: AbstractMatrix{<:Real},
    CI <: AbstractMatrix{<:Real},
    X <: AbstractVector{<:Real},
    U <: AbstractVector{<:Real},
    Z <: AbstractVector{<:Real}}


function assign!(mdl::Model{<:PID}, point::PIDPoint)
    @unpack k_p, k_i, k_d, τ_f = point
    mdl.parameters.k_p[] = k_p
    mdl.parameters.k_i[] = k_i
    mdl.parameters.k_d[] = k_d
    mdl.parameters.τ_f[] = τ_f
end

function assign!(mdl::Model{<:LQR}, point::LQRPoint)
    @unpack C_fbk, C_fwd, C_int, x_trim, u_trim, z_trim = point
    mdl.parameters.C_fbk .= C_fbk
    mdl.parameters.C_fwd .= C_fwd
    mdl.parameters.C_int .= C_int
    mdl.parameters.x_trim .= x_trim
    mdl.parameters.u_trim .= u_trim
    mdl.parameters.z_trim .= z_trim
end


############################# Data Arrays ######################################

const PIDData{N} = PIDParams{<:AbstractArray{<:Real, N}}

const LQRData{N} = LQRParams{CB, CF, CI, X, U, Z} where {
    CB <: AbstractArray{<:AbstractMatrix{<:Real}, N},
    CF <: AbstractArray{<:AbstractMatrix{<:Real}, N},
    CI <: AbstractArray{<:AbstractMatrix{<:Real}, N},
    X <: AbstractArray{<:AbstractVector{<:Real}, N},
    U <: AbstractArray{<:AbstractVector{<:Real}, N},
    Z <: AbstractArray{<:AbstractVector{<:Real}, N}}

function PIDData(params::Array{<:PIDPoint})
    params_nt = StructArrays.components(StructArray(params))
    PIDParams(values(params_nt)...)
end

function LQRData(params::Array{<:LQRPoint})
    params_nt = StructArrays.components(StructArray(params))
    LQRParams(values(params_nt)...)
end

function save_data(params::Union{PIDData{N}, LQRData{N}}, bounds::LookupBounds{N},
                    fname::String = joinpath(@__DIR__, "test.h5")) where {N}

    fid = h5open(fname, "w")

    create_group(fid, "params")
    foreach(propertynames(params)) do name
        array = getproperty(params, name)
        fid["params"][string(name)] = stack(array)
    end

    fid["bounds"] = stack(bounds) #2xN matrix

    close(fid)

end

save_data(params::Array{<:PIDPoint}, args...) = save_data(PIDData(params), args...)
save_data(params::Array{<:LQRPoint}, args...) = save_data(LQRData(params), args...)


function load_data_pid(fname::String)

    fid = h5open(fname, "r")

    #read fieldnames as ordered in PIDParams and into an instance
    params_stacked = map(name -> read(fid["params"][string(name)]), fieldnames(PIDParams))
    bounds_stacked = read(fid["bounds"])

    close(fid)

    #arrange bounds back into a Tuple of 2-Tuples
    bounds = mapslices(x->tuple(x...), bounds_stacked, dims = 1) |> vec |> Tuple
    N = length(bounds) #number of interpolation dimensions

    #PID parameters are scalars, so these are already N-dimensional arrays
    params_tuple = params_stacked

    return (params = PIDParams(params_tuple...), bounds = bounds)

end


function load_data_lqr(fname::String)

    fid = h5open(fname, "r")

    #read fieldnames as ordered in LQRParams and into an instance
    params_stacked = map(name -> read(fid["params"][string(name)]), fieldnames(LQRParams))
    bounds_stacked = read(fid["bounds"])

    close(fid)

    #arrange bounds back into a Tuple of 2-Tuples
    bounds = mapslices(x->tuple(x...), bounds_stacked, dims = 1) |> vec |> Tuple
    N = length(bounds) #number of interpolation dimensions

    #generate Tuple of N-dimensional arrays of either SVectors (for x_trim,
    #u_trim and z_trim) or SMatrices (for C_fbk, C_fwd and C_int)
    params_tuple = map(params_stacked) do p_stacked
        if ndims(p_stacked) == N+1 #vector parameter
            return map(SVector{size(p_stacked)[1]}, eachslice(p_stacked; dims = Tuple(2:N+1)))
        elseif ndims(p_stacked) == N+2 #matrix parameter
            return map(SMatrix{size(p_stacked)[1:2]...}, eachslice(p_stacked; dims = Tuple(3:N+2)))
        else
            error("Number of interpolation dimensions was determined as $N, "*
                "so stacked arrays must be either either $(N+1)-dimensional "*
                "for vector parameters or $(N+2)-dimensional for matrix parameters. "*
                "Stacked array is $(ndims(p_stacked))-dimensional for $p_name")
        end
    end

    return (params = LQRParams(params_tuple...), bounds = bounds)

end


############################### Lookups ########################################

const PIDLookup = PIDParams{T} where {T <: AbstractInterpolation}

const LQRLookup = LQRParams{CB, CF, CI, X, U, Z} where {
    CB <: AbstractInterpolation,
    CF <: AbstractInterpolation,
    CI <: AbstractInterpolation,
    X <: AbstractInterpolation,
    U <: AbstractInterpolation,
    Z <: AbstractInterpolation}


function build_lookup_pid(params::PIDData{N}, bounds::LookupBounds{N}) where {N}

    @assert allequal(size.(params))
    itp_lengths = size(params.k_p)

    #define interpolation mode and ranges, handling singleton dimensions
    itp_args = map(bounds, itp_lengths) do b, l
        r = range(b..., length = l)
        (mode, scaling) = length(r) > 1 ? (BSpline(Linear()), r) : (NoInterp(), 1:1)
        return (mode = mode, scaling = scaling)
    end |> collect |> StructArray

    @unpack mode, scaling = itp_args
    interps = [extrapolate(scale(interpolate(getproperty(params, p), tuple(mode...)), scaling...), Flat())
                for p in propertynames(params)]

    return PIDParams(interps...)

end

function build_lookup_lqr(params::LQRData{N}, bounds::LookupBounds{N}) where {N}

    #lengths of N interpolation dimensions must be consistent among params
    sizes = map(n -> size(getproperty(params, n)), propertynames(params))
    @assert allequal(sizes)
    itp_lengths = sizes[1]

    #define interpolation mode and ranges, handling singleton dimensions
    itp_args = map(bounds, itp_lengths) do b, l
        r = range(b..., length = l)
        (mode, scaling) = length(r) > 1 ? (BSpline(Linear()), r) : (NoInterp(), 1:1)
        return (mode = mode, scaling = scaling)
    end |> collect |> StructArray

    @unpack mode, scaling = itp_args
    interps = [extrapolate(scale(interpolate(getproperty(params, p), tuple(mode...)), scaling...), Flat())
                for p in propertynames(params)]

    return LQRParams(interps...)

end


build_lookup_pid(fname::String) = build_lookup_pid(load_data_pid(fname)...)
build_lookup_lqr(fname::String) = build_lookup_lqr(load_data_lqr(fname)...)


function (lookup::PIDLookup)(args::Vararg{Real, N}) where {N}
    @unpack k_p, k_i, k_d, τ_f = lookup
    PIDParams(; k_p = k_p(args...),
                k_i = k_i(args...),
                k_d = k_d(args...),
                τ_f = τ_f(args...)
                )
end

function (lookup::LQRLookup)(args::Vararg{Real, N}) where {N}
    @unpack C_fbk, C_fwd, C_int, x_trim, u_trim, z_trim = lookup
    LQRParams(;
        C_fbk = C_fbk(args...),
        C_fwd = C_fwd(args...),
        C_int = C_int(args...),
        x_trim = x_trim(args...),
        u_trim = u_trim(args...),
        z_trim = z_trim(args...),
        )
end


end #submodule

################################################################################
############################ PID Optimization ##################################
################################################################################

module PIDOpt

using StaticArrays, UnPack, NLopt, ControlSystems
using RobustAndOptimalControl: hinfnorm2
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