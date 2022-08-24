module Essentials

using ComponentArrays, StaticArrays, UnPack, LinearAlgebra
using Flight.Systems
using Flight.Plotting

import ControlSystems #avoids name clash with ControlSystems.StateSpace into scope

export StateSpace, PICompensator


################################################################################
########################### StateSpace Model ###################################

const tV = AbstractVector{<:Float64}
const tM = AbstractMatrix{<:Float64}

struct StateSpace{  LX, LU, LY, #state, input and output vector lengths
                tX <: tV, tU <: tV, tY <: tV,
                tA <: tM, tB <: tM, tC <: tM, tD <: tM} <: Component

    ẋ0::tX; x0::tX; u0::tU; y0::tY; #reference values (for linearized systems)
    A::tA; B::tB; C::tC; D::tD; #state-space matrices
    x_cache::tX; y_cache::tY; y_cache_out::tY;
    Δx_cache::tX; Δu_cache::tU

    function StateSpace(ẋ0, x0, u0, y0, A, B, C, D)

        lengths = map(length, (x0, u0, y0))
        types = map(typeof, (x0, u0, y0, A, B, C, D))

        vectors = map(copy, (ẋ0, x0, u0, y0))
        matrices = map(copy, (A, B, C, D))
        caches = map(copy, (x0, y0, y0, x0, u0))

        new{lengths..., types...}(vectors..., matrices..., caches...)

    end

end

StateSpace(; ẋ0, x0, u0, y0, A, B, C, D) = StateSpace(ẋ0, x0, u0, y0, A, B, C, D)

ControlSystems.ss(cmp::StateSpace) = ControlSystems.ss(cmp.A, cmp.B, cmp.C, cmp.D)

Systems.init(::SystemX, cmp::StateSpace) = copy(cmp.x0)
Systems.init(::SystemU, cmp::StateSpace) = copy(cmp.u0)
Systems.init(::SystemY, cmp::StateSpace) = SVector{length(cmp.y0)}(cmp.y0)

function Systems.f_ode!(sys::System{<:StateSpace{LX, LU, LY}}) where {LX, LU, LY}

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

function Base.filter(cmp::StateSpace; x = (), u = (), y = ())

    x_ind = (!isempty(x) ? x : keys(cmp.x0))
    u_ind = (!isempty(u) ? u : keys(cmp.u0))
    y_ind = (!isempty(y) ? y : keys(cmp.y0))

    ẋ0 = cmp.ẋ0[x_ind]
    x0 = cmp.x0[x_ind]
    u0 = cmp.u0[u_ind]
    y0 = cmp.y0[y_ind]
    A = cmp.A[x_ind, x_ind]
    B = cmp.B[x_ind, u_ind]
    C = cmp.C[y_ind, x_ind]
    D = cmp.D[y_ind, u_ind]

    return StateSpace(; ẋ0, x0, u0, y0, A, B, C, D)

end

#convert a mutable ComponentVector into a (non-allocating) immutable
#ComponentVector with an underlying StaticVector, preserving the original
#ComponentVector's axes
function StaticArrays.SVector{L}(x::ComponentVector) where {L}
    ComponentVector(SVector{L}(getdata(x)), getaxes(x))
end


####################### Proportional-Integral Compensator ######################
################################################################################

struct PICompensator{N} <: Component
    k_p::SVector{N,Float64} #proportional gain
    k_i::SVector{N,Float64} #integral gain
    k_l::SVector{N,Float64} #integrator leak factor
    bounds::NTuple{2,MVector{N,Float64}} #output bounds
end

function PICompensator{N}(; k_p::Real = 1.0, k_i::Real = 0.1, k_l::Real = 0.0,
                            bounds::NTuple{2,Real} = (-1.0, 1.0)) where {N}
    s2v = (x)->fill(x,N)
    PICompensator{N}(s2v(k_p), s2v(k_i), s2v(k_l), (s2v(bounds[1]), s2v(bounds[2])))
end

Base.@kwdef struct PICompensatorU{N}
    input::MVector{N,Float64} = zeros(N)
    reset::MVector{N,Bool} = zeros(Bool, N)
    sat_enable::MVector{N,Bool} = ones(Bool, N)
end

Base.@kwdef struct PICompensatorY{N}
    reset::SVector{N,Bool} = zeros(SVector{N, Bool}) #reset input
    input::SVector{N,Float64} = zeros(SVector{N}) #input signal
    state::SVector{N,Float64} = zeros(SVector{N}) #integrator state
    out_p::SVector{N,Float64} = zeros(SVector{N}) #proportional term
    out_i::SVector{N,Float64} = zeros(SVector{N}) #integral term
    out_free::SVector{N,Float64} = zeros(SVector{N}) #total output, free
    out::SVector{N,Float64} = zeros(SVector{N}) #total output
    sat_status::SVector{N,Bool} = zeros(SVector{N, Bool}) #saturation status
end

Systems.init(::SystemX, ::PICompensator{N}) where {N} = zeros(N)
Systems.init(::SystemY, ::PICompensator{N}) where {N} = PICompensatorY{N}()
Systems.init(::SystemU, ::PICompensator{N}) where {N} = PICompensatorU{N}()

function Systems.f_ode!(sys::System{<:PICompensator{N}}) where {N}

    @unpack k_p, k_i, k_l, bounds = sys.params
    @unpack sat_enable = sys.u

    state = SVector{N, Float64}(sys.x)
    reset = SVector{N, Bool}(sys.u.reset)
    input = SVector{N, Float64}(sys.u.input)

    out_p = k_p .* input
    out_i = k_i .* state
    out_free = out_p + out_i #raw output
    out_clamped = clamp.(out_free, bounds[1], bounds[2]) #clamped output
    out = out_free .* .!sat_enable .+ out_clamped .* sat_enable
    sat_status = out .!= out_free #saturated?
    sys.ẋ .= (input - k_l .* state) .* .!sat_status .* .!reset

    sys.y = PICompensatorY(; reset, input, state, out_p, out_i,
                            out_free, out, sat_status)

end

function Systems.f_step!(sys::System{<:PICompensator{N}}) where {N}

    # @show sys.y.input
    # @show sys.y.out
    # @show sys.y.sat_status

    x = SVector{N,Float64}(sys.x)
    x_new = x .* .!sys.u.reset
    x_mod = any(x .!= x_new)

    sys.x .= x_new
    return x_mod

end


# ############################## Plotting ########################################

function Plotting.make_plots(th::TimeHistory{<:PICompensatorY}; kwargs...)

    pd = OrderedDict{Symbol, Plots.Plot}()

    splt_u = plot(th.input; title = "Input",
        ylabel = L"$u$", kwargs...)

    splt_s = plot(th.state; title = "Integrator State",
        ylabel = L"$s$", kwargs...)

    splt_y_p = plot(th.out_p; title = "Proportional Term",
        ylabel = L"$y_p$", kwargs...)

    splt_y_i = plot(th.out_i; title = "Integral Term",
        ylabel = L"$y_i$", kwargs...)

    splt_y_free = plot(th.out_free; title = "Free Output",
        ylabel = L"$y_{free}$", kwargs...)

    splt_y = plot(th.out; title = "Actual Output",
        ylabel = L"$y$", kwargs...)

    splt_sat = plot(th.sat_status; title = "Saturation Status",
        ylabel = L"$S$", kwargs...)

    pd[:us] = plot(splt_u, splt_s, splt_sat;
        plot_title = "Input & Integrator State",
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


end #module