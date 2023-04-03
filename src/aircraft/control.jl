module Control

using ComponentArrays, StaticArrays, UnPack, LinearAlgebra
using CImGui, CImGui.CSyntax
using ControlSystemsBase: ControlSystemsBase, ss

using Flight.FlightCore

export LinearStateSpace, PIContinuous

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

ControlSystemsBase.ss(cmp::LinearStateSpace) = ControlSystemsBase.ss(cmp.A, cmp.B, cmp.C, cmp.D)

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
    input::MVector{N,Float64} = zeros(N)
    hold::MVector{N,Bool} = zeros(Bool, N)
    reset::MVector{N,Bool} = zeros(Bool, N)
    sat_enable::MVector{N,Bool} = ones(Bool, N)
end

Base.@kwdef struct PIContinuousY{N}
    hold::SVector{N,Bool} = zeros(SVector{N, Bool}) #hold integrator state
    reset::SVector{N,Bool} = zeros(SVector{N, Bool}) #reset integrator state
    input::SVector{N,Float64} = zeros(SVector{N}) #input signal
    state::SVector{N,Float64} = zeros(SVector{N}) #integrator state
    out_p::SVector{N,Float64} = zeros(SVector{N}) #proportional term
    out_i::SVector{N,Float64} = zeros(SVector{N}) #integral term
    out_free::SVector{N,Float64} = zeros(SVector{N}) #total output, free
    out::SVector{N,Float64} = zeros(SVector{N}) #total output
    sat_enable::SVector{N,Bool} = zeros(SVector{N, Bool}) #saturation status
    sat_status::SVector{N,Int64} = zeros(SVector{N, Int64}) #saturation status
    int_status::SVector{N,Bool} = zeros(SVector{N, Bool}) #integrator accumulating
end

Systems.init(::SystemX, ::PIContinuous{N}) where {N} = zeros(N)
Systems.init(::SystemY, ::PIContinuous{N}) where {N} = PIContinuousY{N}()
Systems.init(::SystemU, ::PIContinuous{N}) where {N} = PIContinuousU{N}()

function Systems.f_ode!(sys::System{<:PIContinuous{N}}) where {N}

    @unpack k_p, k_i, k_l, bounds = sys.params
    @unpack sat_enable = sys.u

    state = SVector{N, Float64}(sys.x)

    hold = SVector(sys.u.hold)
    reset = SVector(sys.u.reset)
    input = SVector(sys.u.input)
    sat_enable = SVector(sys.u.sat_enable)

    out_p = k_p .* input
    out_i = k_i .* state
    out_free = out_p + out_i #raw output
    out_clamped = clamp.(out_free, bounds[1], bounds[2]) #clamped output
    out = (out_free .* .!sat_enable) .+ (out_clamped .* sat_enable)

    sat_upper = (out_free .>= bounds[2]) .* sat_enable
    sat_lower = (out_free .<= bounds[1]) .* sat_enable
    sat_status = sat_upper - sat_lower
    int_status = sign.(input .* sat_status) .<= 0 #enable integrator?

    sys.ẋ .= (input .* int_status - k_l .* state) .* .!reset .* .!hold

    sys.y = PIContinuousY(; hold, reset, input, state, out_p, out_i,
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

    splt_u = plot(th.input; title = "Input",
        ylabel = L"$u$", kwargs...)

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

#################################### GUI #######################################


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
                reset[i] = reset_ref[]
                sat_enable[i] = sat_ref[]
                CImGui.TreePop()
            end
        end
        CImGui.TreePop()
    end

    CImGui.End()

end

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
        @unpack hold, reset, input, sat_enable, sat_status, out_free, out = y
        CImGui.Text("Hold = $hold")
        CImGui.Text("Reset = $reset")
        CImGui.Text("Input = $input")
        CImGui.Text("Unbounded Output = $out_free")
        CImGui.Text("Output = $out")
        CImGui.Text("Saturation Enabled = $sat_enable")
        CImGui.Text("Saturation Status = $sat_status")
        CImGui.TreePop()
    end

    CImGui.End()

end #function

end #module