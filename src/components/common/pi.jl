
######################### PICompensator ########################################
################################################################################
struct PICompensator{N} <: Component
    k_p::SVector{N,Float64} #proportional gain
    k_i::SVector{N,Float64} #integral gain
    k_l::SVector{N,Float64} #integrator leak factor
end

# PICompensator{N}(; k_p = 5.0, k_i = 400.0, k_l = 0.2) where {N} = PICompensator{N}(k_p, k_i, k_l)

function PICompensator{N}(k_p::Real, k_i::Real, k_l::Real) where {N}
    PICompensator{N}(fill(k_p, N), fill(k_i, N), fill(k_l, N))
end

Base.@kwdef struct PICompensatorU{N}
    reset::MVector{N,Bool} = zeros(Bool, N)
    sat_enable::MVector{N,Bool} = zeros(N)
    sat_bounds::MVector{N,Float64} = zeros(N)
    input::MVector{N,Float64} = zeros(N)
end

Base.@kwdef struct PICompensatorY{N}
    reset::SVector{N,Bool} = zeros(SVector{N, Bool}) #reset input
    input::SVector{N,Float64} = zeros(SVector{N}) #input signal
    e_p::SVector{N,Float64} = zeros(SVector{N}) #proportional output
    e_i::SVector{N,Float64} = zeros(SVector{N}) #integral output
    e_raw::SVector{N,Float64} = zeros(SVector{N}) #total output, raw
    output::SVector{N,Float64} = zeros(SVector{N}) #total output, clipped
    sat::SVector{N,Bool} = zeros(SVector{N, Bool}) #saturation flag
end


################################### UP TO HERE #############################

init(::SystemX, ::PICompensator{N}) where {N} = zeros(N) #v friction integrator states
init(::SystemY, ::PICompensator{N}) where {N} = PICompensatorY{N}()
init(::SystemU, ::PICompensator{N}) where {N} = PICompensatorU{N}()

f_ode!(sys::System{<:PICompensator{1}}, v_in::Real) = f_ode!(sys, SVector{1, Float64}(v_in))

function f_ode!(sys::System{<:PICompensator{N}}, v_in::AbstractVector{<:Real}) where {N}

    @unpack k_p, k_i, k_l = sys.params

    v = SVector{N, Float64}(v_in)
    s = SVector{N, Float64}(sys.x)
    reset = SVector{N, Bool}(sys.u.reset)

    α_p = -k_p .* v
    α_i = -k_i .* s
    α_raw = α_p + α_i #raw μ scaling
    α = clamp.(α_raw, -1, 1) #clipped μ scaling
    sat = abs.(α_raw) .> abs.(α) #saturated?
    sys.ẋ .= (v - k_l .* s) .* .!sat .* .!reset #if not, integrator accumulates

    sys.y = PICompensatorY(; reset, v, s, α_p, α_i, α_raw, α, sat)

end

function f_step!(sys::System{<:PICompensator{N}}) where {N}

    x = sys.x
    reset = sys.u.reset

    #not vectorized to avoid allocations
    x_mod = false
    for i in 1:N
        x_tmp = x[i]
        x[i] *= !reset[i]
        x_mod = x_mod || (x[i] != x_tmp)
    end
    return x_mod

end

############################## Plotting ########################################

function make_plots(th::TimeHistory{<:PICompensatorY}; kwargs...)

    pd = OrderedDict{Symbol, Plots.Plot}()

    splt_v = plot(th.v; title = "Velocity",
        ylabel = L"$v \ (m/s)$", kwargs...)

    splt_s = plot(th.s; title = "Velocity Integral",
        ylabel = L"$s \ (m)$", kwargs...)

    splt_α_p = plot(th.α_p; title = "Proportional Term",
        ylabel = L"$\alpha_p$", kwargs...)

    splt_α_i = plot(th.α_i; title = "Integral Term",
        ylabel = L"$\alpha_i$", kwargs...)

    splt_α_raw = plot(th.α_raw; title = "Raw Output",
        ylabel = L"$\alpha_{raw}$", kwargs...)

    splt_α = plot(th.α; title = "Clipped Output",
        ylabel = L"$\alpha$", kwargs...)

    splt_sat = plot(th.sat; title = "Saturation",
        ylabel = L"$S$", kwargs...)

    pd[:vs] = plot(splt_v, splt_s, splt_sat;
        plot_title = "Contact Point Kinematics",
        layout = (1,3),
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    pd[:pi] = plot(splt_α_p, splt_α_i, splt_sat;
        plot_title = "Proportional and Integral Terms",
        layout = (1,3),
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    pd[:output] = plot(splt_α_raw, splt_α, splt_sat;
        plot_title = "PICompensator Output",
        layout = (1,3),
        kwargs..., plot_titlefontsize = 20) #override titlefontsize after kwargs

    return pd

end
