# module Linearization

#this module should probably be put together with Trim in c172r/tools.jl

using ComponentArrays
using UnPack
using FiniteDiff: finite_difference_jacobian! as jacobian!

using Flight

#flaps and mixture are omitted from the input vector, because they will often be
#saturated at one of their range limits (for example, when operating with fully
#retracted or fully extended flaps, full lean or full rich mixture), and
#linearization will not go well in such cases. flaps and mixture are best kept
#as configuration parameters

#when trimming, we enforced ω_eb_b_dot = v_eOb_b_dot = 0. as parameters, we
#imposed θ_lb_dot and ψ_lb_dot. here, linearization uses ψ_nb, θ_nb, φ_nb. the
#only expected discrepancy is between ψ_lb_dot and ψ_nb_dot due to the z
#component of the NED transport rate, which is absent in ψ_lb_dot

const XLinearTemplate = ComponentVector(
    ψ_lb = 0.0, #heading (lb)
    θ_lb = 0.0, #inclination (lb)
    φ_lb = 0.0, #bank (lb)
    ϕ = 0.0, #latitude
    λ = 0.0, #longitude
    h_e = 0.0, #ellipsoidal altitude
    p_eb = 0.0, #roll rate (ω_eb_b)
    q_eb = 0.0, #pitch rate (ω_eb_b)
    r_eb = 0.0, #yaw rate (ω_eb_b)
    v_xb = 0.0, #x-body ground-referenced velocity (v_eOb_b)
    v_yb = 0.0, #y-body ground-referenced velocity (v_eOb_b)
    v_zb = 0.0, #z-body ground-referenced velocity (v_eOb_b)
    α_filt = 0.0, #filtered AoA
    β_filt = 0.0, #filtered AoS
    ω_eng = 0.0, #engine speed
    m_fuel = 0.0, #fuel load
)

const AxesXLinear = getaxes(XLinearTemplate)
const XLinear{T, D} = ComponentVector{T, D, typeof(AxesXLinear)}

Base.@kwdef struct ULinear
    throttle::Float64 = 0.0
    yoke_y::Float64 = 0.0 #elevator control
    yoke_x::Float64 = 0.0 #aileron control
    pedals::Float64 = 0.0 #rudder control
end

const AxesULinear = (ComponentArrays.Axis(fieldnames(ULinear)),)
const ULinearVector{T, D} = ComponentVector{T, D, typeof(AxesULinear)}

Base.@kwdef struct YLinear
    ψ_lb::Float64 = 0.0 #heading (lb)
    θ_lb::Float64 = 0.0 #inclination (lb)
    φ_lb::Float64 = 0.0 #bank (lb)
    ϕ::Float64 = 0.0 #latitude
    λ::Float64 = 0.0 #longitude
    h_e::Float64 = 0.0 #ellipsoidal altitude
    p_eb::Float64 = 0.0 #roll rate (ω_eb_b)
    q_eb::Float64 = 0.0 #pitch rate (ω_eb_b)
    r_eb::Float64 = 0.0 #yaw rate (ω_eb_b)
    V::Float64 = 0.0 #true airspeed
    α::Float64 = 0.0 #AoA
    β::Float64 = 0.0 #AoS
    f_xb::Float64 = 0.0 #x-body specific force (f_iOb_b)
    f_yb::Float64 = 0.0 #y-body specific force (f_iOb_b)
    f_zb::Float64 = 0.0 #z-body specific force (f_iOb_b)
    ω_eng::Float64 = 0.0 #engine speed
    m_fuel::Float64 = 0.0 #fuel load
end

const AxesYLinear = (ComponentArrays.Axis(fieldnames(YLinear)),)
const YLinearVector{T, D} = ComponentVector{T, D, typeof(AxesYLinear)}


function update_aircraft!(ac::System{<:Cessna172R{NED}}, x::XLinear)

    @unpack x = ψ_lb, θ_lb, φ_lb, ϕ, λ, h_e, p_eb, q_eb, r_eb, v_xb, v_yb, v_zb, α_filt, β_filt, ω_eng, m_fuel

    e_nb = e_lb = REuler(ψ_lb, θ_lb, φ_lb) #arbitrarily set LTF = NED

    ω_eb_b = SVector(p_eb, q_eb, r_eb)
    v_eOb_b = SVector(v_xb, v_yb, v_zb)

    x.pos.e_nb .= SVector(e_nb.ψ, e_nb.θ, e_nb.φ)
    x.pos.ϕ = ϕ
    x.pos.λ = λ
    x.pos.h_e = h_e
    x.vel.ω_eb_b .= ω_eb_b
    x.vel.v_eOb_b .= v_eOb_b

    #airframe-specific
    ac.x.airframe.aero.α_filt .= x.α_filt
    ac.x.airframe.aero.β_filt .= x.β_filt
    ac.x.airframe.pwp.engine.ω = x.ω_eng
    ac.x.airframe.fuel .= x.m_fuel

end

function update_aircraft!(ac::System{<:Cessna172R}, u::ULinearVector)

    @unpack throttle, yoke_y, yoke_x, pedals = u
    @pack! ac.u.avionics = throttle, yoke_y, yoke_x, pedals

end


function update_x!(x::XLinear, ac::System{<:Cessna172R{NED}})

    @unpack e_nb, ϕ, λ, h_e, ω_eb_b, v_eOb_b = ac.x.kinematics

    e_lb = e_nb #arbitrary
    ψ_lb, θ_lb, φ_lb = e_lb.ψ, e_lb.θ, e_lb.φ
    p_eb, q_eb, r_eb = ω_eb_b
    v_xb, v_yb, v_zb = v_eOb_b

    #airframe-specific
    α_filt = ac.x.airframe.aero.α_filt
    β_filt = ac.x.airframe.aero.β_filt
    ω_eng = ac.x.airframe.pwp.engine.ω
    m_fuel = ac.x.airframe.fuel.m

    @pack! x = ψ_lb, θ_lb, φ_lb, ϕ, λ, h_e, p_eb, q_eb, r_eb, v_xb, v_yb, v_zb, α_filt, β_filt, ω_eng, m_fuel

end


function update_ẋ!(ẋ::XLinear, ac::System{<:Cessna172R})

    @unpack α_eb_b, v̇_eOb_b = ac.y.rigidbody
    @unpack q_nb, v_eOb_n, n_e, h_e = ac.y.kinematics

    q_lb = q_nb
    e_lb = REuler(q_lb)
    Ob = Geographic(n_e, h_e)
    ϕ_λ = LatLon(n_e)

    ω_eb_b = SVector(p_eb, q_eb, r_eb)
    ω_en_n = Kinematics.get_ω_en_n(v_eOb_n, Ob)
    ω_el_n = Kinematics.get_ω_el_n(v_eOb_n, Ob)
    ω_el_b = q_nb'(ω_el_n)
    ω_lb_b = ω_eb_b - ω_el_b

    ė_lb = Attitude.dt(e_lb, ω_lb_b)
    ϕ_λ_dot = Geodesy.dt(ϕ_λ, ω_en_n)

    ẋ.ψ_lb = ė_lb.ψ
    ẋ.θ_lb = ė_lb.θ
    ẋ.φ_lb = ė_lb.φ
    ẋ.ϕ_dot = ϕ_λ_dot.ϕ
    ẋ.λ_dot = ϕ_λ_dot.λ
    ẋ.h_e = -v_eOb_n[3]
    ẋ.p_eb, ẋ.q_eb, ẋ.r_eb = α_eb_b
    ẋ.v_xb, ẋ.v_yb, ẋ.v_zb = v̇_eOb_b

    #airframe-specific
    ẋ.α_filt = ac.ẋ.airframe.aero.α_filt
    ẋ.β_filt = ac.ẋ.airframe.aero.β_filt
    ẋ.ω_eng = ac.ẋ.airframe.pwp.engine.ω
    ẋ.m_fuel = ac.ẋ.airframe.fuel.m

function update_u!(u::ULinearVector, ac::System{<:Cessna172R})

    @unpack throttle, yoke_y, yoke_x, pedals = ac.u.avionics
    @pack! u = throttle, yoke_y, yoke_x, pedals

end

function update_y!(y::YLinearVector, ac::System{<:Cessna172R})

    @unpack q_nb, n_e, h_e, ω_eb_b = ac.y.kinematics
    @unpack α, β = ac.y.airframe.aero
    @unpack ψ, θ, φ = REuler(q_nb)
    @unpack ϕ, λ = LatLon(n_e)

    p_eb, q_eb, r_eb = ω_eb_b
    f_xb, f_yb, f_zb = ac.y.rigidbody.f_Ob_b
    V = ac.y.airflow.TAS
    ω_eng = ac.y.airframe.pwp.engine.ω
    m_fuel = ac.y.airframe.fuel.m

    @pack! y = ψ_lb, θ_lb, φ_lb, ϕ, λ, h_e, p_eb, q_eb, r_eb, V, α, β, f_xb, f_yb, f_zb, ω_eng, m_fuel

end

function test()

    ac = System(Cessna172R(NED()))
    env = System(SimpleEnvironment())
    params = C172R.Trim.Parameters()

    # state = C172R.Trim.State() #optional initial trim guess
    (exit_flag, trim_state) = C172R.Trim.trim!(; ac, env, params) #ac is now trimmed

    #save the trimmed aircraft ẋ, x, u and y for later
    ẋ_ref = similar(XLinearTemplate); update_ẋ(ẋ_ref, ac)
    x_ref = similar(XLinearTemplate); update_x(x_ref, ac)
    u_ref = similar(UTemplate); assign!(u_ref, ac) #get the reference value from the trimmed aircraft
    y_ref = similar(YTemplate); assign!(y_ref, ac) #idem

    # @show exit_flag
    # @show trim_state
    # @show ẋ_ref
    # @show x_ref
    # @show u_ref
    # @show y_ref

    #function wrapper around f_ode!() mutating ẋ and y. define a cache in case
    #we get generic y or u
    f_nonlinear! = let ac = ac, env = env, params = params, state = trim_state,
                       u_axes = getaxes(UTemplate), y_axes = getaxes(YTemplate)

        function (ẋ, y, x, u)

            #we need these because finite_difference_jacobian! sometimes
            #provides generic Vectors instead of ComponentVectors. creating
            #these does not allocate, because the underlying data is already
            #provided by u and y. updating u_cv and y_cv also updates u and y,
            #while satisfying the interface requirements of the assign!
            #functions
            u_cv = ComponentVector(u, u_axes)
            y_cv = ComponentVector(y, y_axes)

            #make sure any input or state not set by x and u is at its reference
            #trim value. this reverts changes to the aircraft done by functions
            #sharing the same aircraft instance
            C172R.Trim.assign!(ac, env, params, state)

            # ac.x .= x
            assign!(ac, x_cv)
            assign!(ac, u_cv)

            f_ode!(ac, env)

            ẋ .= ac.ẋ
            assign!(y_cv, ac) #this also updates y

        end

    end

    ẋ_tmp = similar(x_ref)
    y_tmp = similar(y_ref)

    # @btime $f_nonlinear!($ẋ_tmp, $y_tmp, $x_ref, $u_ref)
    f_nonlinear!(ẋ_tmp, y_tmp, x_ref, u_ref)

    @assert ẋ_tmp ≈ ẋ_ref #sanity check
    @assert y_tmp ≈ y_ref #sanity check


    f_A! = let u = u_ref, y = similar(y_ref) #y is discarded
        (ẋ, x) -> f_nonlinear!(ẋ, y, x, u)
    end

    f_B! = let x = x_ref, y = similar(y_ref) #y is discarded
        (ẋ, u) -> f_nonlinear!(ẋ, y, x, u)
    end

    f_C! = let u = u_ref, ẋ = similar(ẋ_ref) #ẋ is discarded
        (y, x) -> f_nonlinear!(ẋ, y, x, u)
    end

    f_D! = let x = x_ref, ẋ = similar(ẋ_ref) #ẋ is discarded
        (y, u) -> f_nonlinear!(ẋ, y, x, u)
    end

    #preallocate
    A = x_ref * x_ref'
    B = x_ref * u_ref'
    C = y_ref * x_ref'
    D = y_ref * u_ref'

    jacobian!(A, f_A!, x_ref)
    jacobian!(B, f_B!, u_ref)
    jacobian!(C, f_C!, x_ref)
    jacobian!(D, f_D!, u_ref)

    # A_θ = A[:kinematics,:][:pos,:][:e_nb,:][2,:][:kinematics][:pos]

    println("It may be convenient to restrict the linearization functionality to NED kinematics, ",
    "then define also an XTemplate with only the continous state variables that matter, and",
    "ordered conveniently to partition the state in longitudinal and lateral directional dynamics")

    #despues podria crear tambien un LinearAircraft{A, B, C, D, X, U, Y} <:
    #SystemDescriptor.
    #ahi me guardo ẋ_ref, y_ref, etc
    #esto seria generico, no tendria por que ser de C172R

    #maybe define a LinearModel storing ẋ_ref, x_ref, y_ref, u_ref, and the state space
    #matrices (or a state space system directly). the linear model is given by:
    """
    Δẋ = A*Δx + B*Δu
    Δy = C*Δx + D*Δu

    #if we receive ẋ, x, u and y, we need ẋ_ref, x_ref, u_ref and y_ref to compute the
    increments. this is needed for simulation, not for controller design

    """



    return (A, B, C, D)




end



# end #module