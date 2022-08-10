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

UTemplate = ComponentVector(
    throttle = 0.0,
    yoke_y = 0.0, #elevator control
    yoke_x = 0.0, #aileron control
    pedals = 0.0, #rudder control
    )

YTemplate = ComponentVector(
    ψ = 0.0, θ = 0.0, φ = 0.0, #heading, inclination, bank
    ϕ = 0.0, λ = 0.0, h_e = 0.0, #latitude, longitude, altitude
    p = 0.0, q = 0.0, r = 0.0, #angular rates
    TAS = 0.0, α = 0.0, β = 0.0, #airspeed, AoA, AoS
    f_x = 0.0, f_y = 0.0, f_z = 0.0, #specific force
    ω_eng = 0.0, m_fuel = 0.0 #engine speed, fuel load
)

const U{T, D} = ComponentVector{T, D, typeof(getaxes(UTemplate))} where {T, D}
const Y{T, D} = ComponentVector{T, D, typeof(getaxes(YTemplate))} where {T, D}

function assign!(u::U, ac::System{<:Cessna172R})

    @unpack throttle, yoke_y, yoke_x, pedals = ac.u.avionics
    @pack! u = throttle, yoke_y, yoke_x, pedals

end

function assign!(ac::System{<:Cessna172R}, u::U)

    @unpack throttle, yoke_y, yoke_x, pedals = u
    @pack! ac.u.avionics = throttle, yoke_y, yoke_x, pedals

end

function assign!(y::Y, ac::System{<:Cessna172R})

    @unpack q_nb, n_e, h_e, ω_lb_b = ac.y.kinematics
    @unpack α, β = ac.y.airframe.aero
    @unpack ψ, θ, φ = REuler(q_nb)
    @unpack ϕ, λ = LatLon(n_e)

    p, q, r = ω_lb_b
    f_x, f_y, f_z = ac.y.rigidbody.f_Ob_b
    TAS = ac.y.airflow.TAS
    ω_eng = ac.y.airframe.pwp.engine.ω
    m_fuel = ac.y.airframe.fuel.m

    @pack! y = ψ, θ, φ, ϕ, λ, h_e, p, q, r, TAS, α, β, f_x, f_y, f_z, ω_eng, m_fuel

end

function test()

    ac = System(Cessna172R(NED()))
    env = System(SimpleEnvironment())
    params = C172R.Trim.Parameters()

    # state = C172R.Trim.State() #optional initial trim guess
    (exit_flag, trim_state) = C172R.Trim.trim!(; ac, env, params) #ac is now trimmed

    #save the trimmed aircraft ẋ, x, u and y for later
    ẋ_ref = copy(ac.ẋ) #save the trimmed aircraft state vector derivative for later
    x_ref = copy(ac.x) #
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

            ac.x .= x
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