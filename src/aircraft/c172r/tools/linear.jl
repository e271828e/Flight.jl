module Linear

using ComponentArrays
using UnPack
using FiniteDiff: finite_difference_jacobian! as jacobian!

using Flight.FlightCore.Systems

using Flight.FlightPhysics.Attitude
using Flight.FlightPhysics.Geodesy
using Flight.FlightPhysics.Kinematics

using Flight.FlightComponents.Environment
using Flight.FlightComponents.Control: LinearStateSpace

using ..C172R
using ..Trim

export linearize!

############################### Linearization ##################################

#flaps and mixture are omitted from the input vector and considered parameters,
#note that they will often be set at one of their range limits (for example,
#when operating with fully retracted or fully extended flaps, full lean or full
#rich mixture), and the resulting linearization would not be valid in such cases
#anyway

const UTemplate = ComponentVector(
    throttle = 0.0,
    elevator = 0.0, #yoke forward
    aileron = 0.0, #yoke right
    rudder = 0.0, #right pedal forward
    )

const YTemplate = ComponentVector(
    ψ = 0.0, θ = 0.0, φ = 0.0, #heading, inclination, bank (body/NED)
    ϕ = 0.0, λ = 0.0, h = 0.0, #latitude, longitude, ellipsoidal altitude
    p = 0.0, q = 0.0, r = 0.0, #angular rates (ω_eb_b)
    TAS = 0.0, α = 0.0, β = 0.0, #airspeed, AoA, AoS
    f_x = 0.0, f_y = 0.0, f_z = 0.0, #specific force at G (f_iG_b)
    ω_eng = 0.0, m_fuel = 0.0 #engine speed, fuel load
)

const U{T, D} = ComponentVector{T, D, typeof(getaxes(UTemplate))} where {T, D}
const Y{T, D} = ComponentVector{T, D, typeof(getaxes(YTemplate))} where {T, D}

function assign!(u::U, ac::System{<:Cessna172R})

    @unpack throttle, elevator_trim, aileron_trim, rudder_trim = ac.u.avionics
    u.throttle = throttle
    u.elevator = elevator_trim
    u.aileron = aileron_trim
    u.rudder = rudder_trim

end

function assign!(ac::System{<:Cessna172R}, u::U)

    @unpack throttle, elevator, aileron, rudder = u
    ac.u.avionics.throttle = throttle
    ac.u.avionics.elevator_trim = elevator
    ac.u.avionics.aileron_trim = aileron
    ac.u.avionics.rudder_trim = rudder

end

function assign!(y::Y, ac::System{<:Cessna172R})

    @unpack q_nb, n_e, h_e, ω_eb_b = ac.y.kinematics
    @unpack α, β = ac.y.airframe.aero
    @unpack ψ, θ, φ = REuler(q_nb)
    @unpack ϕ, λ = LatLon(n_e)

    h = h_e
    p, q, r = ω_eb_b
    f_x, f_y, f_z = ac.y.rigidbody.f_G_b
    TAS = ac.y.air.TAS
    ω_eng = ac.y.airframe.pwp.engine.ω
    m_fuel = ac.y.airframe.fuel.m

    @pack! y = ψ, θ, φ, ϕ, λ, h, p, q, r, TAS, α, β, f_x, f_y, f_z, ω_eng, m_fuel

end

linearize(ac::System{<:Cessna172R{NED}}; kwargs...) = LinearStateSpace(ac; kwargs...)

function LinearStateSpace( ac::System{<:Cessna172R{NED}};
    env::System{<:AbstractEnvironment} = System(SimpleEnvironment()),
    trim_params::Trim.Parameters = Trim.Parameters(),
    trim_state::Trim.State = Trim.State())

    (_, trim_state) = Trim.trim!(ac, env, trim_params, trim_state)

    #save the trimmed aircraft's ẋ, x, u and y for later
    ẋ0_full = copy(ac.ẋ)
    x0_full = copy(ac.x)
    u0 = similar(UTemplate); assign!(u0, ac) #get reference value from trimmed aircraft
    y0 = similar(YTemplate); assign!(y0, ac) #idem

    #function wrapper around f_ode!(), mutates ẋ and y.
    f_nonlinear! = let ac = ac, env = env,
                       params = trim_params, trim_state = trim_state,
                       u_axes = getaxes(u0), y_axes = getaxes(y0)

        function (ẋ, y, x, u)
            # cast y and u into ComponentVectors in case we get generic Vectors
            # from FiniteDiff. these do not allocate, because the underlying
            #data is already in u and y
            u_cv = ComponentVector(u, u_axes)
            y_cv = ComponentVector(y, y_axes)

            #make sure any input or state not set by x and u is at its reference
            #trim value. this reverts any potential changes to the aircraft done
            #by functions sharing the same aircraft instance
            C172R.Trim.assign!(ac, env, trim_params, trim_state)

            assign!(ac, u_cv)
            ac.x .= x
            f_ode!(ac, env)

            ẋ .= ac.ẋ
            assign!(y_cv, ac) #this also updates y (shares its data with y_cv)

        end

    end

    # ẋ_tmp = similar(x0_full)
    # y_tmp = similar(y0)
    # # @btime $f_nonlinear!($ẋ_tmp, $y_tmp, $x0_full, $u0)
    # f_nonlinear!(ẋ_tmp, y_tmp, x0_full, u0)

    # @assert ẋ_tmp ≈ ẋ0_full #sanity check
    # @assert y_tmp ≈ y0 #sanity check

    (A_full, B_full, C_full, D_full) = ss_matrices(f_nonlinear!;
                                            ẋ0 = ẋ0_full, y0, x0 = x0_full, u0)

    #once we're done, make the ac is returned to its trimmed status
    C172R.Trim.assign!(ac, env, trim_params, trim_state)

    #so far we have worked with the nonlinear aircraft model's full state
    #vector, whose layout is given by the hierarchical structure of the
    #aircraft's subsystems, and contains many states irrelevant for an in-flight
    #linearized model. we now replace this complete nonlinear state vector with
    #a reduced one, more suitable for the linear model, containing only relevant
    #states arranged in a flat hierarchy of scalars

    x_labels = (
        ψ = "kinematics.pos.ψ_nb",
        θ = "kinematics.pos.θ_nb",
        φ = "kinematics.pos.φ_nb",
        ϕ = "kinematics.pos.ϕ",
        λ = "kinematics.pos.λ",
        h = "kinematics.pos.h_e",
        p = "kinematics.vel.ω_eb_b[1]",
        q = "kinematics.vel.ω_eb_b[2]",
        r = "kinematics.vel.ω_eb_b[3]",
        v_x = "kinematics.vel.v_eOb_b[1]",
        v_y = "kinematics.vel.v_eOb_b[2]",
        v_z = "kinematics.vel.v_eOb_b[3]",
        α_filt = "airframe.aero.α_filt",
        β_filt = "airframe.aero.β_filt",
        ω_eng = "airframe.pwp.engine.ω",
        fuel = "airframe.fuel[1]",
    )

    #find the indices for the variables in x_labels
    x_indices = [ComponentArrays.label2index(x0_full, s)[1] for s in x_labels]
    #create a new ComponentArrays.Axis with the selected labels
    x_axis = Axis(keys(x_labels))

    #extract the required elements from ẋ0_full and x0_full and rebuild them
    #with the new x_axis
    ẋ0 = ComponentVector(ẋ0_full[x_indices], x_axis)
    x0 = ComponentVector(x0_full[x_indices], x_axis)

    #extract the required rows and columns from A, the required rows from B, and
    #the required columns from C. then rebuild them with the new x_axis and the
    #previous u_axis and y_axis
    A = ComponentMatrix(A_full[x_indices, x_indices], x_axis, x_axis)
    B = ComponentMatrix(B_full[x_indices, :], x_axis, getaxes(u0)[1])
    C = ComponentMatrix(C_full[:, x_indices], getaxes(y0)[1], x_axis)
    D = D_full

    return LinearStateSpace(ẋ0, x0, u0, y0, A, B, C, D)

end


function ss_matrices(f_nonlinear!::Function; ẋ0, y0, x0, u0)

    f_A! = let u = u0, y = similar(y0) #y is discarded
        (ẋ, x) -> f_nonlinear!(ẋ, y, x, u)
    end

    f_B! = let x = x0, y = similar(y0) #y is discarded
        (ẋ, u) -> f_nonlinear!(ẋ, y, x, u)
    end

    f_C! = let u = u0, ẋ = similar(ẋ0) #ẋ is discarded
        (y, x) -> f_nonlinear!(ẋ, y, x, u)
    end

    f_D! = let x = x0, ẋ = similar(ẋ0) #ẋ is discarded
        (y, u) -> f_nonlinear!(ẋ, y, x, u)
    end

    #preallocate
    A = x0 * x0'
    B = x0 * u0'
    C = y0 * x0'
    D = y0 * u0'

    jacobian!(A, f_A!, x0)
    jacobian!(B, f_B!, u0)
    jacobian!(C, f_C!, x0)
    jacobian!(D, f_D!, u0)

    return (A, B, C, D)

end


end #module