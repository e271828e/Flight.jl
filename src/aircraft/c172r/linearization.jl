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

const UTemplate = ComponentVector(
    throttle = 0.0,
    yoke_y = 0.0, #elevator control
    yoke_x = 0.0, #aileron control
    pedals = 0.0, #rudder control
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

    @unpack throttle, yoke_y, yoke_x, pedals = ac.u.avionics
    @pack! u = throttle, yoke_y, yoke_x, pedals

end

function assign!(ac::System{<:Cessna172R}, u::U)

    @unpack throttle, yoke_y, yoke_x, pedals = u
    @pack! ac.u.avionics = throttle, yoke_y, yoke_x, pedals

end

function assign!(y::Y, ac::System{<:Cessna172R})

    @unpack q_nb, n_e, h_e, ω_eb_b = ac.y.kinematics
    @unpack α, β = ac.y.airframe.aero
    @unpack ψ, θ, φ = REuler(q_nb)
    @unpack ϕ, λ = LatLon(n_e)

    h = h_e
    p, q, r = ω_eb_b
    f_x, f_y, f_z = ac.y.rigidbody.f_G_b
    TAS = ac.y.airflow.TAS
    ω_eng = ac.y.airframe.pwp.engine.ω
    m_fuel = ac.y.airframe.fuel.m

    @pack! y = ψ, θ, φ, ϕ, λ, h, p, q, r, TAS, α, β, f_x, f_y, f_z, ω_eng, m_fuel

end

function linearize!(; ac::System{<:Cessna172R{NED}},
    env::System{<:AbstractEnvironment} = System(SimpleEnvironment()),
    params::C172R.Trim.Parameters = C172R.Trim.Parameters(),
    state::C172R.Trim.State = C172R.Trim.State())

    #trim the aircraft
    (exit_flag, trim_state) = C172R.Trim.trim!(; ac, env, params, state)

    #save the trimmed aircraft ẋ, x, u and y for later
    ẋ_ref = copy(ac.ẋ) #save the trimmed aircraft state vector derivative for later
    x_ref = copy(ac.x) #
    u_ref = similar(UTemplate); assign!(u_ref, ac) #get the reference value from the trimmed aircraft
    y_ref = similar(YTemplate); assign!(y_ref, ac) #idem

    #function wrapper around f_ode!() mutating ẋ and y.
    f_nonlinear! = let ac = ac, env = env, params = params, state = trim_state,
                       u_axes = getaxes(u_ref), y_axes = getaxes(y_ref)

        function (ẋ, y, x, u)
            # cast y and u into ComponentVectors in case we get generic Vectors
            # from FiniteDiff. these do not allocate, because the underlying
            #data is already in u and y
            u_cv = ComponentVector(u, u_axes)
            y_cv = ComponentVector(y, y_axes)

            #make sure any input or state not set by x and u is at its reference
            #trim value. this reverts any potential changes to the aircraft done
            #by functions sharing the same aircraft instance
            C172R.Trim.assign!(ac, env, params, state)

            assign!(ac, u_cv)
            ac.x .= x
            f_ode!(ac, env)
            # @show ac.y.rigidbody.f_G_b

            ẋ .= ac.ẋ
            assign!(y_cv, ac) #this also updates y (shares its data with y_cv)

        end

    end

    # fA, fB, fC, fD = get_functions(f_nonlinear!; ẋ_ref, y_ref, x_ref, u_ref )

    # u_tmp = copy(u_ref)
    # y_tmp = copy(y_ref)

    # Δyoke_y = 0.01
    # fD(y_tmp, u_ref)
    # @show y_tmp.f_x
    # u_tmp.yoke_y += Δyoke_y
    # fD(y_tmp, u_tmp)
    # @show y_tmp.f_x
    # @show (y_tmp - y_ref).f_x /Δyoke_y

    # C172R.Trim.assign!(ac, env, params, trim_state)
    # @show f_x_ref = ac.y.rigidbody.f_G_b[1]
    # ac.u.avionics.yoke_y += Δyoke_y
    # f_ode!(ac, env)
    # @show f_x = ac.y.rigidbody.f_G_b[1]
    # @show (f_x - f_x_ref) / Δyoke_y

    # return

    (A, B, C, D) = ss_matrices(f_nonlinear!; ẋ_ref, y_ref, x_ref, u_ref )

    x_lin_labels = (
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
        m_fuel = "airframe.fuel[1]",
    )

    #first need to these is a reduced and reordered set of numeric indices,
    #forming a reduced state vector, suitable for the linearized model. we need
    #these to extract the desired rows and columns from A, the corresponding
    #rows from B, and the corresponding columns from C. only after we do this,
    #we create a new XAxis and rebuild the affected matrices (A, B, C)
    x_lin_indices = [ComponentArrays.label2index(x_ref, s)[1] for s in x_lin_labels]
    x_axis_lin = Axis(keys(x_lin_labels))

    x_ref_lin = ComponentVector(x_ref[x_lin_indices], x_axis_lin)
    A_lin = ComponentMatrix(A[x_lin_indices, x_lin_indices], x_axis_lin, x_axis_lin)
    B_lin = ComponentMatrix(B[x_lin_indices, :], x_axis_lin, getaxes(u_ref)[1])
    C_lin = ComponentMatrix(C[:, x_lin_indices], getaxes(y_ref)[1], x_axis_lin)
    D_lin = D

    return A_lin, B_lin, C_lin, D_lin
end


function ss_matrices(f_nonlinear!::Function; ẋ_ref, y_ref, x_ref, u_ref)

    # ẋ_tmp = similar(x_ref)
    # y_tmp = similar(y_ref)

    # # @btime $f_nonlinear!($ẋ_tmp, $y_tmp, $x_ref, $u_ref)
    # f_nonlinear!(ẋ_tmp, y_tmp, x_ref, u_ref)

    # @assert ẋ_tmp ≈ ẋ_ref #sanity check
    # @assert y_tmp ≈ y_ref #sanity check


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

    return (A, B, C, D)

end
# function get_functions(f_nonlinear!::Function; ẋ_ref, y_ref, x_ref, u_ref)

#     f_A! = let u = u_ref, y = similar(y_ref) #y is discarded
#         (ẋ, x) -> f_nonlinear!(ẋ, y, x, u)
#     end

#     f_B! = let x = x_ref, y = similar(y_ref) #y is discarded
#         (ẋ, u) -> f_nonlinear!(ẋ, y, x, u)
#     end

#     f_C! = let u = u_ref, ẋ = similar(ẋ_ref) #ẋ is discarded
#         (y, x) -> f_nonlinear!(ẋ, y, x, u)
#     end

#     f_D! = let x = x_ref, ẋ = similar(ẋ_ref) #ẋ is discarded
#         (y, u) -> f_nonlinear!(ẋ, y, x, u)
#     end

#     return f_A!, f_B!, f_C!, f_D!

# end

# 2) realmente necesito definirme structs para U y para Y del
#    System{LinearDynamics}? No: U es mutable de por si, asi que no necesito
#    definir una struct para el. Y ahora que he descubierto que puedo usar
#    ComponentVectors con SVectors, tampoco necesito hacerlo para Y, porque un
#    ComponentVector creado a partir de un SVector no genera allocations. Asi que
#    me basta con los templates que tenia antes

"""

    #despues podria crear tambien un LinearAircraft{A, B, C, D, X, U, Y} <:
    #SystemDescriptor.
    #ahi me guardo ẋ_ref, y_ref, etc
    #esto seria generico, no tendria por que ser de C172R

    #maybe define a LinearModel storing ẋ_ref, x_ref, y_ref, u_ref, and the state space
    #matrices (or a state space system directly). the linear model is given by:
    # Δẋ = A*Δx + B*Δu
    # Δy = C*Δx + D*Δu

    #if we receive ẋ, x, u and y, we need ẋ_ref, x_ref, u_ref and y_ref to compute the
    # increments. this is needed for simulation, not for controller design


    """




# end #module