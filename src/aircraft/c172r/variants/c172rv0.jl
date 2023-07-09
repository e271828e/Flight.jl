module C172Rv0

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays
using NLopt
using FiniteDiff: finite_difference_jacobian! as jacobian!

using Flight.FlightCore.Systems
using Flight.FlightCore.IODevices
using Flight.FlightCore.Joysticks

using Flight.FlightPhysics.Attitude
using Flight.FlightPhysics.Geodesy
using Flight.FlightPhysics.Kinematics
using Flight.FlightPhysics.Environment

using Flight.FlightComponents.Piston
using Flight.FlightComponents.Aircraft
using Flight.FlightComponents.World

using ..Airframe

export Cessna172Rv0


################################################################################
############################### Cessna172Rv0 #####################################

#Barebones Cessna172R variant, with NoAvionics and manual control over the
#airframe's MechanicalActuation
const Cessna172Rv0{K, F} = AircraftTemplate{K, F, NoAvionics} where {K, F <: C172RAirframe}
Cessna172Rv0(kinematics = LTF()) = AircraftTemplate(kinematics, C172RAirframe(), NoAvionics())

############################ Joystick Mappings #################################

#redirect input assignments directly to the airframe
function IODevices.assign!(sys::System{<:Cessna172Rv0}, joystick::Joystick,
                           mapping::InputMapping)
    IODevices.assign!(sys.airframe, joystick, mapping)
end


############################### Trimming #######################################
################################################################################

#for trimming, the incremental actuation inputs aileron, elevator and rudder are
#set to zero. when a simulation is run from a trimmed state, leaving the
#controller at a neutral position, which means zero aileron, elevator and
#rudder, means that the control surfaces will be left at their trimmed values
#aileron_trim, elevator_trim and rudder_trim

#first 2 are aircraft-agnostic
const TrimStateTemplate = ComponentVector(
    α_a = 0.0, #angle of attack, aerodynamic axes
    φ_nb = 0.0, #bank angle
    n_eng = 0.5, #normalized engine speed (ω/ω_rated)
    throttle = 0.5,
    aileron = 0.0,
    elevator = 0.0,
    rudder = 0.0, #rudder↑ -> aero.u.r↓ -> right yaw
)

const TrimState{T, D} = ComponentVector{T, D, typeof(getaxes(TrimStateTemplate))} where {T, D}

function TrimState(; α_a = 0.1, φ_nb = 0.0, n_eng = 0.75,
    throttle = 0.47, aileron = 0.014, elevator = -0.0015, rudder = 0.02)

    x = copy(TrimStateTemplate)
    @pack! x = α_a, φ_nb, n_eng, throttle, aileron, elevator, rudder
    return x

end

struct TrimParameters
    Ob::Geographic{NVector, Ellipsoidal} #position
    ψ_nb::Float64 #geographic heading
    TAS::Float64 #true airspeed
    γ_wOb_n::Float64 #wind-relative flight path angle
    ψ_lb_dot::Float64 #LTF-relative turn rate
    θ_lb_dot::Float64 #LTF-relative pitch rate
    β_a::Float64 #sideslip angle measured in the aerodynamic reference frame
    fuel::Float64 #fuel load, 0 to 1
    mixture::Float64 #engine mixture control, 0 to 1
    flaps::Float64 #flap setting, 0 to 1
end

function TrimParameters(;
    loc::Abstract2DLocation = LatLon(), h::Altitude = HOrth(1000),
    ψ_nb = 0.0, TAS = 40.0, γ_wOb_n = 0.0, ψ_lb_dot = 0.0, θ_lb_dot = 0.0,
    β_a = 0.0, fuel = 0.5, mixture = 0.5, flaps = 0.0)

    Ob = Geographic(loc, h)
    TrimParameters(Ob, ψ_nb, TAS, γ_wOb_n, ψ_lb_dot, θ_lb_dot, β_a, fuel, mixture, flaps)
end

#given the body-axes wind-relative velocity, the wind-relative flight path angle
#and the bank angle, the pitch angle is unambiguously determined
function θ_constraint(; v_wOb_b, γ_wOb_n, φ_nb)
    TAS = norm(v_wOb_b)
    a = v_wOb_b[1] / TAS
    b = (v_wOb_b[2] * sin(φ_nb) + v_wOb_b[3] * cos(φ_nb)) / TAS
    sγ = sin(γ_wOb_n)

    return atan((a*b + sγ*√(a^2 + b^2 - sγ^2))/(a^2 - sγ^2))
    # return asin((a*sγ + b*√(a^2 + b^2 - sγ^2))/(a^2 + b^2)) #equivalent

end

function Kinematics.Initializer(trim_state::TrimState,
                                trim_params::TrimParameters,
                                env::System{<:AbstractEnvironment})

    @unpack TAS, β_a, γ_wOb_n, ψ_nb, ψ_lb_dot, θ_lb_dot, Ob = trim_params
    @unpack α_a, φ_nb = trim_state

    v_wOb_a = Atmosphere.get_velocity_vector(TAS, α_a, β_a)
    v_wOb_b = C172Rv0.Airframe.f_ba.q(v_wOb_a) #wind-relative aircraft velocity, body frame

    θ_nb = θ_constraint(; v_wOb_b, γ_wOb_n, φ_nb)
    e_nb = REuler(ψ_nb, θ_nb, φ_nb)
    q_nb = RQuat(e_nb)

    e_lb = e_nb #initialize LTF arbitrarily to NED
    ė_lb = SVector(ψ_lb_dot, θ_lb_dot, 0.0)
    ω_lb_b = Attitude.ω(e_lb, ė_lb)

    loc = NVector(Ob)
    h = HEllip(Ob)

    v_wOb_n = q_nb(v_wOb_b) #wind-relative aircraft velocity, NED frame
    v_ew_n = AtmosphericData(env.atm, Ob).wind.v_ew_n
    v_eOb_n = v_ew_n + v_wOb_n

    Kinematics.Initializer(; q_nb, loc, h, ω_lb_b, v_eOb_n, Δx = 0.0, Δy = 0.0)

end

#assigns trim state and parameters to the aircraft system, and then updates it
#by calling f_ode!
function assign!(ac::System{<:Cessna172Rv0},
                env::System{<:AbstractEnvironment},
                trim_params::TrimParameters,
                trim_state::TrimState)

    @unpack TAS, β_a, fuel, flaps, mixture = trim_params
    @unpack n_eng, α_a, throttle, aileron, elevator, rudder = trim_state

    init_kinematics!(ac, Kinematics.Initializer(trim_state, trim_params, env))

    ω_eng = n_eng * ac.airframe.pwp.engine.params.ω_rated

    ac.x.airframe.aero.α_filt = α_a #ensures zero state derivative
    ac.x.airframe.aero.β_filt = β_a #ensures zero state derivative
    ac.x.airframe.pwp.engine.ω = ω_eng
    ac.x.airframe.fuel .= fuel

    ac.u.airframe.act.throttle = throttle
    ac.u.airframe.act.aileron_trim = aileron
    ac.u.airframe.act.elevator_trim = elevator
    ac.u.airframe.act.rudder_trim = rudder
    ac.u.airframe.act.flaps = flaps
    ac.u.airframe.act.mixture = mixture

    #incremental control inputs to zero
    ac.u.airframe.act.elevator = 0
    ac.u.airframe.act.aileron = 0
    ac.u.airframe.act.rudder = 0

    #engine must be running, no way to trim otherwise
    ac.s.airframe.pwp.engine.state = Piston.eng_running
    @assert ac.x.airframe.pwp.engine.ω > ac.airframe.pwp.engine.params.ω_idle

    #engine idle compensator: as long as the engine remains at normal
    #operational speeds, well above its nominal idle speed, the idle controller
    #compensator's output will be saturated at its lower bound by proportional
    #error. its integrator will be disabled, its state will not change nor have
    #any effect on the engine. we can simply set it to zero
    ac.x.airframe.pwp.engine.idle .= 0.0

    #engine friction compensator: with the engine running at normal operational
    #speeds, the engine's friction constraint compensator will be saturated, so
    #its integrator will be disabled and its state will not change. furthermore,
    #with the engine running friction is ignored. we can simply set it to zero.
    ac.x.airframe.pwp.engine.frc .= 0.0

    f_ode!(ac, env)

    #check assumptions concerning airframe systems states & derivatives
    wow = reduce(|, SVector{3,Bool}(leg.strut.wow for leg in ac.y.airframe.ldg))
    @assert wow === false
    @assert ac.ẋ.airframe.pwp.engine.idle[1] .== 0
    @assert ac.ẋ.airframe.pwp.engine.frc[1] .== 0

end

function cost(ac::System{<:Cessna172Rv0})

    v_nd_dot = SVector{3}(ac.ẋ.kinematics.vel.v_eOb_b) / norm(ac.y.kinematics.common.v_eOb_b)
    ω_dot = SVector{3}(ac.ẋ.kinematics.vel.ω_eb_b) #ω should already of order 1
    n_eng_dot = ac.ẋ.airframe.pwp.engine.ω / ac.airframe.pwp.engine.params.ω_rated

    sum(v_nd_dot.^2) + sum(ω_dot.^2) + n_eng_dot^2

end

function get_target_function(ac::System{<:Cessna172Rv0},
    env::System{<:AbstractEnvironment}, trim_params::TrimParameters = TrimParameters())

    let ac = ac, env = env, trim_params = trim_params
        function (x::TrimState)
            assign!(ac, env, trim_params, x)
            return cost(ac)
        end
    end

end

function Aircraft.trim!( ac::System{<:Cessna172Rv0};
                env::System{<:AbstractEnvironment} = System(SimpleEnvironment()),
                trim_params::TrimParameters = TrimParameters())

    f_target = get_target_function(ac, env, trim_params)

    #wrapper around f_target with the interface required by NLopt
    trim_state = TrimState() #could define this initial condition as an optional input
    ax = getaxes(trim_state)
    function f_opt(x::Vector{Float64}, ::Vector{Float64})
        s = ComponentVector(x, ax)
        return f_target(s)
    end

    n = length(trim_state)
    x0 = zeros(n); lower_bounds = similar(x0); upper_bounds = similar(x0); initial_step = similar(x0)

    x0[:] .= trim_state

    lower_bounds[:] .= TrimState(
        α_a = -π/12,
        φ_nb = -π/3,
        n_eng = 0.4,
        throttle = 0,
        aileron = -1,
        elevator = -1,
        rudder = -1)

    upper_bounds[:] .= TrimState(
        α_a = ac.airframe.aero.params.α_stall[2], #critical AoA is 0.28 < 0.36
        φ_nb = π/3,
        n_eng = 1.1,
        throttle = 1,
        aileron = 1,
        elevator = 1,
        rudder = 1)

    initial_step[:] .= 0.05 #safe value for all optimization variables

    # @show initial_cost = f_opt(x0, x0)
    # @btime $f_opt($x0, $x0)

    #any of these three algorithms works
    # opt = Opt(:LN_NELDERMEAD, length(x0))
    opt = Opt(:LN_BOBYQA, length(x0))
    # opt = Opt(:GN_CRS2_LM, length(x0))
    opt.min_objective = f_opt
    opt.maxeval = 100000
    opt.stopval = 1e-14
    opt.lower_bounds = lower_bounds
    opt.upper_bounds = upper_bounds
    opt.initial_step = initial_step

    # @btime optimize($opt, $x0)

    (minf, minx, exit_flag) = optimize(opt, x0)

    success = (exit_flag === :STOPVAL_REACHED)
    if !success
        println("Warning: Optimization failed with exit_flag $exit_flag")
    end
    trim_state_opt = ComponentVector(minx, ax)
    assign!(ac, env, trim_params, trim_state_opt)
    return (success = success, result = trim_state_opt)


end

function Aircraft.trim!(
    world::System{<:SimpleWorld{<:Cessna172Rv0, <:AbstractEnvironment}};
    trim_params::TrimParameters = TrimParameters())

    trim!(world.ac; env = world.env, trim_params = trim_params)

end


################################################################################
############################### Linearization ##################################

#flaps and mixture are omitted from the control input vector and considered
#parameters instead,

const LinearUTemplate = ComponentVector(
    throttle = 0.0,
    aileron = 0.0,
    elevator = 0.0,
    rudder = 0.0,
    )

const LinearXTemplate = ComponentVector(
    ψ = 0.0, θ = 0.0, φ = 0.0, #heading, inclination, bank (body/NED)
    ϕ = 0.0, λ = 0.0, h = 0.0, #latitude, longitude, ellipsoidal altitude
    p = 0.0, q = 0.0, r = 0.0, #angular rates (ω_eb_b)
    v_x = 0.0, v_y = 0.0, v_z = 0.0, #Ob/ECEF velocity, body axes
    α_filt = 0.0, β_filt = 0.0, #filtered airflow angles
    ω_eng = 0.0, fuel = 0.0, #engine speed, fuel fraction
)

const LinearYTemplate = ComponentVector(
    ψ = 0.0, θ = 0.0, φ = 0.0, #heading, inclination, bank (body/NED)
    ϕ = 0.0, λ = 0.0, h = 0.0, #latitude, longitude, ellipsoidal altitude
    p = 0.0, q = 0.0, r = 0.0, #angular rates (ω_eb_b)
    TAS = 0.0, α = 0.0, β = 0.0, #airspeed, AoA, AoS
    f_x = 0.0, f_y = 0.0, f_z = 0.0, #specific force at G (f_iG_b)
    ω_eng = 0.0, m_fuel = 0.0 #engine speed, fuel mass
)

const LinearU{T, D} = ComponentVector{T, D, typeof(getaxes(LinearUTemplate))} where {T, D}
const LinearX{T, D} = ComponentVector{T, D, typeof(getaxes(LinearXTemplate))} where {T, D}
const LinearY{T, D} = ComponentVector{T, D, typeof(getaxes(LinearYTemplate))} where {T, D}

#access labels for LinearX components within the overall Cessna172Rv0{NED} state vector
const x_labels = (
        "kinematics.pos.ψ_nb", "kinematics.pos.θ_nb", "kinematics.pos.φ_nb",
        "kinematics.pos.ϕ", "kinematics.pos.λ", "kinematics.pos.h_e",
        "kinematics.vel.ω_eb_b[1]", "kinematics.vel.ω_eb_b[2]", "kinematics.vel.ω_eb_b[3]",
        "kinematics.vel.v_eOb_b[1]", "kinematics.vel.v_eOb_b[2]", "kinematics.vel.v_eOb_b[3]",
        "airframe.aero.α_filt", "airframe.aero.β_filt",
        "airframe.pwp.engine.ω", "airframe.fuel[1]",
    )

function assign!(u::LinearU, ac::System{<:Cessna172Rv0})

    @unpack throttle, aileron, elevator, rudder = ac.u.airframe.act
    @pack! u = throttle, aileron, elevator, rudder

end

function assign!(ac::System{<:Cessna172Rv0}, u::LinearU)

    @unpack throttle, aileron, elevator, rudder = u
    @pack! ac.u.airframe.act = throttle, aileron, elevator, rudder

end

function assign!(y::LinearY, ac::System{<:Cessna172Rv0})

    @unpack q_nb, n_e, h_e, ω_eb_b = ac.y.kinematics
    @unpack α, β = ac.y.airframe.aero
    @unpack ψ, θ, φ = REuler(q_nb)
    @unpack ϕ, λ = LatLon(n_e)

    h = h_e
    p, q, r = ω_eb_b
    f_x, f_y, f_z = ac.y.rigidbody.f_G_b
    TAS = ac.y.air.TAS
    ω_eng = ac.y.airframe.pwp.engine.ω
    m_fuel = ac.y.airframe.fuel.m_avail

    @pack! y = ψ, θ, φ, ϕ, λ, h, p, q, r, TAS, α, β, f_x, f_y, f_z, ω_eng, m_fuel

end

function Aircraft.linearize!(ac::System{<:Cessna172Rv0{NED}};
    env::System{<:AbstractEnvironment} = System(SimpleEnvironment()),
    trim_params::TrimParameters = TrimParameters())

    (_, trim_state) = trim!(ac; env, trim_params)

    #save the trimmed aircraft's ẋ, x, u and y for later
    ẋ0_full = copy(ac.ẋ)
    x0_full = copy(ac.x)
    u0 = similar(LinearUTemplate); assign!(u0, ac) #get reference value from trimmed aircraft
    y0 = similar(LinearYTemplate); assign!(y0, ac) #idem

    #function wrapper around f_ode!(), mutates ẋ and y.
    f_nonlinear! = let ac = ac, env = env,
                       trim_params = trim_params, trim_state = trim_state,
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
            assign!(ac, env, trim_params, trim_state)

            assign!(ac, u_cv)
            ac.x .= x
            f_ode!(ac, env)

            ẋ .= ac.ẋ
            assign!(y_cv, ac) #this also updates y (shares its data with y_cv)

        end

    end

    (A_full, B_full, C_full, D_full) = ss_matrices(f_nonlinear!;
                                            ẋ0 = ẋ0_full, y0, x0 = x0_full, u0)

    #once we're done, ensure the aircraft is restored to its trimmed status, so
    #the response can be compared with that of its linear counterpart
    assign!(ac, env, trim_params, trim_state)

    #find the indices for the components in the reduced LinearX state vector
    #within the overall aircraft state vector
    x_indices = [ComponentArrays.label2index(x0_full, s)[1] for s in x_labels]

    x_axis = getaxes(LinearXTemplate)[1]

    #extract the required elements from ẋ0_full and x0_full and rebuild them
    #with the LinearX axis
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