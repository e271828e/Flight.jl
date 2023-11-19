module Aircraft

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays
using FiniteDiff: finite_difference_jacobian! as jacobian!

using Flight.FlightCore
using Flight.FlightPhysics

using ..Control

export AbstractAirframe, EmptyAirframe
export AbstractAvionics, NoAvionics
export AbstractTrimParameters, AbstractTrimState
export init_kinematics!, trim!, linearize!

################################################################################
########################### AbstractAirframe ###################################

abstract type AbstractAirframe <: SystemDefinition end
RigidBody.MassTrait(::System{<:AbstractAirframe}) = HasMass()
RigidBody.AngMomTrait(::System{<:AbstractAirframe}) = HasAngularMomentum()
RigidBody.WrenchTrait(::System{<:AbstractAirframe}) = GetsExternalWrench()

################################ EmptyAirframe #################################

@kwdef struct EmptyAirframe <: AbstractAirframe
    mass_distribution::RigidBodyDistribution = RigidBodyDistribution(1, SA[1.0 0 0; 0 1.0 0; 0 0 1.0])
end

RigidBody.AngMomTrait(::System{EmptyAirframe}) = HasNoAngularMomentum()
RigidBody.WrenchTrait(::System{EmptyAirframe}) = GetsNoExternalWrench()
RigidBody.get_mp_Ob(sys::System{EmptyAirframe}) = MassProperties(sys.constants.mass_distribution)

################################################################################
############################## Aircraft Physics ################################

@kwdef struct Physics{F <: AbstractAirframe,
                      K <: AbstractKinematicDescriptor,
                      T <: AbstractTerrain} <: SystemDefinition
    airframe::F = EmptyAirframe()
    kinematics::K = LTF()
    terrain::T = HorizontalTerrain() #shared with other aircraft instances
    atmosphere::LocalAtmosphere = LocalAtmosphere() #externally controlled
end

struct PhysicsY{F, K}
    airframe::F
    kinematics::K
    rigidbody::RigidBodyData
    air::AirData
end

Systems.init(::SystemY, ac::Physics) = PhysicsY(
    init_y(ac.airframe),
    init_y(ac.kinematics),
    RigidBodyData(),
    AirData())

function init_kinematics!(sys::System{<:Physics}, ic::KinematicInit)
    Kinematics.init!(sys.x.kinematics, ic)
end

###############################################################################
############################# AbstractAvionics #################################

abstract type AbstractAvionics <: SystemDefinition end

################################### NoAvionics #################################

struct NoAvionics <: AbstractAvionics end


###############################################################################
#################### AbstractAirframe update methods ###########################

function Systems.f_ode!(airframe::System{<:AbstractAirframe},
                        kin::KinematicData,
                        air::AirData,
                        trn::AbstractTerrain)
    MethodError(f_ode!, (airframe, avionics, kin, air, trn)) |> throw
end

function Systems.f_ode!(::System{EmptyAirframe},
                        ::KinematicData,
                        ::AirData,
                        ::AbstractTerrain)
    nothing
end

#this method can be extended if required, but in principle Airframe shouldn't
#implement discrete dynamics; discretized algorithms belong in Avionics
function Systems.f_disc!(::System{<:AbstractAirframe}, ::Real)
    return false
end

#f_step! may use the recursive fallback implementation

###############################################################################
#################### AbstractAvionics update methods ###########################

#avionics update methods should only mutate the avionics System

#this method can be extended if required, but in principle avionics shouldn't
#involve continuous dynamics.
function Systems.f_ode!(::System{<:AbstractAvionics},
                        ::System{<:Physics})
    nothing
end

function Systems.f_disc!(avionics::System{<:AbstractAvionics},
                        physics::System{<:Physics},
                        Δt::Real)
    MethodError(f_disc!, (avionics, physics, Δt)) |> throw
end

function Systems.f_disc!(::System{NoAvionics},
                        ::System{<:Physics},
                        ::Real)
    return false
end

#f_step! can use the recursive fallback implementation


################################################################################
###################### AircraftPhysics Update methods ##########################

function Systems.f_ode!(physics::System{<:Physics})

    @unpack ẋ, x, subsystems, constants = physics
    @unpack kinematics, airframe, atmosphere = subsystems
    @unpack terrain = constants

    f_ode!(kinematics)
    f_ode!(atmosphere) #currently does nothing
    kin_data = KinematicData(kinematics)
    atm_data = LocalAtmosphericData(atmosphere)

    air_data = AirData(kin_data, atm_data)

    #update airframe
    f_ode!(airframe, kin_data, air_data, terrain)

    #get inputs for rigid body dynamics
    mp_Ob = get_mp_Ob(airframe)
    wr_b = get_wr_b(airframe)
    hr_b = get_hr_b(airframe)

    #update velocity derivatives and rigid body data
    rb_data = f_rigidbody!(kinematics.ẋ.vel, kin_data, mp_Ob, wr_b, hr_b)

    physics.y = PhysicsY(airframe.y, kinematics.y, rb_data, air_data)

    return nothing

end

#f_step! will use the recursive fallback implementation

#within Physics, only the airframe may be modified by f_disc! (and it generally
#shouldn't)
function Systems.f_disc!(physics::System{<:Physics}, Δt::Real)

    @unpack kinematics, airframe = physics
    @unpack rigidbody, air = physics.y

    x_mod = false
    x_mod |= f_disc!(airframe, Δt)

    #airframe might have modified its outputs, so we need to reassemble our y
    physics.y = PhysicsY(airframe.y, kinematics.y, rigidbody, air)

    return x_mod
end

#these are meant to map avionics outputs to airframe control inputs, they are
#called both within the aircraft's f_ode! and f_disc! before the physics update
function assign!(physics::System{<:Physics}, avionics::System{<:AbstractAvionics})
    assign!(physics.airframe, avionics)
end

function assign!(airframe::System{<:AbstractAirframe},
                avionics::System{<:AbstractAvionics})
    MethodError(assign!, (airframe, avionics)) |> throw
end

assign!(::System{<:AbstractAirframe}, ::System{NoAvionics}) = nothing


################################################################################
################################## Template ####################################

@kwdef struct Template{P <: Physics, A <: AbstractAvionics} <: SystemDefinition
    physics::P = Physics()
    avionics::A = NoAvionics()
end

struct TemplateY{P <: PhysicsY, A}
    physics::P
    avionics::A
end

Systems.init(::SystemY, ac::Template) = TemplateY(init_y(ac.physics), init_y(ac.avionics))

init_kinematics!(ac::System{<:Template}, ic::KinematicInit) = init_kinematics!(ac.physics, ic)

function Systems.f_ode!(ac::System{<:Template})

    @unpack physics, avionics = ac.subsystems

    f_ode!(avionics, physics)
    assign!(physics, avionics)
    f_ode!(physics)

    ac.y = TemplateY(physics.y, avionics.y)

    return nothing

end

function Systems.f_disc!(ac::System{<:Template}, Δt::Real)

    @unpack physics, avionics = ac.subsystems

    x_mod = false
    x_mod |= f_disc!(avionics, physics, Δt)
    assign!(physics, avionics)
    x_mod |= f_disc!(physics, Δt)

    ac.y = TemplateY(physics.y, avionics.y)

    return x_mod
end

#f_step! can use the recursive fallback implementation


############################# XPlaneConnect ####################################

Visualization.set_position!(xpc::XPCDevice, y::TemplateY) = Visualization.set_position!(xpc, y.physics)

function Visualization.set_position!(xpc::XPCDevice, y::PhysicsY)

    aircraft = 0

    @unpack ϕ_λ, e_nb, h_o = y.kinematics

    lat = rad2deg(ϕ_λ.ϕ)
    lon = rad2deg(ϕ_λ.λ)

    psi = rad2deg(e_nb.ψ)
    theta = rad2deg(e_nb.θ)
    phi = rad2deg(e_nb.φ)

    Visualization.set_position!(xpc; lat, lon, h_o, psi, theta, phi, aircraft)

end


################################### Trimming ###################################
################################################################################

abstract type AbstractTrimParameters end
const AbstractTrimState{N} = FieldVector{N, Float64}

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

trim!( ac::System{<:Template}, args...; kwargs...) = trim!(ac.physics, args...; kwargs...)

function trim!( physics::System{<:Physics}, args...)
    MethodError(trim!, (physics, args...)) |> throw
end

function assign!(::System{<:Physics}, ::AbstractTrimParameters, ::AbstractTrimState)
    error("An assign! method must be defined by each Aircraft.Physics subtype")
end


################################################################################
################################ Linearization #################################

ẋ_linear(physics::System{<:Physics})::FieldVector = throw(MethodError(ẋ_linear!, (physics,)))
x_linear(physics::System{<:Physics})::FieldVector = throw(MethodError(x_linear!, (physics,)))
u_linear(physics::System{<:Physics})::FieldVector = throw(MethodError(u_linear!, (physics,)))
y_linear(physics::System{<:Physics})::FieldVector = throw(MethodError(y_linear!, (physics,)))

assign_x!(physics::System{<:Physics}, x::AbstractVector{Float64}) = throw(MethodError(assign_x!, (physics, x)))
assign_u!(physics::System{<:Physics}, u::AbstractVector{Float64}) = throw(MethodError(assign_u!, (physics, u)))

linearize!(ac::System{<:Aircraft.Template}, args...) = linearize!(ac.physics, args...)

function Aircraft.linearize!( physics::System{<:Aircraft.Physics},
                            trim_params::AbstractTrimParameters)

    (_, trim_state) = trim!(physics, trim_params)

    ẋ0 = ẋ_linear(physics)::FieldVector
    x0 = x_linear(physics)::FieldVector
    u0 = u_linear(physics)::FieldVector
    y0 = y_linear(physics)::FieldVector

    #f_main will not be returned for use in another scope, so we don't need to
    #capture physics with a let block, because they are guaranteed not be
    #reassigned within the scope of linearize!
    function f_main(x, u)

        assign_x!(physics, x)
        assign_u!(physics, u)
        f_ode!(physics)

        return (ẋ = ẋ_linear(physics), y = y_linear(physics))

    end

    (A, B, C, D) = ss_matrices(f_main, x0, u0)

    #restore the System to its trimmed condition
    assign!(physics, trim_params, trim_state)

    #now we need to rebuild vectors and matrices for the LinearStateSpace as
    #ComponentArrays, because we want matrix components to remain labelled,
    #which cannot be achieved with FieldVectors

    x_axis = Axis(propertynames(x0))
    u_axis = Axis(propertynames(u0))
    y_axis = Axis(propertynames(y0))

    ẋ0_cv = ComponentVector(ẋ0, x_axis)
    x0_cv = ComponentVector(x0, x_axis)
    u0_cv = ComponentVector(u0, u_axis)
    y0_cv = ComponentVector(y0, y_axis)

    A_cv = ComponentMatrix(A, x_axis, x_axis)
    B_cv = ComponentMatrix(B, x_axis, u_axis)
    C_cv = ComponentMatrix(C, y_axis, x_axis)
    D_cv = ComponentMatrix(D, y_axis, u_axis)

    return LinearStateSpace(ẋ0_cv, x0_cv, u0_cv, y0_cv, A_cv, B_cv, C_cv, D_cv)

end

#given a trim point (x0, u0) for the nonlinear system (ẋ, y) = f(x,u), compute
#the state space matrices (A, B, C, D) for the system's linearization around
#(x0, u0), given by: Δẋ = AΔx + BΔu; Δy = CΔx + DΔu
function ss_matrices(f_main::Function, x0::AbstractVector{Float64},
                                       u0::AbstractVector{Float64})

    f_ẋ(x, u) = f_main(x, u).ẋ
    f_y(x, u) = f_main(x, u).y
    y0 = f_y(x0, u0)

    #none of these closures will be returned for use in another scope, so we
    #don't need to capture x0 and u0 with a let block, because they are
    #guaranteed not be reassigned within the scope of ss_matrices
    f_A!(ẋ, x) = (ẋ .= f_ẋ(x, u0))
    f_B!(ẋ, u) = (ẋ .= f_ẋ(x0, u))
    f_C!(y, x) = (y .= f_y(x, u0))
    f_D!(y, u) = (y .= f_y(x0, u))

    #preallocate mutable arrays
    A = (x0 * x0') |> Matrix
    B = (x0 * u0') |> Matrix
    C = (y0 * x0') |> Matrix
    D = (y0 * u0') |> Matrix

    jacobian!(A, f_A!, Vector(x0))
    jacobian!(B, f_B!, Vector(u0))
    jacobian!(C, f_C!, Vector(x0))
    jacobian!(D, f_D!, Vector(u0))

    return (A, B, C, D)

end


function Control.LinearStateSpace(
            ac::System{<:Aircraft.Template}, args...; kwargs...)
    LinearStateSpace(ac.physics, args...; kwargs...)
end

############################### Plotting #######################################

function Plotting.make_plots(th::TimeHistory{<:PhysicsY}; kwargs...)

    return OrderedDict(
        :airframe => make_plots(th.airframe; kwargs...),
        :kinematics => make_plots(th.kinematics; kwargs...),
        :rigidbody => make_plots(th.rigidbody; kwargs...),
        :air => make_plots(th.air; kwargs...),
    )

end

function Plotting.make_plots(th::TimeHistory{<:TemplateY}; kwargs...)

    return OrderedDict(
        :physics => make_plots(th.physics; kwargs...),
        :avionics => make_plots(th.avionics; kwargs...),
    )

end

################################### GUI ########################################


function GUI.draw!(sys::System{<:Template}, label::String = "Aircraft")

    @unpack y = sys

    CImGui.Begin(label)

    show_physics = @cstatic check=false @c CImGui.Checkbox("Physics", &check)
    show_avionics = @cstatic check=false @c CImGui.Checkbox("Avionics", &check)

    show_physics && GUI.draw!(sys.physics, sys.avionics)
    show_avionics && GUI.draw!(sys.avionics, sys.physics)

    CImGui.End()

end

function GUI.draw!(physics::System{<:Physics},
                   avionics::System{<:AbstractAvionics},
                   label::String = "Aircraft Physics")

    @unpack airframe, atmosphere = physics.subsystems
    @unpack terrain = physics.constants
    @unpack kinematics, rigidbody, air = physics.y

    CImGui.Begin(label)

    show_airframe = @cstatic check=false @c CImGui.Checkbox("Airframe", &check)
    show_atmosphere = @cstatic check=false @c CImGui.Checkbox("Local Atmosphere", &check)
    show_terrain = @cstatic check=false @c CImGui.Checkbox("Terrain", &check)
    show_dyn = @cstatic check=false @c CImGui.Checkbox("Dynamics", &check)
    show_kin = @cstatic check=false @c CImGui.Checkbox("Kinematics", &check)
    show_air = @cstatic check=false @c CImGui.Checkbox("Air", &check)

    show_airframe && GUI.draw!(physics.airframe, avionics)
    show_atmosphere && GUI.draw!(physics.atmosphere)
    show_terrain && GUI.draw(terrain)
    show_dyn && GUI.draw(rigidbody, "Dynamics")
    show_kin && GUI.draw(kinematics, "Kinematics")
    show_air && GUI.draw(air, "Air")

    CImGui.End()

end

end #module