module AircraftBase

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays
using FiniteDiff: finite_difference_jacobian! as jacobian!

using Flight.FlightCore
using Flight.FlightPhysics
using Flight.FlightComponents

export AbstractPlatform, EmptyPlatform
export AbstractAvionics, NoAvionics
export AbstractTrimParameters, AbstractTrimState
export trim!, linearize!

################################################################################
########################### AbstractPlatform ###################################

abstract type AbstractPlatform <: SystemDefinition end
Dynamics.MassTrait(::System{<:AbstractPlatform}) = HasMass()
Dynamics.AngularMomentumTrait(::System{<:AbstractPlatform}) = HasAngularMomentum()
Dynamics.ExternalWrenchTrait(::System{<:AbstractPlatform}) = GetsExternalWrench()

################################ EmptyPlatform #################################

@kwdef struct EmptyPlatform <: AbstractPlatform
    mass_distribution::RigidBodyDistribution = RigidBodyDistribution(1, SA[1.0 0 0; 0 1.0 0; 0 0 1.0])
end

Dynamics.AngularMomentumTrait(::System{EmptyPlatform}) = HasNoAngularMomentum()
Dynamics.ExternalWrenchTrait(::System{EmptyPlatform}) = GetsNoExternalWrench()
Dynamics.get_mp_Ob(sys::System{EmptyPlatform}) = MassProperties(sys.constants.mass_distribution)

################################################################################
############################## Aircraft Physics ################################

@kwdef struct Physics{F <: AbstractPlatform,
                      K <: AbstractKinematicDescriptor,
                      T <: AbstractTerrain} <: SystemDefinition
    platform::F = EmptyPlatform()
    kinematics::K = LTF()
    terrain::T = HorizontalTerrain() #shared with other aircraft instances
    atmosphere::LocalAtmosphere = LocalAtmosphere() #externally controlled
end

struct PhysicsY{F, K}
    platform::F
    kinematics::K
    dynamics::DynamicsData
    air::AirData
end

Systems.Y(ac::Physics) = PhysicsY(
    Systems.Y(ac.platform),
    Systems.Y(ac.kinematics),
    DynamicsData(),
    AirData())

function Systems.init!(sys::System{<:Physics}, ic::KinematicInit)
    Systems.init!(sys.kinematics, ic)
    f_ode!(sys) #update state derivatives and outputs
end

###############################################################################
############################# AbstractAvionics #################################

abstract type AbstractAvionics <: SystemDefinition end

################################### NoAvionics #################################

struct NoAvionics <: AbstractAvionics end


###############################################################################
#################### AbstractPlatform update methods ###########################

function Systems.f_ode!(platform::System{<:AbstractPlatform},
                        kin::KinematicData,
                        air::AirData,
                        trn::AbstractTerrain)
    MethodError(f_ode!, (platform, avionics, kin, air, trn)) |> throw
end

function Systems.f_ode!(::System{EmptyPlatform},
                        ::KinematicData,
                        ::AirData,
                        ::AbstractTerrain)
    nothing
end

#this method can be extended if required, but in principle Platform shouldn't
#implement discrete dynamics; discretized algorithms belong in Avionics
function Systems.f_disc!(::System{<:AbstractPlatform}, ::Real)
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
    @unpack kinematics, platform, atmosphere = subsystems
    @unpack terrain = constants

    f_ode!(kinematics)
    f_ode!(atmosphere) #currently does nothing
    kin_data = KinematicData(kinematics)
    atm_data = LocalAtmosphericData(atmosphere)

    air_data = AirData(kin_data, atm_data)

    #update platform
    f_ode!(platform, kin_data, air_data, terrain)

    #get inputs for rigid body dynamics
    mp_Ob = get_mp_Ob(platform)
    wr_b = get_wr_b(platform)
    hr_b = get_hr_b(platform)

    #update velocity derivatives and rigid body data
    dyn_data = Dynamics.update!(kinematics.ẋ.vel, kin_data, mp_Ob, wr_b, hr_b)

    physics.y = PhysicsY(platform.y, kinematics.y, dyn_data, air_data)

    return nothing

end

#f_step! will use the recursive fallback implementation

#within Physics, only the platform may be modified by f_disc! (and it generally
#shouldn't)
function Systems.f_disc!(physics::System{<:Physics}, Δt::Real)

    @unpack kinematics, platform = physics
    @unpack dynamics, air = physics.y

    x_mod = false
    x_mod |= f_disc!(platform, Δt)

    #platform might have modified its outputs, so we need to reassemble our y
    physics.y = PhysicsY(platform.y, kinematics.y, dynamics, air)

    return x_mod
end

#these are meant to map avionics outputs to platform control inputs, they are
#called both within the aircraft's f_ode! and f_disc! before the physics update
function assign!(physics::System{<:Physics}, avionics::System{<:AbstractAvionics})
    assign!(physics.platform, avionics)
end

function assign!(platform::System{<:AbstractPlatform},
                avionics::System{<:AbstractAvionics})
    MethodError(assign!, (platform, avionics)) |> throw
end

assign!(::System{<:AbstractPlatform}, ::System{NoAvionics}) = nothing


################################################################################
################################## Aircraft ####################################

@kwdef struct Aircraft{P <: Physics, A <: AbstractAvionics} <: SystemDefinition
    physics::P = Physics()
    avionics::A = NoAvionics()
end

struct AircraftY{P <: PhysicsY, A}
    physics::P
    avionics::A
end

Systems.Y(ac::Aircraft) = AircraftY(Systems.Y(ac.physics), Systems.Y(ac.avionics))

function Systems.init!(ac::System{<:Aircraft}, ic::KinematicInit)
    Systems.init!(ac.physics, ic)
    f_ode!(ac) #update state derivatives and outputs
end

function Systems.f_ode!(ac::System{<:Aircraft})

    @unpack physics, avionics = ac.subsystems

    f_ode!(avionics, physics)
    assign!(physics, avionics)
    f_ode!(physics)

    ac.y = AircraftY(physics.y, avionics.y)

    return nothing

end

function Systems.f_disc!(ac::System{<:Aircraft}, Δt::Real)

    @unpack physics, avionics = ac.subsystems

    x_mod = false
    x_mod |= f_disc!(avionics, physics, Δt)
    assign!(physics, avionics)
    x_mod |= f_disc!(physics, Δt)

    ac.y = AircraftY(physics.y, avionics.y)

    return x_mod
end

#f_step! can use the recursive fallback implementation


############################# XPlaneConnect ####################################

Visualization.set_position!(xpc::XPCDevice, y::AircraftY) = Visualization.set_position!(xpc, y.physics)

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

Systems.init!( ac::System{<:Aircraft}, params::AbstractTrimParameters) = trim!(ac, params)

trim!( ac::System{<:Aircraft}, args...; kwargs...) = trim!(ac.physics, args...; kwargs...)

function trim!( physics::System{<:Physics}, args...)
    MethodError(trim!, (physics, args...)) |> throw
end

function assign!(::System{<:Physics}, ::AbstractTrimParameters, ::AbstractTrimState)
    error("An assign! method must be defined by each AircraftBase.Physics subtype")
end


################################################################################
################################ Linearization #################################

ẋ_linear(physics::System{<:Physics})::FieldVector = throw(MethodError(ẋ_linear!, (physics,)))
x_linear(physics::System{<:Physics})::FieldVector = throw(MethodError(x_linear!, (physics,)))
u_linear(physics::System{<:Physics})::FieldVector = throw(MethodError(u_linear!, (physics,)))
y_linear(physics::System{<:Physics})::FieldVector = throw(MethodError(y_linear!, (physics,)))

assign_x!(physics::System{<:Physics}, x::AbstractVector{Float64}) = throw(MethodError(assign_x!, (physics, x)))
assign_u!(physics::System{<:Physics}, u::AbstractVector{Float64}) = throw(MethodError(assign_u!, (physics, u)))

linearize!(ac::System{<:AircraftBase.Aircraft}, args...) = linearize!(ac.physics, args...)

function AircraftBase.linearize!( physics::System{<:AircraftBase.Physics},
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

    #now we need to rebuild vectors and matrices for the LinearizedSS as
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

    return Control.Continuous.LinearizedSS(ẋ0_cv, x0_cv, u0_cv, y0_cv, A_cv, B_cv, C_cv, D_cv)

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


function Control.Continuous.LinearizedSS(
            ac::System{<:AircraftBase.Aircraft}, args...; kwargs...)
    Control.Continuous.LinearizedSS(ac.physics, args...; kwargs...)
end

############################### Plotting #######################################

function Plotting.make_plots(ts::TimeSeries{<:PhysicsY}; kwargs...)

    return OrderedDict(
        :platform => make_plots(ts.platform; kwargs...),
        :kinematics => make_plots(ts.kinematics; kwargs...),
        :dynamics => make_plots(ts.dynamics; kwargs...),
        :air => make_plots(ts.air; kwargs...),
    )

end

function Plotting.make_plots(ts::TimeSeries{<:AircraftY}; kwargs...)

    return OrderedDict(
        :physics => make_plots(ts.physics; kwargs...),
        :avionics => make_plots(ts.avionics; kwargs...),
    )

end

################################### GUI ########################################


function GUI.draw!(sys::System{<:Aircraft};
                    p_open::Ref{Bool} = Ref(true), label::String = "Aircraft")

    @unpack y = sys

    CImGui.Begin(label, p_open)

    @cstatic c_phy=false c_avs=false begin
        @c CImGui.Checkbox("Physics", &c_phy)
        c_phy && @c GUI.draw!(sys.physics, sys.avionics, &c_phy)
        @c CImGui.Checkbox("Avionics", &c_avs)
        c_avs && @c GUI.draw!(sys.avionics, sys.physics, &c_avs)
    end

    CImGui.End()

end

function GUI.draw!(physics::System{<:Physics},
                   avionics::System{<:AbstractAvionics},
                   p_open::Ref{Bool} = Ref(true),
                   label::String = "Aircraft Physics")

    @unpack platform, atmosphere = physics.subsystems
    @unpack terrain = physics.constants
    @unpack kinematics, dynamics, air = physics.y

    CImGui.Begin(label, p_open)

    @cstatic(c_afm=false, c_atm=false, c_trn=false, c_dyn =false, c_kin=false, c_air=false,
    begin
            @c CImGui.Checkbox("Platform", &c_afm)
            @c CImGui.Checkbox("Atmosphere", &c_atm)
            @c CImGui.Checkbox("Terrain", &c_trn)
            @c CImGui.Checkbox("Dynamics", &c_dyn)
            @c CImGui.Checkbox("Kinematics", &c_kin)
            @c CImGui.Checkbox("Air", &c_air)
            c_afm && @c GUI.draw!(physics.platform, avionics, &c_afm)
            c_atm && @c GUI.draw!(physics.atmosphere, &c_atm)
            c_trn && @c GUI.draw(terrain, &c_trn)
            c_dyn && @c GUI.draw(dynamics, &c_dyn)
            c_kin && @c GUI.draw(kinematics, &c_kin)
            c_air && @c GUI.draw(air, &c_air)
    end)

    CImGui.End()

end

end #module