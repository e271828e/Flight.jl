module Aircraft

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays

using Flight.FlightCore.Systems
using Flight.FlightCore.Plotting
using Flight.FlightCore.GUI
using Flight.FlightCore.XPC

using Flight.FlightPhysics.Attitude
using Flight.FlightPhysics.Geodesy
using Flight.FlightPhysics.Kinematics
using Flight.FlightPhysics.RigidBody
using Flight.FlightPhysics.Environment

export AbstractAirframe, EmptyAirframe
export AbstractAvionics, NoAvionics
export AircraftTemplate
export AbstractTrimParameters
export init_kinematics!, trim!, linearize!

###############################################################################
############################## Airframe #######################################

abstract type AbstractAirframe <: SystemDefinition end
RigidBody.MassTrait(::System{<:AbstractAirframe}) = HasMass()
RigidBody.AngMomTrait(::System{<:AbstractAirframe}) = HasAngularMomentum()
RigidBody.WrenchTrait(::System{<:AbstractAirframe}) = GetsExternalWrench()

function Systems.f_ode!(airframe::System{<:AbstractAirframe},
                        kin::KinematicData, air::AirData,
                        trn::System{<:AbstractTerrain})
    MethodError(f_ode!, (airframe, kin, air, trn)) |> throw
end

########################## EmptyAirframe ###########################

@kwdef struct EmptyAirframe <: AbstractAirframe
    mass_distribution::RigidBodyDistribution = RigidBodyDistribution(1, SA[1.0 0 0; 0 1.0 0; 0 0 1.0])
end

RigidBody.AngMomTrait(::System{EmptyAirframe}) = HasNoAngularMomentum()
RigidBody.WrenchTrait(::System{EmptyAirframe}) = GetsNoExternalWrench()

RigidBody.get_mp_Ob(sys::System{EmptyAirframe}) = MassProperties(sys.params.mass_distribution)

function Systems.f_ode!(::System{EmptyAirframe}, ::KinematicData, ::AirData,
                        ::System{<:AbstractTerrain})
end

###############################################################################
######################### AbstractAvionics ####################################

abstract type AbstractAvionics <: SystemDefinition end

function Systems.f_disc!(avionics::System{<:AbstractAvionics}, Δt::Real,
                        airframe::System{<:AbstractAirframe},
                        kin::KinematicData, rb::RigidBodyData, air::AirData,
                        trn::System{<:AbstractTerrain})
    MethodError(f_disc!, (avionics, Δt, airframe, kin, air, rb, trn)) |> throw
end

function map_controls!(airframe::System{<:AbstractAirframe},
                       avionics::System{<:AbstractAvionics})
    MethodError(map_controls!, (airframe, avionics)) |> throw
end


################################### NoAvionics #################################

struct NoAvionics <: AbstractAvionics end

function Systems.f_disc!(::System{NoAvionics}, ::Real,
                        ::System{<:AbstractAirframe},
                        ::KinematicData, ::RigidBodyData, ::AirData,
                        ::TerrainData)
    return false
end

map_controls!(::System{<:AbstractAirframe}, ::System{<:NoAvionics}) = nothing


###############################################################################
############################## AircraftTemplate ###################################

struct AircraftTemplate{K <: AbstractKinematicDescriptor,
                F <: AbstractAirframe,
                A <: AbstractAvionics} <: SystemDefinition
    kinematics::K
    airframe::F
    avionics::A
end

function AircraftTemplate(kinematics::K = LTF(),
                airframe::F = EmptyAirframe(),
                avionics::A = NoAvionics()) where {K,F,A}
    AircraftTemplate{K,F,A}(kinematics, airframe, avionics)
end

#override the generic NamedTuple to include stuff besides subsystem outputs
@kwdef struct AircraftTemplateY{K, F, A}
    kinematics::K
    airframe::F
    avionics::A
    rigidbody::RigidBodyData
    air::AirData
end

Systems.init(::SystemY, ac::AircraftTemplate) = AircraftTemplateY(
    init_y(ac.kinematics), init_y(ac.airframe), init_y(ac.avionics),
    RigidBodyData(), AirData())

function init_kinematics!(ac::System{<:AircraftTemplate}, ic::KinematicInit)
    Kinematics.init!(ac.x.kinematics, ic)
end


function Systems.f_ode!(sys::System{<:AircraftTemplate}, env::System{<:AbstractEnvironment})

    @unpack ẋ, x, subsystems = sys
    @unpack kinematics, airframe, avionics = subsystems
    @unpack atm, trn = env

    #update kinematics
    f_ode!(kinematics)
    kin_data = KinematicData(kinematics)
    air_data = AirData(kin_data, atm)

    #update airframe components
    f_ode!(airframe, kin_data, air_data, trn)

    #get inputs for rigid body dynamics
    mp_Ob = get_mp_Ob(airframe)
    wr_b = get_wr_b(airframe)
    hr_b = get_hr_b(airframe)

    #update velocity derivatives and rigid body data
    rb_data = f_rigidbody!(kinematics.ẋ.vel, kin_data, mp_Ob, wr_b, hr_b)

    sys.y = AircraftTemplateY(kinematics.y, airframe.y, avionics.y, rb_data, air_data)

    return nothing

end

function Systems.f_disc!(sys::System{<:AircraftTemplate}, Δt::Real, env::System{<:AbstractEnvironment})

    @unpack airframe, avionics, kinematics = sys.subsystems
    y = sys.y

    kin_data = y.kinematics.common
    rb_data = y.rigidbody
    air_data = y.air
    trn_data = TerrainData(env.trn, kin_data.n_e)

    #could use chained | instead, but this is clearer
    x_mod = false
    #in principle, only avionics should have discrete dynamics (it's the only
    #aircraft subsystem in which discretized algorithms should live)
    x_mod |= f_disc!(avionics, Δt, airframe, kin_data, rb_data, air_data, trn_data)
    map_controls!(airframe, avionics)

    #avionics might have modified its outputs, so we need to reassemble everything
    sys.y = AircraftTemplateY(kinematics.y, airframe.y, avionics.y, rb_data, air_data)

    return x_mod
end

function Systems.f_step!(sys::System{<:AircraftTemplate})

    @unpack kinematics, airframe, avionics = sys.subsystems

    #could use chained | instead, but this is clearer
    x_mod = false
    x_mod |= f_step!(kinematics)
    x_mod |= f_step!(airframe)
    x_mod |= f_step!(avionics)

    return x_mod
end


############################# XPlaneConnect ####################################

function XPC.set_position!(xp::XPCDevice, y::AircraftTemplateY)

    aircraft = 0

    @unpack ϕ_λ, e_nb, h_o = y.kinematics

    lat = rad2deg(ϕ_λ.ϕ)
    lon = rad2deg(ϕ_λ.λ)

    psi = rad2deg(e_nb.ψ)
    theta = rad2deg(e_nb.θ)
    phi = rad2deg(e_nb.φ)

    XPC.set_position!(xp; lat, lon, h_o, psi, theta, phi, aircraft)

end


################################# Tools ########################################

abstract type AbstractTrimParameters end

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

function trim!( ac::System, args...; kwargs...)
    MethodError(trim!, (ac, args...)) |> throw
end

function linearize!(ac::System, args...; kwargs...)
    MethodError(trim!, (ac, args...)) |> throw
end

############################### Plotting #######################################

function Plotting.make_plots(th::TimeHistory{<:AircraftTemplateY}; kwargs...)

    return OrderedDict(
        :kinematics => make_plots(th.kinematics; kwargs...),
        :airframe => make_plots(th.airframe; kwargs...),
        :avionics => make_plots(th.avionics; kwargs...),
        :rigidbody => make_plots(th.rigidbody; kwargs...),
        :air => make_plots(th.air; kwargs...),
    )

end

################################### GUI ########################################

function GUI.draw!(sys::System{<:AircraftTemplate}, label::String = "Aircraft")

    @unpack y = sys

    CImGui.Begin(label)

    show_dyn = @cstatic check=false @c CImGui.Checkbox("Dynamics", &check)
    show_kin = @cstatic check=false @c CImGui.Checkbox("Kinematics", &check)
    show_air = @cstatic check=false @c CImGui.Checkbox("Air", &check)
    show_airframe = @cstatic check=false @c CImGui.Checkbox("Airframe", &check)
    show_avionics = @cstatic check=false @c CImGui.Checkbox("Avionics", &check)

    show_dyn && GUI.draw(y.rigidbody, "Dynamics")
    show_kin && GUI.draw(y.kinematics, "Kinematics")
    show_air && GUI.draw(y.air, "Air")
    show_airframe && GUI.draw!(sys.airframe, sys.avionics)
    show_avionics && GUI.draw!(sys.avionics, sys.airframe)

    CImGui.End()

end


end #module