module Aircraft

using LinearAlgebra
using UnPack
using StaticArrays, ComponentArrays
using CImGui, CImGui.CSyntax, CImGui.CSyntax.CStatic

using Flight.Engine.Systems
using Flight.Engine.Plotting
using Flight.Engine.IODevices
using Flight.Engine.XPlane
using Flight.Engine.GUI

using Flight.Physics.Attitude
using Flight.Physics.Geodesy
using Flight.Physics.Kinematics
using Flight.Physics.RigidBody

using Flight.Components.Terrain
using Flight.Components.Environment

export AircraftTemplate, AbstractAirframe, AbstractAerodynamics, AbstractAvionics


###############################################################################
############################## Airframe #######################################

abstract type AbstractAirframe <: Component end
RigidBody.MassTrait(::System{<:AbstractAirframe}) = HasMass()
RigidBody.AngMomTrait(::System{<:AbstractAirframe}) = HasAngularMomentum()
RigidBody.WrenchTrait(::System{<:AbstractAirframe}) = GetsExternalWrench()

########################## EmptyAirframe ###########################

Base.@kwdef struct EmptyAirframe <: AbstractAirframe
    mass_distribution::RigidBodyDistribution = RigidBodyDistribution(1, SA[1.0 0 0; 0 1.0 0; 0 0 1.0])
end

RigidBody.AngMomTrait(::System{EmptyAirframe}) = HasNoAngularMomentum()
RigidBody.WrenchTrait(::System{EmptyAirframe}) = GetsNoExternalWrench()

RigidBody.get_mp_Ob(sys::System{EmptyAirframe}) = MassProperties(sys.params.mass_distribution)

####################### AbstractAerodynamics ##########################

abstract type AbstractAerodynamics <: Component end

RigidBody.MassTrait(::System{<:AbstractAerodynamics}) = HasNoMass()
RigidBody.AngMomTrait(::System{<:AbstractAerodynamics}) = HasNoAngularMomentum()
RigidBody.WrenchTrait(::System{<:AbstractAerodynamics}) = GetsExternalWrench()


###############################################################################
######################### AbstractAvionics ####################################

abstract type AbstractAvionics <: Component end

struct NoAvionics <: AbstractAvionics end

###############################################################################
############################## AircraftTemplate ###################################

struct AircraftTemplate{K <: AbstractKinematicDescriptor,
                F <: AbstractAirframe,
                A <: AbstractAvionics} <: Component
    kinematics::K
    airframe::F
    avionics::A
end

function AircraftTemplate(kinematics::K = LTF(),
                airframe::F = EmptyAirframe(),
                avionics::A = NoAvionics()) where {K,F,A}
    AircraftTemplate{K,F,A}(kinematics, airframe, avionics)
end

#override the default Component update_y! to include stuff besides subsystem
#outputs
Base.@kwdef struct AircraftTemplateY{K, F, A}
    kinematics::K
    airframe::F
    avionics::A
    rigidbody::RigidBodyData
    air::AirData
end

Systems.init(::SystemY, ac::AircraftTemplate) = AircraftTemplateY(
    init_y(ac.kinematics), init_y(ac.airframe), init_y(ac.avionics),
    RigidBodyData(), AirData())

function init!(ac::System{<:AircraftTemplate}, ic::KinematicInit)
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

    #update avionics and airframe components
    f_ode!(avionics, airframe, kin_data, air_data, trn)
    f_ode!(airframe, avionics, kin_data, air_data, trn)

    mp_Ob = get_mp_Ob(airframe)
    wr_b = get_wr_b(airframe)
    hr_b = get_hr_b(airframe)

    #update velocity derivatives
    rb_data = f_rigidbody!(kinematics.ẋ.vel, kin_data, mp_Ob, wr_b, hr_b)

    sys.y = AircraftTemplateY(kinematics.y, airframe.y, avionics.y, rb_data, air_data)

    return nothing

end

function Systems.f_step!(sys::System{<:AircraftTemplate})
    @unpack kinematics, airframe, avionics = sys

    #could use chained | instead, but this is clearer
    x_mod = false
    x_mod |= f_step!(kinematics)
    x_mod |= f_step!(airframe, avionics, kinematics)
    x_mod |= f_step!(avionics, airframe, kinematics)

    return x_mod
end

function Systems.f_disc!(sys::System{<:AircraftTemplate}, Δt)
    @unpack kinematics, airframe, avionics = sys

    #could use chained | instead, but this is clearer
    x_mod = false
    #in principle, only avionics will have discrete dynamics (it's the aircraft
    #subsystem in which discretized algorithms are implemented)
    x_mod |= f_disc!(avionics, airframe, kinematics, Δt)

    return x_mod
end


############################# XPlaneConnect ####################################

function XPlane.set_position!(xp::XPConnect, y::AircraftTemplateY)

    aircraft = 0

    kin = y.kinematics

    ll = LatLon(kin.n_e)
    e_nb = REuler(kin.q_nb)

    lat = rad2deg(ll.ϕ)
    lon = rad2deg(ll.λ)
    h = kin.h_o

    psi = rad2deg(e_nb.ψ)
    theta = rad2deg(e_nb.θ)
    phi = rad2deg(e_nb.φ)

    XPlane.set_position!(xp; lat, lon, h, psi, theta, phi, aircraft)

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

function GUI.draw!(sys::System{<:AircraftTemplate})

    @unpack y = sys

    CImGui.Begin("Aircraft")

    show_dyn = @cstatic check=false @c CImGui.Checkbox("Dynamics", &check)
    show_kin = @cstatic check=false @c CImGui.Checkbox("Kinematics", &check)
    show_air = @cstatic check=false @c CImGui.Checkbox("Air", &check)
    show_airframe = @cstatic check=false @c CImGui.Checkbox("Airframe", &check)
    show_avionics = @cstatic check=false @c CImGui.Checkbox("Avionics", &check)

    show_dyn && GUI.draw!(y.rigidbody, "Dynamics")
    show_kin && GUI.draw!(y.kinematics, "Kinematics")
    show_air && GUI.draw!(y.air, "Air")
    show_airframe && GUI.draw!(sys.airframe) #disallow setting airframe systems inputs directly
    show_avionics && GUI.draw!(sys.avionics)

    CImGui.End()

end


end #module