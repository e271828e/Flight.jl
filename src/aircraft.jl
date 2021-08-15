module Aircraft

using StaticArrays: SVector, SMatrix
using LinearAlgebra
using ComponentArrays
using UnPack

using Flight.Kinematics
using Flight.Airframe
using Flight.Powerplant
using Flight.Airdata
using Flight.System

import Flight.System: X, Y, U, D, f_output! #these we need to extend

export ParametricAircraft


abstract type AbstractTerrainModel end
abstract type AbstractAtmosphericModel end

#TerrainModel does not belong to the Aircraft itself. it must be defined
#separately, and passed as an external data source. the same goes for
#AtmosphereModel. there must be a level above the Aircraft, which will typically
#be the simulation scheduler, that defines both the environmental models and all
#the aircraft participating in the simulation. this may be a block based
#simulation engine or a custom made one.

#the AtmosphericModel contained in the Aircraft should behave as a gateway to
#the actual AtmosphericModel, through which the evolving atmospheric model can
#be queried by the Aircraft for the values at its current location. the
#TerrainModel can be handled similarly, because even if it is not evolving in
#time, it may be an arbitrarily complex terrain database which must serve all
#vehicles. initially, these Model "clients" will be simply references to the
#models themselves, because these will be constant, simple and shared with no
#one else

#an alternative approach, instead of querying the models in a client-server
#architecture, would be simply to provide access references to the environment
#model data and methods and store them in the external data source fields of the
#Aircraft System being simulated. then, the Aircraft can locally evaluate those
#methods for their particular position. in that case, the exchanged information
#is not the product of model evaluation at the Aircraft location, but the data
#required to perform those evaluations within the Aircraft model.

#both scenarios could be realized with an Observer model (or even Reactive
#programming, see Rocket.jl), with any vehicle in the simulation subscribing to
#the AtmosphericModel simulation and the TerrainModel. each time these update,
#they notify all their subscribers so they can also be evolved by the scheduler.
#if this communication occurs through Channels, each Model (Atmospheric,
#Terrain, Vehicle...) can run on separate Tasks, and each of these in a
#different thread. Thread safety.

struct DummyTerrainModel <: AbstractTerrainModel end
struct DummyAtmosphericModel <: AbstractAtmosphericModel end

Base.@kwdef struct Environment
    trn::AbstractTerrainModel = DummyTerrainModel()
    atm::AbstractAtmosphericModel = DummyAtmosphericModel()
end

abstract type AbstractMassModel end

#given some inputs (typically state of the fuel system and external payloads),
#an AbstractMassModel returns a MassData struct (defined in the Dynamics
#module). for now, we can simply define a ConstantMassModel

Base.@kwdef struct ConstantMassModel <: AbstractMassModel
    m::Float64 = 1.0
    J_Ob_b::SMatrix{3, 3, Float64, 9} = SMatrix{3,3,Float64}(I)
    r_ObG_b::SVector{3, Float64} = zeros(SVector{3})
end

get_mass_data(model::ConstantMassModel) = MassData(model.m, model.J_Ob_b, model.r_ObG_b)


################# SKETCH
#ver si esto ofrece alguna ventaja en cuanto a type stability respecto de
#definir pwp como campo de la immutable struct
struct ParametricAircraft{Pwp} <: AbstractSystem
    mass_model::AbstractMassModel
end
ParametricAircraft(mass_model, pwp) = ParametricAircraft{pwp}(mass_model)
function ParametricAircraft()
    pwp = ComponentGroup((left = EThruster(), right = EThruster(), back = EThruster(), front = EThruster()))
    ParametricAircraft(ConstantMassModel(), pwp)
end

X(::ParametricAircraft{Pwp}) where {Pwp} = ComponentVector(kin = X(Kin()), pwp = X(Pwp))
U(::ParametricAircraft{Pwp}) where {Pwp} = ComponentVector(pwp = U(Pwp))
Y(::ParametricAircraft{Pwp}) where {Pwp} = ComponentVector(kin = Y(Kin()), acc = Y(Acc()), pwp = Y(Pwp), air = Y(AirData()))
D(::ParametricAircraft) = Environment()

function f_output!(y, ẋ, x, u, t::Real,
                   data, aircraft::ParametricAircraft{Pwp}) where {Pwp}
    @unpack trn, atm = data

    #update kinematics
    f_pos!(y.kin, ẋ.kin.pos, x.kin)

    mass_data = get_mass_data(aircraft.mass_model)
    # y.air .= get_air_data(). #call air data system here to update air data, passing also as
    # argument data.atmospheric_model

    #update powerplant
    f_output!(y.pwp, ẋ.pwp, x.pwp, u.pwp, t, y.air, Pwp)
    #update landing gear
    # f_output!(y.ldg, ẋ.ldg, x.ldg, u.ldg, t, LdgD(y.kin, d.terrain), Ldg)

    #initialize external Wrench and additional angular momentum
    wr_ext_Ob_b = Wrench()
    h_rot_b = SVector(0.,0.,0.)

    #add powerplant contributions
    wr_ext_Ob_b .+= get_wr_Ob_b(y.pwp, Pwp)
    h_rot_b += get_h_Gc_b(y.pwp, Pwp)

    #update dynamics
    f_vel!(y.acc, ẋ.kin.vel, wr_ext_Ob_b, h_rot_b, mass_data, y.kin)

    return nothing
end


end