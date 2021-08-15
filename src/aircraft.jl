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

export TestAircraft


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

Base.@kwdef struct TestAircraft <: AbstractSystem
    mass_model::ConstantMassModel = ConstantMassModel()
    power_plant::EThruster = EThruster()
    landing_gear::Nothing = nothing
    control_surfaces::Nothing = nothing
end

#AbstractSystem interface
#X and U return initialized state and input vectors. do NOT access Templates directly
const TestAircraftXTemplate = ComponentVector( pv = X(PosVel()), pwp = X(EThruster()))
const TestAircraftUTemplate = ComponentVector( pwp = U(EThruster()))

const TestAircraftX{D} = ComponentVector{Float64, D, typeof(getaxes(TestAircraftXTemplate))} where {D<:AbstractVector{Float64}}
const TestAircraftU{D} = ComponentVector{Float64, D, typeof(getaxes(TestAircraftUTemplate))} where {D<:AbstractVector{Float64}}

const TestAircraftD = Environment

struct TestAircraftY
    kin::KinY
    pwp::EThrusterY
end

X(::TestAircraft) = copy(TestAircraftXTemplate)
U(::TestAircraft) = copy(TestAircraftUTemplate)
Y(::TestAircraft) = copy(TestAircraftYTemplate)
D(::TestAircraft) = TestAircraftD()

function f_output!(ẋ::TestAircraftX, x::TestAircraftX, u::TestAircraftU, t::Real,
                   data::TestAircraftD, aircraft::TestAircraft)

    @unpack trn, atm = data

    y_pv = f_pos!(ẋ.pv.pos, x.pv) #ẋ.pv.pos updated

    mass_data = get_mass_data(aircraft.mass_model)
    air_data = AirDataSensed()

    pwp_data = EThrusterD(air_data)

    y_pwp = f_output!(ẋ.pwp, x.pwp, u.pwp, t, pwp_data, aircraft.power_plant) #ẋ.pwp updated
    wr_ext_Ob_b = y_pwp.wr_Ob_b
    h_rot_b = y_pwp.h_rot_b

    y_ac = f_vel!(ẋ.pv.vel, wr_ext_Ob_b, h_rot_b, mass_data, y_pv) #ẋ.pv.vel updated

    y_kin = KinY(y_pv.pos, y_pv.vel, y_ac)

    return TestAircraftY(y_kin, y_pwp)
end

end