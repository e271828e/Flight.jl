module TestC172Yv2

using Test, UnPack, BenchmarkTools, Sockets, JSON3, Logging

using Flight.FlightCore
using Flight.FlightLib
using Flight.FlightAircraft

using Flight.FlightAircraft.C172Y.C172YControl: ModeControlLon, ModeControlLat,
        AltTrackingState, is_on_gnd
using Flight.FlightAircraft.C172Y.C172YGuidance: Segment

export test_c172y2

y_kin(aircraft::Model{<:Cessna172Yv2}) = aircraft.y.vehicle.kinematics
y_air(aircraft::Model{<:Cessna172Yv2}) = aircraft.y.vehicle.airflow
y_aero(aircraft::Model{<:Cessna172Yv2}) = aircraft.y.vehicle.systems.aero


end #module