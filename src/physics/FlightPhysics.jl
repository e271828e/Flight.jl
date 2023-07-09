module FlightPhysics

using Reexport

include("quaternions.jl"); using .Quaternions
include("attitude.jl"); using .Attitude
include("geodesy.jl"); using .Geodesy
include("kinematics.jl"); using .Kinematics
include("rigidbody.jl"); using .RigidBody
include("atmosphere.jl"); using .Atmosphere
include("terrain.jl"); using .Terrain
include("environment.jl"); using .Environment

end