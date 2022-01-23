module C182T

using LinearAlgebra
using StaticArrays
using ComponentArrays
using UnPack
using Unitful

using Flight.Modeling
using Flight.Plotting
using Flight.Misc
using Flight.Attitude
using Flight.Terrain
using Flight.Airdata
using Flight.Kinematics
using Flight.Dynamics
using Flight.Components
using Flight.Aerodynamics: AbstractAerodynamics
using Flight.Propulsion: EThruster, ElectricMotor, SimpleProp, CW, CCW
using Flight.LandingGear: LandingGearUnit, DirectSteering, DirectBraking, Strut, SimpleDamper
using Flight.Aircraft: AircraftBase, AbstractAircraftID
using Flight.Input: XBoxController, get_axis_value, is_released

import Flight.Modeling: init_x, init_y, init_u, init_d, f_cont!, f_disc!
import Flight.Plotting: plots
import Flight.Components: MassTrait, WrenchTrait, AngularMomentumTrait, get_wr_b, get_mp_b
import Flight.Aircraft: assign_joystick_inputs!

# export BeaverDescriptor

struct ID <: AbstractAircraftID end

struct Payload <: SystemDescriptor
    pilot::MassProperties #mp_b
    copilot::MassProperties #mp_b
    baggage::MassProperties

    function Payload( ; pilot::AbstractMassDistribution = PointMass(75),
                        copilot::AbstractMassDistribution = PointMass(75),
                        baggage::AbstractMassDistribution = PointMass(50))

        pilot_slot = FrameTransform(r = SVector{3}(0.183, -0.356, 0.899))
        copilot_slot = FrameTransform(r = SVector{3}(0.183, 0.356, 0.899))
        baggage_slot = FrameTransform(r = SVector{3}(-1.316, 0, 0.899))

        return new( MassProperties(pilot, pilot_slot),
                    MassProperties(copilot, copilot_slot),
                    MassProperties(baggage, baggage_slot))

    end
end

Base.@kwdef mutable struct PayloadD
    pilot::Bool = true
    copilot::Bool = true
    baggage::Bool = true
end

init_d(::Payload) = PayloadD()

MassTrait(::System{Payload}) = HasMass()
WrenchTrait(::System{Payload}) = HasNoWrench()
AngularMomentumTrait(::System{Payload}) = HasNoAngularMomentum()

function get_mp_b(sys::System{Payload})
    mp_b = MassProperties()
    sys.d.pilot ? mp_b += sys.params.pilot : nothing
    sys.d.copilot ? mp_b += sys.params.copilot : nothing
    sys.d.baggage ? mp_b += sys.params.baggage : nothing
    return mp_b
end


struct Fuel <: SystemDescriptor end

init_x(::Fuel) = ComponentVector(m_left = 0.0, m_right = 0.0) #fuel tank contents

MassTrait(::System{Fuel}) = HasMass()
WrenchTrait(::System{Fuel}) = HasNoWrench()
AngularMomentumTrait(::System{Fuel}) = HasNoAngularMomentum()

function get_mp_b(sys::System{Fuel})
    mp_b = MassProperties()

    frame_left = FrameTransform(r = SVector{3}(0.325, -2.845, 0))
    frame_right = FrameTransform(r = SVector{3}(0.325, 2.845, 0))

    m_left = PointMass(sys.x.m_left)
    m_right = PointMass(sys.x.m_right)

    mp_b += MassProperties(m_left, frame_left)
    mp_b += MassProperties(m_right, frame_right)

    return mp_b
end

#here we would define f_cont! to account for fuel consumption. also could add
#discrete states to define which tank(s) we're drawing from


c = RigidBody(894, SA[1285 0 0; 0 1825 0; 0 0 2667])
t_bc = FrameTransform(r = SVector{3}(0.056, 0, 0.582))
mp_b = MassProperties(c, t_bc)


MTOW = 1406

#aerodynamics (SI)
S = 16.165
b = 10.912
c = 1.494


end #module