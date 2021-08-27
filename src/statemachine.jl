module StateMachine

using Flight.System
import Flight.System: DiscreteSystem, d0, u0, f_disc!

export AbstractStateMachine, NoStateMachine

#a StateMachine is really a DiscreteSystem. therefore, it is described by an
#AbstractComponent (but not an AirframeComponent)
abstract type AbstractStateMachine <: AbstractComponent end

struct NoStateMachine <: AbstractStateMachine end
struct NoStateMachineD <: AbstractD{NoStateMachine} end
struct NoStateMachineU <: AbstractU{NoStateMachine} end
d0(::NoStateMachine) = NoStateMachineD()
u0(::NoStateMachine) = NoStateMachineU()

function DiscreteSystem(stm::AbstractStateMachine, d = d0(stm), u = u0(stm), t = Ref(0.0))
    params = stm #params is the component itself
    subsystems = nothing #no subsystems to define
    DiscreteSystem{map(typeof, (stm, d, u, params, subsystems))...}(d, u, t, params, subsystems)
end


end