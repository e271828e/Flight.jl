module StateMachine

using Flight.System
import Flight.System: DiscreteSystem, get_y0, get_u0, get_d0, f_disc!

export AbstractStateMachine, NoStateMachine

#a StateMachine is really a DiscreteSystem. therefore, it is described by an
#AbstractComponent (but not an AirframeComponent)
abstract type AbstractStateMachine <: AbstractComponent end

struct NoStateMachine <: AbstractStateMachine end
struct NoStateMachineY <: AbstractY{NoStateMachine} end
struct NoStateMachineU <: AbstractU{NoStateMachine} end
struct NoStateMachineD <: AbstractD{NoStateMachine} end
get_y0(::NoStateMachine) = NoStateMachineY()
get_d0(::NoStateMachine) = NoStateMachineD()
get_u0(::NoStateMachine) = NoStateMachineU()

function DiscreteSystem(stm::AbstractStateMachine, y = get_y0(stm),
    u = get_u0(stm), d = get_d0(stm), t = Ref(0.0))
    params = stm #params is the component itself
    subsystems = nothing #no subsystems to define
    DiscreteSystem{map(typeof, (stm, y, u, d, params, subsystems))...}(y, u, d, t, params, subsystems)
end


end