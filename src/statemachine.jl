module StateMachine

using Flight.ModelingTools
import Flight.ModelingTools: System, init_y0, init_u0, init_d0, f_disc!

export AbstractStateMachine, NoStateMachine

#a StateMachine is really a DiscreteSystem. therefore, it is described by an
#AbstractComponent (but not an AirframeComponent)
abstract type AbstractStateMachine <: AbstractComponent end

struct NoStateMachine <: AbstractStateMachine end

function System(stm::AbstractStateMachine, ẋ = init_x0(stm),
    x = init_x0(stm), y = init_y0(stm), u = init_u0(stm), d = init_d0(stm), t = Ref(0.0))
    params = stm #params is the component itself
    subsystems = nothing #no subsystems to define
    System{map(typeof, (stm, x, y, u, d, params, subsystems))...}(ẋ, x, y, u, d, t, params, subsystems)
end


end