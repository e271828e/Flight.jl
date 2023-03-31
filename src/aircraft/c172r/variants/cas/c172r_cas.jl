module C172RCAS

using UnPack
using Printf
using CImGui, CImGui.CSyntax, CImGui.CSyntax.CStatic

using Flight.FlightCore
using Flight.FlightPhysics
using Flight.FlightAircraft

include(normpath(dirname(@__FILE__), "../../airframe.jl")); using .C172RAirframe

export Cessna172RX


################################################################################
############################ AugmentedControls #################################

################################################################################

struct CASLogic <: Component end

@enum FlightPhase begin
    phase_gnd = 0
    phase_air = 1
end

@enum CASState begin
    CAS_disabled = 0
    CAS_standby = 1
    CAS_active = 2
end

Base.@kwdef mutable struct CASLogicU
    enable::Bool = false
end

Base.@kwdef struct CASLogicY
    enable::Bool = false
    flight_phase::FlightPhase = phase_gnd
    state::CASState = CAS_disabled
end

Systems.init(::SystemU, ::CASLogic) = CASLogicU()
Systems.init(::SystemY, ::CASLogic) = CASLogicY()

function compute_outputs(logic::System{CASLogic}, airframe::System{<:Airframe})

    y_ldg = airframe.ldg.y
    nlg_wow = y_ldg.nose.strut.wow
    lmain_wow = y_ldg.left.strut.wow
    rmain_wow = y_ldg.right.strut.wow

    enable = logic.u.enable
    flight_phase = (nlg_wow && lmain_wow && rmain_wow) ? on_air : phase_gnd

    if !enable
        state = CAS_disabled
    else #CAS enabled
        state = (flight_phase == phase_gnd ? CAS_standby : CAS_enabled)
    end

    return CASLogic(; enable, flight_phase, state)

end

#purely periodic system, only updates its outputs here. for now it doesn't need
#a memory state, because all its outputs can be computed on the fly from its
#current inputs
function Systems.f_disc!(logic::System{CASLogic}, ::Real, airframe::System{<:Airframe})
    logic.y = compute_outputs(logic, airframe)
end

################################# RateCAS ######################################

Base.@kwdef struct RateCAS <: Component
    roll::PICompensator{1} = PICompensator{1}()
    pitch::PICompensator{1} = PICompensator{1}()
    yaw::PICompensator{1} = PICompensator{1}()
end

################################# RateCAS ######################################

Base.@kwdef struct AugmentedControls <: AbstractAvionics
    logic::CASLogic = CASLogic()
    rate::RateCAS = RateCAS()
end

#we could reuse ReversibleActuationU here, but noticing that for AugmentedControl aileron,
#elevator and rudder actually mean roll_input, pitch_input, yaw_input. so we may
#be better off redefining them. also, we need the ap_enable input

#with the current sign criteria, positive aileron, elevator and rudder inputs
#yield positive increments to p, q and r, respectively. this means that for
#example, if the output of the pitch rate compensator (proportional plus
#integral pitch rate error, q_dmd - q_actual) is positive, we need a positive
#elevator input to the airframe actuation

Base.@kwdef mutable struct AugmentedControlsU
    eng_start::Bool = false
    eng_stop::Bool = false
    CAS_enable::Bool = false
    throttle::Ranged{Float64, 0, 1} = 0.0
    mixture::Ranged{Float64, 0, 1} = 0.5
    roll_input::Ranged{Float64, -1, 1} = 0.0
    pitch_input::Ranged{Float64, -1, 1} = 0.0
    yaw_input::Ranged{Float64, -1, 1} = 0.0
    aileron_trim::Ranged{Float64, -1, 1} = 0.0
    elevator_trim::Ranged{Float64, -1, 1} = 0.0
    rudder_trim::Ranged{Float64, -1, 1} = 0.0
    flaps::Ranged{Float64, 0, 1} = 0.0
    brake_left::Ranged{Float64, 0, 1} = 0.0
    brake_right::Ranged{Float64, 0, 1} = 0.0
end

const AugmentedCommands = C172RAirframe.ReversibleActuationY

Systems.init(::SystemU, ::AugmentedControls) = FeedthroughActuationU()
function Systems.init(::SystemY, c::AugmentedControls)
    return (logic = init_y(c.logic), rate = init_y(c.rate), cmd = AugmentedCommands())
end


function Systems.f_ode!(sys::System{AugmentedControls},
                        airframe::System{<:Airframe},
                        kin::KinematicData, ::AirData, ::RigidBodyData,
                        ::System{<:AbstractTerrain})

    @unpack logic, rate = sys
    @unpack roll_input, pitch_input, yaw_input = sys.u

    #we can do this in all cases, because if !CAS_active, the compensators will
    #be set to reset anyway
    p_dmd = roll_input
    p = kin.common.ω_lb_b[1]

    rate.roll.u.input .= p_dmd - p
    rate.pitch.u.input .= q_dmd - q

    f_ode!(rate) #update rate ẋ and y

    if logic.y.state == CAS_active
        aileron = rate.y.roll.out[1]
        elevator = rate.y.pitch.out[1]
        rudder = rate.y.yaw.out[1]
    else #standby or disabled, direct controls
        aileron = roll_input
        elevator = pitch_input
        rudder = yaw_input
    end

    cmd = AugmentedCommands(;
        eng_start, eng_stop, throttle, mixture, aileron, elevator, rudder,
        aileron_trim, moar
    )

    sys.y = (logic = logic.y, rate = rate.y, cmd = cmd)

end

function Systems.f_disc!(sys::System{AugmentedControls}, Δt::Real,
                        airframe::System{<:Airframe},
                        ::KinematicData, ::RigidBodyData, ::AirData,
                        ::System{<:AbstractTerrain})

    @unpack logic, rate = sys

    logic.u.enable = sys.u.CAS_enable
    f_disc!(logic, Δt, airframe)

    if logic.y.state != CAS_active #reset
        rate.roll.u.reset .= true
        rate.pitch.u.reset .= true
        rate.yaw.u.reset .= true
    else
        rate.roll.u.reset .= false
        rate.pitch.u.reset .= false
        rate.yaw.u.reset .= false
    end

    #we need to update the outputs to wrap the updated logic.y
    sys.y = (logic = logic.y, rate = rate.y, cmd = sys.y.cmd)

    return false

end

#f_step! can safely use the fallback method


# function Aircraft.map_controls!(airframe::System{<:Airframe}, avionics::System{AugmentedControls})

#     @unpack sm, ap, act = avionics

#     @unpack

#     # @unpack throttle, aileron_trim, aileron, elevator_trim, elevator,
#     #         rudder_trim, rudder, brake_left, brake_right, flaps, mixture,
#     #         eng_start, eng_stop = act.y

#     @unpack aero, pwp, ldg = airframe

#     pwp.u.engine.start = eng_start
#     pwp.u.engine.stop = eng_stop
#     pwp.u.engine.thr = throttle
#     pwp.u.engine.mix = mixture
#     ldg.u.nose.steering[] = (rudder_trim + rudder) #rudder↑ (right pedal forward) -> nose wheel steering right
#     ldg.u.left.braking[] = brake_left
#     ldg.u.right.braking[] = brake_right
#     aero.u.e = (elevator_trim + elevator) #elevator↑ (stick forward) -> e↑ -> pitch down
#     aero.u.a = (aileron_trim + aileron) #aileron↑ (stick right) -> a↑ -> roll right
#     aero.u.r = -(rudder_trim + rudder) #rudder↑ (right pedal forward) -> r↓ -> yaw right
#     aero.u.f = flaps #flaps↑ -> δf↑

#     return nothing
# end

#opcion 1:
#definimos un AugmentedControlsU identico a RevControlsU pero sustituyendo elevator por
#pitch_input, aileron por roll_input y rudder por yaw_input

#internamente hacemos que el origen de elevator de RevAvionicsU se determine de
#una de dos maneras, wow o no wow, y en funcion de ello le asignamos pitch_input
#o pitch_cas_output

#opcion 2:
#prescindimos de FeedthroughActuation, y consideramos solo AugmentedControlsU / AutopilotU
#mas adelante podemos meter actuators



end #module