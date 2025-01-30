module AircraftBase

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays
using FiniteDiff: finite_difference_jacobian! as jacobian!

using Flight.FlightCore
using Flight.FlightLib

export AbstractComponents, NoComponents
export AbstractAvionics, NoAvionics
export AbstractTrimParameters, AbstractTrimState
export trim!, linearize!


################################################################################
########################### AbstractComponents ###################################

abstract type AbstractComponents <: SystemDefinition end

#AbstractComponents subtypes will generally be too specific for the fallback
function Systems.f_ode!(components::System{<:AbstractComponents},
                        kin::KinData,
                        air::AirData,
                        trn::AbstractTerrain)
    MethodError(f_ode!, (components, kin, air, trn)) |> throw
end

#for efficiency, we disallow using the fallback method. it would traverse the
#whole Components System hierarchy, and without good reason, because in
#principle Components shouldn't implement discrete dynamics; discretized
#algorithms belong in Avionics. can still be overridden by subtypes if required
Systems.f_disc!(::NoScheduling, ::System{<:AbstractComponents}) = nothing


################################ NoComponents #################################

@kwdef struct NoComponents <: AbstractComponents
    mass_distribution::RigidBodyDistribution = RigidBodyDistribution(1, SA[1.0 0 0; 0 1.0 0; 0 0 1.0])
end

Dynamics.get_hr_b(::System{NoComponents}) = zeros(SVector{3})
Dynamics.get_wr_b(::System{NoComponents}) = Wrench()
Dynamics.get_mp_b(sys::System{NoComponents}) = MassProperties(sys.constants.mass_distribution)

function Systems.f_ode!(::System{NoComponents},
                        ::KinData,
                        ::AirData,
                        ::AbstractTerrain)
    nothing
end


################################################################################
############################## Vehicle ################################

@kwdef struct Vehicle{F <: AbstractComponents,
                      K <: AbstractKinematicDescriptor,
                      T <: AbstractTerrain} <: SystemDefinition
    components::F = NoComponents()
    kinematics::K = LTF()
    dynamics::RigidBodyDynamics = RigidBodyDynamics()
    terrain::T = HorizontalTerrain() #shared with other aircraft instances
    atmosphere::LocalAtmosphere = LocalAtmosphere() #externally controlled
end

struct VehicleY{F, K}
    components::F
    kinematics::K
    mass::MassProperties #complete vehicle body mass properties
    actions::Actions
    accelerations::Accelerations
    air::AirData
end

Systems.Y(ac::Vehicle) = VehicleY(
    Systems.Y(ac.components),
    Systems.Y(ac.kinematics),
    MassProperties(),
    Actions(),
    Accelerations(),
    AirData())

function Systems.init!(sys::System{<:Vehicle}, ic::KinInit)
    @unpack kinematics, dynamics = sys.subsystems
    Systems.init!(kinematics, ic)
    dynamics.x .= kinematics.u #ESSENTIAL!
    f_ode!(sys) #update ẋ and y
end


###############################################################################
############################# AbstractAvionics #################################

abstract type AbstractAvionics <: SystemDefinition end

Systems.init!( ::System{<:AbstractAvionics}, ic::KinInit) = nothing

################################### NoAvionics #################################

struct NoAvionics <: AbstractAvionics end

################################################################################
######################## Vehicle/Avionics update methods #######################

#in order to be able to trim the standalone Vehicle without any auxiliary
#Avionics, the Vehicle's update methods must not require Avionics as an
#argument. to make this possible, any changes to Vehicle should be done through
#the joint assign! methods. accordingly, Avionics update methods should only
#mutate the Avionics System itself, not the Vehicle

# function Systems.f_ode!(::System{<:AbstractAvionics},
#                         ::System{<:Vehicle})
#     nothing
# end

# function Systems.f_disc!(sch::NoScheduling, avionics::System{<:AbstractAvionics},
#                         vehicle::System{<:Vehicle})
#     MethodError(f_disc!, (sch, avionics, vehicle)) |> throw
# end

# function Systems.f_disc!(::NoScheduling, ::System{NoAvionics}, ::System{<:Vehicle})
# end

#f_step! can use the recursive fallback implementation

function Systems.f_ode!(vehicle::System{<:Vehicle})

    @unpack ẋ, x, subsystems, constants = vehicle
    @unpack kinematics, dynamics, components, atmosphere = subsystems
    @unpack terrain = constants

    kinematics.u .= dynamics.x
    f_ode!(kinematics) #update ẋ and y before extracting kinematics data
    f_ode!(atmosphere) #update ẋ and y before extracting atmospheric data
    kin_data = KinData(kinematics)
    atm_data = AtmData(atmosphere)
    air_data = AirData(kin_data, atm_data)

    #update components' ẋ and y
    f_ode!(components, kin_data, air_data, terrain)

    #components is the only subsystem that has mass and can receive actions
    mp_Ob = get_mp_b(components)
    wr_ext_Ob = get_wr_b(components)
    hr_b = get_hr_b(components)
    actions = Actions(components; mp_Ob, wr_ext_Ob, hr_b, kin_data)

    #update velocity derivatives and rigid body data
    f_ode!(dynamics, mp_Ob, kin_data, actions)
    accelerations = dynamics.y

    vehicle.y = VehicleY(components.y, kinematics.y, mp_Ob, actions, accelerations, air_data)

    return nothing

end

#f_step! will use the recursive fallback implementation

#within Vehicle, only the components may be modified by f_disc!
function Systems.f_disc!(::NoScheduling, vehicle::System{<:Vehicle})

    @unpack components = vehicle.subsystems
    @unpack kinematics, mass, actions, accelerations, air = vehicle.y

    f_disc!(components)

    # components might have modified its outputs, so we need to reassemble y
    vehicle.y = VehicleY(components.y, kinematics, mass, actions, accelerations, air)

end

#these are meant to map avionics outputs to the vehicle, and in particular to
#components inputs. they are called both within the aircraft's f_ode! and
#f_disc! before the vehicle update
function assign!(vehicle::System{<:Vehicle}, avionics::System{<:AbstractAvionics})
    assign!(vehicle.components, avionics)
end

function assign!(components::System{<:AbstractComponents},
                avionics::System{<:AbstractAvionics})
    MethodError(assign!, (components, avionics)) |> throw
end

assign!(::System{<:AbstractComponents}, ::System{NoAvionics}) = nothing


################################################################################
################################## Aircraft ####################################

@kwdef struct Aircraft{P <: Vehicle, A <: AbstractAvionics} <: SystemDefinition
    vehicle::P = Vehicle()
    avionics::A = NoAvionics()
end

struct AircraftY{P <: VehicleY, A}
    vehicle::P
    avionics::A
end

Systems.Y(ac::Aircraft) = AircraftY(Systems.Y(ac.vehicle), Systems.Y(ac.avionics))

function Systems.init!(ac::System{<:Aircraft}, ic::KinInit)
    Systems.init!(ac.vehicle, ic)
    Systems.init!(ac.avionics, ic)
    f_ode!(ac) #update state derivatives and outputs
end

function Systems.f_ode!(ac::System{<:Aircraft})

    @unpack vehicle, avionics = ac.subsystems

    f_ode!(avionics, vehicle)
    assign!(vehicle, avionics)
    f_ode!(vehicle)

    ac.y = AircraftY(vehicle.y, avionics.y)
end

function Systems.f_disc!(::NoScheduling, ac::System{<:Aircraft})

    @unpack vehicle, avionics = ac.subsystems

    f_disc!(avionics, vehicle)
    assign!(vehicle, avionics)
    f_disc!(vehicle)

    ac.y = AircraftY(vehicle.y, avionics.y)
end

#f_step! can use the recursive fallback implementation

Kinematics.KinData(ac::System{<:Aircraft}) = KinData(ac.y.vehicle.kinematics)


################################# XPCClient ####################################

#extract attitude and position for X-Plane Connect client
function Systems.extract_output(ac::System{<:Aircraft}, ::Type{XPCPosition}, ::IOMapping)
    return XPCPosition(KinData(ac))
end


################################### Trimming ###################################
################################################################################

abstract type AbstractTrimParameters end
const AbstractTrimState{N} = FieldVector{N, Float64}

#given the body-axes wind-relative velocity, the wind-relative flight path angle
#and the bank angle, the pitch angle is unambiguously determined
function θ_constraint(; v_wOb_b, γ_wOb_n, φ_nb)
    TAS = norm(v_wOb_b)
    a = v_wOb_b[1] / TAS
    b = (v_wOb_b[2] * sin(φ_nb) + v_wOb_b[3] * cos(φ_nb)) / TAS
    sγ = sin(γ_wOb_n)

    return atan((a*b + sγ*√(a^2 + b^2 - sγ^2))/(a^2 - sγ^2))
    # return asin((a*sγ + b*√(a^2 + b^2 - sγ^2))/(a^2 + b^2)) #equivalent

end

Systems.init!( ac::System{<:Aircraft}, params::AbstractTrimParameters) = trim!(ac, params)

function trim!( ac::System{<:Aircraft}, params::AbstractTrimParameters)
    result = trim!(ac.vehicle, params) #compute vehicle trim state
    trim!(ac.avionics, ac.vehicle) #make avionics consistent with vehicle trim state
    ac.y = AircraftY(ac.vehicle.y, ac.avionics.y)
    return result
end

function trim!( vehicle::System{<:Vehicle}, args...)
    MethodError(trim!, (vehicle, args...)) |> throw
end

function assign!(::System{<:Vehicle}, ::AbstractTrimParameters, ::AbstractTrimState)
    error("An assign! method must be defined by each AircraftBase.Vehicle subtype")
end

function trim!( ::System{NoAvionics}, ::System{<:Vehicle})
end


################################################################################
################################ Linearization #################################

ẋ_linear(vehicle::System{<:Vehicle})::FieldVector = throw(MethodError(ẋ_linear!, (vehicle,)))
x_linear(vehicle::System{<:Vehicle})::FieldVector = throw(MethodError(x_linear!, (vehicle,)))
u_linear(vehicle::System{<:Vehicle})::FieldVector = throw(MethodError(u_linear!, (vehicle,)))
y_linear(vehicle::System{<:Vehicle})::FieldVector = throw(MethodError(y_linear!, (vehicle,)))

assign_x!(vehicle::System{<:Vehicle}, x::AbstractVector{Float64}) = throw(MethodError(assign_x!, (vehicle, x)))
assign_u!(vehicle::System{<:Vehicle}, u::AbstractVector{Float64}) = throw(MethodError(assign_u!, (vehicle, u)))

linearize!(ac::System{<:Aircraft}, args...) = linearize!(ac.vehicle, args...)

function linearize!( vehicle::System{<:Vehicle}, trim_params::AbstractTrimParameters)

    (_, trim_state) = trim!(vehicle, trim_params)

    ẋ0 = ẋ_linear(vehicle)::FieldVector
    x0 = x_linear(vehicle)::FieldVector
    u0 = u_linear(vehicle)::FieldVector
    y0 = y_linear(vehicle)::FieldVector

    #f_main will not be returned for use in another scope, so we don't need to
    #capture vehicle with a let block, because they are guaranteed not be
    #reassigned within the scope of linearize!
    function f_main(x, u)

        assign_x!(vehicle, x)
        assign_u!(vehicle, u)
        f_ode!(vehicle)

        return (ẋ = ẋ_linear(vehicle), y = y_linear(vehicle))

    end

    (A, B, C, D) = ss_matrices(f_main, x0, u0)

    #restore the System to its trimmed condition
    assign!(vehicle, trim_params, trim_state)

    #now we need to rebuild vectors and matrices for the LinearizedSS as
    #ComponentArrays, because we want matrix components to remain labelled,
    #which cannot be achieved with FieldVectors

    x_axis = Axis(propertynames(x0))
    u_axis = Axis(propertynames(u0))
    y_axis = Axis(propertynames(y0))

    ẋ0_cv = ComponentVector(ẋ0, x_axis)
    x0_cv = ComponentVector(x0, x_axis)
    u0_cv = ComponentVector(u0, u_axis)
    y0_cv = ComponentVector(y0, y_axis)

    A_cv = ComponentMatrix(A, x_axis, x_axis)
    B_cv = ComponentMatrix(B, x_axis, u_axis)
    C_cv = ComponentMatrix(C, y_axis, x_axis)
    D_cv = ComponentMatrix(D, y_axis, u_axis)

    return Control.Continuous.LinearizedSS(ẋ0_cv, x0_cv, u0_cv, y0_cv, A_cv, B_cv, C_cv, D_cv)

end

#given a trim point (x0, u0) for the nonlinear system (ẋ, y) = f(x,u), compute
#the state space matrices (A, B, C, D) for the system's linearization around
#(x0, u0), given by: Δẋ = AΔx + BΔu; Δy = CΔx + DΔu
function ss_matrices(f_main::Function, x0::AbstractVector{Float64},
                                       u0::AbstractVector{Float64})

    f_ẋ(x, u) = f_main(x, u).ẋ
    f_y(x, u) = f_main(x, u).y
    y0 = f_y(x0, u0)

    #none of these closures will be returned for use in another scope, so we
    #don't need to capture x0 and u0 with a let block, because they are
    #guaranteed not be reassigned within the scope of ss_matrices
    f_A!(ẋ, x) = (ẋ .= f_ẋ(x, u0))
    f_B!(ẋ, u) = (ẋ .= f_ẋ(x0, u))
    f_C!(y, x) = (y .= f_y(x, u0))
    f_D!(y, u) = (y .= f_y(x0, u))

    #preallocate mutable arrays
    A = (x0 * x0') |> Matrix
    B = (x0 * u0') |> Matrix
    C = (y0 * x0') |> Matrix
    D = (y0 * u0') |> Matrix

    jacobian!(A, f_A!, Vector(x0))
    jacobian!(B, f_B!, Vector(u0))
    jacobian!(C, f_C!, Vector(x0))
    jacobian!(D, f_D!, Vector(u0))

    return (A, B, C, D)

end


function Control.Continuous.LinearizedSS(
            ac::System{<:Aircraft}, args...; kwargs...)
    Control.Continuous.LinearizedSS(ac.vehicle, args...; kwargs...)
end

############################### Plotting #######################################

function Plotting.make_plots(ts::TimeSeries{<:VehicleY}; kwargs...)

    return OrderedDict(
        :components => make_plots(ts.components; kwargs...),
        :kinematics => make_plots(ts.kinematics; kwargs...),
        :accelerations => make_plots(ts.accelerations; kwargs...),
        :air => make_plots(ts.air; kwargs...),
    )

end

function Plotting.make_plots(ts::TimeSeries{<:AircraftY}; kwargs...)

    return OrderedDict(
        :vehicle => make_plots(ts.vehicle; kwargs...),
        :avionics => make_plots(ts.avionics; kwargs...),
    )

end

################################### GUI ########################################


function GUI.draw!(sys::System{<:Aircraft};
                    p_open::Ref{Bool} = Ref(true), label::String = "Aircraft")

    @unpack y = sys

    CImGui.Begin(label, p_open)

    @cstatic c_phy=false c_avs=false begin
        @c CImGui.Checkbox("Vehicle", &c_phy)
        c_phy && @c GUI.draw!(sys.vehicle, sys.avionics, &c_phy)
        @c CImGui.Checkbox("Avionics", &c_avs)
        c_avs && @c GUI.draw!(sys.avionics, sys.vehicle, &c_avs)
    end

    CImGui.End()

end

function GUI.draw!(vehicle::System{<:Vehicle},
                   avionics::System{<:AbstractAvionics},
                   p_open::Ref{Bool} = Ref(true),
                   label::String = "Vehicle")

    @unpack components, atmosphere = vehicle.subsystems
    @unpack terrain = vehicle.constants
    @unpack mass, air, actions, accelerations, kinematics = vehicle.y

    CImGui.Begin(label, p_open)

    @cstatic(c_afm=false, c_atm=false, c_trn=false, c_mas = false, c_act = false, c_acc =false, c_kin=false, c_air=false,
    begin
            @c CImGui.Checkbox("Components", &c_afm)
            @c CImGui.Checkbox("Atmosphere", &c_atm)
            @c CImGui.Checkbox("Terrain", &c_trn)
            @c CImGui.Checkbox("Air", &c_air)
            @c CImGui.Checkbox("Mass", &c_mas)
            @c CImGui.Checkbox("Actions", &c_act)
            @c CImGui.Checkbox("Accelerations", &c_acc)
            @c CImGui.Checkbox("Kinematics", &c_kin)
            c_afm && @c GUI.draw!(components, avionics, &c_afm)
            c_atm && @c GUI.draw!(atmosphere, &c_atm)
            c_trn && @c GUI.draw(terrain, &c_trn)
            c_air && @c GUI.draw(air, &c_air)
            c_mas && @c GUI.draw(mass, &c_mas)
            c_act && @c GUI.draw(actions, &c_act)
            c_acc && @c GUI.draw(accelerations, &c_acc)
            c_kin && @c GUI.draw(kinematics.data, &c_kin)
    end)

    CImGui.End()

end

end #module