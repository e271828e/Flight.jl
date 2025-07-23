module AircraftBase

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays
using FiniteDiff: finite_difference_jacobian! as jacobian!
using DataStructures: OrderedDict

using Flight.FlightCore
using Flight.FlightLib

export AbstractSystems, NoSystems
export AbstractSystemsInitializer, NoSystemsInitializer
export AbstractAvionics, NoAvionics
export AbstractTrimParameters, AbstractTrimState
export trim!, linearize!



################################################################################
############################ Vehicle Systems ###################################

abstract type AbstractSystems <: ModelDefinition end

################################## NoSystems ###################################

@kwdef struct NoSystems <: AbstractSystems
    mass_distribution::RigidBodyDistribution = RigidBodyDistribution(1, SA[1.0 0 0; 0 1.0 0; 0 0 1.0])
end

Dynamics.get_hr_b(::Model{NoSystems}) = zeros(SVector{3})
Dynamics.get_wr_b(::Model{NoSystems}) = Wrench()
Dynamics.get_mp_b(mdl::Model{NoSystems}) = MassProperties(mdl.mass_distribution)

@no_updates NoSystems

############################# Initialization ###################################

abstract type AbstractSystemsInitializer end

struct NoSystemsInitializer <: AbstractSystemsInitializer end

Modeling.init!(::Model{<:AbstractSystems}, init::NoSystemsInitializer) = nothing


################################################################################
############################## Vehicle ################################

@kwdef struct Vehicle{S <: AbstractSystems,
                      K <: AbstractKinematicDescriptor } <: ModelDefinition
    systems::S = NoSystems()
    kinematics::K = WA()
    dynamics::VehicleDynamics = VehicleDynamics()
end

struct VehicleY{S, K}
    systems::S
    kinematics::K
    dynamics::DynamicsData
    airflow::AirflowData
end

Modeling.Y(ac::Vehicle) = VehicleY(
    Modeling.Y(ac.systems),
    KinData(),
    DynamicsData(),
    AirflowData())


################################ Initialization ################################

struct VehicleInitializer{S <: AbstractSystemsInitializer}
    kin::KinInit
    sys::S
end

function Modeling.init!( mdl::Model{<:Vehicle},
                        init::VehicleInitializer,
                        atmosphere::Model{<:AbstractAtmosphere} = Model(SimpleAtmosphere()),
                        terrain::Model{<:AbstractTerrain} = Model(HorizontalTerrain()))

    @unpack kinematics, dynamics, systems = mdl.submodels
    Modeling.init!(kinematics, init.kin)
    Modeling.init!(systems, init.sys)
    dynamics.x .= kinematics.u #essential
    f_ode!(mdl, atmosphere, terrain) #update vehicle's ẋ and y
end


################################### Trimming ###################################

abstract type AbstractTrimParameters end
const AbstractTrimState{N} = FieldVector{N, Float64}

function assign!(vehicle::Model{<:Vehicle},
                params::AbstractTrimParameters,
                state::AbstractTrimState,
                args...)
    MethodError(assign!, (vehicle, params, state, args...)) |> throw
end

#to be implemented by each Vehicle subtype
function Modeling.init!( mdl::Model{<:Vehicle},
                        condition::AbstractTrimParameters, args...)
    MethodError(Modeling.init!, (mdl, condition, args...)) |> throw
end

function trim!( ac::Model, condition::AbstractTrimParameters, args...)
    MethodError(trim!, (ac, condition, args...)) |> throw
end

#trim constraint: given the body-axes wind-relative velocity, the wind-relative
#flight path angle and the bank angle, the pitch angle is unambiguously
#determined
function θ_constraint(; v_wb_b, γ_wb_n, φ_nb)
    TAS = norm(v_wb_b)
    a = v_wb_b[1] / TAS
    b = (v_wb_b[2] * sin(φ_nb) + v_wb_b[3] * cos(φ_nb)) / TAS
    sγ = sin(γ_wb_n)

    return atan((a*b + sγ*√(a^2 + b^2 - sγ^2))/(a^2 - sγ^2))
    # return asin((a*sγ + b*√(a^2 + b^2 - sγ^2))/(a^2 + b^2)) #equivalent
end



###############################################################################
############################# AbstractAvionics #################################

abstract type AbstractAvionics <: ModelDefinition end

################################### NoAvionics #################################

struct NoAvionics <: AbstractAvionics end
@no_updates NoAvionics

Modeling.init!(::Model{NoAvionics}, args...) = nothing

################################################################################
######################## Vehicle/Avionics update methods #######################

#we should be able to trim the standalone Vehicle without any auxiliary
#Avionics. therefore, the Vehicle's update methods must not require Avionics as
#an argument. to make this possible, any changes to Vehicle should be done
#through the joint assign! methods. accordingly, Avionics update methods should
#only mutate the Avionics Model itself, not the Vehicle

function Modeling.f_ode!(vehicle::Model{<:Vehicle},
                        atmosphere::Model{<:AbstractAtmosphere},
                        terrain::Model{<:AbstractTerrain})

    @unpack ẋ, x, submodels, constants = vehicle
    @unpack kinematics, dynamics, systems = submodels

    kinematics.u .= dynamics.x
    f_ode!(kinematics) #update ẋ and y before extracting kinematics data
    kin_data = KinData(kinematics)
    airflow_data = AirflowData(atmosphere, kin_data)

    #update vehicle systems
    f_ode!(systems, kin_data, airflow_data, terrain)

    #update vehicle dynamics
    f_ode!(dynamics, systems, kin_data)

    vehicle.y = VehicleY(systems.y, kinematics.y, dynamics.y, airflow_data)
    nothing

end

function Modeling.f_step!(vehicle::Model{<:Vehicle},
                         atm::Model{<:AbstractAtmosphere},
                         trn::Model{<:AbstractTerrain})

    @unpack systems, kinematics, dynamics = vehicle.submodels

    f_step!(kinematics)
    f_step!(systems, atm, trn)

end

#within Vehicle, only vehicle systems may be modified by f_disc!
function Modeling.f_disc!(::NoScheduling, vehicle::Model{<:Vehicle},
                         atmosphere::Model{<:AbstractAtmosphere},
                         terrain::Model{<:AbstractTerrain})

    @unpack systems, kinematics, dynamics = vehicle.submodels

    f_disc!(systems, atmosphere, terrain)

    # systems.y might have changed, so we should update vehicle.y
    vehicle.y = VehicleY(systems.y, kinematics.y, dynamics.y, vehicle.y.airflow)
    nothing

end

#these map avionics outputs to the vehicle, and in particular to
#systems inputs. they are called both within the aircraft's f_ode! and
#f_disc! before the vehicle update
function assign!(vehicle::Model{<:Vehicle}, avionics::Model{<:AbstractAvionics})
    assign!(vehicle.systems, avionics)
end

function assign!(systems::Model{<:AbstractSystems},
                avionics::Model{<:AbstractAvionics})
    MethodError(assign!, (systems, avionics)) |> throw
end

assign!(::Model{<:AbstractSystems}, ::Model{NoAvionics}) = nothing


################################################################################
################################## Aircraft ####################################

@kwdef struct Aircraft{V <: Vehicle, A <: AbstractAvionics} <: ModelDefinition
    vehicle::V = Vehicle()
    avionics::A = NoAvionics()
end

function Modeling.f_ode!(ac::Model{<:Aircraft},
                        atmosphere::Model{<:AbstractAtmosphere},
                        terrain::Model{<:AbstractTerrain})

    @unpack vehicle, avionics = ac.submodels
    f_ode!(avionics, vehicle)
    assign!(vehicle, avionics)
    f_ode!(vehicle, atmosphere, terrain)
    update_output!(ac)
end

function Modeling.f_disc!(::NoScheduling,
                        ac::Model{<:Aircraft},
                        atmosphere::Model{<:AbstractAtmosphere},
                        terrain::Model{<:AbstractTerrain})

    @unpack vehicle, avionics = ac.submodels
    f_disc!(avionics, vehicle)
    assign!(vehicle, avionics)
    f_disc!(vehicle, atmosphere, terrain)
    update_output!(ac)
end

function Modeling.f_step!(ac::Model{<:Aircraft},
                         atmosphere::Model{<:AbstractAtmosphere},
                         terrain::Model{<:AbstractTerrain})

    @unpack vehicle, avionics = ac.submodels
    f_step!(vehicle, atmosphere, terrain)
end

#the Vehicle's initialization methods (kinematics and trimming) accept
#atmosphere and terrain Models as optional arguments. we pass them if provided;
#otherwise, they will be instantiated ad hoc by the Vehicle's methods
function Modeling.init!( aircraft::Model{<:Aircraft},
                        condition::Union{<:VehicleInitializer, <:AbstractTrimParameters},
                        args...)

    @unpack vehicle, avionics = aircraft.submodels
    Modeling.init!(vehicle, condition, args...)
    Modeling.init!(avionics, vehicle) #avionics init only relies on vehicle
    update_output!(aircraft)
end

Kinematics.KinData(ac::Model{<:Aircraft}) = KinData(ac.vehicle.kinematics)


################################################################################
############################### XPlane12Control #################################

function IODevices.extract_output(ac::Model{<:Aircraft}, ::XPlane12ControlMapping)
    return Network.xpmsg_set_pose(XPlanePose(KinData(ac))) #UDP message
end


################################################################################
################################ Linearization #################################

ẋ_linear(vehicle::Model{<:Vehicle})::FieldVector = throw(MethodError(ẋ_linear, (vehicle,)))
x_linear(vehicle::Model{<:Vehicle})::FieldVector = throw(MethodError(x_linear, (vehicle,)))
u_linear(vehicle::Model{<:Vehicle})::FieldVector = throw(MethodError(u_linear, (vehicle,)))
y_linear(vehicle::Model{<:Vehicle})::FieldVector = throw(MethodError(y_linear, (vehicle,)))

assign_x!(vehicle::Model{<:Vehicle}, x::AbstractVector{Float64}) = throw(MethodError(assign_x!, (vehicle, x)))
assign_u!(vehicle::Model{<:Vehicle}, u::AbstractVector{Float64}) = throw(MethodError(assign_u!, (vehicle, u)))

linearize!(ac::Model{<:Aircraft}, args...) = linearize!(ac.vehicle, args...)

function linearize!( vehicle::Model{<:Vehicle},
                    trim_params::AbstractTrimParameters)

    #The velocity states in the linearized model must be aerodynamic so that
    #they can be readily used for flight control design. Since the velocity
    #states in the nonlinear model are Earth-relative, we should always set
    #wind velocity to zero for linearization
    atmosphere = Model(SimpleAtmosphere(; wind = NoWind()))
    terrain = Model(HorizontalTerrain())

    (_, trim_state) = Modeling.init!(vehicle, trim_params, atmosphere, terrain)

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
        f_ode!(vehicle, atmosphere, terrain)

        return (ẋ = ẋ_linear(vehicle), y = y_linear(vehicle))

    end

    (A, B, C, D) = ss_matrices(f_main, x0, u0)

    #restore the Model to its trimmed condition
    assign!(vehicle, trim_params, trim_state, atmosphere, terrain)

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

    #none of these closures will be returned for use in another scope or
    #reassigned within ss_matrices, so there is no need to capture x0 and u0
    #with a let block
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
            ac::Model{<:Aircraft}, args...; kwargs...)
    Control.Continuous.LinearizedSS(ac.vehicle, args...; kwargs...)
end

############################### Plotting #######################################

function Plotting.make_plots(ts::TimeSeries{<:VehicleY}; kwargs...)

    return OrderedDict(
        :systems => make_plots(ts.systems; kwargs...),
        :kinematics => make_plots(ts.kinematics; kwargs...),
        :dynamics => make_plots(ts.dynamics; kwargs...),
        :air => make_plots(ts.air; kwargs...),
    )

end


################################### GUI ########################################


function GUI.draw!(ac::Model{<:Aircraft}, p_open::Ref{Bool} = Ref(true),
                    label::String = "Aircraft")

    @unpack vehicle, avionics = ac.submodels
    CImGui.Begin(label, p_open)

    @cstatic c_phy=false c_avs=false begin
        @c CImGui.Checkbox("Vehicle", &c_phy)
        c_phy && @c GUI.draw!(vehicle, avionics, &c_phy)
        @c CImGui.Checkbox("Avionics", &c_avs)
        c_avs && @c GUI.draw!(avionics, vehicle, &c_avs)
    end

    CImGui.End()

end

function GUI.draw!(vehicle::Model{<:Vehicle},
                   avionics::Model{<:AbstractAvionics},
                   p_open::Ref{Bool} = Ref(true),
                   label::String = "Vehicle")

    @unpack systems = vehicle.submodels
    @unpack kinematics, dynamics, airflow = vehicle.y

    CImGui.Begin(label, p_open)

    @cstatic(c_afm=false, c_kin =false, c_dyn=false, c_air=false,
    begin
            @c CImGui.Checkbox("Systems", &c_afm)
            @c CImGui.Checkbox("Kinematics", &c_kin)
            @c CImGui.Checkbox("Dynamics", &c_dyn)
            @c CImGui.Checkbox("Airflow", &c_air)
            c_afm && @c GUI.draw!(systems, avionics, &c_afm)
            c_kin && @c GUI.draw(kinematics, &c_kin)
            c_dyn && @c GUI.draw(dynamics, &c_dyn)
            c_air && @c GUI.draw(airflow, &c_air)
    end)

    CImGui.End()

end

end #module