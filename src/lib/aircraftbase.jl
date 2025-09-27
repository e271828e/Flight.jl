module AircraftBase

using LinearAlgebra, UnPack, StaticArrays, ComponentArrays
using FiniteDiff: finite_difference_jacobian! as jacobian!
using DataStructures: OrderedDict
using CImGui: Begin, End, BeginTable, EndTable, TableNextColumn, TableNextRow,
              Text, CollapsingHeader, Separator, SameLine

using Flight.FlightCore
using Flight.FlightLib

export AbstractVehicleSystems, NoVehicleSystems
export AbstractVehicleSystemsInitializer, NoVehicleSystemsInitializer
export AbstractAvionics, NoAvionics
export AbstractTrimParameters, AbstractTrimState
export trim!, linearize!



################################################################################
############################ Vehicle Systems ###################################

abstract type AbstractVehicleSystems <: ModelDefinition end

################################## NoVehicleSystems ###################################

@kwdef struct NoVehicleSystems <: AbstractVehicleSystems
    mass_distribution::RigidBodyDistribution = RigidBodyDistribution(1, SMatrix{3,3,Float64}(I))
end

Dynamics.get_hr_b(::Model{NoVehicleSystems}) = zeros(SVector{3})
Dynamics.get_wr_b(::Model{NoVehicleSystems}) = Wrench()
Dynamics.get_mp_b(mdl::Model{NoVehicleSystems}) = MassProperties(mdl.mass_distribution)

@no_updates NoVehicleSystems

############################# Initialization ###################################

abstract type AbstractVehicleSystemsInitializer end

struct NoVehicleSystemsInitializer <: AbstractVehicleSystemsInitializer end

Modeling.init!(::Model{<:AbstractVehicleSystems}, init::NoVehicleSystemsInitializer) = nothing


################################################################################
############################## Vehicle ################################

@kwdef struct Vehicle{S <: AbstractVehicleSystems,
                      K <: AbstractKinematicDescriptor } <: ModelDefinition
    systems::S = NoVehicleSystems()
    kinematics::K = WA()
    dynamics::VehicleDynamics = VehicleDynamics()
end

struct VehicleY{S}
    systems::S
    kinematics::KinData
    dynamics::DynamicsData
    airflow::AirData
end

Modeling.Y(vehicle::Vehicle) = VehicleY(
    Modeling.Y(vehicle.systems),
    KinData(),
    DynamicsData(),
    AirData())


################################ Initialization ################################

struct VehicleInitializer{S <: AbstractVehicleSystemsInitializer}
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

function trim!( aircraft::Model, params::AbstractTrimParameters, args...)
    MethodError(trim!, (aircraft, params, args...)) |> throw
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

    @unpack ẋ, x, submodels, parameters = vehicle
    @unpack kinematics, dynamics, systems = submodels

    kinematics.u .= dynamics.x
    f_ode!(kinematics) #update ẋ and y before extracting kinematic data

    kin_data = KinData(kinematics)
    airflow_data = AirData(atmosphere, kin_data)

    #update vehicle systems
    f_ode!(systems, terrain, kin_data, airflow_data)

    #update vehicle dynamics
    mp_Σ_b = get_mp_b(systems)
    wr_Σ_b = get_wr_b(systems)
    ho_Σ_b = get_hr_b(systems)
    @unpack q_eb, r_eb_e = kin_data
    @pack! dynamics.u = mp_Σ_b, wr_Σ_b, ho_Σ_b, q_eb, r_eb_e

    f_ode!(dynamics)

    vehicle.y = VehicleY(systems.y, kin_data, dynamics.y, airflow_data)
    nothing

end

function Modeling.f_step!(vehicle::Model{<:Vehicle},
                         atmosphere::Model{<:AbstractAtmosphere},
                         terrain::Model{<:AbstractTerrain})

    @unpack systems, kinematics, dynamics = vehicle.submodels

    f_step!(kinematics)
    f_step!(systems, atmosphere, terrain)

end

#within Vehicle, only vehicle systems may be modified by f_periodic!
function Modeling.f_periodic!(::NoScheduling, vehicle::Model{<:Vehicle},
                         atmosphere::Model{<:AbstractAtmosphere},
                         terrain::Model{<:AbstractTerrain})

    @unpack systems, kinematics, dynamics = vehicle.submodels

    f_periodic!(systems, atmosphere, terrain)

    # systems.y might have changed, so we should update vehicle.y
    vehicle.y = VehicleY(systems.y, kinematics.y, dynamics.y, vehicle.y.airflow)
    nothing

end

#these map avionics outputs to the vehicle, and in particular to
#systems inputs. they are called both within the aircraft's f_ode! and
#f_periodic! before the vehicle update
function assign!(vehicle::Model{<:Vehicle}, avionics::Model{<:AbstractAvionics})
    assign!(vehicle.systems, avionics)
end

function assign!(systems::Model{<:AbstractVehicleSystems},
                avionics::Model{<:AbstractAvionics})
    MethodError(assign!, (systems, avionics)) |> throw
end

assign!(::Model{<:AbstractVehicleSystems}, ::Model{NoAvionics}) = nothing


################################################################################
################################## Aircraft ####################################

@kwdef struct Aircraft{V <: Vehicle, A <: AbstractAvionics} <: ModelDefinition
    vehicle::V = Vehicle()
    avionics::A = NoAvionics()
end

function Modeling.f_ode!(aircraft::Model{<:Aircraft},
                        atmosphere::Model{<:AbstractAtmosphere},
                        terrain::Model{<:AbstractTerrain})

    @unpack vehicle, avionics = aircraft
    f_ode!(avionics, vehicle)
    assign!(vehicle, avionics)
    f_ode!(vehicle, atmosphere, terrain)
    f_output!(aircraft)
end

function Modeling.f_periodic!(::NoScheduling,
                        aircraft::Model{<:Aircraft},
                        atmosphere::Model{<:AbstractAtmosphere},
                        terrain::Model{<:AbstractTerrain})

    @unpack vehicle, avionics = aircraft
    f_periodic!(avionics, vehicle)
    assign!(vehicle, avionics)
    f_periodic!(vehicle, atmosphere, terrain)
    f_output!(aircraft)
end

function Modeling.f_step!(aircraft::Model{<:Aircraft},
                         atmosphere::Model{<:AbstractAtmosphere},
                         terrain::Model{<:AbstractTerrain})

    @unpack vehicle, avionics = aircraft
    f_step!(avionics, vehicle)
    assign!(vehicle, avionics)
    f_step!(vehicle, atmosphere, terrain)
end

#the Vehicle's initialization methods (kinematics and trimming) accept
#atmosphere and terrain Models as optional arguments. we pass them if provided;
#otherwise, they will be instantiated ad hoc by the Vehicle's methods
function Modeling.init!( aircraft::Model{<:Aircraft},
                        condition::Union{<:VehicleInitializer, <:AbstractTrimParameters},
                        args...)

    @unpack vehicle, avionics = aircraft
    Modeling.init!(vehicle, condition, args...)
    Modeling.init!(avionics, vehicle) #avionics init only relies on vehicle
    f_output!(aircraft)
end

Kinematics.KinData(aircraft::Model{<:Aircraft}) = KinData(aircraft.vehicle.kinematics)


################################################################################
############################### XPlane12Control #################################

function IODevices.extract_output(aircraft::Model{<:Aircraft}, ::XPlane12ControlMapping)
    return Network.xpmsg_set_pose(XPlanePose(KinData(aircraft))) #UDP message
end


################################################################################
################################ Linearization #################################

ẋ_linear(vehicle::Model{<:Vehicle})::FieldVector = throw(MethodError(ẋ_linear, (vehicle,)))
x_linear(vehicle::Model{<:Vehicle})::FieldVector = throw(MethodError(x_linear, (vehicle,)))
u_linear(vehicle::Model{<:Vehicle})::FieldVector = throw(MethodError(u_linear, (vehicle,)))
y_linear(vehicle::Model{<:Vehicle})::FieldVector = throw(MethodError(y_linear, (vehicle,)))

assign_x!(vehicle::Model{<:Vehicle}, x::AbstractVector{Float64}) = throw(MethodError(assign_x!, (vehicle, x)))
assign_u!(vehicle::Model{<:Vehicle}, u::AbstractVector{Float64}) = throw(MethodError(assign_u!, (vehicle, u)))

linearize!(aircraft::Model{<:Aircraft}, args...) = linearize!(aircraft.vehicle, args...)

function linearize!( vehicle::Model{<:Vehicle}, trim_params::AbstractTrimParameters)

    #The velocity states in the linearized model must be aerodynamic so that
    #they can be readily used for flight control design. Since the velocity
    #states in the nonlinear model are Earth-relative, we should always set
    #wind velocity to zero for linearization
    atmosphere = Model(SimpleAtmosphere(; wind = NoWind()))
    terrain = Model(HorizontalTerrain())

    (_, trim_state) = Modeling.init!(vehicle, trim_params, atmosphere, terrain)

    x0 = x_linear(vehicle)::FieldVector
    u0 = u_linear(vehicle)::FieldVector

    #f and g will not be returned for use in another scope, so we don't need to
    #capture vehicle with a let block, because they are guaranteed not be
    #reassigned within the scope of linearize!
    function f(x, u)
        assign_x!(vehicle, x)
        assign_u!(vehicle, u)
        f_ode!(vehicle, atmosphere, terrain)
        return ẋ_linear(vehicle)
    end

    function g(x, u)
        assign_x!(vehicle, x)
        assign_u!(vehicle, u)
        f_ode!(vehicle, atmosphere, terrain)
        y_linear(vehicle)
    end

    #restore the Model to its trimmed condition
    assign!(vehicle, trim_params, trim_state, atmosphere, terrain)

    return Control.Continuous.LinearizedSS(f, g, x0, u0)

end


############################### Plotting #######################################

function Plotting.make_plots(ts::TimeSeries{<:VehicleY}; kwargs...)

    return OrderedDict(
        :systems => make_plots(ts.systems; kwargs...),
        :kinematics => make_plots(ts.kinematics; kwargs...),
        :dynamics => make_plots(ts.dynamics; kwargs...),
        :airflow => make_plots(ts.airflow; kwargs...),
    )

end


################################### GUI ########################################


function GUI.draw!(aircraft::Model{<:Aircraft}, p_open::Ref{Bool} = Ref(true),
                    label::String = "Aircraft")

    @unpack vehicle, avionics = aircraft
    CImGui.Begin(label, p_open)

    @cstatic c_veh=false c_avs=false begin
        @c CImGui.Checkbox("Vehicle", &c_veh)
        c_veh && @c GUI.draw!(vehicle, avionics, &c_veh)
        @c CImGui.Checkbox("Avionics", &c_avs)
        c_avs && @c GUI.draw!(avionics, vehicle, &c_avs)
    end

    CImGui.End()

end


function GUI.draw!(vehicle::Model{<:Vehicle},
                   avionics::Model{<:AbstractAvionics},
                   p_open::Ref{Bool} = Ref(true))

    @unpack systems = vehicle.submodels
    @unpack kinematics, dynamics, airflow = vehicle.y

    CImGui.Begin("Vehicle", p_open)

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


function GUI.draw( vehicle::VehicleY)

    @unpack kinematics, dynamics, airflow = vehicle

    @unpack e_nb, ω_wb_b, n_e, ϕ_λ, h_e, h_o, v_gnd, χ_gnd, γ_gnd, v_eb_n = kinematics
    @unpack CAS, EAS, TAS, M, T, p, ρ, Δp, v_wb_b = airflow
    @unpack ψ, θ, φ = e_nb
    @unpack ϕ, λ = ϕ_λ

    α, β = Atmosphere.get_airflow_angles(v_wb_b)
    clm = -v_eb_n[3]

    if BeginTable("Flight Data", 2, CImGui.ImGuiTableFlags_SizingStretchSame | CImGui.ImGuiTableFlags_BordersInner)
        TableNextRow()
            TableNextColumn()

            Text("Airspeed (Calibrated)"); SameLine(240)
            Text(@sprintf("%.3f m/s | %.3f kts", CAS, Atmosphere.SI2kts(CAS)))
            Text("Airspeed (Equivalent)"); SameLine(240)
            Text(@sprintf("%.3f m/s | %.3f kts", EAS, Atmosphere.SI2kts(EAS)))
            Text("Airspeed (True)"); SameLine(240)
            Text(@sprintf("%.3f m/s | %.3f kts", TAS, Atmosphere.SI2kts(TAS)))
            Text("Angle of Attack"); SameLine(240)
            Text(@sprintf("%.3f deg", rad2deg(α)))
            Text("Sideslip Angle"); SameLine(240)
            Text(@sprintf("%.3f deg", rad2deg(β)))

            Separator()

            Text("Dynamic Pressure"); SameLine(240)
            Text(@sprintf("%.3f Pa", Δp))
            Text("Impact Pressure"); SameLine(240)
            Text(@sprintf("%.3f Pa", Δp))
            Text("Mach"); SameLine(240)
            Text(@sprintf("%.3f", M))

            Separator()

            Text("Static Temperature"); SameLine(240)
            Text(@sprintf("%.3f K", T))
            Text("Static Pressure"); SameLine(240)
            Text(@sprintf("%.3f Pa", p))
            Text("Density"); SameLine(240)
            Text(@sprintf("%.3f kg/m^3", ρ))

            Separator()

            Text("Specific Force (x)"); SameLine(240)
            Text(@sprintf("%.3f g", dynamics.f_c_c[1]/Dynamics.g₀))
            Text("Specific Force (y)"); SameLine(240)
            Text(@sprintf("%.3f g", dynamics.f_c_c[2]/Dynamics.g₀))
            Text("Specific Force (z)"); SameLine(240)
            Text(@sprintf("%.3f g", dynamics.f_c_c[3]/Dynamics.g₀))


        TableNextColumn()

            Text("Roll Rate"); SameLine(240)
            Text(@sprintf("%.3f deg/s", rad2deg(ω_wb_b[1])))
            Text("Pitch Rate"); SameLine(240)
            Text(@sprintf("%.3f deg/s", rad2deg(ω_wb_b[2])))
            Text("Yaw Rate"); SameLine(240)
            Text(@sprintf("%.3f deg/s", rad2deg(ω_wb_b[3])))

            Separator()

            Text("Heading"); SameLine(240);
            Text(@sprintf("%.3f deg", rad2deg(ψ)))
            Text("Inclination"); SameLine(240)
            Text(@sprintf("%.3f deg", rad2deg(θ)))
            Text("Bank"); SameLine(240)
            Text(@sprintf("%.3f deg", rad2deg(φ)))

            Separator()

            Text("Latitude"); SameLine(240)
            Text(@sprintf("%.6f deg", rad2deg(ϕ)))
            Text("Longitude"); SameLine(240)
            Text(@sprintf("%.6f deg", rad2deg(λ)))
            Text("Altitude (Ellipsoidal)"); SameLine(240)
            Text(@sprintf("%.3f m | %.3f ft", Float64(h_e), Float64(h_e)/0.3048))
            Text("Altitude (Orthometric)"); SameLine(240)
            Text(@sprintf("%.3f m | %.3f ft", Float64(h_o), Float64(h_o)/0.3048))
            # Text("Height Over Ground"); SameLine(240)
            # Text(@sprintf("%.3f m | %.3f ft", hog, hog/0.3048))

            Separator()

            Text("Ground Speed"); SameLine(240)
            Text(@sprintf("%.3f m/s | %.3f kts", v_gnd, Atmosphere.SI2kts(v_gnd)))
            Text("Course Angle"); SameLine(240)
            Text(@sprintf("%.3f deg", rad2deg(χ_gnd)))
            Text("Flight Path Angle"); SameLine(240)
            Text(@sprintf("%.3f deg", rad2deg(γ_gnd)))
            Text("Climb Rate"); SameLine(240)
            Text(@sprintf("%.3f m/s", clm))


        EndTable()
    end

end


end #module