module Robot2D

using LinearAlgebra, StaticArrays, ComponentArrays
using ControlSystems, RobustAndOptimalControl, ComponentArrays, LinearAlgebra
using CImGui: Begin, End, Text, IsItemActive, Checkbox, CollapsingHeader,
            SameLine, PushItemWidth, PopItemWidth, BeginTable, EndTable,
            TableNextRow, TableNextColumn, AlignTextToFramePadding, Separator,
            ImVec2

using Flight.FlightCore
using Flight.FlightLib
using Flight.FlightLib.Linearization: delete_vars
using Flight.FlightLib.Control.Discrete: LQR, PID, LQROutput, PIDOutput
using Flight.FlightLib.Control.Discrete: LQRDataPoint


const g = 9.80665 #m/s^2, standard gravity

################################################################################
################################## Vehicle #####################################

#body 1: main body, comprising the vehicle chassis and the DC motor's case and stator
#body 2: rolling body, comprising the wheels, axle, and the DC motor's rotor.

@kwdef struct Vehicle <: ModelDefinition
    L::Float64 = 0.15 #distance from main body's origin to its CoM (m)
    R::Float64 = 0.05 #wheel radius (m)
    m_1::Float64 = 1.0 #mass of main body (kg)
    m_2::Float64 = 0.1 #mass of rolling body (kg)
    J_1::Float64 = 1/12 * m_1 * (2L)^2 #main body's moment of inertia with respect to its CoM (kg*m^2)
    J_2::Float64 = 1/2 * m_2 * R^2 #rolling body's moment of inertia with respect to its CoM (kg*m^2)
    k_m::Float64 = 0.32 #motor's torque constant (N*m)
    b_m::Float64 = 0.0189 #motor's effective damping coefficient (N*m*s/rad)
    J_m::Float64 = 0.0014 #motor's effective moment of inertia (kg*m^2)
end

@kwdef struct VehicleY
    ω::Float64 = 0.0 #angular velocity of main body (rad/s)
    v::Float64 = 0.0 #horizontal velocity of the system's origin (m/s)
    θ::Float64 = 0.0 #angle of main body with respect to vertical (rad)
    η::Float64 = 0.0 #horizontal position of the system's origin (m)
    u_m::Float64 = 0.0 #motor control input (normalized)
    τ_m::Float64 = 0.0 #motor torque (N*m)
    ω_dot::Float64 = 0.0 #angular acceleration of main body (rad/s^2)
    v_dot::Float64 = 0.0 #linear acceleration of the system's origin (m/s^2)
end

Modeling.X(::Vehicle) = ComponentVector( ω = 0.0, v = 0.0, θ = 0.0, η = 0.0)
Modeling.U(::Vehicle) = Ref(Ranged(0.0, -1., 1.))
Modeling.Y(::Vehicle) = VehicleY()

function Modeling.f_ode!(mdl::Model{Vehicle})

    (; ẋ, x, u, parameters) = mdl
    (; ω, v, θ, η) = x
    (; L, R, m_1, m_2, J_1, J_2, k_m, b_m, J_m) = parameters

    u_m = Float64(u[])
    ω_m = v / R - ω
    τ_ss = k_m * u_m - b_m * ω_m #steady-state motor torque

    sθ = sin(θ)
    cθ = cos(θ)

    M_11 = m_1 * L^2 + J_1 + J_m
    M_22 = m_1 + m_2 + (J_2 + J_m) / R^2
    M_12 = m_1 * L * cθ - J_m / R

    M = @SMatrix[
        M_11    M_12
        M_12    M_22
    ]

    b = @SVector[
        -τ_ss + m_1 * L * g * sθ
        τ_ss / R + m_1 * L * ω^2 * sθ
    ]

    ω_dot, v_dot = M\b
    ω_m_dot = v_dot/R - ω_dot

    θ_dot = ω
    η_dot = v

    τ_m = τ_ss - J_m * ω_m_dot

    ẋ.ω = ω_dot
    ẋ.v = v_dot
    ẋ.θ = θ_dot
    ẋ.η = η_dot

    mdl.y = VehicleY(; ω, v, θ, η, u_m, τ_m, ω_dot, v_dot)

end

@no_step Vehicle
@no_periodic Vehicle


function GUI.draw!(mdl::Model{<:Vehicle}, p_open::Ref{Bool} = Ref(true))

    (; u, y, parameters) = mdl
    (; ω, v, θ, η) = y
    (; R, L) = parameters

    Begin("Vehicle", p_open)

    if BeginTable("Text Data", 2, CImGui.ImGuiTableFlags_SizingStretchProp)# | CImGui.ImGuiTableFlags_Resizable)# | CImGui.ImGuiTableFlags_BordersInner)

        TableNextRow()
        TableNextColumn()
        AlignTextToFramePadding()
        Text("Motor Command")
        TableNextColumn()
        PushItemWidth(-10)
        u[] = GUI.safe_slider("##Motor Command", u[], "%.6f")
        PopItemWidth()

        TableNextRow()
        TableNextColumn()
        Text("Angular Rate")
        TableNextColumn()
        Text(@sprintf("%.6f deg/s", rad2deg(ω)))

        TableNextRow()
        TableNextColumn()
        Text("Angle")
        TableNextColumn()
        Text(@sprintf("%.6f deg", rad2deg(θ)))

        TableNextRow()
        TableNextColumn()
        Text("Velocity")
        TableNextColumn()
        Text(@sprintf("%.6f m/s", v))

        TableNextRow()
        TableNextColumn()
        Text("Position")
        TableNextColumn()
        Text(@sprintf("%.6f m", η))

        EndTable()
    end

    # 1. Define canvas dimensions and reserve space
    window_sz = CImGui.GetWindowSize()
    canvas_p0 = CImGui.GetCursorScreenPos() #Top-left of drawing area
    canvas_sz = ImVec2(0.95window_sz.x, 0.8window_sz.x)
    CImGui.Dummy(canvas_sz) #Reserve space in the layout

    # 2. Get the DrawList for the current window
    draw_list = CImGui.GetWindowDrawList()

    # 3. Calculate geometry
    # Colors (ABGR packed format for ImGui: 0xAABBGGRR)
    col_bg     = 0xFF202020 # Dark Gray Background
    col_wheel  = 0xFFFFFFFF # White
    col_chassis= 0xFF00FFFF # Yellow
    col_ground = 0xFF808080 # Gray
    col_spoke  = 0xFF0000FF # Red
    col_com  = 0xFF000000 # Black

    #for ImGui, θ=0 is upwards, CW is positive, y-axis points downwards
    scale = 0.6canvas_sz.y / 2L #total chassis length spans 60% of canvas height
    R_px = R * scale # Wheel radius in pixels (approx)
    L_px = L * scale # Distance to CoM

    #origin position
    cx = canvas_p0.x + canvas_sz.x * 0.5
    cy = canvas_p0.y + canvas_sz.y * 0.8  # Wheel sits near the bottom

    #CoM position
    gx = cx + L_px * sin(θ)
    gy = cy - L_px * cos(θ)

    #chassis end point
    ex = cx + 2L_px * sin(θ)
    ey = cy - 2L_px * cos(θ)

    #wheel radius line end point
    φ = η/R
    wx = cx + R_px * sin(φ)
    wy = cy - R_px * cos(φ)

    # 4. Draw commands
    # Background rectangle
    CImGui.AddRectFilled(draw_list, canvas_p0,
        ImVec2(canvas_p0.x + canvas_sz.x, canvas_p0.y + canvas_sz.y), col_bg)

    # Ground
    CImGui.AddLine(draw_list, ImVec2(canvas_p0.x, cy + R_px),
        ImVec2(canvas_p0.x + canvas_sz.x, cy + R_px), col_ground, 2.0)

    # Wheel Body
    CImGui.AddCircle(draw_list, ImVec2(cx, cy), R_px, col_wheel, 36, 2.0)

    # Wheel Radius
    CImGui.AddLine(draw_list, ImVec2(cx, cy), ImVec2(wx, wy), col_spoke, 2.0)

    # Chassis
    CImGui.AddLine(draw_list, ImVec2(cx, cy), ImVec2(ex, ey), col_chassis, 8.0)

    # CoM (Small Circle at tip)
    CImGui.AddCircleFilled(draw_list, ImVec2(gx, gy), 3.0, col_com)

    End()

end


################################################################################
################################ Initialization ################################

@kwdef struct InitParameters <: FieldVector{5, Float64}
    u_m::Float64 = 0.0
    ω_dot::Float64 = 0.0
    ω::Float64 = 0.0
    θ::Float64 = 0.0
    η::Float64 = 0.0
end

function Modeling.init!( mdl::Model{Vehicle}, ip::InitParameters = InitParameters())

    (; u_m, ω_dot, ω, θ, η) = ip
    (; x, u, parameters) = mdl
    (; L, R, m_1, m_2, J_1, J_2, k_m, b_m, J_m) = parameters

    sθ = sin(θ)
    cθ = cos(θ)

    M_11 = m_1 * L^2 + J_1 + J_m
    M_22 = m_1 + m_2 + (J_2 + J_m) / R^2
    M_12 = m_1 * L * cθ - J_m / R

    A = @SMatrix[
        1       M_12
        -1/R    M_22
    ]

    b = @SVector[
        m_1 * L * g * sθ - M_11 * ω_dot
        m_1 * L * ω^2 * sθ - M_12 * ω_dot
    ]

    τ_ss, _ = A\b

    ω_m = (k_m * u_m - τ_ss) / b_m
    v = (ω + ω_m) * R

    (x.ω, x.v, x.θ, x.η)  = (ω, v, θ, η)
    u[] = u_m

    f_ode!(mdl)

end

################################################################################
################################ State Space ###################################

@kwdef struct XStateSpace <: FieldVector{4, Float64}
    ω::Float64 = 0.0
    v::Float64 = 0.0
    θ::Float64 = 0.0
    η::Float64 = 0.0
end

@kwdef struct UStateSpace <: FieldVector{1, Float64}
    m::Float64 = 0.0
end

@kwdef struct YStateSpace <: FieldVector{6, Float64}
    ω::Float64 = 0.0
    v::Float64 = 0.0
    θ::Float64 = 0.0
    η::Float64 = 0.0
    u_m::Float64 = 0.0
    τ_m::Float64 = 0.0
end

#get state space system's state vector derivative from Model
function get_ẋ_ss(mdl::Model{Vehicle})
    (; ω, v, θ, η) = mdl.ẋ
    XStateSpace(; ω, v, θ, η)
end

#get state space system's state vector from Model
function get_x_ss(mdl::Model{Vehicle})
    (; ω, v, θ, η) = mdl.x
    XStateSpace(; ω, v, θ, η)
end

#get state space system's input vector from Model
function get_u_ss(mdl::Model{Vehicle})
    UStateSpace(mdl.u[])
end

#get state space system's output vector from Model
function get_y_ss(mdl::Model{Vehicle})
    (; ω, v, θ, η, u_m, τ_m) = mdl.y
    YStateSpace(; ω, v, θ, η, u_m, τ_m)
end

#assign state space system's state vector to Model
function assign_x_ss!(mdl::Model{Vehicle}, x_ss::AbstractVector{<:Real})
    (; ω, v, θ, η) = XStateSpace(x_ss)
    mdl.x.ω, mdl.x.v, mdl.x.θ, mdl.x.η = ω, v, θ, η
end

#assign state space system's input vector to Model
function assign_u_ss!(mdl::Model{Vehicle}, u_ss::AbstractVector{<:Real})
    mdl.u[] = UStateSpace(u_ss[1])[1]
end

#build state space system's update function
function get_f_ss(mdl::Model{Vehicle})
    let mdl = mdl
        function (x::AbstractVector{<:Real}, u::AbstractVector{<:Real})
            assign_x_ss!(mdl, x)
            assign_u_ss!(mdl, u)
            f_ode!(mdl)
            get_ẋ_ss(mdl)
        end
    end
end

#build state space system's output function
function get_h_ss(mdl::Model{Vehicle})
    let mdl = mdl
        function (x::AbstractVector{<:Real}, u::AbstractVector{<:Real})
            assign_x_ss!(mdl, x)
            assign_u_ss!(mdl, u)
            f_ode!(mdl)
            get_y_ss(mdl)
        end
    end
end


################################################################################
################################ Linearization #################################

function Linearization.linearize(mdl::Model{Vehicle}, ip::InitParameters = InitParameters())

    init!(mdl, ip)

    #state space system's functions
    f = get_f_ss(mdl)
    h = get_h_ss(mdl)

    #linearization point
    x0 = get_x_ss(mdl)
    u0 = get_u_ss(mdl)

    lss = linearize(f, h, x0, u0)

    #restore Vehicle state
    init!(mdl, ip)

    return lss

end


################################################################################
################################# Controller ###################################

#state vector for velocity LQR
@kwdef struct XController <: FieldVector{3, Float64}
    ω::Float64 = 0.0
    v::Float64 = 0.0
    θ::Float64 = 0.0
end

@enum ControlMode mode_m = 0 mode_v = 1 mode_η = 2

@kwdef struct Controller <: ModelDefinition
    v2m::LQR{3,1,1,3,1} = LQR{3,1,1}() #velocity controller
    η2v::PID = PID() #position controller
    v_lim::Ref{Float64} = Ref(0.0)
end

@kwdef mutable struct ControllerU
    mode::ControlMode = mode_v
    m_ref::Ranged{Float64, -1., 1.} = 0.0 #motor command reference
    v_ref::Float64 = 0.0 #velocity reference
    η_ref::Float64 = 0.0 #position reference
end

@kwdef struct ControllerY
    mode::ControlMode = mode_v
    m_ref::Ranged{Float64,-1.,1.} = 0.0 #motor command reference
    v_ref::Float64 = 0.0 #velocity reference
    η_ref::Float64 = 0.0 #position reference
    m_cmd::Ranged{Float64,-1.,1.} = 0.0 #motor command output
    v2m::LQROutput{3,1,1,3,1} = LQROutput{3,1,1}()
    η2v::PIDOutput = PIDOutput()
end

Modeling.U(::Controller) = ControllerU()
Modeling.Y(::Controller) = ControllerY()


function Modeling.f_periodic!(::NoScheduling, mdl::Model{<:Controller}, vehicle::Model{<:Vehicle})

    (; v2m, η2v) = mdl.submodels
    (; v_lim) = mdl.parameters
    (; mode, m_ref, v_ref, η_ref) = mdl.u
    (; ω, v, θ, η) = vehicle.y

    #external motor command
    m_cmd = m_ref

    #position mode active
    if mode === mode_η
        η2v.u.input = η_ref - η #assign position error
        f_periodic!(η2v) #update position control submodel
        v_ref = η2v.y.output #override external velocity reference
    end

    #position or velocity mode active
    if (mode === mode_v) || (mode === mode_η)
        v2m.u.x .= XController(; ω, v, θ) #state vector value
        v2m.u.z .= v #command vector value
        v2m.u.z_ref .= clamp(v_ref, -v_lim[], v_lim[]) #command vector reference
        f_periodic!(v2m) #update velocity control submodel
        m_cmd = v2m.y.output[1] #override external motor command reference
    end

    mdl.y = ControllerY(; mode, m_ref, v_ref, η_ref, m_cmd, v2m = v2m.y, η2v = η2v.y)

end


function Modeling.init!(mdl::Model{<:Controller}, vehicle::Model{<:Vehicle})

    (; v2m, η2v) = mdl.submodels
    (; k_m, b_m, R) = vehicle.parameters

    v_max = k_m * R / b_m #maximum steady-state velocity
    mdl.parameters.v_lim[] = 0.4v_max #maximum velocity reference

    #compute and assign velocity controller gains
    Control.Discrete.assign!(v2m, design_v2m_lqr(vehicle))

    #set velocity controller output bounds
    v2m.parameters.bound_lo .= -1 #minimum motor command
    v2m.parameters.bound_hi .= 1 #maximum motor command

    #assign position controller gains
    η2v.parameters.k_p[] = 0.6
    η2v.parameters.k_i[] = 0.0
    η2v.parameters.k_d[] = 0.0

    #set position controller output bounds
    η2v.parameters.bound_lo[] = -mdl.parameters.v_lim[]
    η2v.parameters.bound_hi[] = mdl.parameters.v_lim[]

    #reset model inputs
    mdl.u.mode = mode_v
    mdl.u.m_ref = 0
    mdl.u.v_ref = 0
    mdl.u.η_ref = 0

    #reset submodels
    foreach(values(mdl.submodels)) do ss
        init!(ss)
    end


end


function design_v2m_lqr(mdl::Model{<:Vehicle})

    lss = linearize(mdl)
    lss = delete_vars(lss, :η) #build reduced design model

    x_trim = lss.x0
    n_x = length(x_trim)
    x_labels = collect(keys(x_trim))
    @assert tuple(x_labels...) === propertynames(Robot2D.XController())

    u_labels = [:m, ]
    u_trim = lss.u0[u_labels]

    z_labels = [:v, ]
    z_trim = lss.y0[z_labels]

    A = lss.A
    B = lss.B
    C = lss.C[z_labels, :]
    D = lss.D[z_labels, :]

    C_int = C[z_labels, :]
    D_int = D[z_labels, :]
    n_int = size(C_int, 1)

    A_aug = [A zeros(size(A, 1), n_int); C_int zeros(n_int, n_int)]
    B_aug = [B; D_int]
    C_aug = [C zeros(size(C, 1), n_int)]
    D_aug = D

    P_aug = ss(A_aug, B_aug, C_aug, D_aug)

    x_aug_labels = push!(copy(x_labels), :ξ_v)
    Q_diag = ComponentVector(zeros(length(x_aug_labels)), Axis(x_aug_labels))
    R_diag = ComponentVector(zeros(length(u_labels)), Axis(u_labels))

    Q_diag.ω = 1e-3
    Q_diag.v = 1e-2
    Q_diag.θ = 0
    Q_diag.ξ_v = 5e-2

    R_diag.m = 1e-1

    Q = diagm(Q_diag)
    R = diagm(R_diag)

    #compute augmented feedback matrix
    K_aug = lqr(P_aug, Q, R)

    L = [A B; C D]
    M = inv(L)
    M_12 = M[1:n_x, n_x+1:end]
    M_22 = M[n_x+1:end, n_x+1:end]

    #extract system state and integrator blocks from the feedback matrix
    K_x = K_aug[:, 1:n_x]
    K_ξ = K_aug[:, n_x+1:end]

    K_fbk = K_x
    K_fwd = M_22 + K_x * M_12
    K_int = K_ξ

    return LQRDataPoint(; K_fbk, K_fwd, K_int, x_trim, u_trim, z_trim)

end

@no_ode Controller
@no_step Controller


function GUI.draw!(mdl::Model{<:Controller}, vehicle::Model{<:Vehicle}, p_open::Ref{Bool} = Ref(true))

    (; u, y) = mdl
    (; v2m, η2v) = mdl.submodels
    (; v_lim) = mdl.parameters
    (; v, η, u_m) = vehicle.y


    Begin("Controller", p_open)
    if BeginTable("Controller", 3, CImGui.ImGuiTableFlags_SizingStretchProp)# | CImGui.ImGuiTableFlags_Resizable)# | CImGui.ImGuiTableFlags_BordersInner)

        TableNextRow()
        TableNextColumn()
        AlignTextToFramePadding()
        Text("Control Mode")
        TableNextColumn()
        mode_button("Direct", mode_m, u.mode, y.mode)
        IsItemActive() && (u.mode = mode_m; u.m_ref = u_m)
        SameLine()
        mode_button("Velocity", mode_v, u.mode, y.mode)
        IsItemActive() && (u.mode = mode_v; u.v_ref = v)
        SameLine()
        mode_button("Position", mode_η, u.mode, y.mode)
        IsItemActive() && (u.mode = mode_η; u.η_ref = η)

        TableNextRow()
        TableNextColumn()
        AlignTextToFramePadding()
        Text("Motor Command")
        TableNextColumn()
        PushItemWidth(-10)
        u.m_ref = safe_slider("##Motor Reference", u.m_ref, "%.6f")
        PopItemWidth()
        TableNextColumn()
        Text(@sprintf("%.6f", u_m))

        TableNextRow()
        TableNextColumn()
        AlignTextToFramePadding()
        Text("Velocity (m/s)")
        TableNextColumn()
        PushItemWidth(-10)
        u.v_ref = safe_slider("##Velocity Reference", u.v_ref, -v_lim[], v_lim[], "%.6f")
        PopItemWidth()
        TableNextColumn()
        Text(@sprintf("%.6f", v))

        TableNextRow()
        TableNextColumn()
        AlignTextToFramePadding()
        Text("Position (m)")
        TableNextColumn()
        PushItemWidth(-10)
        u.η_ref = safe_slider("##Position Reference", u.η_ref, -5, 5, "%.6f")
        PopItemWidth()
        TableNextColumn()
        Text(@sprintf("%.6f", η))

        EndTable()

        Separator()

        CollapsingHeader("Velocity Controller") && GUI.draw(v2m)
        CollapsingHeader("Position Controller") && GUI.draw(η2v)
    end

    End()

end


################################################################################
################################# Robot ########################################

@kwdef struct Robot{C} <: ModelDefinition
    vehicle::Vehicle = Vehicle()
    controller::C = Controller()
end

@kwdef struct LostBalance <: Exception
    msg::String = ""
end

function Modeling.f_ode!( mdl::Model{<:Robot})

    f_ode!(mdl.vehicle)
    f_output!(mdl)

end

function Modeling.f_periodic!(::NoScheduling, mdl::Model{<:Robot})

    (; vehicle, controller) = mdl.submodels
    f_periodic!(controller, vehicle)
    vehicle.u[] = controller.y.m_cmd
    f_output!(mdl)

end

function Modeling.f_step!(mdl::Model{<:Robot})

    θ_max = deg2rad(45)
    (; θ) = mdl.vehicle.y
    (abs(θ) > θ_max) && throw(LostBalance(
        "Chassis tilt θ = $(rad2deg(θ)) deg exceeded allowable limit (θ = $(rad2deg(θ_max)))" *
        "at t = $(mdl.t[]) s"))

end

function Modeling.init!( mdl::Model{<:Robot}, ip::InitParameters = InitParameters())

    (; vehicle, controller) = mdl.submodels
    init!(vehicle, ip)
    init!(controller, vehicle) #run controller design and assign gains
    f_output!(mdl) #update parent model output

end

function GUI.draw!(mdl::Model{<:Robot}, p_open::Ref{Bool} = Ref(true))

    (; vehicle, controller) = mdl

    Begin("Robot", p_open)
    @cstatic c_veh=false c_ctl=false begin
        @c Checkbox("Vehicle", &c_veh)
        c_veh && @c GUI.draw!(vehicle, &c_veh)
        @c Checkbox("Controller", &c_ctl)
        c_ctl && @c GUI.draw!(controller, vehicle, &c_ctl)
    end
    End()

end

end #module
