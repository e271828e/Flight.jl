#Alternative to sketch.jl using derivative binding instead of the explicit
#computer/integrator component split. Same illustrative, non-committed syntax.
#
#The feature: a component may declare that some (or all) fields of its ẋ are
#bound to designated output ports of its own, instead of computing them in an
#f method. Bindings are applied by the framework after the output sweep, when
#the component's own outputs are guaranteed fresh. Mixing is allowed: bind some
#fields, compute the rest in f. Here, neither component needs an f at all, and
#each derivative-producing computation runs exactly once.
#
#Compared with sketch.jl: two components instead of four, eight connections
#instead of thirteen, and everything derivable from pose alone migrates from
#stage 2 to stage 1 (pose is now WA's own state, not an input), so ne's g_s2
#only waits on the wrench chain.

######################## NewtonEuler (merged) ##################################

struct NewtonEuler <: AbstractComponent end

init_x(::NewtonEuler) = (ω_eb_b = zeros(SVector{3}), v_eb_b = zeros(SVector{3}))
init_z(::NewtonEuler) = nothing

g_s1(::NewtonEuler, x, m) = x   #ω_eb_b, v_eb_b published as stage-1 outputs

function g_s2(::NewtonEuler, x, u, m)

    (; mp_Σ_b, wr_Σ_b, ho_Σ_b, q_eb, r_eb_e) = u   #velocities no longer inputs:
    (; ω_eb_b, v_eb_b) = x                          #read own state directly

    ω_ie_e = SVector{3}(0, 0, ω_ie) #use WGS84 constant
    ω_ie_b = q_eb'(ω_ie_e)

    #frame transform from c (CoM) to b (body)
    r_bc_b = mp_Σ_b.r_OG
    t_cb = FrameTransform(r = -r_bc_b) #pure translation

    ###################### angular & linear momentum equations #################

    #translate data to frame c
    mp_Σ_c = t_cb(mp_Σ_b)
    wr_Σ_c = t_cb(wr_Σ_b)
    ho_Σ_c = ho_Σ_b

    F_Σ_c = wr_Σ_c.F; τ_Σ_c = wr_Σ_c.τ
    m_Σ = mp_Σ_c.m; J_Σ_c = mp_Σ_c.J

    ω_ec_c = ω_eb_b
    v_ec_c = v_eb_b + ω_ec_c × r_bc_b

    ω_ie_c = ω_ie_b
    ω_ic_c = ω_ie_c + ω_ec_c

    #compute geographic position of Oc
    r_bc_e = q_eb(r_bc_b)
    r_ec_e = r_eb_e + r_bc_e
    Oc = Cartesian(r_ec_e)

    #define auxiliary local-level frame l with Ol = Oc
    q_el = ltf(Oc)
    q_cl = q_eb' ∘ q_el

    #compute gravity at c
    g_c_c = q_cl(SVector{3}(0, 0, gravity(Oc)))

    #solve dynamic equations at c
    hc_Σ_c = J_Σ_c * ω_ic_c + ho_Σ_c
    ω̇_ec_c = J_Σ_c \ (τ_Σ_c - J_Σ_c * (ω_ie_c × ω_ec_c) - ω_ic_c × hc_Σ_c)
    v̇_ec_c = 1/m_Σ * F_Σ_c + g_c_c - (ω_ec_c + 2ω_ie_c) × v_ec_c

    #translate derivatives back to b
    ω̇_eb_b = ω̇_ec_c
    v̇_eb_b = v̇_ec_c - ω̇_ec_c × r_bc_b

    return NewtonEulerOutput(; ω̇_eb_b, v̇_eb_b)
    #plus the remaining acceleration-like quantities (a_eb_b, f_c_c, ...),
    #published as ports for accelerometer-style consumers

end

#no f method — Newton-Euler is solved exactly once, in g_s2:
ẋ_bindings(::NewtonEuler) = (ω_eb_b = :ω̇_eb_b, v_eb_b = :v̇_eb_b)


########################### WA (merged) ########################################

struct WA <: AbstractComponent end

init_x(::WA) = (q_wb = RQuat(), q_ew = RQuat(), h_e = HEllip())
init_z(::WA) = nothing

function g_s1(::WA, x, m)   #everything derivable from pose alone: now stage 1

    (; q_wb, q_ew, h_e) = x

    ψ_nw = get_ψ_nw(q_ew)
    q_nw = Rz(ψ_nw)
    q_nb = q_nw ∘ q_wb
    q_eb = q_ew ∘ q_wb
    e_nb = REuler(q_nb)

    n_e = NVector(q_ew)
    ϕ_λ = LatLon(n_e)
    h_o = HOrth(h_e, n_e)
    Ob = Geographic(n_e, h_e)
    r_eb_e = Cartesian(Ob)

    return KinPose(; q_wb, q_ew, q_nb, q_eb, e_nb, n_e, ϕ_λ, h_e, h_o, r_eb_e)

end

function g_s2(::WA, x, u, m)

    (; ω_eb_b, v_eb_b) = u          #from NewtonEuler's stage-1 outputs
    (; q_wb, q_ew, h_e) = x

    #cheap pose recomputation (shared intermediates across groups are either
    #recomputed, as here, or hoisted; see framework_design.md §5.2)
    ψ_nw = get_ψ_nw(q_ew)
    q_nw = Rz(ψ_nw)
    q_nb = q_nw ∘ q_wb

    n_e = NVector(q_ew)
    Ob = Geographic(n_e, h_e)

    v_eb_n = q_nb(v_eb_b)
    ω_ew_n = get_ω_ew_n(v_eb_n, Ob)
    ω_ew_w = q_nw'(ω_ew_n)
    ω_ew_b = q_wb'(ω_ew_w)
    ω_wb_b = ω_eb_b - ω_ew_b

    v_gnd, χ_gnd, γ_gnd = get_vχγ(v_eb_n)

    q̇_wb = Attitude.dt(q_wb, ω_wb_b)
    q̇_ew = Attitude.dt(q_ew, ω_ew_w)
    ḣ_e = -v_eb_n[3]

    return KinVel(; v_eb_n, ω_wb_b, ω_ew_b, v_gnd, χ_gnd, γ_gnd, q̇_wb, q̇_ew, ḣ_e)

end

ẋ_bindings(::WA) = (q_wb = :q̇_wb, q_ew = :q̇_ew, h_e = :ḣ_e)

#projection follows the state: quaternion renormalization lives here
project(::WA, x) = (q_wb = normalize(x.q_wb), q_ew = normalize(x.q_ew), h_e = x.h_e)


############################## Assembly ########################################

asm = Assembly()
add!(asm, :ne,  NewtonEuler())
add!(asm, :wa,  WA())
add!(asm, :sys, Systems())

#connect!(assembly, destination path, origin path)
connect!(asm, (:wa, :u, :ω_eb_b), (:ne, :y, :ω_eb_b))      #stage 1 → stage 2
connect!(asm, (:wa, :u, :v_eb_b), (:ne, :y, :v_eb_b))      #stage 1 → stage 2

connect!(asm, (:ne, :u, :q_eb),   (:wa, :y, :q_eb))        #stage 1 → stage 2
connect!(asm, (:ne, :u, :r_eb_e), (:wa, :y, :r_eb_e))      #stage 1 → stage 2

connect!(asm, (:ne, :u, :wr_Σ_b), (:sys, :y, :wr_Σ_b))     #reduce-port (+)
connect!(asm, (:ne, :u, :mp_Σ_b), (:sys, :y, :mp_Σ_b))     #reduce-port (+)
connect!(asm, (:ne, :u, :ho_Σ_b), (:sys, :y, :ho_Σ_b))     #reduce-port (+)

connect!(asm, (:sys, :u, :v_eb_n), (:wa, :y, :v_eb_n))

#derived schedule:
#  s1: ne.g_s1, wa.g_s1, sys stage-1        (any order)
#  s2: wa.g_s2 → sys.g_s2 → ne.g_s2         (topological)
#  post-sweep: ẋ bindings, projection on wa
