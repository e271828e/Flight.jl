struct NewtonEuler <: AbstractComponent end

init_x(::NewtonEuler) = nothing
init_z(::NewtonEuler) = nothing

function g_s2(::NewtonEuler, x, u, m)

    (; mp_Σ_b, wr_Σ_b, ho_Σ_b, q_eb, r_eb_e, ω_eb_b, v_eb_b) = u

    ω_ie_e = SVector{3, Float64}(0, 0, ω_ie) #use WGS84 constant
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
    q_be = q_eb'
    q_ce = q_be
    q_cl = q_ce ∘ q_el

    #compute gravity at c
    g_c_l = SVector{3,Float64}(0, 0, gravity(Oc)) #gravity at c, l coordinates\
    g_c_c = q_cl(g_c_l) #gravity at c, c coordinates

    #solve dynamic equations at c
    hc_Σ_c = J_Σ_c * ω_ic_c + ho_Σ_c
    ω̇_ec_c = J_Σ_c \ (τ_Σ_c - J_Σ_c * (ω_ie_c × ω_ec_c) - ω_ic_c × hc_Σ_c)
    v̇_ec_c = 1/m_Σ * F_Σ_c + g_c_c - (ω_ec_c + 2ω_ie_c) × v_ec_c

    #translate derivatives back to b
    ω̇_eb_b = ω̇_ec_c
    v̇_eb_b = v̇_ec_c - ω̇_ec_c × r_bc_b

    return NewtonEulerOutput(ω̇_eb_b = ω̇_eb_b, v̇_eb_b = v̇_eb_b)
    #NewtonEulerOutput would also include the rest of acceleration-like quantities

end


######################

struct NewtonEulerIntegrator <: AbstractComponent end

init_x(::NewtonEulerIntegrator) = (ω_eb_b = zeros(SVector{3}), v_eb_b = zeros(SVector{3}))
init_z(::NewtonEulerIntegrator) = nothing

f(::NewtonEulerIntegrator, x, u, m, t) = (ω_eb_b = u.ω̇_eb_b, v_eb_b = u.v̇_eb_b)
g_s1(::NewtonEulerIntegrator, x, m) = x


#######################

struct WA <: AbstractComponent end

init_x(::WA) = nothing
init_z(::WA) = nothing

function g_s2(::WA, x, u, m)

    (; ω_eb_b, v_eb_b, q_wb, q_ew, h_e) = u

    ψ_nw = get_ψ_nw(q_ew)
    q_nw = Rz(ψ_nw)
    q_nb = q_nw ∘ q_wb
    q_eb = q_ew ∘ q_wb
    q_en = q_eb ∘ q_nb'
    e_nb = REuler(q_nb)

    n_e = NVector(q_ew)
    ϕ_λ = LatLon(n_e)
    h_o = HOrth(h_e, n_e)

    v_eb_n = q_nb(v_eb_b)
    Ob = Geographic(n_e, h_e)
    r_eb_e = Cartesian(Ob)
    ω_ew_n = get_ω_ew_n(v_eb_n, Ob)

    ω_ew_w = q_nw'(ω_ew_n)
    ω_ew_b = q_wb'(ω_ew_w)
    ω_wb_b = ω_eb_b - ω_ew_b

    v_gnd, χ_gnd, γ_gnd = get_vχγ(v_eb_n)

    q̇_wb = Attitude.dt(q_wb, ω_wb_b)
    q̇_ew = Attitude.dt(q_ew, ω_ew_w)
    ḣ_e = -v_eb_n[3]

    return KinOutput( q̇_wb = q̇_wb, q̇_ew = q̇_ew, ḣ_e = ḣ_e)
    #KinOutput would also include all the other velocity, position and attitude quantities

end


struct WAIntegrator <: AbstractComponent end

init_x(::WAIntegrator) = (q_wb = RQuat(), q_ew = RQuat(), h_e = HEllip())
init_z(::WAIntegrator) = nothing

f(::WAIntegrator, x, u, m, t) = (q_wb = u.q̇_wb, q_ew = u.q̇_ew, h_e = u.ḣ_e)
g_s1(::WAIntegrator, x, m) = x


#####################################

asm = Assembly()
add!(asm, :ne, NewtonEuler())
add!(asm, :ne_int, NewtonEulerIntegrator())
add!(asm, :wa, WA())
add!(asm, :wa_int, WAIntegrator())
add!(asm, :sys, Systems())

#connect!(assembly, destination path, origin path)
connect!(asm, (:wa_int, :u, :q̇_wb), (:wa, :y, :q̇_wb))
connect!(asm, (:wa_int, :u, :q̇_ew), (:wa, :y, :q̇_ew))
connect!(asm, (:wa_int, :u, :ḣ_e), (:wa, :y, :ḣ_e))

connect!(asm, (:wa, :u, :ω_eb_b), (:ne_int, :y, :ω_eb_b))
connect!(asm, (:wa, :u, :v_eb_b), (:ne_int, :y, :v_eb_b))

connect!(asm, (:wa, :u, :q_wb), (:wa_int, :y, :q_wb))
connect!(asm, (:wa, :u, :q_ew), (:wa_int, :y, :q_ew))
connect!(asm, (:wa, :u, :h_e), (:wa_int, :y, :h_e))

connect!(asm, (:ne_int, :u, :ω̇_eb_b), (:ne, :y, :ω̇_eb_b))
connect!(asm, (:ne_int, :u, :v̇_eb_b), (:ne, :y, :v̇_eb_b))

connect!(asm, (:ne, :u, :q_eb), (:wa, :y, :q_eb))
connect!(asm, (:ne, :u, :r_eb_e), (:wa, :y, :r_eb_e))

connect!(asm, (:ne, :u, :wr_Σ_b), (:sys, :y, :wr_Σ_b))

connect!(asm, (:ne, :u, :ω_eb_b), (:ne_int, :y, :ω_eb_b))
connect!(asm, (:ne, :u, :v_eb_b), (:ne_int, :y, :v_eb_b))

connect!(asm, (:sys, :u, :v_eb_n), (:wa, :y, :v_eb_n)) #
