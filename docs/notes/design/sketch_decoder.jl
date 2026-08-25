#Model-side sketch: the rigid-body kinematics/dynamics core spelled against the
#settled design (framework_design.md v0.19) — two-stage interfaces (§5.2), the
#declaration layer (§13), explicit summing junctions (§6) and type-based
#assemblies (§13.5-§13.8). Illustrative, non-committed syntax; not runnable.
#
#UNCHECKED against the current spec: unlike condition_demo.jl, this sketch was
#never remapped after the framework_spec.md renumberings, so every § below is
#presumed stale — cited section numbers are the v0.19-era design doc's.
#
#Interfaces (continuous component), in evaluation order. Every function takes
#(comp, args): one NamedTuple bundle of views, destructured by name — the
#author names only what the body reads (§5.2, the bundle law). Maximal legal
#sets, ambient t (and ws where declared) riding unnamed in the suffix:
#
#   y_x      = h_x(comp, args)      # ⊆ (; x, m, t, ws) — no-feedthrough stage
#   y_xu     = h_xu(comp, args)     # + u, y_x — all wired inputs
#   ẋ        = f(comp, args)        # ⊆ (; x, m, y, u, t, ws); y complete, fresh
#
#at step boundaries (none used in this sketch):
#
#   fired    = guard(comp, args)    # ⊆ (; x, m, y, u, t)
#   (x⁺, m⁺) = handler(comp, args)  #then: project → re-run the output stages
#   x⁺       = project(comp, x)     #manifold projection; positional (§5.2)
#
#A declared output that names a state field no stage produces is auto-published
#at stage-1 position (§5.2) — no identity-decoder boilerplate anywhere below.
#Every derivative-producing computation runs exactly once, in h_xu; f is a copy
#from the fresh table. The explicit computer/integrator split (a stateless
#derivative-computing component wired into a trivial state holder — the earlier
#sketch.jl form, now retired) remains expressible without framework support and
#is the idiom of choice when the factoring earns reuse (§5.2, §9.4).
#
#Port types are spelled at their Float64 nominal faces (§13.2). The migration's
#Tier-1 parametrization pass (§9.2) makes the underlying structs parametric so
#the activation leaf walk can re-scalar them (Wrench{Float64} → Wrench{Dual});
#the declarations themselves never mention the activation scalar (v0.20).

######################## NewtonEuler ###########################################

struct NewtonEuler <: AbstractComponent end

init_x(::NewtonEuler) = (ω_eb_b = zeros(SVector{3}), v_eb_b = zeros(SVector{3}))

#input contract: bare types at their Float64 faces (§13.2)
input_types(::NewtonEuler) = (
    mp_Σ_b = MassProperties{Float64},   #total mass properties, body frame
    wr_Σ_b = Wrench{Float64},           #total external wrench, body frame
    ho_Σ_b = SVector{3, Float64},       #total internal angular momentum, body frame
    q_eb   = RQuat{Float64},
    r_eb_e = SVector{3, Float64},
)

#output contract = the public interface, at concrete nominal types (§13.2).
#ω_eb_b / v_eb_b name state fields no stage produces → auto-published (stage 1).
#The accelerations are genuine interface — accelerometer-style consumers read
#them — so they live in outputs, not locals.
output_types(::NewtonEuler) = (
    ω_eb_b = SVector{3, Float64},
    v_eb_b = SVector{3, Float64},
    ω̇_eb_b = SVector{3, Float64},
    v̇_eb_b = SVector{3, Float64},
    #plus the remaining acceleration-like quantities (a_eb_b, f_c_c, ...)
)

function h_xu(::NewtonEuler, (; x, u))

    (; mp_Σ_b, wr_Σ_b, ho_Σ_b, q_eb, r_eb_e) = u
    (; ω_eb_b, v_eb_b) = x          #own state, direct view (§5.2)

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

    return (; ω̇_eb_b, v̇_eb_b)

end

#Newton-Euler is solved exactly once, in h_xu; f copies from the fresh table:
f(::NewtonEuler, (; y)) = (ω_eb_b = y.ω̇_eb_b, v_eb_b = y.v̇_eb_b)


########################### WA #################################################

struct WA <: AbstractComponent end

init_x(::WA) = (q_wb = RQuat(), q_ew = RQuat(), h_e = HEllip())

input_types(::WA) = (ω_eb_b = SVector{3, Float64}, v_eb_b = SVector{3, Float64})

#q_wb / q_ew / h_e: state fields no stage produces → auto-published (stage 1);
#the rest of stage-1 membership (pose-derived vs velocity-dependent) is derived
#by the build, invisible in the contract (§13.2 — no stage tags)
output_types(::WA) = (
    q_wb = RQuat{Float64}, q_ew = RQuat{Float64}, q_nw = RQuat{Float64},
    q_nb = RQuat{Float64}, q_eb = RQuat{Float64}, e_nb = REuler{Float64},
    n_e = NVector{Float64}, ϕ_λ = LatLon{Float64}, h_e = HEllip{Float64},
    h_o = HOrth{Float64}, Ob = Geographic{Float64}, r_eb_e = SVector{3, Float64},
    v_eb_n = SVector{3, Float64}, ω_wb_b = SVector{3, Float64},
    ω_ew_b = SVector{3, Float64},
    v_gnd = Float64, χ_gnd = Float64, γ_gnd = Float64,
)

#consumed only by WA's own f — component-local cross-stage cells, not
#interface (§13.2). q̇'s are quaternion-backing derivatives (4 scalars).
local_types(::WA) = (
    q̇_wb = SVector{4, Float64}, q̇_ew = SVector{4, Float64}, ḣ_e = Float64,
)

function h_x(::WA, (; x))   #everything derivable from pose alone: stage 1

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

    return (; q_nw, q_nb, q_eb, e_nb, n_e, ϕ_λ, h_o, Ob, r_eb_e)

end

function h_xu(::WA, (; x, u, y_x))

    (; ω_eb_b, v_eb_b) = u              #from NewtonEuler's auto-published state
    (; q_nw, q_nb, Ob) = y_x           #own stage-1 results — nothing recomputed
    (; q_wb, q_ew) = x                  #own state, direct view

    v_eb_n = q_nb(v_eb_b)
    ω_ew_n = get_ω_ew_n(v_eb_n, Ob)
    ω_ew_w = q_nw'(ω_ew_n)
    ω_ew_b = q_wb'(ω_ew_w)
    ω_wb_b = ω_eb_b - ω_ew_b

    v_gnd, χ_gnd, γ_gnd = get_vχγ(v_eb_n)

    q̇_wb = Attitude.dt(q_wb, ω_wb_b)
    q̇_ew = Attitude.dt(q_ew, ω_ew_w)
    ḣ_e = -v_eb_n[3]

    return (; v_eb_n, ω_wb_b, ω_ew_b, v_gnd, χ_gnd, γ_gnd, q̇_wb, q̇_ew, ḣ_e)

end

f(::WA, (; y)) = (q_wb = y.q̇_wb, q_ew = y.q̇_ew, h_e = y.ḣ_e)

#projection follows the state: quaternion renormalization lives here
project(::WA, x) = (q_wb = normalize(x.q_wb), q_ew = normalize(x.q_ew), h_e = x.h_e)


########################### Systems ############################################
#Toy force contributors — declarations shown, bodies elided. Contributors are
#ragged by design (§6): Aero has wrench but no mass, only PWP has angular
#momentum; nobody publishes zero-filled identity noise.

struct Aero <: AbstractComponent end

input_types(::Aero) = (v_eb_n = SVector{3, Float64},)
output_types(::Aero) = (wr_b = Wrench{Float64},)
#h_xu body elided

struct PWP <: AbstractComponent end

init_x(::PWP) = (ω_prop = 0.0,)
input_types(::PWP) = (throttle = Float64,)
output_types(::PWP) = (wr_b = Wrench{Float64}, ho_b = SVector{3, Float64})
#h_xu, f bodies elided

struct MassStore <: AbstractComponent end   #placeholder for airframe/fuel mass

output_types(::MassStore) = (mp_b = MassProperties{Float64},)
#h_x body elided (no inputs, no feedthrough: stage 1)

#Named site-specific summing junction (§6): the contributor set documented in
#the contract; an ordinary component, no framework privileges. The generic
#positional SumJunction{V, N} is the library tool for configuration-variable
#sites.
struct WrenchSum <: AbstractComponent end

input_types(::WrenchSum) = (aero = Wrench{Float64}, pwp = Wrench{Float64})
output_types(::WrenchSum) = (; Σ = Wrench{Float64})
h_xu(::WrenchSum, (; u)) = (; Σ = u.aero + u.pwp)

#Assembly = plain struct; component-typed fields are the children (field names
#= path segments), connections is the kind marker (§13.5). Where two children
#contribute, a junction sums them; where one does, the export is a direct
#wire — the hierarchical aggregation idiom (§6).
struct Systems <: AbstractComponent
    aero::Aero
    pwp::PWP
    mass::MassStore
    wr_sum::WrenchSum
end

connections(::Systems) = (
    "aero/wr_b" => "wr_sum/aero",
    "pwp/wr_b"  => "wr_sum/pwp",
)

exports(::Systems) = (
    "wr_Σ_b"  => "wr_sum/Σ",            #aggregate: a first-class signal
    "mp_Σ_b"  => "mass/mp_b",           #single contributor: direct export
    "ho_Σ_b"  => "pwp/ho_b",            #single contributor: direct export
    "v_eb_n"  => "aero/v_eb_n",         #input face: fed by the parent's wire
    "throttle" => "pwp/throttle",       #input face
)


########################### Vehicle (root) #####################################
#A parametric child slot (sys::S, today's Cessna172X{K, A} shape) is the
#substitutability idiom (§13.5); concrete here for brevity.

struct Vehicle <: AbstractComponent
    ne::NewtonEuler
    wa::WA
    sys::Systems
end

connections(::Vehicle) = (
    #"src/port" => "dst/port", child-port → child-port (§13.6)
    "ne/ω_eb_b"   => "wa/ω_eb_b",       #stage 1 → stage 2
    "ne/v_eb_b"   => "wa/v_eb_b",       #stage 1 → stage 2
    "wa/q_eb"     => "ne/q_eb",         #stage 1 → stage 2
    "wa/r_eb_e"   => "ne/r_eb_e",       #stage 1 → stage 2
    "wa/v_eb_n"   => "sys/v_eb_n",
    "sys/wr_Σ_b"  => "ne/wr_Σ_b",
    "sys/mp_Σ_b"  => "ne/mp_Σ_b",
    "sys/ho_Σ_b"  => "ne/ho_Σ_b",
)

#At the root, an exported input face IS the write surface: "throttle" becomes
#a root slot with no extra vocabulary (§13.6). Face names are contract tokens
#(no slash); bulk re-export with prefixing is the faces(...) helper (§13.8).
exports(::Vehicle) = (
    "throttle" => "sys/throttle",       #input face → root slot
    "pose"     => "wa/Ob",              #output faces: the observation contract
    "v_gnd"    => "wa/v_gnd",
)

#derived schedule (§5.2):
#  s1: auto-publications (ne, wa state) + wa.h_x, mass.h_x        (any order)
#  s2: wa.h_xu → {aero.h_xu, pwp.h_xu} → wr_sum.h_xu → ne.h_xu (topological)
#  post-sweep: all f against the complete signal table; projection on wa
