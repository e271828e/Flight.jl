#Runnable demo of the §14 condition algebra (framework_spec.md):
#fragment/at/merge/override lazy trees, flattening with provenance, the
#duplicate-leaf error, and §14.9 mounting — a TrimProblem relocated whole by
#at(prefix, problem). Written as the inspection aid for the mounting
#semantics: every step prints the tree or the flattened entry list.
#
#Self-contained, no framework, no packages:
#
#   julia docs/notes/design/condition_demo.jl
#
#The node structs and the resolution pass are miniature but semantically
#faithful implementations of §14.2/§14.3/§14.6; the one stand-in is slot
#resolution, which here uses the demo rule «root slot = mount prefix with
#slashes → dots, then the face name» ("throttle" at "wing" → "wing.throttle",
#§14.9's own example) in place of the Build's export-chain walk.
#
#NOT COVERED BY THE CHECK BATTERY: tools/check_refs.jl and tools/linkify.jl
#scan the COMPANIONS .md roster only, so every § here is hand-verified.
#Re-verify by hand after any spec renumbering. Last verified 2026-08-16.

module ConditionDemo

######################## The algebra (§14.2, §14.6) ############################

#x is the fused state store — continuous and discrete leaves alike (D-173);
#there is no separate z payload, z being the shift operator only
struct Fragment{X<:NamedTuple, M<:NamedTuple, S<:NamedTuple}
    x::X; m::M; slots::S            #self-vocabulary payloads; no paths
end
fragment(; x = (;), m = (;), slots = (;)) = Fragment(x, m, slots)

struct Scoped{N}                    #at(prefix, node): stores, never applies
    prefix::String
    node::N
end
at(prefix::AbstractString, node) = Scoped(String(prefix), node)

struct Merged{T<:Tuple}             #merge(ns...): collects; order = diagnostics only
    nodes::T
end
merge(nodes...) = Merged(nodes)

struct Override{B, P}               #override(base, patch): ordered, asymmetric (§14.6)
    base::B
    patch::P
end
override(base, patch) = Override(base, patch)
override(base, patch, rest...) = override(Override(base, patch), rest...)

######################## Resolution: flatten → validate (§14.3) ################

struct Entry
    path::String                    #component path from the resolution root
    store::Symbol                   #:x, :m, or :slot
    field::Symbol                   #state/mode field, or face name for :slot
    value::Any
    provenance::String
end

entrykey(e::Entry) = (e.path, e.store, e.field)

struct DuplicateLeafError <: Exception
    a::Entry
    b::Entry
end

function Base.showerror(io::IO, e::DuplicateLeafError)
    print(io, "DuplicateLeafError: two entries for ", loc_string(e.a),
          " in one layer:\n    1. ", e.a.provenance, " = ", repr(e.a.value),
          "\n    2. ", e.b.provenance, " = ", repr(e.b.value))
end

function check_collisions(entries)
    seen = Dict{Any, Entry}()
    for e in entries
        haskey(seen, entrykey(e)) && throw(DuplicateLeafError(seen[entrykey(e)], e))
        seen[entrykey(e)] = e
    end
    entries
end

#flattening is the only place path strings are concatenated (§14.3)
join_path(prefix, p) = isempty(prefix) ? p : prefix * "/" * p

function flatten!(out, f::Fragment, path, prov)
    prov *= "fragment"
    for store in (:x, :m)
        for (field, value) in pairs(getfield(f, store))
            push!(out, Entry(path, store, field, value, prov))
        end
    end
    for (face, value) in pairs(f.slots)     #face vocabulary of the authoring level
        push!(out, Entry(path, :slot, face, value, prov))
    end
end

flatten!(out, s::Scoped, path, prov) =
    flatten!(out, s.node, join_path(path, s.prefix), prov * "at(\"$(s.prefix)\") → ")

function flatten!(out, m::Merged, path, prov)
    for (i, n) in enumerate(m.nodes)
        flatten!(out, n, path, prov * "merge[$i] → ")
    end
end

#each override side is a layer: collision-checked independently; a leaf present
#in both takes the patch's value, provenance recording both sources (§14.6)
function flatten!(out, o::Override, path, prov)
    base = Entry[]; flatten!(base, o.base, path, prov * "base → ")
    check_collisions(base)
    patch = Entry[]; flatten!(patch, o.patch, path, prov * "patch → ")
    check_collisions(patch)
    bykey = Dict(entrykey(p) => p for p in patch)
    overridden = Set{Any}()
    for b in base
        if haskey(bykey, entrykey(b))
            p = bykey[entrykey(b)]
            push!(out, Entry(p.path, p.store, p.field, p.value,
                  p.provenance * "  [overrode " * b.provenance * " = " * repr(b.value) * "]"))
            push!(overridden, entrykey(b))
        else
            push!(out, b)
        end
    end
    for p in patch
        entrykey(p) in overridden || push!(out, p)
    end
end

function resolve(node)      #the §14.3 flatten + the validation batch's collision check
    out = Entry[]
    flatten!(out, node, "", "")
    check_collisions(out)
end

#demo stand-in for the Build's export-chain walk from the mount point (§14.9);
#the real resolution errors if the face never surfaces at the root
root_slot(e::Entry) = isempty(e.path) ? String(e.field) :
    replace(e.path, "/" => ".") * "." * String(e.field)

######################## Printing ##############################################

function label(f::Fragment)
    parts = String[]
    for store in (:x, :m, :slots)
        nt = getfield(f, store)
        isempty(nt) || push!(parts, "$store = $nt")
    end
    "fragment(" * join(parts, ", ") * ")"
end
label(s::Scoped) = "at \"$(s.prefix)\""
label(::Merged) = "merge"
label(::Override) = "override (base, patch — patch wins)"

children(::Fragment) = ()
children(s::Scoped) = (s.node,)
children(m::Merged) = m.nodes
children(o::Override) = (o.base, o.patch)

function print_tree(node, prefix = "", islast = true, isroot = true)
    if isroot
        println(label(node))
    else
        println(prefix, islast ? "└─ " : "├─ ", label(node))
    end
    ch = children(node)
    childprefix = isroot ? "" : prefix * (islast ? "   " : "│  ")
    for (i, c) in enumerate(ch)
        print_tree(c, childprefix, i == length(ch), false)
    end
end

loc_string(e::Entry) = e.store === :slot ?
    "slot \"" * root_slot(e) * "\"" :
    (isempty(e.path) ? "" : e.path * ":") * "$(e.store).$(e.field)"

function print_entries(entries)
    locs = [loc_string(e) for e in entries]
    vals = [repr(e.value) for e in entries]
    w1 = maximum(length.(locs)); w2 = maximum(length.(vals))
    for (e, l, v) in zip(entries, locs, vals)
        println("  ", rpad(l, w1), " = ", rpad(v, w2), "   ⟵ ", e.provenance)
    end
end

function banner(s)
    println(); println("═"^78); println("  ", s); println("═"^78)
end

######################## Toy aircraft (dispatch targets only) ##################

struct PistonEngine; ω_rated::Float64; end
struct PWP; engine::PistonEngine; end
struct Aero; end
struct C172XSystems; pwp::PWP; aero::Aero; end
struct Cessna172X; sys::C172XSystems; end
struct SimpleWorld; lead::Cessna172X; wing::Cessna172X; end

#fragment functions (§14.2): shipped beside the component, dispatched on it,
#composed by pull — each level speaks its own fields, children and faces
condition(eng::PistonEngine; n_eng) =
    fragment(x = (ω = n_eng * eng.ω_rated,), m = (phase = :running,))

condition(sys::C172XSystems; n_eng, α_a, β_a) = merge(
    at("pwp/engine", condition(sys.pwp.engine; n_eng)),
    at("aero",       fragment(x = (α_filt = α_a, β_filt = β_a))))

#aircraft-shipped baseline (§14.6): full slot coverage, one authoritative home
ready_for_taxi(ac::Cessna172X) = merge(
    at("sys", condition(ac.sys; n_eng = 0.25, α_a = 0.0, β_a = 0.0)),
    fragment(slots = (throttle = 0.0, elevator = 0.0, mixture = 0.5)))

######################## TrimProblem (§14.7, §14.9) ############################

struct TrimProblem{G, C, R, F}
    guess::G; lower::G; upper::G    #same-shaped Float64 NamedTuples
    condition::C                    #condition-valued function over the decisions
    reads::R                        #inert selector data
    residuals::F                    #path-free residual math over the gathered reads
end

struct Deriv; path::String; field::Symbol; end
struct Output; path::String; field::Symbol; end
deriv(path, field) = Deriv(path, field)
output(path, field) = Output(path, field)

#at lifts to problems in five lines (§14.9's lift, less the tolerances field
#this toy problem does not model):
at(prefix::AbstractString, p::TrimProblem) = TrimProblem(
    p.guess, p.lower, p.upper,                 #path-free: pass through
    d -> at(prefix, p.condition(d)),           #post-compose: wrap each returned tree
    at(prefix, p.reads),                       #reads are inert selector data: same Scoped node
    p.residuals)                               #path-free: pass through

describe(s::Deriv, prefix)  = "ẋ[" * join_path(prefix, s.path) * " : $(s.field)]"
describe(s::Output, prefix) = "y[" * join_path(prefix, s.path) * " : $(s.field)]"

print_reads(s::Scoped, prefix = "") = print_reads(s.node, join_path(prefix, s.prefix))
function print_reads(reads::NamedTuple, prefix = "")
    for (name, sel) in pairs(reads)
        println("  ", rpad(string(name), 6), " = ", describe(sel, prefix))
    end
end

#level flight, toy spelling: α tied to θ inside the condition — the analytic
#elimination of §14.7 surviving as plain user math
trim_problem(ac::Cessna172X) = TrimProblem(
    (throttle = 0.5,  θ = 0.05, n_eng = 0.7),
    (throttle = 0.0,  θ = -0.3, n_eng = 0.4),
    (throttle = 1.0,  θ = 0.3,  n_eng = 1.2),
    d -> merge(
        at("sys", condition(ac.sys; n_eng = d.n_eng, α_a = d.θ, β_a = 0.0)),
        at("kin", fragment(x = (θ = d.θ,))),
        fragment(slots = (throttle = d.throttle,))),
    (θ̇ = deriv("kin", :θ), ω̇ = deriv("sys/pwp/engine", :ω), γ = output("kin", :γ_gnd)),
    v -> (v.θ̇, v.ω̇, v.γ))

######################## The walkthrough #######################################

function main()

    engine = PistonEngine(220.0)
    ac = Cessna172X(C172XSystems(PWP(engine), Aero()))
    world = SimpleWorld(ac, ac)

    banner("0. Instances carry no paths (§8.6)")
    println("  world.lead === world.wing: ", world.lead === world.wing,
            "   — value-identical siblings; only declaration names address them")

    banner("1. A condition is a lazy, inert tree — no path arithmetic at composition")
    taxi = ready_for_taxi(ac)
    print_tree(taxi)

    banner("2. Resolution: flatten (the one place paths concatenate) → validate")
    print_entries(resolve(taxi))

    banner("3. Duplicate leaf in one layer: error with both provenance chains")
    bad = merge(taxi, at("sys/pwp/engine", fragment(x = (ω = 0.0,))))
    try
        resolve(bad)
    catch e
        println("  ", sprint(showerror, e))
    end

    banner("4. override(base, patch): ordered layering — collision is the intent")
    warm = override(taxi, fragment(slots = (throttle = 0.3, mixture = 0.8)))
    print_entries(resolve(warm))

    banner("5. A TrimProblem is an implicitly specified condition (§14.9)")
    p = trim_problem(ac)
    println("  guess = ", p.guess, "\n")
    println("  condition(guess):")
    print_tree(p.condition(p.guess))
    println("\n  reads (aircraft-relative):")
    print_reads(p.reads)

    banner("6. at(\"wing\", problem) relocates it whole — five lines, nothing else")
    p_wing = at("wing", p)
    println("  condition(guess), now mounted:")
    print_entries(resolve(p_wing.condition(p_wing.guess)))
    println("\n  reads, now mounted (the Scoped wrapper, entered at resolution):")
    print_reads(p_wing.reads)
    println("\n  note the slot: face \"throttle\", authored at the aircraft,")
    println("  resolves from the mount point → root slot \"wing.throttle\"")

    banner("7. Commit = override(baseline, at(mount, condition(d*))) (§14.9)")
    baseline = merge(
        at("lead", ready_for_taxi(world.lead)),
        at("wing", ready_for_taxi(world.wing)),
        fragment(slots = (wind_N = 0.0,)))          #environment: world-level face
    d_star = (throttle = 0.42, θ = 0.031, n_eng = 0.83)
    commit = override(baseline, at("wing", p.condition(d_star)))
    print_entries(resolve(commit))
    println("\n  lead untouched; wing's solved leaves patched with dual provenance;")
    println("  method nesting became value layering — the world wrapper dissolved")

end

main()

end #module
