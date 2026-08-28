# --- the condition algebra and `init!` (§14.1–§14.6; increment 18) ---------------
# The inert lazy tree, its resolution against a build, root-input totality, and
# the service that puts all three together. The fixtures live at top level for the
# README's local-scope reason.

# Every store a condition can address, and two root inputs in declaration order:
# a continuous `x` (the plant's `q`), a discrete `s` (the integrator's `acc`), a
# mode set (the trigger's), and the `u`, `e` root inputs whose ordering is what
# `UninitializedInputs` has to report.
tri() = Group((; plant = Plant(), ctl = DiscreteIntegrator(3.0), trig = Trigger(0.5));
              wires = ("plant/y" => "trig/sig",),
              inputs = ("u" => "plant/u", "e" => "ctl/e"))

# A named assembly one level down, its input face fed from the root's: the
# per-level `at` addressing D-207's total face graph makes resolvable.
nested() = Group((; loop = SampledLoop()); inputs = ("in" => "loop/ref",),
                 outputs = ("loop/y" => "y",))

# A workspace declarer, for the "never workspace" half of §14.1's rule.
scratchy() = Group((; sm = Smoother(0.5));
                   inputs = ("a" => "sm/a", "b" => "sm/b"))

# An offset pair, one stage shape each, both at `Relative(2, 1)` so neither is
# due at boundary zero: `hold` samples the root input through `h_su`, `off`
# publishes its own `s` through `h_s` and accumulates in `g`. What D-205 rules
# on is exactly what these two read at `t₀`.
offset_pair() = Group((; hold = ZOH(), off = DiscreteIntegrator(1.0));
                      wires = ("hold/out" => "off/e",),
                      inputs = ("in" => "hold/in",),
                      outputs = ("hold/out" => "held", "off/u" => "acc"),
                      rates = (; hold = Relative(2, 1),
                                 off = Relative(2, 1)))

@testset "composition is inert and lazy: no path arithmetic, no validation (§14.2)" begin
    # A deep tree over a path that resolves against nothing. Constructing it
    # performs no lookup, concatenates no string and checks no field.
    n = combine(at("nowhere", fragment(x = (q = 1.0,))),
                at("children", combine(at("plant", fragment(x = (q = 2.0,))),
                                       fragment(s = (acc = 3.0,)))))
    @test n isa Combined
    @test n.nodes[1] isa Scoped && n.nodes[1].prefix == "nowhere"
    @test n.nodes[1].node isa Fragment
    @test n.nodes[2].node isa Combined                  # `at` stores, never applies
    @test n.nodes[2].node.nodes[1].prefix == "plant"    # unconcatenated with "children"
    # Every node is isbits but for the prefix strings (§14.2): rebuilding the
    # tree per iteration is stack-only construction.
    @test isbits(fragment(x = (q = 1.0,), inputs = (u = 2.0,)))
    @test isbits(combine(fragment(), fragment()))
    # It fails at resolution, where the build is finally in hand.
    @test failure(() -> resolve_condition(n, build(tri()))) isa BuildError
end

@testset "a `combine` collision names both provenance chains and the layering combinator (§14.2)" begin
    b = build(tri())
    e = failure(() -> resolve_condition(combine(at("plant", condition(Plant(); y = 1.0)),
                                      at("plant",
                                         fragment(x = (q = SVector(2.0, 0.0),)))), b))
    @test e isa BuildError
    @test occursin("DuplicateConditionLeaf", e.msg)
    @test occursin("`x.q` at `plant`", e.msg)
    @test occursin("combine[1] → at(\"plant\")", e.msg)
    @test occursin("combine[2] → at(\"plant\")", e.msg)
    @test occursin("collision-intolerant by design — use `override(base, patch)`", e.msg)
end

@testset "`override` layers: the patch wins, untouched leaves pass through (§14.6)" begin
    b = build(tri())
    base = combine(at("plant", fragment(x = (q = SVector(1.0, 2.0),))),
                   fragment(inputs = (u = 1.0, e = 2.0)))
    input(p, f) = only(v for (face, _, v) in p.inputs if face === f)

    p = resolve_condition(override(base, fragment(inputs = (u = 9.0,))), b)
    @test input(p, :u) === 9.0                     # the patch wins on the shared leaf
    @test input(p, :e) === 2.0                     # untouched leaves pass through
    @test only(v for (_, v) in p.xs) === SVector(1.0, 2.0)
    @test length(p.inputs) == 2                     # the overridden leaf is replaced, not doubled

    # Layering is variadic, and the last layer wins.
    p3 = resolve_condition(override(base, fragment(inputs = (u = 9.0,)),
                                    fragment(inputs = (u = 7.0,))), b)
    @test input(p3, :u) === 7.0

    # Provenance keeps both sources: the patch's own chain, and the base's
    # beside it — surfaced here through a violation on the overridden leaf.
    e = failure(() -> resolve_condition(override(fragment(inputs = (u = 1.0,)),
                                       fragment(inputs = (u = "high",))), b))
    @test occursin("override[patch 1] → fragment(inputs).u", e.msg)
    @test occursin("overrode override[base] → fragment(inputs).u", e.msg)

    # A collision *within* one layer is still an error (§14.6).
    e = failure(() -> resolve_condition(override(combine(fragment(inputs = (u = 1.0,)),
                                               fragment(inputs = (u = 2.0,))),
                                       fragment(inputs = (u = 3.0,))), b))
    @test occursin("DuplicateConditionLeaf", e.msg)

    # §14.6's central use case: a full-coverage baseline authored at the root,
    # under a patch a component's own `condition` method ships against its own
    # face. Two spellings of one root input are one leaf, so they layer.
    p4 = resolve_condition(override(fragment(inputs = (u = 1.0, e = 2.0)),
                          at("plant", fragment(inputs = (u = 9.0,)))), b)
    @test input(p4, :u) === 9.0 && input(p4, :e) === 2.0
    @test length(p4.inputs) == 2

    # The same two spellings under `combine` still collide — and with layering
    # no longer reaching this branch, its directive is advice that works.
    e = failure(() -> resolve_condition(combine(fragment(inputs = (u = 1.0,)),
                                      at("plant",
                                         fragment(inputs = (u = 9.0,)))), b))
    @test occursin("DuplicateConditionLeaf", e.msg) && occursin("root input `u`", e.msg)
    @test occursin("use `override(base, patch)`", e.msg)
end

@testset "blending a node with a bare NamedTuple is a directive error method (§14.2)" begin
    # Raised at composition time, before any resolution pass or provenance
    # chain exists — which is why it carries its own kind. No build in hand.
    for f in (() -> combine(fragment(), (q = 1.0,)),        # node × NamedTuple
              () -> combine((q = 1.0,), fragment()),        # and the other order
              () -> combine(fragment(), fragment(), (q = 1.0,)),   # at any arity
              () -> at("plant", (q = 1.0,)),
              () -> override(fragment(), (q = 1.0,)),
              () -> fragment(x = 3.0))                      # and a non-NamedTuple payload
        e = failure(f)
        @test e isa BuildError && occursin("ConditionNodeMisuse", e.msg)
    end
    @test occursin("wrap the NamedTuple in `fragment(…)`",
                   failure(() -> combine(fragment(), (q = 1.0,))).msg)
    # And the service entry point itself: a bare NamedTuple where a condition
    # belongs gets the directive, never a `MethodError`.
    e = failure(() -> init!(Simulation(tri(); h = 1//10), (u = 1.0, e = 2.0)))
    @test e isa BuildError && occursin("ConditionNodeMisuse", e.msg)
end

@testset "resolution collects every violation into one throw (§14.3, §13.1)" begin
    b = build(tri())
    bad = combine(at("nope", fragment(x = (q = 1.0,))),                  # unknown path
                  at("plant", fragment(x = (nope = 1.0,))),     # undeclared field
                  at("ctl", fragment(s = (acc = "high",))),     # unconvertible
                  at("trig", fragment(inputs = (sig = 1.0,))),  # never a root input
                  fragment(inputs = (u = 1.0,)),
                  at("plant", fragment(inputs = (u = 2.0,))))   # one root input, twice
    e = failure(() -> resolve_condition(bad, b))
    @test e isa BuildError
    @test occursin("5 violations", e.msg)                  # the full list, one throw
    @test occursin("no component of this build", e.msg)
    @test occursin("`init_x` at `plant` declares `q`", e.msg)
    @test occursin("does not convert", e.msg)
    @test occursin("wired internally", e.msg)
    @test occursin("root input `u` is written twice", e.msg)

    # An assembly path owns no state, and saying so beats "no such path".
    e = failure(() -> resolve_condition(at("loop", fragment(x = (q = 1.0,))), build(Vehicle())))
    @test occursin("is an assembly", e.msg)

    # A tier's own state letter: `s` on a continuous component is not a typo
    # the resolver should guess at.
    e = failure(() -> resolve_condition(at("plant", fragment(s = (q = 1.0,))), b))
    @test occursin("continuous component and declares no `init_s`", e.msg)
end

@testset "a condition names state, modes and root inputs — never outputs, never workspace (§14.1)" begin
    e = failure(() -> resolve_condition(at("plant", fragment(x = (y = 1.0,))), build(tri())))
    @test occursin("`y` is an output port", e.msg)
    e = failure(() -> resolve_condition(at("sm", fragment(s = (tmp = 1.0,))), build(scratchy())))
    @test occursin("`tmp` is a workspace entry", e.msg)
end

@testset "input faces resolve through the export chain to a root input (§14.2)" begin
    b = build(tri())
    # The authoring level names a face of its own contract; resolution walks
    # the chain and lands on the root input the obligation ends at.
    p = resolve_condition(combine(at("plant", fragment(inputs = (u = 1.0,))),
                        fragment(inputs = (e = 2.0,))), b)
    @test Set(p.faces) == Set([:u, :e])
    # An internally wired input reaches no root input: writing it would be
    # meaningless, the first sweep overwriting it. Unexported stays unpokeable.
    @test occursin("reaches no root input",
                   failure(() -> resolve_condition(at("trig",
                                            fragment(inputs = (sig = 1.0,))), b)).msg)
    @test occursin("declares no input face",
                   failure(() -> resolve_condition(at("plant",
                                            fragment(inputs = (nope = 1.0,))), b)).msg)
    # The discrimination an `inputs` payload makes on its own: a prefix naming no
    # level of this build reads "no component", not the face-typo message the
    # real component earns above — an empty face list alone cannot tell them
    # apart, assemblies leaving no row in the flat list (§14.3).
    @test occursin("no component of this build",
                   failure(() -> resolve_condition(at("nope",
                                                       fragment(inputs = (dead = 1.0,))),
                                                    b)).msg)
    @test occursin("no root input face",
                   failure(() -> resolve_condition(fragment(inputs = (nope = 1.0,)), b)).msg)
end

@testset "an `at` prefix stopping at an assembly resolves its faces (§14.2, D-207)" begin
    b = build(nested())
    # The face graph is total, so a prefix may stop at *any* child's faces: the
    # sub-assembly's own `ref` is looked up at that level and followed to the
    # root input the chain ends at. The plan is the one the root spelling
    # produces, entry for entry — two spellings of one root input.
    p = resolve_condition(at("loop", fragment(inputs = (ref = 1.0,))), b)
    q = resolve_condition(fragment(inputs = (in = 1.0,)), b)
    @test p.inputs == q.inputs && p.faces == q.faces == [:in]

    # A component-fed face is still unpokeable, at an assembly prefix as at a
    # primitive's: `Vehicle` feeds the loop's `ref` from its own `trim`.
    @test occursin("reaches no root input",
                   failure(() -> resolve_condition(at("loop", fragment(inputs = (ref = 1.0,))),
                                         build(Vehicle()))).msg)

    # A face the level does not declare, with the level's face list in hand.
    e = failure(() -> resolve_condition(at("loop", fragment(inputs = (nope = 1.0,))), b))
    @test occursin("declares no input face `nope`", e.msg) && occursin("`ref`", e.msg)

    # State at an assembly prefix stays refused: assemblies own no state, and
    # only the `inputs` payload gained a level to resolve at.
    e = failure(() -> resolve_condition(at("loop", fragment(x = (q = 1.0,))), b))
    @test occursin("is an assembly", e.msg) && occursin("components and root inputs", e.msg)
end

@testset "root-input totality is checked pre-write, and a rejection changes nothing (§14.6, D-068)" begin
    sim = Simulation(tri(); h = 1//10)
    init!(sim, fragment(inputs = (u = 1.0, e = 2.0)))
    snap, q = latest(sim), state(sim, "plant").q
    acc, lc = state(sim, "ctl").acc, lifecycle(sim)

    e = failure(() -> init!(sim, combine(at("plant",
                                            fragment(x = (q = SVector(5.0, 5.0),))),
                                         fragment(inputs = (u = 3.0,)))))
    @test e isa BuildError && occursin("UninitializedInputs", e.msg)
    @test occursin("`e`", e.msg) && !occursin("`u`", e.msg)   # only the uncovered face
    # All-or-nothing: the plan's x write and its root-input write both stayed home.
    @test state(sim, "plant").q === q
    @test state(sim, "ctl").acc === acc
    @test port(sim, "", :u) === 1.0
    @test lifecycle(sim) === lc && latest(sim) === snap

    # Every uncovered face, in declaration order (§14.6).
    fresh = Simulation(tri(); h = 1//10)
    @test occursin("`u`, `e`", failure(() -> init!(fresh)).msg)
    @test lifecycle(fresh) === :built
end

@testset "`init!(sim)` is legal exactly where the build has no root inputs (§14.6)" begin
    free = Simulation(single(Sawtooth(1.0)); h = 1//10)
    init!(free)                                    # nothing to cover: total by construction
    @test lifecycle(free) === :initialized
    @test isempty(free.build.flat.root_inputs)
end

@testset "a sparse overlay lands on the declared defaults (§14.1, §14.3)" begin
    sim = Simulation(tri(); h = 1//10)
    init!(sim, combine(at("ctl", fragment(s = (acc = 4.0,))),
                       at("trig", fragment(m = (count = 7,))),
                       fragment(inputs = (u = 1.0, e = 0.0))))
    # Authored fields land; the fields no fragment named hold their declared
    # defaults, the overlay being `merge(defaults, overlay)` and nothing more.
    @test port(sim, "ctl", :u) === 4.0              # `s`, published before `g`
    @test state(sim, "plant").q === SVector(0.0, 0.0)
    @test modes(sim, "trig") === (state = :armed, count = 7)
    @test port(sim, "", :u) === 1.0
end

@testset "fresh run: `init!` restarts from the declared defaults, not the trajectory (D-063)" begin
    sim = Simulation(tri(); h = 1//10)
    init!(sim, fragment(inputs = (u = 4.0, e = 1.0)))
    run!(sim; t_end = 2.0)
    @test lifecycle(sim) === :stopped
    @test state(sim, "plant").q !== SVector(0.0, 0.0)
    @test state(sim, "ctl").acc != 0.0
    @test modes(sim, "trig").state === :fired        # the run moved every store

    # The overlay base is always the declared defaults, never the last
    # trajectory's final state: bitwise the `init_*` values again.
    init!(sim, fragment(inputs = (u = 0.0, e = 0.0)))
    @test state(sim, "plant").q === SVector(0.0, 0.0)
    @test port(sim, "ctl", :u) === 0.0
    @test modes(sim, "trig") === (state = :armed, count = 0)
end

@testset "boundary zero publishes every output stage, due or not (§14.5, D-205)" begin
    # Neither component is due at boundary zero — first tick at Φ·Δt_base — and
    # both publish there all the same, evaluated from the authored world: the
    # ZOH from the condition's root-input value, the integrator from its authored `s`.
    # The probe's synthesized 0.0 reaches no published cell.
    total = combine(at("off", fragment(s = (acc = 7.0,))),
                    fragment(inputs = (in = 3.0,)))
    sim = Simulation(offset_pair(); h = 1//100)
    init!(sim, total)
    @test port(sim, "", :held) === 3.0             # u(t₀), not the probe's 0.0
    @test port(sim, "", :acc) === 7.0              # the authored `s`, published
    @test port(latest(sim), "", :held) === 3.0     # and the t₀ snapshot carries both
    @test port(latest(sim), "", :acc) === 7.0

    # `g` stays gated by Φ (§10.5): the evaluation above is establishment, not a
    # scheduled sample, so the authored `s` survives boundary zero untouched and
    # the first sample `off` *consumes* is its own tick's, at Φ·Δt_base.
    @test state(sim, "off").acc === 7.0
    step!(sim)                                     # base tick 1 = Φ·Δt_base: the first tick
    @test state(sim, "off").acc ≈ 7.0 + 0.02 * 3.0    # one period, one sample

    # Cold and warm are the same world: the authored condition determines the
    # `t₀` table outright, so a restart after a run reproduces it exactly.
    run!(sim; t_end = 0.5)
    @test state(sim, "off").acc > 7.0     # the run moved the store
    init!(sim, total)
    @test port(sim, "", :held) === 3.0 && port(sim, "", :acc) === 7.0
    @test state(sim, "off").acc === 7.0
end

@testset "an authored state fires its guard at t₀ (§14.5, §10.6)" begin
    # A condition landing a predicate in holding territory — authored, not
    # staged: boundary zero establishes every prior as not-holding, so the
    # event fires visibly at t₀ rather than one step later.
    m = tri()
    sim = Simulation(m; h = 1//10)
    init!(sim, combine(at("plant", condition(m.children.plant; y = 1.0)),
                       fragment(inputs = (u = 0.0, e = 0.0))))
    @test modes(sim, "trig") === (state = :fired, count = 1)
    @test port(sim, "trig", :on) === true
    @test latest(sim).t == 0.0 && sim.exec.clock.step == 0
end

@testset "the fragment-function idiom composes by pull across two levels (§14.2)" begin
    veh = Vehicle(; k = 2.0)
    sim = Simulation(veh; h = 1//50)
    init!(sim, condition(veh; ref = 1.0, y = 0.3, v = -0.2, cmd = 2.0))
    # Deep paths are compiled derivatives of the nesting `at` recorded: the
    # owner wrote "loop", the loop wrote "plant", and nothing wrote "loop/plant".
    @test state(sim, "loop/plant").q === SVector(0.3, -0.2)
    @test port(sim, "loop", :cmd) === 2.0
    @test port(sim, "", :ref) === 1.0

    # And the baseline-plus-tweak spelling the same function supports (§14.6).
    init!(sim, override(condition(veh; ref = 1.0), fragment(inputs = (ref = 5.0,))))
    @test port(sim, "", :ref) === 5.0
    @test state(sim, "loop/plant").q === SVector(0.0, 0.0)   # the baseline's own defaults
end

# --- the specialized application register (§14.3, §14.4, D-066) -----------------
# The other register over the same checks: a plan compiled from a tree's shape,
# applied to every later tree of that shape. The fixtures stay at top level for
# the README's local-scope reason.

# A discrete component with two `s` fields. The composite store merge needs one
# store whose fields two layers can author separately, and no coverage component
# declares one — a single-field store lets `override` hide the property behind
# its own last-wins.
struct Ledger <: AbstractComponent end

init_s(::Ledger) = (a = 0.0, b = 0.0)
output_types(::Ledger) = (total = Float64,)

h_s(::Ledger, (; s)) = (total = s.a + s.b,)
g(::Ledger, (; s)) = (a = s.a, b = s.b)

# One shape over `tri()`'s four homes, its values the parameters: the `at`
# literals sit at one source location, which is what makes the `===` prefix
# sweep a pointer compare rather than a string compare (§14.4).
tri_tree(q, acc, state, u, e) =
    combine(at("plant", fragment(x = (q = q,))),
            at("ctl", fragment(s = (acc = acc,))),
            at("trig", fragment(m = (state = state,))),
            fragment(inputs = (u = u, e = e)))

# The composite a service compiles from: a baseline authoring one field of the
# store and the per-iteration layer authoring the other (§14.3's fork).
ledger_tree(a, b) = override(at("led", fragment(s = (a = a,))),
                             at("led", fragment(s = (b = b,))))

# Everything an `apply!` can land in, for the comparisons below.
landed(sim) = (copy(sim.exec.xbuf),
               [s === nothing ? nothing : s[] for s in sim.exec.sstores],
               [m === nothing ? nothing : m[] for m in sim.exec.mstores],
               [port(sim, "", f) for f in sim.build.flat.root_inputs])

@testset "a shape-compiled plan lands what the dynamic walk lands (§14.4, D-066)" begin
    specialized = Simulation(tri(); h = 1//10)
    dynamic = Simulation(tri(); h = 1//10)
    # Compiled from one tree, applied to a second of the same shape with other
    # values everywhere: the plan holds lenses, not the values it was shown.
    plan = compile_plan(tri_tree(SVector(1.0, 2.0), 3.0, :fired, 4.0, 5.0), specialized.build)
    later = tri_tree(SVector(9.0, 8.0), 7.0, :armed, 6.0, 5.5)

    apply!(specialized.exec, plan, later)
    apply!(dynamic.exec, resolve_condition(later, dynamic.build))
    @test landed(specialized) == landed(dynamic)

    # And what it landed is the second tree's values, in all four homes.
    @test specialized.exec.xbuf == [9.0, 8.0]
    @test state(specialized, "ctl") === (acc = 7.0,)
    @test modes(specialized, "trig") === (state = :armed, count = 0)
    @test port(specialized, "", :u) === 6.0 && port(specialized, "", :e) === 5.5

    # The tree type is the plan's own type parameter, and the positions the
    # flattening recorded are the step tuples down to the authored values.
    @test plan isa SpecializedPlan{typeof(later)}
    @test only(plan.xs).authored isa Authored{(:nodes, 1, :node, :x, :q),SVector{2,Float64}}
end

@testset "the specialized register writes without allocating (§14.4, §7.5)" begin
    sim = Simulation(tri(); h = 1//10)
    plan = compile_plan(tri_tree(SVector(1.0, 2.0), 3.0, :fired, 4.0, 5.0), sim.build)
    tree = tri_tree(SVector(9.0, 8.0), 7.0, :armed, 6.0, 5.5)
    # The whole path: the prefix sweep, the flat-buffer write, both stores as
    # whole values, and the two root-input scatters. The tree is handed in
    # already built — an `at` node holds a `String` and is therefore not isbits,
    # so its construction is the caller's cost, not the register's.
    @test (@ballocated apply!($(sim.exec), $plan, $tree)) == 0
end

@testset "the store merge is the composite's, not one layer's (§14.3)" begin
    sim = Simulation(Group((; led = Ledger())); h = 1//10)
    # The plan is compiled from the whole composite tree a service builds, so
    # the two layers' fields meet in one `merge(defaults, overlay)`. A plan over
    # the patch alone would write `(a = 0.0, b = 4.0)` — the baseline's `a`
    # replaced by the declared default, which is the trap the composite avoids.
    plan = compile_plan(ledger_tree(1.0, 2.0), sim.build)
    apply!(sim.exec, plan, ledger_tree(3.0, 4.0))
    @test state(sim, "led") === (a = 3.0, b = 4.0)
    @test (@ballocated apply!($(sim.exec), $plan, $(ledger_tree(3.0, 4.0)))) == 0
end

@testset "shape drift is a structured error, and nothing is written (§14.4, §9.5)" begin
    sim = Simulation(tri(); h = 1//10)
    plan = compile_plan(tri_tree(SVector(1.0, 2.0), 3.0, :fired, 4.0, 5.0), sim.build)
    before = landed(sim)

    # A tree of another type never reaches a write: the shape is proven by
    # dispatch, and the fallback method names both types.
    e = failure(() -> apply!(sim.exec, plan, at("plant", fragment(x = (q = SVector(1.0, 2.0),)))))
    @test e isa BuildError
    @test occursin("ConditionShapeDrift", e.msg)
    @test occursin(string(typeof(tri_tree(SVector(1.0, 2.0), 3.0, :fired, 4.0, 5.0))), e.msg)
    @test occursin("Scoped{Fragment", e.msg)          # the observed tree's own type
    @test landed(sim) == before

    # The prefixes are runtime fields the type cannot carry, so the `===` sweep
    # closes the remainder — and it runs before any write. The sweep is a
    # pointer compare, which is honest because `at` stores the string it was
    # given (`String(::String)` returns its argument) and equal literals are one
    # object: the all-literal case §14.4 assumes, where the compares fold away.
    @test at("plant", fragment()).prefix === at("plant", fragment()).prefix
    drifted = combine(at("plant", fragment(x = (q = SVector(9.0, 8.0),))),
                      at("plant", fragment(s = (acc = 7.0,))),   # was "ctl"
                      at("trig", fragment(m = (state = :armed,))),
                      fragment(inputs = (u = 6.0, e = 5.5)))
    e2 = failure(() -> apply!(sim.exec, plan, drifted))
    @test e2 isa BuildError
    @test occursin("ConditionShapeDrift", e2.msg)
    @test occursin("(.nodes[2].prefix)", e2.msg)       # the position, in tree-step spelling
    @test occursin("\"ctl\"", e2.msg) && occursin("\"plant\"", e2.msg)
    @test landed(sim) == before
end

@testset "the converters are baked per leaf, at the activation (§14.3)" begin
    sim = Simulation(tri(), D8; h = 1//10)
    b = sim.build

    # A plain `Float64` leaf against a seeded activation: the zero-partial
    # embedding, which is semantically exact for a value held at the operating
    # point and in no other case.
    held = at("plant", fragment(x = (q = SVector(1.0, 2.0),)))
    apply!(sim.exec, compile_plan(held, b, D8), held)
    @test sim.exec.xbuf == D8[1.0, 2.0]
    @test all(iszero, ForwardDiff.partials(sim.exec.xbuf[1]))

    # A leaf already at the activation's scalar — decision-descended — takes the
    # type's own methods, partials flowing through untouched.
    d = ForwardDiff.Dual{Nothing}(2.5, ntuple(i -> i == 1 ? 1.0 : 0.0, 8)...)
    seeded = at("plant", fragment(x = (q = SVector(d, zero(d)),)))
    apply!(sim.exec, compile_plan(seeded, b, D8), seeded)
    @test ForwardDiff.value(sim.exec.xbuf[1]) === 2.5
    @test ForwardDiff.partials(sim.exec.xbuf[1])[1] === 1.0

    # And the one case no converter covers: a discrete `s` is frozen at a
    # non-nominal activation (§9.4), so a decision variable authored into it is
    # refused at resolution, with the clause that says why.
    e = failure(() -> compile_plan(at("ctl", fragment(s = (acc = d,))), b, D8))
    @test e isa BuildError
    @test occursin("ConditionResolution: `s.acc` at `ctl` takes Float64", e.msg)
    @test occursin("a decision variable descends into neither a frozen discrete `s` nor a " *
                   "pinned leaf", e.msg)
    # The nominal activation's own refusals are unchanged: no clause where the
    # value is simply the wrong kind of thing.
    e0 = failure(() -> compile_plan(at("ctl", fragment(s = (acc = :nope,))), b))
    @test occursin("which does not convert (§14.3)", e0.msg)
end
