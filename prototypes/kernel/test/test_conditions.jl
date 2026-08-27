# --- the condition algebra and `init!` (§14.1–§14.6; increment 18) ---------------
# The inert lazy tree, its resolution against a build, slot totality, and the
# service that puts all three together. The fixtures live at top level for the
# README's local-scope reason.

# Every store a condition can address, and two root slots in declaration order:
# a continuous `x` (the plant's `q`), a discrete `s` (the integrator's `acc`), a
# mode set (the trigger's), and the `u`, `e` slots whose ordering is what
# `UninitializedSlots` has to report.
tri() = Group((; plant = Plant(), ctl = DiscreteIntegrator(3.0), trig = Trigger(0.5)),
              ("children/plant/y" => "children/trig/sig",),
              ("u" => "children/plant/u", "e" => "children/ctl/e"), ())

# A workspace declarer, for the "never workspace" half of §14.1's rule.
scratchy() = Group((; sm = Smoother(0.5)), (),
                   ("a" => "children/sm/a", "b" => "children/sm/b"), ())

# An offset pair, one stage shape each, both at `Relative(2, 1)` so neither is
# due at boundary zero: `hold` samples the root slot through `h_su`, `off`
# publishes its own `s` through `h_s` and accumulates in `g`. What D-205 rules
# on is exactly what these two read at `t₀`.
offset_pair() = Group((; hold = ZOH(), off = DiscreteIntegrator(1.0)),
                      ("children/hold/out" => "children/off/e",),
                      ("in" => "children/hold/in",),
                      ("children/hold/out" => "held", "children/off/u" => "acc"),
                      (; var"children/hold" = Relative(2, 1),
                         var"children/off" = Relative(2, 1)))

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
    @test isbits(fragment(x = (q = 1.0,), slots = (u = 2.0,)))
    @test isbits(combine(fragment(), fragment()))
    # It fails at resolution, where the build is finally in hand.
    @test failure(() -> resolve(n, build(tri()))) isa BuildError
end

@testset "a `combine` collision names both provenance chains and the layering combinator (§14.2)" begin
    b = build(tri())
    e = failure(() -> resolve(combine(at("children/plant", condition(Plant(); y = 1.0)),
                                      at("children/plant",
                                         fragment(x = (q = SVector(2.0, 0.0),)))), b))
    @test e isa BuildError
    @test occursin("DuplicateConditionLeaf", e.msg)
    @test occursin("`x.q` at `children/plant`", e.msg)
    @test occursin("combine[1] → at(\"children/plant\")", e.msg)
    @test occursin("combine[2] → at(\"children/plant\")", e.msg)
    @test occursin("collision-intolerant by design — use `override(base, patch)`", e.msg)
end

@testset "`override` layers: the patch wins, untouched leaves pass through (§14.6)" begin
    b = build(tri())
    base = combine(at("children/plant", fragment(x = (q = SVector(1.0, 2.0),))),
                   fragment(slots = (u = 1.0, e = 2.0)))
    slot(p, f) = only(v for (face, _, v) in p.slots if face === f)

    p = resolve(override(base, fragment(slots = (u = 9.0,))), b)
    @test slot(p, :u) === 9.0                      # the patch wins on the shared leaf
    @test slot(p, :e) === 2.0                      # untouched leaves pass through
    @test only(v for (_, v) in p.xs) === SVector(1.0, 2.0)
    @test length(p.slots) == 2                     # the overridden leaf is replaced, not doubled

    # Layering is variadic, and the last layer wins.
    p3 = resolve(override(base, fragment(slots = (u = 9.0,)), fragment(slots = (u = 7.0,))), b)
    @test slot(p3, :u) === 7.0

    # Provenance keeps both sources: the patch's own chain, and the base's
    # beside it — surfaced here through a violation on the overridden leaf.
    e = failure(() -> resolve(override(fragment(slots = (u = 1.0,)),
                                       fragment(slots = (u = "high",))), b))
    @test occursin("override[patch 1] → fragment(slots).u", e.msg)
    @test occursin("overrode override[base] → fragment(slots).u", e.msg)

    # A collision *within* one layer is still an error (§14.6).
    e = failure(() -> resolve(override(combine(fragment(slots = (u = 1.0,)),
                                               fragment(slots = (u = 2.0,))),
                                       fragment(slots = (u = 3.0,))), b))
    @test occursin("DuplicateConditionLeaf", e.msg)

    # §14.6's central use case: a full-coverage baseline authored at the root,
    # under a patch a component's own `condition` method ships against its own
    # face. Two spellings of one root slot are one leaf, so they layer.
    p4 = resolve(override(fragment(slots = (u = 1.0, e = 2.0)),
                          at("children/plant", fragment(slots = (u = 9.0,)))), b)
    @test slot(p4, :u) === 9.0 && slot(p4, :e) === 2.0
    @test length(p4.slots) == 2

    # The same two spellings under `combine` still collide — and with layering
    # no longer reaching this branch, its directive is advice that works.
    e = failure(() -> resolve(combine(fragment(slots = (u = 1.0,)),
                                      at("children/plant",
                                         fragment(slots = (u = 9.0,)))), b))
    @test occursin("DuplicateConditionLeaf", e.msg) && occursin("root slot `u`", e.msg)
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
                  at("children/plant", fragment(x = (nope = 1.0,))),     # undeclared field
                  at("children/ctl", fragment(s = (acc = "high",))),     # unconvertible
                  at("children/trig", fragment(slots = (sig = 1.0,))),   # never a root slot
                  fragment(slots = (u = 1.0,)),
                  at("children/plant", fragment(slots = (u = 2.0,))))    # one root slot, twice
    e = failure(() -> resolve(bad, b))
    @test e isa BuildError
    @test occursin("5 violations", e.msg)                  # the full list, one throw
    @test occursin("no component of this build", e.msg)
    @test occursin("`init_x` at `children/plant` declares `q`", e.msg)
    @test occursin("does not convert", e.msg)
    @test occursin("wired internally", e.msg)
    @test occursin("root slot `u` is written twice", e.msg)

    # An assembly path owns no state, and saying so beats "no such path".
    e = failure(() -> resolve(at("loop", fragment(x = (q = 1.0,))), build(Vehicle())))
    @test occursin("is an assembly", e.msg)

    # A tier's own state letter: `s` on a continuous component is not a typo
    # the resolver should guess at.
    e = failure(() -> resolve(at("children/plant", fragment(s = (q = 1.0,))), b))
    @test occursin("continuous component and declares no `init_s`", e.msg)
end

@testset "a condition names state, modes and root slots — never outputs, never workspace (§14.1)" begin
    e = failure(() -> resolve(at("children/plant", fragment(x = (y = 1.0,))), build(tri())))
    @test occursin("`y` is an output port", e.msg)
    e = failure(() -> resolve(at("children/sm", fragment(s = (tmp = 1.0,))), build(scratchy())))
    @test occursin("`tmp` is a workspace entry", e.msg)
end

@testset "slot faces resolve through the export chain to a root slot (§14.2)" begin
    b = build(tri())
    # The authoring level names a face of its own contract; resolution walks
    # the chain and lands on the root slot the obligation ends at.
    p = resolve(combine(at("children/plant", fragment(slots = (u = 1.0,))),
                        fragment(slots = (e = 2.0,))), b)
    @test Set(p.faces) == Set([:u, :e])
    # An internally wired input has no slot: writing it would be meaningless,
    # the first sweep overwriting it. Unexported stays unpokeable.
    @test occursin("reaches no root slot",
                   failure(() -> resolve(at("children/trig",
                                            fragment(slots = (sig = 1.0,))), b)).msg)
    @test occursin("declares no input face",
                   failure(() -> resolve(at("children/plant",
                                            fragment(slots = (nope = 1.0,))), b)).msg)
    @test occursin("no root input face",
                   failure(() -> resolve(fragment(slots = (nope = 1.0,)), b)).msg)
end

@testset "slot totality is checked pre-write, and a rejection changes nothing (§14.6, D-068)" begin
    sim = Simulation(tri(); h = 1//10)
    init!(sim, fragment(slots = (u = 1.0, e = 2.0)))
    snap, q = latest(sim), state(sim, "children/plant").q
    acc, lc = state(sim, "children/ctl").acc, lifecycle(sim)

    e = failure(() -> init!(sim, combine(at("children/plant",
                                            fragment(x = (q = SVector(5.0, 5.0),))),
                                         fragment(slots = (u = 3.0,)))))
    @test e isa BuildError && occursin("UninitializedSlots", e.msg)
    @test occursin("`e`", e.msg) && !occursin("`u`", e.msg)   # only the uncovered face
    # All-or-nothing: the plan's x write and its slot write both stayed home.
    @test state(sim, "children/plant").q === q
    @test state(sim, "children/ctl").acc === acc
    @test port(sim, "", :u) === 1.0
    @test lifecycle(sim) === lc && latest(sim) === snap

    # Every uncovered face, in declaration order (§14.6).
    fresh = Simulation(tri(); h = 1//10)
    @test occursin("`u`, `e`", failure(() -> init!(fresh)).msg)
    @test lifecycle(fresh) === :built
end

@testset "`init!(sim)` is legal exactly where the build has no root slots (§14.6)" begin
    free = Simulation(single(Sawtooth(1.0)); h = 1//10)
    init!(free)                                    # nothing to cover: total by construction
    @test lifecycle(free) === :initialized
    @test isempty(free.flat.slots)
end

@testset "a sparse overlay lands on the declared defaults (§14.1, §14.3)" begin
    sim = Simulation(tri(); h = 1//10)
    init!(sim, combine(at("children/ctl", fragment(s = (acc = 4.0,))),
                       at("children/trig", fragment(m = (count = 7,))),
                       fragment(slots = (u = 1.0, e = 0.0))))
    # Authored fields land; the fields no fragment named hold their declared
    # defaults, the overlay being `merge(defaults, overlay)` and nothing more.
    @test port(sim, "children/ctl", :u) === 4.0              # `s`, published before `g`
    @test state(sim, "children/plant").q === SVector(0.0, 0.0)
    @test modes(sim, "children/trig") === (state = :armed, count = 7)
    @test port(sim, "", :u) === 1.0
end

@testset "fresh run: `init!` restarts from the declared defaults, not the trajectory (D-063)" begin
    sim = Simulation(tri(); h = 1//10)
    init!(sim, fragment(slots = (u = 4.0, e = 1.0)))
    run!(sim; t_end = 2.0)
    @test lifecycle(sim) === :stopped
    @test state(sim, "children/plant").q !== SVector(0.0, 0.0)
    @test state(sim, "children/ctl").acc != 0.0
    @test modes(sim, "children/trig").state === :fired        # the run moved every store

    # The overlay base is always the declared defaults, never the last
    # trajectory's final state: bitwise the `init_*` values again.
    init!(sim, fragment(slots = (u = 0.0, e = 0.0)))
    @test state(sim, "children/plant").q === SVector(0.0, 0.0)
    @test port(sim, "children/ctl", :u) === 0.0
    @test modes(sim, "children/trig") === (state = :armed, count = 0)
end

@testset "boundary zero publishes every output stage, due or not (§14.5, D-205)" begin
    # Neither component is due at boundary zero — first tick at Φ·Δt_base — and
    # both publish there all the same, evaluated from the authored world: the
    # ZOH from the condition's slot value, the integrator from its authored `s`.
    # The probe's synthesized 0.0 reaches no published cell.
    total = combine(at("children/off", fragment(s = (acc = 7.0,))),
                    fragment(slots = (in = 3.0,)))
    sim = Simulation(offset_pair(); h = 1//100)
    init!(sim, total)
    @test port(sim, "", :held) === 3.0             # u(t₀), not the probe's 0.0
    @test port(sim, "", :acc) === 7.0              # the authored `s`, published
    @test port(latest(sim), "", :held) === 3.0     # and the t₀ snapshot carries both
    @test port(latest(sim), "", :acc) === 7.0

    # `g` stays gated by Φ (§10.5): the evaluation above is establishment, not a
    # scheduled sample, so the authored `s` survives boundary zero untouched and
    # the first sample `off` *consumes* is its own tick's, at Φ·Δt_base.
    @test state(sim, "children/off").acc === 7.0
    step!(sim)                                     # base tick 1 = Φ·Δt_base: the first tick
    @test state(sim, "children/off").acc ≈ 7.0 + 0.02 * 3.0    # one period, one sample

    # Cold and warm are the same world: the authored condition determines the
    # `t₀` table outright, so a restart after a run reproduces it exactly.
    run!(sim; t_end = 0.5)
    @test state(sim, "children/off").acc > 7.0     # the run moved the store
    init!(sim, total)
    @test port(sim, "", :held) === 3.0 && port(sim, "", :acc) === 7.0
    @test state(sim, "children/off").acc === 7.0
end

@testset "an authored state fires its guard at t₀ (§14.5, §10.6)" begin
    # A condition landing a predicate in holding territory — authored, not
    # staged: boundary zero establishes every prior as not-holding, so the
    # event fires visibly at t₀ rather than one step later.
    m = tri()
    sim = Simulation(m; h = 1//10)
    init!(sim, combine(at("children/plant", condition(m.children.plant; y = 1.0)),
                       fragment(slots = (u = 0.0, e = 0.0))))
    @test modes(sim, "children/trig") === (state = :fired, count = 1)
    @test port(sim, "children/trig", :on) === true
    @test latest(sim).t == 0.0 && sim.clock.step == 0
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
    init!(sim, override(condition(veh; ref = 1.0), fragment(slots = (ref = 5.0,))))
    @test port(sim, "", :ref) === 5.0
    @test state(sim, "loop/plant").q === SVector(0.0, 0.0)   # the baseline's own defaults
end
