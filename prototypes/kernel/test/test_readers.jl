# --- the read-selector family, the compiled reader and `capture` (§14.4, §14.1;
# increment 21) ------------------------------------------------------------------
# The five deferred reads, their resolution against a build in §13.1's
# collecting register, the gather twin of `apply!` over an executor, and the
# service that reads the committed world back as a condition. The fixtures live
# at top level for the README's local-scope reason.

# Every home a read can come from, and nothing that fires: a continuous `x`
# (the plant's `q`) with its derivative, a discrete `s` (the integrator's
# `acc`), a mode store no event transitions (`ModedSource`, §8.2), two root
# inputs and one root-exported output face.
readable() = Group((; plant = Plant(), ctl = DiscreteIntegrator(3.0), src = ModedSource());
                   inputs = ("u" => "plant/u", "e" => "ctl/e"),
                   outputs = ("plant/y" => "y",))

# The captured world, authored: every store off its declared default, and `e`
# at zero so the integrator's `g` is stationary — boundary zero's outgoing
# transition (§14.5) is a mover like any handler, and a bit-for-bit round trip
# is a claim about the *establishment*, not about a `g` that would run again.
readable_condition(q = SVector(0.3, -0.2), acc = 4.0) =
    combine(at("plant", fragment(x = (q = q,))),
            at("ctl", fragment(s = (acc = acc,))),
            at("src", fragment(m = (phase = :running,))),
            fragment(inputs = (u = 1.5, e = 0.0)))

# The read set the two activations share.
readable_reads() = reads(q = get_state("plant", :q), v = get_state("plant", :q, 2),
                         acc = get_state("ctl", :acc), q̇ = get_deriv("plant", :q),
                         a = get_deriv("plant", :q, 2), y = get_output("plant", :y),
                         u = get_input(:u), face = get_face(:y))

# The stores a capture has to reproduce, read straight out of an executor.
world(sim) = (copy(sim.exec.xbuf),
              [s === nothing ? nothing : s[] for s in sim.exec.sstores],
              [m === nothing ? nothing : m[] for m in sim.exec.mstores],
              [port(sim, "", f) for f in sim.build.flat.root_inputs],
              sim.exec.clock.t)

@testset "the five selectors read what they name, at either activation (§14.4)" begin
    for T in (Float64, D8)
        sim = Simulation(readable(), T; h = 1//10)
        init!(sim, readable_condition())
        evaluate!(sim.exec)                     # `ẋ` is integrator scratch: fill it first
        r = compile_reads(readable_reads(), sim.build, T)
        v = gather(r, sim.exec)

        @test keys(v) === (:q, :v, :acc, :q̇, :a, :y, :u, :face)
        @test v.q == SVector{2,T}(0.3, -0.2)     # the whole leaf, out of `xbuf`
        @test v.v === v.q[2]                     # `i` indexes the read value
        @test v.acc === 4.0                      # the discrete store, pinned Float64
        @test v.q̇[1] === v.q[2]                  # `f`'s own output, out of `ẋbuf`
        @test v.a === v.q̇[2]
        @test v.y === v.q[1]                     # the stage-1 port's cell
        @test v.u === T(1.5)                     # the root input cell
        @test v.face === v.y                     # the exported face is its producer's cell

        # The leaf types are the activation's, so a `Dual` world reads `Dual`s
        # and the frozen discrete store stays pinned (§9.4, D-166).
        @test v.q isa SVector{2,T} && v.q̇ isa SVector{2,T} && v.y isa T
        @test v.acc isa Float64
    end
end

@testset "the reader is the gather twin: allocation-free over an executor (§14.4, §7.5)" begin
    sim = Simulation(readable(); h = 1//10)
    init!(sim, readable_condition())
    evaluate!(sim.exec)
    r, ex = compile_reads(readable_reads(), sim.build), sim.exec
    gather(r, ex)
    @test @ballocated(gather($r, $ex)) == 0
    @test @inferred(gather(r, ex)) isa NamedTuple
    @test gather(compile_reads(reads(), sim.build), ex) === (;)   # the empty set reads nothing
end

@testset "resolution collects every violation into one refusal (§14.4, §13.1)" begin
    b = build(readable())
    e = failure(() -> compile_reads(reads(a = get_state("plnt", :q),
                                          b = get_output("plant", :thrust),
                                          c = get_deriv("ctl", :acc),
                                          d = get_face(:nope)), b))
    @test e isa BuildError && occursin("4 violations", e.msg)
    @test occursin("`plnt` is no component of this build", e.msg)         # the offender, plainly
    @test occursin("declares no output port `thrust`", e.msg) &&
          occursin("`y`, `power`", e.msg)                                 # the list in hand
    @test occursin("get_deriv(\"ctl\", :acc)", e.msg) &&
          occursin("a discrete `s` has no derivative", e.msg)
    @test occursin("`nope` is no root-exported output face", e.msg) && occursin("`y`", e.msg)
    for label in ("`a`", "`b`", "`c`", "`d`")                             # each read, by label
        @test occursin(label, e.msg)
    end

    # An assembly path, a root input read as a face, an index on a scalar leaf,
    # and a state field the component does not declare.
    e = failure(() -> compile_reads(reads(a = get_output("", :y), b = get_face(:u),
                                          c = get_output("plant", :y, 1),
                                          d = get_state("plant", :ω)), b))
    @test occursin("the root component is an assembly", e.msg) &&
          occursin("read with `get_face`", e.msg)
    @test occursin("`u` is a root *input* face", e.msg) && occursin("get_input", e.msg)
    @test occursin("a scalar has no index", e.msg)
    @test occursin("declares no state field `ω`", e.msg) && occursin("`q`", e.msg)

    # The read set is a type, not a NamedTuple: the bare spelling is refused
    # with a directive, not a `MethodError` (§14.2's rule, one register over).
    e = failure(() -> compile_reads((q = get_state("plant", :q),), b))
    @test e isa BuildError && occursin("is not a read set", e.msg) && occursin("reads(", e.msg)
    @test failure(() -> reads(q = 2.0)) isa BuildError                     # nor is 2.0 a selector
end

@testset "the source rule: a snapshot-bound reader may not name a store selector (§14.4)" begin
    sim = Simulation(readable(); h = 1//10)
    e = failure(() -> attach!(sim, Pad("t"), Readout(q = get_state("plant", :q))))
    @test e isa BuildError && occursin("ReadBindingUnresolved", e.msg)
    @test occursin("a *store* selector", e.msg) && occursin("get_state(\"plant\", :q)", e.msg)
    e = failure(() -> attach!(sim, Pad("t"), Readout(y = get_output("plant", :y, 1))))
    @test e isa BuildError && occursin("a binding read is a whole cell", e.msg)
    @test isempty(sim.plane.roster)              # every rejection left the roster untouched
end

@testset "a reader and a plan belong to one activation, by dispatch (§9.4, §14.4)" begin
    # The offsets, store types and cell addresses a compiled product bakes are
    # one activation's, so a `Float64` product against a `Dual` executor would
    # read and write another cell's slot in silence. The scalar rides in each
    # product's type, the pairing is dispatch, and the mismatch is refused with
    # both scalars named — before anything is touched.
    nominal = Simulation(readable(); h = 1//10)
    seeded = Simulation(readable(), D8; h = 1//10)
    init!(seeded, readable_condition())
    before = world(seeded)

    e = failure(() -> gather(compile_reads(readable_reads(), nominal.build), seeded.exec))
    @test e isa BuildError && startswith(e.msg, "ActivationMismatch:")
    @test occursin("compiled at Float64", e.msg) && occursin("Dual{Nothing, Float64, 8}", e.msg)

    c = readable_condition()
    e2 = failure(() -> apply!(seeded.exec, resolve_condition(c, nominal.build)))
    @test e2 isa BuildError && startswith(e2.msg, "ActivationMismatch:")
    e3 = failure(() -> apply!(seeded.exec, compile_plan(c, nominal.build), c))
    @test e3 isa BuildError && startswith(e3.msg, "ActivationMismatch:")

    @test world(seeded) == before                # every refusal left the executor alone
end

@testset "`capture` reads the committed world back as a total condition (§14.1, §14.10)" begin
    sim = Simulation(readable(); h = 1//10)
    twin = Simulation(readable(); h = 1//10)
    init!(sim, readable_condition())

    (c, t) = capture(sim)
    @test t === 0.0                              # the condition is time-free; `t` rides beside
    @test c isa ConditionNode
    init!(twin, c; t₀ = t)
    @test world(twin) == world(sim)              # x, every `s` and `m`, root inputs, clock

    # It is total by construction (§14.6): no baseline underneath, and the
    # authored values are what a re-application establishes — the defaults
    # would show as `phase = :idle` and `acc = 0.0`.
    @test resolve_condition(c, twin.build).faces == twin.build.flat.root_inputs
    @test modes(twin, "src") === (phase = :running,)
    @test state(twin, "ctl").acc === 4.0

    # And after a trajectory: the same pair, taken at the run's end, is the
    # warm restart's baseline — clock included, which is what `t₀ = t` is for.
    run!(sim; t_end = 0.5)
    @test lifecycle(sim) === :stopped && state(sim, "plant").q != SVector(0.3, -0.2)
    (c2, t2) = capture(sim)
    @test t2 === 0.5
    init!(twin, c2; t₀ = t2)
    @test world(twin) == world(sim)
    @test twin.exec.clock.t === 0.5 && twin.exec.clock.t₀ === 0.5
end

@testset "`capture` is legal in `initialized` and `stopped`, and nowhere else (§14)" begin
    sim = Simulation(readable(); h = 1//10)
    e = failure(() -> capture(sim))
    @test e isa BuildError && occursin("ServiceLifecycle", e.msg)
    @test occursin("never been initialized", e.msg)      # `built`: no committed stores yet
    init!(sim, readable_condition())
    @test capture(sim) isa Tuple                         # `initialized`
    run!(sim; t_end = 0.1)
    @test capture(sim) isa Tuple                         # `stopped`

    # `running` is the §11.3 freeze: the loop owns the stores between drains.
    # Both ends of the run are test-controlled, exactly as in test_lifecycle.
    live = Simulation(armed(); h = 1//100, t_end = 3.0e5, stop_on = ("stop",))
    init!(live, fragment(inputs = (in = 0.0,)))
    task = Threads.@spawn run!(live)
    while lifecycle(live) !== :running
        yield()
    end
    err = failure(() -> capture(live))
    stage!(live, "in" => 1.0)
    wait(task)
    @test err isa BuildError && occursin("ServiceLifecycle", err.msg)
    @test occursin("stopped-sim operation and the simulation is running", err.msg)
end
