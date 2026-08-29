# --- the input trace (§11.5, increment 23) --------------------------------------
# The header captured at `init!` and one sparse record per drained batch behind
# it: what replay (§12.7) consumes, and the primary record the log is derived
# from (D-038). The fixtures live at top level for the README's local-scope
# reason.

# Three root inputs, so a record's *position* against the writer's schema is a
# fact worth asserting rather than a coincidence at width one.
three_root_inputs() = Group((; s = Sum(sa = 1.0, sb = 1.0), g = Gain(2.0));
                            inputs = ("a" => "s/a", "b" => "s/b", "c" => "g/e"))

# A model whose boundary-zero sequence moves both discrete state homes (§14.5):
# the trigger's guard holds in the authored state and fires at `t₀`, and the
# integrator's due `g` runs there too. What the header must hold is what stood
# *before* either did.
boundary_movers() = Group((; t = Trigger(0.5), d = DiscreteIntegrator(1.0));
                          inputs = ("sig" => "t/sig", "e" => "d/e"))

# One record's face, resolved the way a consumer resolves it: through the
# writer's schema in the header (§11.5).
recorded_faces(trc, b::TraceBatch) =
    Symbol[last(trc.header.schemas[b.writer])[i] for (i, _) in b.entries]

@testset "one sparse record per drained batch, against the writer's schema (§11.5, D-176)" begin
    sim = Simulation(three_root_inputs(); h = 1//10)
    init!(sim, fragment(inputs = (a = 0.0, b = 0.0, c = 0.0)))
    @test isempty(trace(sim).batches) && trace(sim).frames == 0

    stage!(sim, "b" => 1.0)
    step!(sim)
    trc = trace(sim)
    b = only(trc.batches)
    @test b.frame == 1                     # the drain precedes the step increment
    @test b.writer == 1                    # the harness register, sole writer here
    @test first(trc.header.schemas[b.writer]) == "harness"
    (pos, v) = only(b.entries)             # sparse: the touched position alone
    @test pos == 2 && v === 1.0            # `b` is position 2 of {a, b, c}
    @test recorded_faces(trc, b) == [:b]
    @test trc.frames == 1

    # A frame nobody staged into records nothing and still advances the count:
    # `frames` is the recording's length, not its batch count.
    step!(sim)
    trc = trace(sim)
    @test length(trc.batches) == 1 && trc.frames == 2

    # A merged batch is recorded once, coalesced (§11.4): the drain is where the
    # conversion happens because the drained tuple is the coalesced truth.
    stage!(sim, "a" => 1.0)
    stage!(sim, "c" => 2.0, "a" => 3.0)
    step!(sim)
    b = last(trace(sim).batches)
    @test b.frame == 3 && recorded_faces(trace(sim), b) == [:a, :c]
    @test [v for (_, v) in b.entries] == [3.0, 2.0]

    # The quiet frame stays free: the register's frame stamp and its counter are
    # field writes, and nothing is recorded where nothing was drained (§11.1).
    @test @ballocated(drain!($sim)) == 0
end

@testset "the header is the pre-sequence state, resolved (§11.5, §14.5, D-038)" begin
    sim = Simulation(boundary_movers(); h = 1//10)
    init!(sim, fragment(inputs = (sig = 1.0, e = 2.0)))
    h = trace(sim).header
    # `flat.paths` is ["t", "d"]: the trigger's modes and the integrator's state.
    @test h.m == Any[(state = :armed, count = 0), nothing]
    @test h.s == Any[nothing, (acc = 0.0,)]
    # …while boundary zero has already fired the guard and run the due `g`.
    @test modes(sim, "t") == (state = :fired, count = 1)
    @test state(sim, "d") == (acc = 0.2,)

    # The root inputs ride along as resolved values (§11.5): neither face is
    # ever staged, so no batch would carry them and replay would have nothing.
    @test h.root_inputs == Pair{Symbol,Any}[:sig => 1.0, :e => 2.0]

    # The deployment block, captured at the same instant as the stores (§11.5):
    # the effective termination pair is the one `init!` knows, the constructor's.
    d = h.deployment
    @test d.t₀ === 0.0 && d.h === 0.1 && d.n == 1 && d.method === :RK4
    @test d.firing_budget == 4 && d.localization_budget == 8
    @test d.t_end === nothing && isempty(d.stop_on)
    @test h.layout.paths == ["t", "d"] && h.layout.root_faces == [:sig, :e]
end

@testset "the kill switch, and the clearing at `init!` (§11.5, D-029)" begin
    off = Simulation(three_root_inputs(); h = 1//10, trace = false)
    init!(off, fragment(inputs = (a = 0.0, b = 0.0, c = 0.0)))
    err = failure(() -> trace(off))
    d = only(err.diagnostics)
    @test err isa BuildError && d isa ArgumentInvalid
    @test d.call === :trace && d.reason === :disabled
    stage!(off, "a" => 1.0)
    step!(off)
    @test isempty(off.trace.batches) && off.trace.frames == 0   # nothing recorded

    # Before the first `init!` there is no header, so there is no recording, and
    # the refusal is the one an advance entry gives (§12.6).
    sim = Simulation(three_root_inputs(); h = 1//10)
    err = failure(() -> trace(sim))
    d = only(err.diagnostics)
    @test err isa BuildError && d isa MissingInit
    @test d.op === :trace && d.status === :built

    init!(sim, fragment(inputs = (a = 0.0, b = 0.0, c = 0.0)))
    stage!(sim, "a" => 1.0)
    step!(sim)
    @test length(trace(sim).batches) == 1 && trace(sim).frames == 1
    # A warm restart is a new trajectory: the trace is cleared and re-captured.
    init!(sim, fragment(inputs = (a = 5.0, b = 0.0, c = 0.0)))
    trc = trace(sim)
    @test isempty(trc.batches) && trc.frames == 0
    @test trc.header.root_inputs == Pair{Symbol,Any}[:a => 5.0, :b => 0.0, :c => 0.0]
end

@testset "a roster change appends the writer set; earlier records keep their schema (§11.5)" begin
    sim = Simulation(three_root_inputs(); h = 1//10)
    init!(sim, fragment(inputs = (a = 0.0, b = 0.0, c = 0.0)))
    stage!(sim, "b" => 1.0)
    step!(sim)                                   # frame 1, against the whole surface

    hd = attach!(sim, Pad("d"), Enumerated("a"))  # a stopped-sim roster change
    stage!(sim, "b" => 2.0)                      # the harness surface is {b, c} now
    stage!(hd, "a" => 5.0)
    step!(sim)                                   # frame 2, in the drain's own order

    trc = trace(sim)
    # The list only grows: the first set stays, the new one is appended whole.
    @test length(trc.header.schemas) == 3
    @test trc.header.schemas[1] == ("harness" => [:a, :b, :c])
    @test trc.header.schemas[2] == ("device 1 (Pad)" => [:a])
    @test trc.header.schemas[3] == ("harness" => [:b, :c])
    @test sim.trace.live_writers == 2:3           # the current writers, named

    (b1, b2, b3) = trc.batches
    # `b` is position 2 of the old schema and position 1 of the new one, and both
    # records resolve to the same face — which is why the list may not be rewritten.
    @test b1.frame == 1 && b1.writer == 1 && b1.entries == Pair{Int,Any}[2 => 1.0]
    @test b2.frame == 2 && b2.writer == 2 && b2.entries == Pair{Int,Any}[1 => 5.0]
    @test b3.frame == 2 && b3.writer == 3 && b3.entries == Pair{Int,Any}[1 => 2.0]
    @test recorded_faces(trc, b1) == [:b] && recorded_faces(trc, b3) == [:b]
    @test recorded_faces(trc, b2) == [:a]

    # A detach appends again, and the batches recorded under the wider set stand.
    detach!(sim, sim.plane.roster[1].dev)
    stage!(sim, "a" => 7.0)
    step!(sim)
    trc = trace(sim)
    @test length(trc.header.schemas) == 4 && sim.trace.live_writers == 4:4
    b4 = last(trc.batches)
    @test b4.frame == 3 && b4.writer == 4 && recorded_faces(trc, b4) == [:a]
    @test trc.batches[1:3] == [b1, b2, b3]        # the earlier records, untouched
end

# --- replay's entry pass (§12.7, increment 23 stage 2) ---------------------------
# Validation and normalization, tested through `_compile_feed` alone: `replay!`
# and its loop are stage 3, and the pass is what a trace has to survive before a
# single write happens.

# The same three root inputs under one extra component: same faces, a different
# `Build` — §12.7's structural line, on the error side of it.
extra_component() = Group((; s = Sum(sa = 1.0, sb = 1.0), g = Gain(2.0), k = TickCounter());
                          inputs = ("a" => "s/a", "b" => "s/b", "c" => "g/e"))

# A short recorded session over `three_root_inputs()`: two frames, one batch each.
function recorded_session()
    sim = Simulation(three_root_inputs(); h = 1//10)
    init!(sim, fragment(inputs = (a = 0.0, b = 0.0, c = 0.0)))
    stage!(sim, "b" => 1.0)
    step!(sim)
    stage!(sim, "a" => 2.0, "c" => 3.0)
    step!(sim)
    trace(sim)
end

# A fresh, initialized target of the recording's own build.
function replay_target()
    sim = Simulation(three_root_inputs(); h = 1//10)
    init!(sim, fragment(inputs = (a = 0.0, b = 0.0, c = 0.0)))
    sim
end

@testset "the scalar is dispatch, not a comparison (§12.7)" begin
    trc = recorded_session()
    err = failure(() -> _compile_feed(Simulation(three_root_inputs(), D8; h = 1//10), trc))
    d = only(err.diagnostics)
    @test err isa BuildError && d isa ReplayHeaderMismatch
    @test d.what === :scalar && d.expected === Float64 && d.found === D8
    # …and the matching pair compiles, which is what makes the refusal dispatch
    # rather than a comparison anyone could forget to write.
    @test _compile_feed(replay_target(), trc) isa ReplayFeed
end

@testset "the header is compared against the build and the deployment binding (§12.7)" begin
    trc = recorded_session()

    # An extra component is a different store layout: the path list and the
    # cell-size list both move, and both are `:store`.
    err = failure(() -> _compile_feed(Simulation(extra_component(); h = 1//10), trc))
    @test err isa BuildError && all(d isa ReplayHeaderMismatch for d in err.diagnostics)
    paths = only(d for d in err.diagnostics if d.name === :paths)
    @test paths.what === :store && paths.expected == ["s", "g"] && paths.found == ["s", "g", "k"]
    @test any(d -> d.what === :store && d.name === :sizes, err.diagnostics)

    # The seven trajectory-determining parameters, one assertion each where the
    # keyword is independently settable. `h` moves `Δt_base` with it — the two
    # are one grid — so the collection carries both.
    err = failure(() -> _compile_feed(Simulation(three_root_inputs(); h = 1//20), trc))
    d = only(x for x in err.diagnostics if x.name === :h)
    @test d.what === :deployment && d.expected === 0.1 && d.found === 0.05
    @test any(x -> x.what === :deployment && x.name === :Δt_base, err.diagnostics)

    err = failure(() -> _compile_feed(
        Simulation(three_root_inputs(); h = 1//10, localization_budget = 4), trc))
    d = only(err.diagnostics)
    @test d.what === :deployment && d.name === :localization_budget
    @test d.expected == 8 && d.found == 4

    err = failure(() -> _compile_feed(
        Simulation(three_root_inputs(); h = 1//10, firing_budget = 2), trc))
    d = only(err.diagnostics)
    @test d.what === :deployment && d.name === :firing_budget
    @test d.expected == 4 && d.found == 2

    # `t_end` and `stop_on` are a recorded fact of the recorded session, never a
    # constraint on this one, and `t₀` is applied rather than compared (§12.7).
    @test _compile_feed(Simulation(three_root_inputs(); h = 1//10, t_end = 3.0,
                                   stop_on = ()), trc) isa ReplayFeed
end

@testset "a recorded schema is validated against the target's own faces (§11.5, §12.7)" begin
    trc = recorded_session()
    # The header rebuilt with one foreign name in the harness schema: the
    # positional records are meaningless without the schema, so the schema is
    # what has to agree with this model.
    bent = _reschema(trc.header, [("harness" => [:a, :zzz, :c])])
    err = failure(() -> _compile_feed(replay_target(), Trace(bent, trc.batches, trc.frames)))
    d = only(err.diagnostics)
    @test err isa BuildError && d isa ReplaySchemaMismatch
    @test d.writer == "harness" && d.unknown == [:zzz]
    @test d.schema == [:a, :zzz, :c] && d.faces == [:a, :b, :c]

    # The header's own disagreements come first and alone: a trace that is both
    # schema-bent and entry-bent reports the schema, because a record resolved
    # through a contradicted schema would report noise (§12.7's two stages).
    bad = Trace(bent, [TraceBatch(1, 1, Pair{Int,Any}[9 => 1.0])], 1)
    @test kinds(failure(() -> _compile_feed(replay_target(), bad))) == [ReplaySchemaMismatch]
end

@testset "every position resolves through its schema, and the misses collect (§12.7)" begin
    trc = recorded_session()
    h = trc.header
    # Three bad records in one trace: a position past the schema's end, a
    # position below its start, and a batch naming a writer the schema list does
    # not have — which is the same failure with the tag spelled positionally.
    bad = Trace(h, [TraceBatch(1, 1, Pair{Int,Any}[7 => 1.0]),
                    TraceBatch(2, 1, Pair{Int,Any}[0 => 1.0, 2 => 5.0]),
                    TraceBatch(2, 9, Pair{Int,Any}[1 => 1.0])], 2)
    err = failure(() -> _compile_feed(replay_target(), bad))
    @test err isa BuildError && all(d isa ReplayUnknownFace for d in err.diagnostics)
    @test [(d.face, d.frame, d.writer) for d in err.diagnostics] ==
          [(7, 1, "harness"), (0, 2, "harness"), (1, 2, "writer #9")]
    @test all(d -> d.faces == [:a, :b, :c], err.diagnostics)

    # An unconvertible recorded value is not a face problem: the face is known
    # and the *value* is what the target's compiled scatter cannot take.
    err = failure(() -> _compile_feed(replay_target(),
                                      Trace(h, [TraceBatch(1, 1, Pair{Int,Any}[1 => "x"])], 1)))
    d = only(err.diagnostics)
    @test d isa ReplayHeaderMismatch && d.what === :root_input
    @test d.name === :a && d.expected === Float64 && d.found == "x"
end

@testset "a valid trace normalizes to compiled scatters, in drain order (§12.7, D-101)" begin
    trc = recorded_session()
    tgt = replay_target()
    feed = _compile_feed(tgt, trc)
    @test [r.frame for r in feed.records] == [1, 2]
    @test feed.next == 1 && feed.last_frame == trc.frames == 2
    # The recorded batch rides beside its thunk: a replay re-records, and what it
    # re-records is the recording's own value (§12.7).
    @test [r.record for r in feed.records] == trc.batches

    # The thunks *are* the recording, applied to this build's cells: sparse, so
    # an untouched face keeps whatever the target's own header put there.
    cells() = (port(tgt, "", :a), port(tgt, "", :b), port(tgt, "", :c))
    @test cells() == (0.0, 0.0, 0.0)
    feed.records[1].thunk()
    @test cells() == (0.0, 1.0, 0.0)
    feed.records[2].thunk()
    @test cells() == (2.0, 1.0, 3.0)

    # A superseded schema entry is compiled against the target's layout like any
    # other, so a recording that outlived a roster change replays whole.
    sim = Simulation(three_root_inputs(); h = 1//10)
    init!(sim, fragment(inputs = (a = 0.0, b = 0.0, c = 0.0)))
    stage!(sim, "b" => 1.0)
    step!(sim)
    attach!(sim, Pad("d"), Enumerated("a"))
    stage!(sim, "b" => 2.0)
    step!(sim)
    grown = _compile_feed(replay_target(), trace(sim))
    @test [r.record.writer for r in grown.records] == [1, 3]     # both schema entries live
end
