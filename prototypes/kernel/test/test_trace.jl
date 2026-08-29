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

@testset "a frame ordinal outside the recording's own length is refused (§12.7)" begin
    h = recorded_session().header
    # The drain visits each frame of `1:frames` exactly once, so a record stamped
    # outside that range names a frame that never comes round: under an exactly
    # keyed cursor it would silently never apply. Both misses collect, named by
    # the writer whose schema the record was written under.
    bad = Trace(h, [TraceBatch(0, 1, Pair{Int,Any}[1 => 9.0]),
                    TraceBatch(99, 1, Pair{Int,Any}[1 => 9.0])], 2)
    err = failure(() -> _compile_feed(replay_target(), bad))
    @test err isa BuildError && all(d isa ReplayHeaderMismatch for d in err.diagnostics)
    @test all(d -> d.what === :frame && d.name === :harness && d.expected == 1:2,
              err.diagnostics)
    @test [d.found for d in err.diagnostics] == [0, 99]

    # With no schema to name the writer either, the tag is spelled positionally —
    # and the frame is the one fault reported, the record being unreachable.
    d = only(failure(() -> _compile_feed(replay_target(),
                Trace(h, [TraceBatch(0, 9, Pair{Int,Any}[1 => 9.0])], 2))).diagnostics)
    @test d isa ReplayHeaderMismatch && d.what === :frame && d.name === Symbol("writer #9")
end

@testset "a valid trace normalizes to compiled scatters, in drain order (§12.7, D-101)" begin
    trc = recorded_session()
    tgt = replay_target()
    feed = _compile_feed(tgt, trc)
    @test [r.frame for r in feed.records] == [1, 2]
    @test feed.next == 1 && length(feed.records) == 2
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

# --- replay as the ordinary loop (§12.7, increment 23 stage 3) -------------------
# D-101's two substitutions, and the properties they exist to buy: bit-identity
# against the identical build, partial replay and the reproduction it opens,
# continuation, the discard of live staging, and the what-if register.

# One model reaching all three state homes with a localized event in the middle:
# the reference-fed plant of `feedback_model` (continuous state, a feedback
# wire), a discrete integrator on a second root input (discrete state), and an
# autonomous bouncer whose reset lands *between* frame tops — so a replay has to
# reproduce `t*` boundaries (§10.4) and mode state, not just a grid of frames.
# `k` is the parameter the what-if register moves.
replay_model(k = 4.0) =
    Group((plant = Plant(; ω = 2.0, ζ = 0.1), ctl = Gain(k), sum = Sum(),
           acc = DiscreteIntegrator(1.0), b = Bouncer(1.0, 0.32));
          wires = ("ctl/out" => "plant/u", "sum/e" => "ctl/e", "plant/y" => "sum/b"),
          inputs = ("ref" => "sum/a", "rate" => "acc/e"),
          outputs = ("b/q" => "bq", "acc/u" => "u"))

# A device that stages once from its own task and departs (§12.4(6)'s voluntary
# exit): deterministic in *what* it stages — which is what a replay test can
# assert — never in when, the frame a wall-clock task reaches being exactly the
# scheduler-determined timing §11.1 indicts.
mutable struct Nudge <: AbstractDevice
    face::String
    v::Float64
end
loop(d::Nudge, h) = (stage!(h, d.face => d.v); nothing)

# Every cell of a published table, in a build-independent order: two sessions
# are two `Simulation`s, so their layouts are equal *values* rather than one
# object, and the sorted address list is what makes "the same table" mean the
# same thing on both sides. Bit-identity is asserted with `==` — that is the
# claim (§12.7), D-163's tolerance rule being about numerical agreement.
snap_cells(s::Snapshot) =
    Any[gather(s.store, s.layout.addr[k]) for k in sort!(collect(keys(s.layout.addr)))]

# Two sessions' logs, boundary for boundary: the `t` stamps and every cell.
same_trajectory(a, b) =
    length(a) == length(b) && all(x.t == y.t && x.frame == y.frame &&
                                  snap_cells(x) == snap_cells(y) for (x, y) in zip(a, b))

# The frame-top publication of frame `k`: a frame with a `t*` boundary publishes
# more than once under the same ordinal, and the frame's own boundary is last.
at_frame(snaps, k::Int) = last(s for s in snaps if s.frame == k)

# A recorded session over `replay_model()`, staged from the harness across
# frames: eight frames, batches at 1 and 4, two localized resets inside.
function recorded_run()
    sim = Simulation(replay_model(); h = 1//10)
    init!(sim, fragment(inputs = (ref = 1.0, rate = 0.0)))
    stage!(sim, "ref" => 2.0)
    step!(sim; frames = 3)
    stage!(sim, "rate" => 1.0)
    step!(sim; frames = 5)
    (sim, trace(sim))
end

# A fresh target of that build, initialized from a *different* condition — so
# nothing but the header can put it on the recorded trajectory.
function replay_twin(k = 4.0; kw...)
    sim = Simulation(replay_model(k); h = 1//10, kw...)
    init!(sim, fragment(inputs = (ref = 0.0, rate = 0.0)))
    sim
end

@testset "a replay reproduces the recorded trajectory bitwise (§12.7, D-101)" begin
    (sim, trc) = recorded_run()
    sim2 = replay_twin()
    replay!(sim2, trc)

    # §12.7: the replay ends `initialized`, never `stopped` — boundary-consistent
    # and ready to advance, which is what makes inspect/step/continue real.
    @test lifecycle(sim2) === :initialized && termination(sim2) === nothing
    @test sim2.exec.clock.step == trc.frames
    @test same_trajectory(logged(sim2), logged(sim))
    # …including the localized boundaries the recording never stored: `t*` is
    # derived from state (§10.4), so reproducing the state reproduces the timing.
    @test length(logged(sim)) == trc.frames + 3     # boundary zero + two resets
    # the live stores agree too, the snapshots deliberately not carrying them
    @test state(sim2, "plant").q === state(sim, "plant").q
    @test state(sim2, "acc") === state(sim, "acc") && modes(sim2, "b") === modes(sim, "b")
    # and the replay re-records: the header inherited, the batches re-drained
    @test trace(sim2).frames == trc.frames
    @test trace(sim2).batches == trc.batches

    # …as an *equal value*, never the recording's own objects: the two traces are
    # two values, so a continuation's growth cannot reach the `Trace` in hand.
    rec = trace(sim2)
    @test all(rec.batches[i] !== trc.batches[i] for i in eachindex(trc.batches))
    @test rec.header.x !== trc.header.x && rec.header.x == trc.header.x
    @test rec.header.s !== trc.header.s && rec.header.m !== trc.header.m
    @test rec.header.root_inputs !== trc.header.root_inputs
    @test rec.header.root_inputs == trc.header.root_inputs
    @test rec.header.schemas !== trc.header.schemas

    # §12.7's own register — `Simulation(world)` then `replay!`, no `init!` at
    # all: `replay!` *is* a door into `initialized`, so it owes one to nothing.
    raw = Simulation(replay_model(); h = 1//10)
    @test lifecycle(raw) === :built
    replay!(raw, trc)
    @test lifecycle(raw) === :initialized
    @test same_trajectory(logged(raw), logged(sim))
end

@testset "a device's recorded batches replay on a deviceless twin (§12.7)" begin
    sim = Simulation(replay_model(); h = 1//10, t_end = 2.0)
    attach!(sim, Nudge("rate", 3.0), Enumerated("rate"))
    init!(sim, fragment(inputs = (ref = 1.0, rate = 0.0)))
    stage!(sim, "ref" => 2.0)      # the harness surface is {ref}: the device holds {rate}
    run!(sim)
    trc = trace(sim)
    @test lifecycle(sim) === :stopped

    # "No devices or mappings present" (§12.7): the recorded batches carry the
    # device's writes, and the schemas resolve them against this build's faces.
    sim2 = replay_twin()
    replay!(sim2, trc)
    @test lifecycle(sim2) === :initialized
    @test same_trajectory(logged(sim2), logged(sim))
    @test state(sim2, "acc") === state(sim, "acc")
    @test trace(sim2).batches == trc.batches
end

@testset "partial replay halts at a frame top, and the next frame reproduces (§12.7, §13.4)" begin
    (sim, trc) = recorded_run()
    sim2 = replay_twin()
    replay!(sim2, trc; to_boundary = 5)
    @test lifecycle(sim2) === :initialized      # the pointer's register: ready to advance
    @test sim2.exec.clock.step == 5 * sim2.n
    @test trace(sim2).frames == 5
    @test same_trajectory(logged(sim2), [s for s in logged(sim) if s.frame ≤ 5])

    # §13.4's workflow minus the error: step the next frame under whatever
    # instrumentation is wanted, and it is the recording's frame 6 bitwise.
    @test step!(sim2) == 1
    @test snap_cells(latest(sim2)) == snap_cells(at_frame(logged(sim), 6))
    @test latest(sim2).t == at_frame(logged(sim), 6).t
end

@testset "a continuation is a live session from the replayed boundary (§12.7)" begin
    (sim, trc) = recorded_run()
    sim2 = replay_twin()
    replay!(sim2, trc)
    run!(sim2; t_end = 1.4)                     # `run!` after `replay!`
    @test lifecycle(sim2) === :stopped && termination(sim2).source === EndTimeReached()
    @test sim2.exec.clock.step == 14            # it proceeded from frame 8, not from zero
    # The session leaves behind a complete, valid trace of *itself*, with the
    # recording as a bit-identical prefix and no special stitching (§12.7).
    cont = trace(sim2)
    @test cont.frames == 14 > trc.frames
    @test cont.batches[1:length(trc.batches)] == trc.batches
    @test cont.header.root_inputs == trc.header.root_inputs   # the header inherited
    # the recording's schema entries stand, this session's appended behind them
    @test cont.header.schemas[1:length(trc.header.schemas)] == trc.header.schemas
    @test sim2.trace.live_writers == (length(trc.header.schemas) + 1):length(cont.header.schemas)
end

# Every discard this writer's account carried, rendered: off a published
# status's `recent` where the frame-top fold caught it (§11.8), and off the
# run's-end sweep's warning where the device's stage landed past the final frame
# top (§12.4) — a replay ends `initialized`, so the sweep has no termination
# record to file its residue in and the log is where it surfaces.
discard_reports(sim, logs, who::String) =
    vcat(String[string(d) for s in logged(sim) for w in s.status.writers if w.who == who
                for d in w.recent if d isa ReplayDiscardedStaging],
         String[string(l.message) for l in logs
                if occursin("ReplayDiscardedStaging from $who, past the final",
                            string(l.message))])

@testset "live staging met by a replay is discarded, and reported (§12.7, §11.8)" begin
    (sim, trc) = recorded_run()
    sim2 = replay_twin()
    dev = Nudge("rate", 99.0)                   # a value that would move the trajectory
    attach!(sim2, dev, Enumerated("rate"))
    logs, _ = Test.collect_test_logs() do
        replay!(sim2, trc)                      # devices are readers: they init and spawn
    end
    @test lifecycle(sim2) === :initialized

    # The property the discard exists to protect: the trajectory is still the
    # recording's, bit for bit, with a live writer staging into the run.
    @test same_trajectory(logged(sim2), logged(sim))
    @test state(sim2, "acc") === state(sim, "acc")
    @test trace(sim2).batches == trc.batches    # and nothing of the device's is recorded

    # …and the drop is loud: one report, on the device's own cell (§11.8's
    # attribution), naming the faces it cost. Wherever the device's timing put
    # it — the account's totals, or the run's-end sweep — it is accounted once.
    who = "device 1 (Nudge)"
    @test accounted(sim2, logs, who, :replay_discarded, "ReplayDiscardedStaging")
    r = only(discard_reports(sim2, logs, who))
    @test occursin("[:rate]", r)
    seen = [d for s in logged(sim2) for w in s.status.writers if w.who == who
              for d in w.recent if d isa ReplayDiscardedStaging]
    @test all(d -> d.faces == [:rate] && 1 ≤ d.frame ≤ trc.frames, seen)
end

# The harness arm of the same discard, with a synchronous route into it: a
# `needs_calling_task` device runs its loop body *inline* on the calling task
# while the run body is spawned (§11.1), so its `stage!(sim, …)` — the harness
# register's own surface, not the device's claim — lands mid-replay, frame after
# frame, where a spawned device's timing is the scheduler's alone.
mutable struct HarnessPoker <: AbstractDevice
    sim::Any
end
needs_calling_task(::HarnessPoker) = true
loop(d::HarnessPoker, h) = (while running(h); stage!(d.sim, "ref" => 99.0); yield(); end;
                            nothing)

@testset "live staging into the harness is discarded on its own cell (§12.7, §11.8)" begin
    # A long recording, so the inline body is scheduled against a run with frames
    # left to give it: the poker's stage is deterministic, its timing never is.
    sim = Simulation(replay_model(); h = 1//10)
    init!(sim, fragment(inputs = (ref = 1.0, rate = 0.0)))
    stage!(sim, "ref" => 2.0)
    step!(sim; frames = 400)
    trc = trace(sim)

    sim2 = replay_twin()
    dev = HarnessPoker(nothing)
    attach!(sim2, dev, Enumerated("rate"))     # the device claims `rate`; `ref` stays the harness's
    dev.sim = sim2
    logs, _ = Test.collect_test_logs() do
        replay!(sim2, trc)
    end
    @test lifecycle(sim2) === :initialized

    # The same property the device arm buys, through the other writer: the
    # trajectory is the recording's bit for bit, and none of the poker's 99.0 is
    # recorded — the header's `ref` survives to the end.
    @test same_trajectory(logged(sim2), logged(sim))
    @test trace(sim2).batches == trc.batches
    @test port(sim2, "", :ref) == 2.0

    # …and the drop is loud on the *harness* cell (§11.8's attribution), naming
    # the faces it cost — the device's own account untouched, it never staged.
    seen = [d for s in logged(sim2) for w in s.status.writers if w.who == "harness"
              for d in w.recent if d isa ReplayDiscardedStaging]
    @test !isempty(discard_reports(sim2, logs, "harness"))
    @test all(d -> d.faces == [:ref] && 1 ≤ d.frame ≤ trc.frames, seen)
    @test writer_status(latest(sim2), "harness").totals.replay_discarded ≥ length(seen) ≥ 1
    @test writer_status(latest(sim2), "device 1 (HarnessPoker)").totals.replay_discarded == 0
end

@testset "a changed parameter replays deterministically: the what-if register (§12.7)" begin
    (sim, trc) = recorded_run()
    # Same structure, a different gain — *parametric* difference, on the
    # non-error side of §12.7's line: the recorded inputs re-driven through a
    # modified model. Determinism is promised; reproduction is not.
    whatif = replay_twin(9.0)
    replay!(whatif, trc)
    @test lifecycle(whatif) === :initialized
    @test whatif.exec.clock.step == trc.frames
    @test state(whatif, "plant").q != state(sim, "plant").q

    # Deterministic: the same what-if twice is the same trajectory.
    twin = replay_twin(9.0)
    replay!(twin, trc)
    @test same_trajectory(logged(twin), logged(whatif))
end

@testset "the lifecycle and range refusals of `replay!` (§12.7, §12.6)" begin
    (_, trc) = recorded_run()
    # `to_boundary` is a pointer into the recording: past its end there is
    # nothing to replay, and the refusal names the argument and its value.
    for bad in (9, -1, 2.5)
        d = only(failure(() -> replay!(replay_twin(), trc; to_boundary = bad)).diagnostics)
        @test d isa ArgumentInvalid && d.call === :replay! && d.reason === :range
        @test d.argument === :to_boundary && d.value == bad
    end
    @test lifecycle(replay_twin()) === :initialized      # a rejected replay wrote nothing

    # `errored` is terminal (§13.6): never resumable, never re-initialized, and
    # `replay!` is refused there exactly as `init!` is — reproduction is
    # replaying the trace on a *fresh* simulation, which is the arm above.
    crashed = Simulation(fed(Exploder(), "arm"); h = 1//10, t_end = 5.0)
    init!(crashed, fragment(inputs = (in = 0.0,)))
    own = trace(crashed)
    stage!(crashed, "in" => true)
    @test_throws Exploded run!(crashed)
    d = only(failure(() -> replay!(crashed, own)).diagnostics)
    @test d isa ServiceLifecycle && d.op === :replay! && d.status === :errored
    # (`replay!` from `:running` is the same gate one line above it, and reaching
    # it needs the spawned-run register `test_lifecycle.jl` exercises for `init!`
    # and `run!`; it is asserted there, for those two entries, and not here.)
end

@testset "a replay runs under the kill switch, and records nothing (§11.5, §12.7)" begin
    (sim, trc) = recorded_run()
    off = replay_twin(; trace = false)
    replay!(off, trc)                      # the feed is compiled from the `Trace` in hand,
    @test lifecycle(off) === :initialized   # never from the target's own register
    @test same_trajectory(logged(off), logged(sim))
    @test off.trace.header === nothing && isempty(off.trace.batches) && off.trace.frames == 0
    @test only(failure(() -> trace(off)).diagnostics).reason === :disabled
end
