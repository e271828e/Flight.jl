# --- the diagnostic channel and the framework status (§11.8; increments 13–14) --
# report! from the author's loop body, the per-writer single-writer cells with
# their ring-plus-counts bound, the frame-top fold beside the staging drain,
# the published status — delta plus totals, heartbeat and task_state — and the
# run's-end sweep. The devices below live at top level for the README's
# local-scope reason.

# A datum-stream parser: the §11.6 tolerance idiom verbatim — catch its own
# parser error, stage nothing for that datum, report, continue with the next.
mutable struct Parser <: AbstractDevice
    datums::Vector{Any}
    fired::Bool
end
Parser(datums...) = Parser(collect(Any, datums), false)
function loop(d::Parser, h)
    d.fired && return nothing
    d.fired = true
    for datum in d.datums
        try
            datum isa Number || throw(ArgumentError("unparseable: $(repr(datum))"))
            stage!(h, "a" => datum)
        catch e
            e isa ArgumentError || rethrow()     # a bug → wrapper → DeviceCrash
            report!(h, MalformedDatum(e))        # garbage → visible, bounded, alive
        end
    end
    stop!(h)
    nothing
end

# A late reporter: its one report races the run's last frame top, so its
# account may land in the terminal status or in the run's-end sweep.
mutable struct LateReporter <: AbstractDevice
    fired::Bool
end
LateReporter() = LateReporter(false)
function loop(d::LateReporter, h)
    d.fired && return nothing
    d.fired = true
    report!(h, MalformedDatum("late"))
    stop!(h)
    nothing
end

# A snapshot consumer that stays in its loop: the liveness testset's device.
mutable struct Ticker <: AbstractDevice
    n::Int
end
Ticker() = Ticker(0)
function loop(d::Ticker, h)
    while running(h)
        wait_next_snapshot(h)
        (d.n += 1) ≥ 3 && return stop!(h)
    end
    nothing
end

@testset "a bad datum is tolerated: catch, stage nothing, report, continue (§11.6)" begin
    sim = Simulation(two_root_inputs(); h = 1//10)
    dev = Parser(0.7, "garbage", 0.9)
    hp = attach!(sim, dev, Enumerated("a"))
    init!(sim, fragment(inputs = (a = 0.0, b = 0.0)))
    logs, _ = Test.collect_test_logs() do
        run!(sim; t_end = 1000.0)                        # ends by the device's stop
    end
    msgs = [string(l.message) for l in logs]
    # The link survived its truncated datagram: no crash, and the report is
    # accounted device-attributed exactly once — in the terminal status when a
    # frame top folded it, in the run's-end sweep when none remained (§11.8).
    @test !any(occursin("DeviceCrash", m) for m in msgs)
    @test accounted(sim, logs, "device 1 (Parser)", :malformed, "MalformedDatum")
    # The author's cause survives wherever the record landed: in some logged
    # snapshot's recent, or in the sweep's presentation.
    carried = any(d isa MalformedDatum && occursin("unparseable", string(d.cause))
                  for s in logged(sim)
                  for d in writer_status(s, "device 1 (Parser)").recent)
    @test carried ⊻ any(occursin("unparseable", m) for m in msgs)
    # The stream's good datums survived — newest wins within the staged batch —
    # applied by a drain the stop did not beat, or still pending in the cell:
    # exactly one of the two, timing's choice.
    p = @atomic hp.writer.cell.pending
    @test (p === nothing ? port(sim, "", :a) : p[].vals[1]) === 0.9
end

@testset "the ring's bound is the rate limit: 16 retained, excess to the counts (§11.8)" begin
    sim = Simulation(two_root_inputs(); h = 1//10)
    h = attach!(sim, Pad("p"), Enumerated("a"))
    init!(sim, fragment(inputs = (a = 0.0, b = 0.0)))
    for k in 1:20                                # one frame's flood, pending in the cell
        report!(h, MalformedDatum("datum $k"))
    end
    run!(sim; t_end = 0.1)                               # the first frame top folds the cell
    mw = writer_status(latest(sim), "device 1 (Pad)")
    # Earliest-in-frame retained: the first occurrences carry the diagnostic
    # content, the excess becomes exactly a per-kind count beside them.
    @test [d.cause for d in mw.recent] == ["datum $k" for k in 1:DIAG_RING]
    @test mw.suppressed.malformed == 4
    # Nothing is lost by not looking: the totals carry the full account.
    @test mw.totals.malformed == 20
end

@testset "the status: the delta rides one snapshot, totals ride every one (§11.8, §11.2)" begin
    sim = Simulation(two_root_inputs(); h = 1//10)
    h = attach!(sim, Pad("p"), Enumerated("a"))
    init!(sim, fragment(inputs = (a = 0.0, b = 0.0)))
    # Boundary zero's status: the writers in the drain's order — devices in
    # attachment order, the harness register, the loop — every account zero,
    # and no run task to be alive: device tasks are run-scoped observables.
    st0 = latest(sim).status
    @test [w.who for w in st0.writers] == ["device 1 (Pad)", "harness", "loop"]
    @test all(w.totals == KindCounts() && isempty(w.recent) for w in st0.writers)
    @test writer_status(latest(sim), "device 1 (Pad)").task_state === :none
    # The harness's and the loop's records have no task and no heartbeat to
    # judge: never stale, `nothing` for both fields.
    @test writer_status(latest(sim), "loop").heartbeat === nothing
    @test !stale(writer_status(latest(sim), "harness"))
    report!(h, MalformedDatum("one"))            # pending before the run: folded at frame 1's top
    run!(sim; t_end = 0.5)
    snaps = logged(sim)                          # boundary zero, then frames 1..5
    dw(s) = writer_status(s, "device 1 (Pad)")
    # Exactly one snapshot carries the occurrence in `recent` — the first
    # published after the fold — while `totals` is monotone from there on:
    # a 60 Hz reader sees it once, an occasional sampler still reads the
    # complete account, and log decimation loses *which* boundary, never
    # *how many* (§11.8).
    @test [length(dw(s).recent) for s in snaps] == [0, 1, 0, 0, 0, 0]
    @test [dw(s).totals.malformed for s in snaps] == [0, 1, 1, 1, 1, 1]
    @test only(dw(snaps[2]).recent).cause == "one"
end

@testset "liveness: heartbeat and task_state ride the device's record (§11.8, §12.2, §12.4)" begin
    sim = Simulation(two_root_inputs(); h = 1//10)
    dev = Ticker()
    attach!(sim, dev, Enumerated())
    init!(sim, fragment(inputs = (a = 0.0, b = 0.0)))
    run!(sim; t_end = 1000.0)                            # ends by the device's stop, ≥ 3 boundaries in
    @test dev.n ≥ 3
    tw = writer_status(latest(sim), "device 1 (Ticker)")
    # The device consumed boundaries through the handle primitives, each pass
    # storing the heartbeat: the terminal record reads fresh, with a live (or
    # by now returned) task — never `:none`, never the stale silence of a
    # device that was never there (§12.2: a starved or dead device shows as a
    # stale heartbeat with a name on it).
    @test tw.heartbeat > 0 && !stale(tw)
    @test tw.task_state in (:running, :done)
end

@testset "the run's-end sweep: past the final frame top, loud rather than lost (§11.8, §12.4)" begin
    sim = Simulation(two_root_inputs(); h = 1//10)
    dev = LateReporter()
    attach!(sim, dev, Enumerated())
    init!(sim, fragment(inputs = (a = 0.0, b = 0.0)))
    logs, _ = Test.collect_test_logs() do
        run!(sim; t_end = 1000.0)
    end
    # The report races the last frame top: folded into the terminal status, or
    # — no drain remaining — presented by the sweep, the tail's renderer of
    # last resort. Exactly one account either way (§11.8, D-201).
    @test accounted(sim, logs, "device 1 (LateReporter)", :malformed, "MalformedDatum")
    # A fresh trajectory opens a fresh account (§11.8): init! resets the
    # totals, and the stepped frames — deviceless, the reporter never respawns —
    # publish a zeroed record for it.
    init!(sim, fragment(inputs = (a = 0.0, b = 0.0)))
    step!(sim; frames = 2)
    @test writer_status(latest(sim), "device 1 (LateReporter)").totals.malformed == 0
end

@testset "the cell is kind-generic: mixed rings, per-kind suppression (§11.8, §13.2)" begin
    # The closed set rides one union; the counter record is fixed-shape isbits,
    # a type rather than a lookup (§11.8), and + is the fold between records.
    @test isbitstype(KindCounts)
    a = _bump(_bump(KindCounts(), :malformed), :crash)
    b = _bump(KindCounts(), :malformed)
    @test (a + b).malformed == 2 && (a + b).crash == 1
    @test _total(a + b) == 3

    # A raw cell takes any kind of the set; the ring preserves arrival order
    # across kinds, earliest-in-frame retained.
    cell = DiagCell(EMPTY_DIAG)
    _report!(cell, MalformedDatum("m1"))
    _report!(cell, OutOfClaimEntry(:flaps, 1.0, [:a, :b], nothing))
    _report!(cell, ClaimedFaceEntry(:a, "device 1 (Pad)", 2.0, :staging))
    _report!(cell, EntryTypeMismatch(:b, "high", Float64))
    _report!(cell, ChatteringBudget("c", :pop, 0.1, 8, 8))
    _report!(cell, FiringBudget("e", :up, 0.0, 4, 4))
    _report!(cell, DeviceCrash(ErrorException("boom"), false))
    _report!(cell, DeviceJoinTimeout("device 9 (Ghost)", 5.0, 1.0, 10))
    batch = _take!(cell)
    @test length(batch.ring) == 8
    @test batch.ring[1] isa MalformedDatum && batch.ring[8] isa DeviceJoinTimeout
    @test batch.suppressed == KindCounts()

    # Past the bound, suppression counts by kind: the record answers "how many
    # of what", not one blurred integer.
    for k in 1:DIAG_RING
        _report!(cell, MalformedDatum("datum $k"))
    end
    _report!(cell, MalformedDatum("late m"))
    _report!(cell, DeviceCrash(ErrorException("late c"), true))
    _report!(cell, DeviceCrash(ErrorException("later c"), true))
    batch = _take!(cell)
    @test length(batch.ring) == DIAG_RING
    @test all(d isa MalformedDatum for d in batch.ring)
    @test batch.suppressed.malformed == 1 && batch.suppressed.crash == 2
    @test _total(batch.suppressed) == 3
    # The take re-arms the sentinel: the next report allocates afresh on the
    # writer's side, and the drained batch is frozen.
    @test _take!(cell) === EMPTY_DIAG
end

@testset "the heartbeat rides in the cell, stored by the handle primitives (§11.8, §12.2)" begin
    sim = Simulation(two_root_inputs(); h = 1//10)
    h = attach!(sim, Pad("p"), Enumerated("a"))
    # A never-heartbeated cell reads its initial 0.0 — stale against any wall
    # clock, which is what marks a failed init!'s device with no machinery
    # (§12.4): dead is "no recent timestamp", recorded nowhere else.
    @test _heartbeat(h.diag) == 0.0
    before = time()
    running(h)                            # any loop-pass primitive beats —
    hb = _heartbeat(h.diag)               # (false here: no run has started)
    @test hb ≥ before
    stage!(h, "a" => 1.0)                 # and so does each of the others
    @test _heartbeat(h.diag) ≥ hb
end

# --- the kind set and its carrier (§13.1, §13.2, Appendix C; increment 22) ------
# One plausible occurrence per kind, so the constructors, `severity`, `path` and
# `message` are exercised for every one of them, and the severity column of
# Appendix C is checked against the listed warning set rather than trusted.

using InteractiveUtils: subtypes    # the coverage check below

@testset "diagnostic kinds (§13.2, Appendix C, D-214, D-215)" begin
    occurrences = Diagnostic[
        # Stratum A
        UnknownPort(entry = "child_connections at `a`, entry `x => y`", end_ = :destination,
                    path = "a/b", spelling = "b/throtle", port = :throtle,
                    candidates = [:throttle, :mixture]),
        UnknownPort(entry = "input_connections at `a`, entry `:u => ()`", end_ = :connection,
                    path = "a", port = :u),
        UnconnectedInput(path = "a/b", face = :u),
        TwoProducers(path = "a/b", port = :u, incumbent = "a sibling wire",
                     entry = "an interface connection"),
        WireTypeMismatch(path = "a/b", face = :u, declared = Float64, producer_path = "a/c",
                         producer_port = :y, observed = Bool, activation = Float64),
        PathResolution(entry = "`resolve` on `Group`", spelling = "a/b/c",
                       reason = :not_a_terminal),
        PathResolution(entry = "`resolve` on `Group`", spelling = "z/y", reason = :unknown_child,
                       owner = "`a`", segment = "z", candidates = ["b", "c"]),
        PathResolution(entry = "`resolve` on `Group`", spelling = "a/b/c",
                       reason = :reaches_past, level = "a/b", segment = "b", tail = 1),
        PathResolution(entry = "`resolve` on `Group`", spelling = "", reason = :empty_path),
        StoreWithoutUpdate(path = "a/b", store = :init_x),
        EventHalfMissing(path = "a/b", event = :snap, reason = :guard, found = Int),
        EventHalfMissing(path = "a/b", event = :snap, reason = :not_an_event, found = Int),
        ClassUnreadable(path = "a", families = "`init_x`, `init_s`", holds_components = true),
        ClassMixed(path = "a", declarations = [:init_x, :output_types]),
        ContainerMixed(path = "a", field = :kids, types = Any[Int, Float64]),
        DeclarationOnWrongTier(path = "a/b", declaration = :h_x, reason = :tier_form,
                               found = :continuous, announced = :discrete),
        DeclarationOnWrongTier(path = "a/b", declaration = :project, reason = :continuous_only,
                               found = :discrete),
        DeclarationOnWrongTier(path = "a/b", declaration = :project, reason = :no_manifold),
        FaceNameIllegal(path = "a", face = "u/v"),
        FaceNameCollision(path = "a", faces = ["u"], site = :assembly),
        FaceNameCollision(path = "", faces = ["u"], site = :root),
        FaceDirectionConflict(entry = "child_connections at `a`", path = "a/b",
                              spelling = "b/u", found = :input, wanted = :producer),
        UnknownFaceSelection(who = "input_passthrough", path = "a/b", reason = :both_given),
        UnknownFaceSelection(who = "input_passthrough", path = "a/b",
                             reason = :unknown_names, names = ["q"], candidates = ["u", "v"]),
        RatesViolation(path = "a", reason = :declaration_shape),
        RatesViolation(path = "a", reason = :value_vocabulary, key = :b, value = 3),
        RatesViolation(path = "a", reason = :multiplier, key = :b, value = 0),
        RatesViolation(path = "a", reason = :phase, key = :b, value = 5),
        RatesViolation(path = "a", reason = :period, key = :b, value = 0//1),
        RatesViolation(path = "a", reason = :offset, key = :b, value = 3//1),
        RatesViolation(path = "a", reason = :unknown_child, key = :z, candidates = ["b"]),
        RatesViolation(path = "a", reason = :continuous_child, key = :b),
        ChildNameCollision(path = "a", name = "b", reason = :sample_times_sugar,
                           provenance = ["container field `kids`, element `b`"], field = :kids),
        ChildNameCollision(path = "a", name = "b", reason = :sibling_field,
                           provenance = ["container field `kids`, element `b`"], field = :kids),
        ChildNameCollision(path = "a", name = "b", reason = :two_children,
                           provenance = ["field `b`", "container field `kids`, element `b`"]),
        TransparentContainerUnknown(path = "a", field = :kids, component = "Group"),
        TierUnreadable(path = "a/b", declarations = [:init_m]),
        IllegalPortType(path = "a/b", site = :port, name = :y, declared = Nothing),
        # Strata B and C
        AlgebraicCycle(members = ["a/b", "a/c"]),
        ProducedByTwoStages(path = "a/b", ports = [:y]),
        DeclaredNotProduced(path = "a/b", ports = [:y], products = [:z]),
        UndeclaredReturnField(path = "a/b", stage = "h_x", name = :q, candidates = [:y]),
        ConformanceFailure(path = "a/b", what = "h_x", reason = :return_type, shape = :ports,
                           observed = Int),
        ConformanceFailure(path = "a/b", what = "project", reason = :field_set, shape = :state,
                           observed_fields = [:p], declared_fields = [:q]),
        ConformanceFailure(path = "a/b", what = "f", reason = :field_type, shape = :init_x,
                           field = :q, observed = Float64, declared = Bool),
        ConformanceFailure(path = "a/b", what = "h_x", reason = :field_type, shape = :ports,
                           field = :y, observed = Int, declared = Float64,
                           activation = Float64),
        ConformanceFailure(path = "a/b", what = "g", reason = :field_set, shape = :init_s,
                           observed = NamedTuple{(:q,),Tuple{Int}},
                           declared = NamedTuple{(:q,),Tuple{Float64}}),
        ConformanceFailure(path = "a/b", what = "event `e`'s handler", reason = :field_set,
                           shape = :mode, field = :k, declared_fields = [:m]),
        GuardForm(path = "a/b", event = :snap, observed = Int),
        HandlerReturnKey(path = "a/b", event = :snap, key = :s, stores = [:x, :m]),
        # Deployment, periphery and services
        MissingInit(op = "run!", status = :built),
        ServiceLifecycle(op = "attach!", status = :running, legal = [:built, :initialized]),
        ServiceLifecycle(op = "capture", status = :errored),
        ServiceLifecycle(op = "run!", status = :stopped),
        ServiceLifecycle(op = "capture", status = :built, legal = [:initialized, :stopped]),
        StopFaceInvalid(face = :done, reason = :unknown, candidates = [:hit]),
        StopFaceInvalid(face = :done, reason = :root_input),
        StopFaceInvalid(face = :done, reason = :not_bool, declared = Float64),
        DeploymentInvalid(parameter = :firing_budget, reason = :range, value = 0),
        DeploymentInvalid(parameter = :h, reason = :inexact, value = 0.01),
        DeploymentInvalid(parameter = :Δt_base, reason = :not_a_quantity, value = Int),
        DeploymentInvalid(parameter = :h, reason = :missing),
        DeploymentInvalid(parameter = :Δt_base, reason = :unanchored, paths = ["a/b"]),
        DeploymentInvalid(parameter = :Δt_base, reason = :no_constraint),
        DeploymentInvalid(parameter = :Δt_base, reason = :not_harmonic, value = 1//3,
                          related = 1//100),
        DeploymentInvalid(parameter = :Δt_base, reason = :disagrees_with_n, value = 1//50,
                          related = 3, quotient = 2),
        DeploymentInvalid(parameter = :Δt_base, reason = :anchor_period, value = 1//30,
                          related = 1//100, provenance = "`sample_times` at `a`, key `b`",
                          admissible = 1//300),
        DeploymentInvalid(parameter = :Δt_base, reason = :anchor_offset, value = 1//7,
                          related = 1//100, provenance = "`sample_times` at `a`, key `b`",
                          admissible = 1//300),
        AttachUnknownFace(binding = "Enumerated", face = :q, candidates = [:a, :b]),
        AlreadyAttached(device = "Pad", incumbent = "device 1 (Pad)", binding = "Enumerated"),
        CallerTaskConflict(device = "Pad", incumbent = "device 1 (Poller)"),
        ClaimConflict(face = :a, device = "Pad", incumbent = "device 1 (Pad)"),
        EmptyGreedyClaim(device = "device 2 (Pad)", binding = "GreedyPlus"),
        BindingContractMismatch(binding = "NoEnum", reason = :claims_missing),
        BindingContractMismatch(binding = "NoEnum", reason = :reads_missing),
        BindingContractMismatch(binding = "NoEnum", reason = :greedy_without_input),
        BindingContractMismatch(binding = "NoEnum", reason = :neither_side),
        BindingContractMismatch(binding = "NoEnum", reason = :greedy_with_claims),
        BindingContractMismatch(binding = "NoEnum", reason = :claims_without_input),
        BindingContractMismatch(binding = "NoEnum", reason = :reads_without_output),
        BindingContractMismatch(binding = "NoEnum", reason = :reads_not_namedtuple,
                                observed = Int),
        BindingContractMismatch(binding = "NoEnum", reason = :reads_not_selectors,
                                observed = "1"),
        ReadBindingUnresolved(binding = "T", selector = "get_state(\"a\", :x)",
                              reason = :store_selector, path = "a", field = :x),
        ReadBindingUnresolved(binding = "T", selector = "get_output(\"a\", :y[1])",
                              reason = :indexed, path = "a", field = :y),
        ReadBindingUnresolved(binding = "T", selector = "get_output(\"a\", :y)",
                              reason = :unknown_cell, path = "a", field = :y),
        ReadBindingUnresolved(binding = "T", selector = "get_input(:q)",
                              reason = :unknown_root_input, field = :q, candidates = [:a]),
        ReadBindingUnresolved(binding = "T", selector = "get_face(:a)",
                              reason = :root_input_not_output, field = :a),
        ReadBindingUnresolved(binding = "T", selector = "get_face(:q)",
                              reason = :unknown_output_face, field = :q),
        ConditionResolution(path = "a", field = :u, reason = :assembly_path, provenance = "fragment"),
        ConditionResolution(path = "z", field = :u, reason = :unknown_path, provenance = "fragment"),
        ConditionResolution(field = :q, reason = :unexported_face, candidates = [:a],
                            provenance = "fragment"),
        ConditionResolution(path = "a", field = :q, reason = :no_input_face,
                            candidates = [:u], provenance = "fragment"),
        ConditionResolution(path = "a", field = :u, reason = :internally_wired,
                            producer = ("a/b", :y), provenance = "fragment"),
        ConditionResolution(path = "a/b", store = :x, field = :q, reason = :no_store,
                            tier = :discrete, provenance = "fragment"),
        ConditionResolution(path = "a/b", store = :x, field = :q, reason = :undeclared_field,
                            candidates = [:p], role = :output_port, provenance = "fragment"),
        ConditionResolution(path = "a/b", store = :x, field = :p, reason = :unconvertible,
                            declared = Float64, observed = String, value = "x", provenance = "fragment"),
        DuplicateConditionLeaf(path = "a/b", store = :x, field = :p,
                               provenance = ["fragment", "at(\"a\") → fragment"]),
        ConditionNodeMisuse(observed = NamedTuple, in_hand = [:Fragment]),
        ConditionNodeMisuse(observed = Int, reason = :fragment_payload, payload = :x),
        UninitializedInputs(op = "init!", faces = [:a, :b]),
        TapResolution(label = :r, selector = "get_state(\"a\", :x)", reason = :assembly_path,
                      path = "a"),
        TapResolution(label = :r, selector = "get_state(\"z\", :x)", reason = :unknown_path,
                      path = "z"),
        TapResolution(label = :r, selector = "get_state(\"a/b\", :x[1])",
                      reason = :scalar_index, path = "a/b", declared = Float64, index = 1),
        TapResolution(label = :r, selector = "get_state(\"a/b\", :q)", reason = :undeclared,
                      path = "a/b", field = :q, declares = :state_field, candidates = [:p]),
        TapResolution(label = :r, selector = "get_deriv(\"a/b\", :q)",
                      reason = :discrete_deriv, path = "a/b", field = :q),
        TapResolution(label = :r, selector = "get_input(:q)", reason = :unknown_root_input,
                      field = :q, candidates = [:a]),
        TapResolution(label = :r, selector = "get_face(:a)", reason = :root_input_not_face,
                      field = :a),
        TapResolution(label = :r, selector = "get_face(:q)", reason = :unknown_output_face,
                      field = :q, candidates = [:done]),
        TrimProblemInvalid(field = :guess, reason = :not_a_namedtuple, observed = Int),
        TrimProblemInvalid(field = :lower, reason = :key_set, names = [:a],
                           expected = [:a, :b]),
        TrimProblemInvalid(field = :residuals, reason = :key_set, names = [:r],
                           expected = [:s]),
        TrimProblemInvalid(field = :guess, reason = :field_types,
                           bad = Pair{Symbol,Any}[:a => Int]),
        TrimProblemInvalid(field = :residuals, reason = :field_types,
                           bad = Pair{Symbol,Any}[:r => String]),
        TrimProblemInvalid(field = :lower, reason = :inverted_box, key = :a, value = 2.0,
                           bound = 1.0),
        TrimProblemInvalid(field = :tolerances, reason = :nonpositive_tolerance, key = :r,
                           value = 0.0),
        TrimProblemInvalid(field = :reads, reason = :not_a_read_set, observed = NamedTuple),
        TrimProblemInvalid(field = :problem, reason = :not_a_problem, observed = Int),
        TrimCommitEvents(events = [("a/b", :snap)]),
        TrimCommitResiduals(residuals = [(:r, 1.0, 0.5)]),
        ConditionShapeDrift(reason = :tree_type, compiled = Int, observed = Float64),
        ConditionShapeDrift(reason = :prefix, compiled = "a", observed = "b",
                            position = (:x, 1)),
        ArgumentInvalid(call = :Period, reason = :inexact, value = 0.02),
        ArgumentInvalid(call = :Hz, reason = :inexact, value = 0.5),
        ArgumentInvalid(call = :Absolute, reason = :inexact, value = 0.5),
        ArgumentInvalid(call = :Absolute, reason = :not_a_quantity, value = 1),
        ArgumentInvalid(call = :step!, reason = :both_given),
        ArgumentInvalid(call = :step!, reason = :range, argument = :frames, value = 0),
        ArgumentInvalid(call = :step!, reason = :range, argument = :t_plus, value = -1.0),
        ArgumentInvalid(call = :run!, reason = :no_clock_bound),
        ArgumentInvalid(call = :trim!, reason = :non_nominal, value = "Simulation{Dual}"),
        ArgumentInvalid(call = :selector, reason = :index_not_integer, value = 1.5),
        ArgumentInvalid(call = :TableBinding, reason = :entry_shape, entry = :a),
        ArgumentInvalid(call = :TableBinding, reason = :no_face, entry = :a),
        ArgumentInvalid(call = :TableBinding, reason = :face_name, entry = :a, value = 1),
        ArgumentInvalid(call = :TableBinding, reason = :vocabulary, entry = :a,
                        argument = :gain, vocabulary = [:face, :deadzone, :expo]),
        ArgumentInvalid(call = :TableBinding, reason = :deadzone, entry = :a, value = 1.5),
        ArgumentInvalid(call = :TableBinding, reason = :expo, entry = :a, value = 2.0),
        ReadSetMisuse(observed = Int, reason = :not_a_selector, label = :r),
        ReadSetMisuse(observed = NamedTuple, reason = :not_a_read_set),
        NotAttached(device = "Pad", roster = ["device 1 (Pad)"]),
        # the runtime stream's eight, re-parented (§11.8)
        MalformedDatum(ArgumentError("bad")),
        OutOfClaimEntry(:a, 1.0, [:b], "device 1 (Pad)"),
        ClaimedFaceEntry(:a, "device 1 (Pad)", 1.0, :staging),
        EntryTypeMismatch(:a, "x", Float64),
        ChatteringBudget("a/b", :snap, 1.0, 8, 9),
        FiringBudget("a/b", :snap, 1.0, 4, 5),
        DeviceCrash(ArgumentError("bad"), false),
        DeviceJoinTimeout("device 1 (Pad)", 5.0, 1.0, 10),
    ]

    # Appendix C's severity column, as the list it is: every other kind is an
    # error, so a kind added on the wrong side of the line fails here.
    warning_kinds = Set{DataType}([EmptyGreedyClaim, TrimCommitEvents, TrimCommitResiduals,
                                  MalformedDatum, OutOfClaimEntry, ClaimedFaceEntry,
                                  EntryTypeMismatch, ChatteringBudget, FiringBudget,
                                  DeviceCrash, DeviceJoinTimeout])
    for d in occurrences
        @test severity(d) === (typeof(d) in warning_kinds ? :warning : :error)
        @test path(d) isa String
        m = message(d)
        @test m isa String && !isempty(m)
    end

    # Every kind of the prototype's set has an occurrence above: the coverage
    # check is over `Diagnostic`'s own subtypes, so adding a kind without an
    # occurrence fails here rather than going unrendered.
    @test Set(typeof.(occurrences)) == Set(subtypes(Diagnostic))

end

# --- rendering: the suite's one deliberate exception -----------------------------
# Everywhere else a test matches a diagnostic's *kind and payload*, never its
# text (§13.2). Here, and only here, the rendered string is the claim: the
# carrier's compiler-style layout, and the didactic register `message` is for —
# state the fix, show the list in hand. Nothing outside this testset may
# `occursin` on a rendered diagnostic; a wording change is free everywhere else
# and lands here.
@testset "rendering: the carrier compiler-style, the register didactic (§13.1, §13.2)" begin
    # Two kinds × two paths: groups in first-appearance order, paths sorted
    # within a group, the kind name leading each line, the count line above.
    e = BuildError(Diagnostic[UnconnectedInput(path = "b", face = :u),
                              FaceNameIllegal(path = "b", face = "p/q"),
                              UnconnectedInput(path = "a", face = :v),
                              FaceNameIllegal(path = "a", face = "r/s")])
    @test kinds(e) == [UnconnectedInput, FaceNameIllegal]
    lines = split(sprint(showerror, e), '\n')
    @test lines[1] == "BuildError: 4 diagnostics"
    @test startswith(lines[2], "  UnconnectedInput: `a`.v")
    @test startswith(lines[3], "  UnconnectedInput: `b`.u")
    @test startswith(lines[4], "  FaceNameIllegal: ") && occursin("`r/s`", lines[4])
    @test startswith(lines[5], "  FaceNameIllegal: ") && occursin("`p/q`", lines[5])

    # A fail-fast site's single diagnostic renders on one line, no count.
    @test sprint(showerror, BuildError(UnconnectedInput(path = "a", face = :v))) ==
          "BuildError: UnconnectedInput: " * message(UnconnectedInput(path = "a", face = :v))

    # The did-you-mean list is carried, not ranked (README): the candidates the
    # site had in hand are printed, and no edit distance orders them.
    m = message(UnknownPort(entry = "wires", end_ = :destination, path = "a/b",
                            spelling = "a/b", port = :throtle,
                            candidates = [:throttle, :brake]))
    @test occursin("names no `throtle`", m) && occursin("throttle, brake", m)

    # The remedy register: the shortfall, then the fix, with the list in hand.
    m = message(UninitializedInputs(op = "init!", faces = [:u, :e]))
    @test occursin("`init!`", m) && occursin("`u`, `e`", m)
    @test occursin("nothing was written", m)

    # A `logged` warning renders like a carrier line — the kind name, then the
    # message — because it never gets a `showerror` to lead it (D-214).
    d = TrimCommitResiduals(residuals = [(:torque, 1.0, 0.5)])
    @test logged(d) == "TrimCommitResiduals: " * message(d)
    @test startswith(logged(d), "TrimCommitResiduals: ") && severity(d) === :warning
end
