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
    sim = Simulation(two_slots(); h = 1//10)
    dev = Parser(0.7, "garbage", 0.9)
    hp = attach!(sim, dev, Enumerated("a"))
    init!(sim)
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
    sim = Simulation(two_slots(); h = 1//10)
    h = attach!(sim, Pad("p"), Enumerated("a"))
    init!(sim)
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
    sim = Simulation(two_slots(); h = 1//10)
    h = attach!(sim, Pad("p"), Enumerated("a"))
    init!(sim)
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
    sim = Simulation(two_slots(); h = 1//10)
    dev = Ticker()
    attach!(sim, dev, Enumerated())
    init!(sim)
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
    sim = Simulation(two_slots(); h = 1//10)
    dev = LateReporter()
    attach!(sim, dev, Enumerated())
    init!(sim)
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
    init!(sim)
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
    _report!(cell, ChatteringBudget("children/c", :pop, 0.1, 8, 8))
    _report!(cell, FiringBudget("children/e", :up, 0.0, 4, 4))
    _report!(cell, DeviceCrash(ErrorException("boom"), false))
    batch = _take!(cell)
    @test length(batch.ring) == 7
    @test batch.ring[1] isa MalformedDatum && batch.ring[7] isa DeviceCrash
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
    sim = Simulation(two_slots(); h = 1//10)
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
