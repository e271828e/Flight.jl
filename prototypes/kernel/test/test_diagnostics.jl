# --- the bad-datum channel and the diagnostic cell (§11.6, §11.8; increment 13) --
# report! from the author's loop body, the per-device single-writer cell with
# its ring-plus-counts bound, the frame-top sentinel-swap drain beside the
# staging drain, and the per-run totals behind diagnostics(sim). The devices
# below live at top level for the README's local-scope reason.

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

# A late reporter: its one report lands after the last frame top the run
# drains, so only the run-end drain can account for it.
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

@testset "a bad datum is tolerated: catch, stage nothing, report, continue (§11.6)" begin
    sim = Simulation(two_slots(); h = 1//10)
    dev = Parser(0.7, "garbage", 0.9)
    attach!(sim, dev, Enumerated("a"))
    init!(sim)
    logs, _ = Test.collect_test_logs() do
        run!(sim, 1000.0)                        # ends by the device's stop
    end
    msgs = [string(l.message) for l in logs]
    # The link survived its truncated datagram: no crash, and the report
    # surfaced device-attributed, carrying the author's cause.
    @test !any(occursin("DeviceCrash", m) for m in msgs)
    @test count(occursin("MalformedDatum from device 1 (Parser)", m) for m in msgs) == 1
    @test any(occursin("unparseable", m) for m in msgs)
    @test diagnostics(sim) == ["device 1 (Parser)" => 1]
    # The stream's good datums reached the slot (newest wins within the staged
    # batch, applied at the next run's first drain).
    run!(sim, (sim.clock.step + 1) * sim.h)
    @test port(sim, "", :a) === 0.9
end

@testset "the ring's bound is the rate limit: 16 retained, excess to the count (§11.8)" begin
    sim = Simulation(two_slots(); h = 1//10)
    h = attach!(sim, Pad("p"), Enumerated("a"))
    init!(sim)
    logs, _ = Test.collect_test_logs() do
        for k in 1:20                            # one frame's flood, staged while pending
            report!(h, MalformedDatum("datum $k"))
        end
        run!(sim, 0.1)                           # the first frame top drains the cell
    end
    msgs = [string(l.message) for l in logs]
    reports = [m for m in msgs if occursin("MalformedDatum", m)]
    @test length(reports) == DIAG_RING + 1       # 16 retained + the suppression summary
    # Earliest-in-frame retained: the first occurrences carry the diagnostic
    # content, the excess becomes exactly a count.
    @test any(occursin("datum 16", m) for m in reports)
    @test !any(occursin("datum 17", m) for m in reports)
    @test any(occursin("4 more occurrence", m) for m in reports)
    # Nothing is lost by not looking: the totals carry the full account.
    @test diagnostics(sim) == ["device 1 (Pad)" => 20]
end

@testset "totals are per run, and the run's end drains what remains (§11.8, §12.4)" begin
    sim = Simulation(two_slots(); h = 1//10)
    dev = LateReporter()
    attach!(sim, dev, Enumerated())
    init!(sim)
    # The report lands between the run's last frame top and its stop: only the
    # run-end drain can surface it within this run — and it does, the warning
    # arriving inside run! and the totals counting it.
    @test_logs (:warn, r"MalformedDatum from device 1 \(LateReporter\): late") #=
        =# match_mode=:any run!(sim, 1000.0)
    @test diagnostics(sim) == ["device 1 (LateReporter)" => 1]
    # A fresh run starts a fresh account: totals count since the run began.
    run!(sim, (sim.clock.step + 2) * sim.h)
    @test diagnostics(sim) == ["device 1 (LateReporter)" => 0]
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
