# --- the binding's runtime half (§11.6, §11.4; increment 13) --------------------
# The shipped TableBinding with its generic map_input — the shared pure
# conditioning helper, with an owner — and the handle's `binding` capability.
# The devices below live at top level for the README's local-scope reason.

# A fed plant exporting its output: the root face set the output side reads —
# one root input (`u`), one exported output face (`y`).
outfaced() = Group((; p = Plant()); inputs = ("u" => "p/u",), outputs = ("p/y" => "y",))

# A boundary-driven output device (§11.2): the loop idiom verbatim — wait,
# gather against the handle's compiled reads, map_output through its binding —
# with the wire a vector standing in for the socket.
mutable struct Telemetry <: AbstractDevice
    wire::Vector{Any}
end
Telemetry() = Telemetry(Any[])
function loop(d::Telemetry, h)
    while true
        snap = wait_next_snapshot(h)
        snap === nothing || push!(d.wire, map_output(gather(h, snap), binding(h)))
        running(h) || break
    end
    nothing
end

# The output side's conformance customers (§11.6), at top level (README).
struct NoReads <: AbstractBinding end            # declared side, no enumeration
is_output(::NoReads) = true
struct ReadsUndeclared <: AbstractBinding end    # enumeration written, trait false
is_input(::ReadsUndeclared) = true
claims(::ReadsUndeclared) = ("a",)
reads(::ReadsUndeclared) = (; x = get_input("a"))
struct BadReadsShape <: AbstractBinding end      # the wrong return shape
is_output(::BadReadsShape) = true
reads(::BadReadsShape) = [get_input("a")]
struct BadReadsEntry <: AbstractBinding end      # a non-selector entry
is_output(::BadReadsEntry) = true
reads(::BadReadsEntry) = (; x = 1.0)
struct Duplex <: AbstractBinding end             # both sides, one binding
is_input(::Duplex) = true
is_output(::Duplex) = true
claims(::Duplex) = ("a",)
reads(::Duplex) = (; echo = get_input("a"), e = get_output("s", "e"))

# A poll-once input device: assembles one datum, stages it through the loop
# idiom — map_input against the handle's own binding — then holds its loop
# until a snapshot shows the write applied, and requests the stop: the full
# poll → stage → observe cycle, deterministic by construction.
mutable struct Poller <: AbstractDevice
    datum::NamedTuple
    fired::Bool
end
Poller(datum) = Poller(datum, false)
function loop(d::Poller, h)
    d.fired && return nothing
    d.fired = true
    pairs = map_input(d.datum, binding(h))
    stage!(h, pairs...)
    (f, v) = first(pairs)
    while running(h)
        snap = wait_next_snapshot(h)
        port(snap, "", Symbol(f)) === v && (stop!(h); break)
    end
    nothing
end

@testset "TableBinding construction validates the table's shape (§11.6)" begin
    err = failure(() -> TableBinding(stick = (deadzone = 0.1,)))
    diag = only(err.diagnostics)
    @test err isa BuildError && diag isa ArgumentInvalid &&
          diag.call === :TableBinding && diag.reason === :no_face && diag.entry === :stick
    err = failure(() -> TableBinding(stick = "elevator"))
    diag = only(err.diagnostics)
    @test err isa BuildError && diag isa ArgumentInvalid &&
          diag.reason === :entry_shape && diag.entry === :stick
    err = failure(() -> TableBinding(stick = (face = "a", deadzon = 0.1)))
    diag = only(err.diagnostics)     # the typo, by name
    @test err isa BuildError && diag isa ArgumentInvalid &&
          diag.reason === :vocabulary && diag.entry === :stick && diag.argument === :deadzon
    err = failure(() -> TableBinding(stick = (face = "a", deadzone = 1.0)))
    diag = only(err.diagnostics)
    @test err isa BuildError && diag isa ArgumentInvalid &&
          diag.reason === :deadzone && diag.entry === :stick && diag.value == 1.0
    err = failure(() -> TableBinding(stick = (face = "a", expo = 1.5)))
    diag = only(err.diagnostics)
    @test err isa BuildError && diag isa ArgumentInvalid &&
          diag.reason === :expo && diag.entry === :stick && diag.value == 1.5
end

@testset "the input side is declared and the claim is the table's face set (§11.6)" begin
    b = TableBinding(stick = (face = "a", deadzone = 0.1, expo = 0.5),
                     thr   = (face = "b",))
    @test is_input(b) && !is_output(b) && !is_greedy(b)
    @test sort(claims(b)) == ["a", "b"]
    # Two channels onto one face collapse to one claim, like any enumeration's
    # duplicates (§11.3).
    sim = Simulation(two_root_inputs(); h = 1//10)
    attach!(sim, Pad("p"), TableBinding(x = (face = "a",), y = (face = "a",)))
    @test sim.plane.roster[1].writer.faces == [:a]
end

@testset "the conditioning: deadzone rescales, expo attenuates, endpoints fixed (§11.4)" begin
    b = TableBinding(stick = (face = "a", deadzone = 0.1, expo = 0.5),
                     dzo   = (face = "a", deadzone = 0.1),
                     exo   = (face = "a", expo = 0.5),
                     thr   = (face = "b",))
    val(datum) = last(only(map_input(datum, b)))
    # Inside the band: exactly zero, both polarities.
    @test val((; stick = 0.05)) == 0.0
    @test val((; stick = -0.05)) == 0.0
    # The endpoints stay fixed under both parameters, and beyond-range input
    # clamps to them (the axis convention binds where conditioning is declared).
    @test val((; stick = 1.0)) ≈ 1.0
    @test val((; stick = -1.0)) ≈ -1.0
    @test val((; stick = 1.7)) ≈ 1.0
    # Deadzone alone: the band edge is continuous and the remainder rescales.
    @test val((; dzo = 0.1)) ≈ 0.0 atol = 1e-15
    @test val((; dzo = 0.55)) ≈ 0.5
    # Expo alone: the published blend a = (1-e)·a + e·a³, midrange attenuated,
    # odd-symmetric.
    @test val((; exo = 0.5)) ≈ 0.5 * 0.5 + 0.5 * 0.5^3
    @test val((; exo = -0.5)) ≈ -(0.5 * 0.5 + 0.5 * 0.5^3)
    @test abs(val((; exo = 0.5))) < 0.5
    # Monotone through the composition.
    levels = [val((; stick = x)) for x in -1.0:0.05:1.0]
    @test issorted(levels)
    # An entry declaring neither passes through untouched — a throttle level
    # outside [-1, 1] semantics, or a press counter (the levels doctrine).
    @test val((; thr = 0.7)) === 0.7
    @test val((; thr = 17)) === 17
end

@testset "map_input is sparse over the datum, and an unknown channel is drift (§11.4, §11.6)" begin
    b = TableBinding(stick = (face = "a",), thr = (face = "b",))
    # Face ⇒ value pairs for exactly what the datum touched: a sparse datum
    # stages a sparse batch, and merge does the rest.
    @test map_input((; thr = 0.7), b) == ("b" => 0.7,)
    @test length(map_input((; stick = 0.2, thr = 0.7), b)) == 2
    # A datum field naming no channel is configuration drift, not a bad datum:
    # the throw is deliberate, and in a loop body it lands as DeviceCrash.
    err = failure(() -> map_input((; wheel = 0.1), b))
    @test err isa ErrorException && occursin("`wheel`", err.msg)
end

@testset "the loop idiom end to end: poll → map_input(binding(handle)) → stage! (§11.6)" begin
    sim = Simulation(two_root_inputs(); h = 1//10)
    dev = Poller((; stick = 0.55, thr = 0.7))
    h = attach!(sim, dev, TableBinding(stick = (face = "a", deadzone = 0.1),
                                       thr   = (face = "b",)))
    @test binding(h) === sim.plane.roster[1].binding
    init!(sim, fragment(inputs = (a = 0.0, b = 0.0)))
    run!(sim; t_end = 1000.0)                # ends by the device's stop, past its observed apply
    @test port(sim, "", :a) ≈ 0.5            # (0.55 − 0.1) / 0.9: conditioned at staging
    @test port(sim, "", :b) === 0.7          # pass-through, bitwise
    # The device-staged trajectory is the directly-staged one: conditioning ran
    # upstream, so the model consumed post-conditioning levels (§11.4).
    ref = Simulation(two_root_inputs(); h = 1//10)
    init!(ref, fragment(inputs = (a = 0.0, b = 0.0)))
    ref_val = last(only(map_input((; stick = 0.55),
                                  TableBinding(stick = (face = "a", deadzone = 0.1)))))
    stage!(ref, "a" => ref_val, "b" => 0.7)
    run!(ref; t_end = sim.exec.clock.step * sim.h)
    @test port(sim, "s", :e) === port(ref, "s", :e)
    # An unknown channel in a real loop body crashes the device by name, the
    # run continuing (§11.6: any non-datum exception propagates to the wrapper).
    sim2 = Simulation(two_root_inputs(); h = 1//10)
    attach!(sim2, Poller((; wheel = 0.1)), TableBinding(stick = (face = "a",)))
    init!(sim2, fragment(inputs = (a = 0.0, b = 0.0)))
    logs, _ = Test.collect_test_logs() do
        run!(sim2; t_end = 0.2)
    end
    @test sim2.exec.clock.step == 2               # the crash is the device's alone
    @test crash_accounted(sim2, logs, "device 1 (Poller)")
end

@testset "the output side completes the conformance check, both directions (§11.6)" begin
    sim = Simulation(two_root_inputs(); h = 1//10)
    err = failure(() -> attach!(sim, Pad("p"), NoReads()))
    diag = only(err.diagnostics)
    @test err isa BuildError && diag isa BindingContractMismatch && diag.reason === :reads_missing
    err = failure(() -> attach!(sim, Pad("p"), ReadsUndeclared()))
    diag = only(err.diagnostics)
    @test err isa BuildError && diag isa BindingContractMismatch &&
          diag.reason === :reads_without_output
    err = failure(() -> attach!(sim, Pad("p"), BadReadsShape()))
    diag = only(err.diagnostics)
    @test err isa BuildError && diag isa BindingContractMismatch &&
          diag.reason === :reads_not_namedtuple
    err = failure(() -> attach!(sim, Pad("p"), BadReadsEntry()))
    diag = only(err.diagnostics)
    @test err isa BuildError && diag isa BindingContractMismatch &&
          diag.reason === :reads_not_selectors
    @test isempty(sim.plane.roster)              # every rejection left the roster untouched
end

@testset "reads resolve at attach: binding drift fails there, never on the wire (§11.2, §14.4)" begin
    sim = Simulation(outfaced(); h = 1//10)
    err = failure(() -> attach!(sim, Pad("t"), Readout(alt = get_output("q", "y"))))
    diag = only(err.diagnostics)
    @test err isa BuildError && diag isa ReadBindingUnresolved && diag.reason === :unknown_cell &&
          diag.selector == "get_output(\"q\", :y)"
    err = failure(() -> attach!(sim, Pad("t"), Readout(v = get_input("nope"))))
    diag = only(err.diagnostics)
    @test err isa BuildError && diag isa ReadBindingUnresolved &&
          diag.reason === :unknown_root_input && diag.candidates == [:u]  # the root-input list, in hand
    err = failure(() -> attach!(sim, Pad("t"), Readout(v = get_face("u"))))
    diag = only(err.diagnostics)
    @test err isa BuildError && diag isa ReadBindingUnresolved && diag.reason === :root_input_not_output
    err = failure(() -> attach!(sim, Pad("t"), Readout(v = get_face("nope"))))
    diag = only(err.diagnostics)
    @test err isa BuildError && diag isa ReadBindingUnresolved && diag.reason === :unknown_output_face
    # A rejected attach consumed no id, and the good one lands as device 1.
    h = attach!(sim, Pad("t"), Readout(alt = get_face("y")))
    @test sim.plane.roster[1].id == 1
    # An output-only binding stakes no claim: the harness keeps every face.
    @test isempty(sim.plane.claimedby)
end

@testset "the compiled gather: wait → gather → map_output on the device task (§11.2, §12.3)" begin
    sim = Simulation(outfaced(); h = 1//10)
    dev = Telemetry()
    h = attach!(sim, dev, Readout(alt = get_face("y"),
                                  raw = get_output("p", "y"),
                                  cmd = get_input("u")))
    init!(sim, fragment(inputs = (u = 0.0,)))
    stage!(sim, "u" => 2.0)
    run!(sim; t_end = 0.5)
    @test !isempty(dev.wire)
    # map_output received the labeled NamedTuple, labels in reads order.
    @test all(keys(nt) == (:alt, :raw, :cmd) for nt in dev.wire)
    # Coherence in every observed snapshot: the exported face and the deep path
    # alias one cell (§11.2), so the two registers agree bitwise, always.
    @test all(nt.alt === nt.raw for nt in dev.wire)
    # The stop wake handed the final world: the last datum is the run's end.
    @test last(dev.wire).raw === port(sim, "p", :y)
    @test last(dev.wire).cmd === 2.0
    # The same compiled gather serves the calling task against any snapshot.
    nt = gather(h, latest(sim))
    @test nt === last(dev.wire)
end

@testset "gather without an output side is a contract misuse, by name (§11.6)" begin
    sim = Simulation(two_root_inputs(); h = 1//10)
    h = attach!(sim, Pad("p"), Enumerated("a"))
    init!(sim, fragment(inputs = (a = 0.0, b = 0.0)))
    err = failure(() -> gather(h, latest(sim)))
    diag = only(err.diagnostics)
    @test err isa BuildError && diag isa DeviceContractMismatch &&
          diag.reason === :no_output_side
end

@testset "a bidirectional binding composes both halves (§11.6)" begin
    sim = Simulation(two_root_inputs(); h = 1//10)
    h = attach!(sim, Pad("p"), Duplex())
    @test sim.plane.roster[1].writer.faces == [:a]     # the input half: the claim staked
    init!(sim, fragment(inputs = (a = 0.0, b = 0.0)))
    stage!(h, "a" => 0.4)
    run!(sim; t_end = 0.2)
    nt = gather(h, latest(sim))                        # the output half: the gather compiled
    @test nt.echo === 0.4                        # the root input read back through get_input
    @test nt.e === port(sim, "s", :e)
end
