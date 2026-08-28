# --- the data plane's core exchange (§11, increment 9) ---------------------------
# Staging → drain inbound, snapshot publication outbound: §11.1's planes 1–3 with
# the roster empty, so the harness register's derived surface is every root face.

# Two root inputs feeding one summing junction: the sparse-batch hazard's shape.
two_root_inputs() = Group((; s = Sum(sa = 1.0, sb = 1.0));
                    inputs = ("a" => "s/a", "b" => "s/b"))

# A chain whose published ports are in lockstep at every boundary — g2 computes
# 2·g1 in the same sweep, so any snapshot mixing two boundaries breaks it.
chain3() = Group((; p = Plant(), g1 = Gain(2.0), g2 = Gain(2.0));
                 wires = ("p/y" => "g1/e",
                          "g1/out" => "g2/e"),
                 inputs = ("u" => "p/u",))

@testset "a staged batch lands at its frame top, and nowhere earlier (§11.1, §11.4)" begin
    sim = Simulation(fed(Plant(), "u"); h = 1//10)
    init!(sim, fragment(inputs = (in = 0.0,)))
    step!(sim; t_plus = 0.3)
    stage!(sim, "in" => 1.0)
    @test port(sim, "", :in) === 0.0    # staging never touches a live root input (§11.1)
    step!(sim; t_plus = 0.5)
    @test port(sim, "", :in) === 1.0         # the next frame's drain applied it

    # The frame's outcome is a pure function of the drained batch: the same
    # write applied directly at the same stopped point is the same trajectory,
    # bitwise.
    ref = Simulation(fed(Plant(), "u"); h = 1//10)
    init!(ref, fragment(inputs = (in = 0.0,)))
    step!(ref; t_plus = 0.3)
    poke!(ref, "in", 1.0)                    # the counterfactual, under the data plane
    step!(ref; t_plus = 0.5)
    @test port(sim, "c", :y) === port(ref, "c", :y)
    @test port(sim, "c", :power) === port(ref, "c", :power)
end

@testset "a staged batch waits for the first frame top; one predating init! clears with it (§11.3, §12.6)" begin
    # The contrast with the authored condition: a root-input value carried into
    # `init!` makes the predicate hold in the authored state and fire at t₀
    # (test_events). A batch staged *before* init! predates the boundary zero
    # it would clobber, so init! clears it (§12.6); staged *after* init! — the
    # pre-run register — it is pending, not applied, and the edge arrives with
    # frame 1's drain, never at boundary zero.
    sim = Simulation(fed(Trigger(0.5), "sig"); h = 1//10)
    stage!(sim, "in" => 1.0)                 # predates boundary zero: cleared by init!
    init!(sim, fragment(inputs = (in = 0.0,)))
    @test (@atomic sim.plane.harness.cell.pending) === nothing
    @test modes(sim, "c").count == 0
    stage!(sim, "in" => 1.0)                 # the pre-run register: init! → stage! → run!
    @test modes(sim, "c").count == 0
    run!(sim; t_end = 0.1)
    @test modes(sim, "c").count == 1
end

@testset "coalescing: merge only — newest wins per face, untouched faces survive (§11.4)" begin
    sim = Simulation(two_root_inputs(); h = 1//10)
    init!(sim, fragment(inputs = (a = 0.0, b = 0.0)))
    stage!(sim, "a" => 1.0)
    stage!(sim, "b" => 2.0)                  # sparse: must not clobber the pending `a`
    stage!(sim, "a" => 3.0)                  # re-staged: the newest level, the per-face ZOH
    run!(sim; t_end = 0.1)
    @test port(sim, "", :a) === 3.0
    @test port(sim, "", :b) === 2.0
    @test port(sim, "s", :e) === 5.0
end

@testset "every check runs at staging, on the writer's side; the drain is pure (§11.4)" begin
    sim = Simulation(two_root_inputs(); h = 1//10)
    init!(sim, fragment(inputs = (a = 0.0, b = 0.0)))
    # Each rejection is written into the harness register's diagnostic cell on
    # the staging task (§11.8) — a stopped-sim staging waits there exactly as
    # its surviving entries wait in the staging cell, both drained at the next
    # run's first frame top.
    stage!(sim, "flaps" => 1.0)
    stage!(sim, "a" => "high")
    # A rejected entry is discarded; the rest of its batch stands.
    stage!(sim, "b" => 4.0, "gear" => 1.0)
    run!(sim; t_end = 0.1)
    @test port(sim, "", :a) === 0.0          # nothing surviving ever named `a`
    @test port(sim, "", :b) === 4.0
    hw = writer_status(latest(sim), "harness")
    @test hw.totals.out_of_claim == 2 && hw.totals.type_mismatch == 1
    @test length(hw.recent) == 3             # the one frame's snapshot carries the delta
    ooc = only(d for d in hw.recent if d isa OutOfClaimEntry && d.face === :flaps)
    @test ooc.value == 1.0 && ooc.incumbent === nothing   # no claim anywhere: no such face
    etm = only(d for d in hw.recent if d isa EntryTypeMismatch)
    @test etm.face === :a && etm.value == "high" && etm.declared === Float64
end

@testset "the shim converts to the activation's root-input types (§11.4)" begin
    sim = Simulation(fed(Plant(), "u"), D8; h = 1//10)
    init!(sim, fragment(inputs = (in = 0.0,)))
    stage!(sim, "in" => 1.0)            # convert to the root input's declared type: D8
    run!(sim; t_end = 0.1)
    @test port(sim, "", :in) === D8(1.0)
end

@testset "publication: one immutable value per frame-top boundary (§11.2)" begin
    sim = Simulation(chain3(); h = 1//10)
    @test latest(sim) === nothing            # nothing is published before init!
    init!(sim, fragment(inputs = (u = 1.0,)))
    snap0 = latest(sim)
    @test snap0.t == 0.0 && snap0.frame == 0 # the boundary-zero snapshot (§14.5)

    step!(sim; t_plus = 0.5)
    snap = latest(sim)
    @test snap.frame == 5 && snap.t == sim.exec.clock.t
    # Boundary-consistent and whole-table: every port bitwise the live table's,
    # the root inputs riding along as the source cells they are (§11.2).
    for (path, name) in (("p", :y), ("p", :power),
                         ("g1", :out), ("g2", :out), ("", :u))
        @test port(snap, path, name) === port(sim, path, name)
    end

    # The binding rule: the run moves on, the published snapshot holds.
    y5 = port(snap, "p", :y)
    step!(sim; t_plus = 0.5)
    @test port(snap, "p", :y) === y5
    @test latest(sim).frame == 10

    # Every frame top publishes, the off-tick boundary included.
    simo = Simulation(chain3(); h = 1//20, n = 2)
    init!(simo, fragment(inputs = (u = 0.0,)))
    run!(simo; t_end = 0.05)                         # one frame, not a base tick
    @test latest(simo).frame == 1
end

@testset "the exchange is wait-free and coherent: no reader ever sees a torn world (§11.2)" begin
    sim = Simulation(chain3(); h = 1//1000)
    init!(sim, fragment(inputs = (u = 1.0,)))
    stop = Threads.Atomic{Bool}(false)
    reader = Threads.@spawn begin
        seen, bad, tprev, mono = 0, 0, -1.0, true
        while true
            s = latest(sim)
            if s !== nothing
                seen += 1
                # In-lockstep within one snapshot: g2 was computed as 2·g1 in
                # the same sweep, so any mix of two boundaries breaks it.
                port(s, "g2", :out) === 2.0 * port(s, "g1", :out) || (bad += 1)
                s.t ≥ tprev || (mono = false)
                tprev = s.t
            end
            stop[] && break                  # one final sample past the stop, so seen ≥ 1
            yield()
        end
        (seen, bad, mono)
    end
    run!(sim; t_end = 2.0)
    stop[] = true
    (seen, bad, mono) = fetch(reader)
    @test seen > 0 && bad == 0 && mono
    @test latest(sim).t == sim.exec.clock.t
end

@testset "staging from another task: the CAS merge loses nothing it shouldn't (§11.4)" begin
    sim = Simulation(two_root_inputs(); h = 1//10)
    init!(sim, fragment(inputs = (a = 0.0, b = 0.0)))
    writer = Threads.@spawn for i in 1:1000
        stage!(sim, "a" => Float64(i))
        yield()
    end
    while !istaskdone(writer)                # drains race the stager, mid-session
        step!(sim; frames = 1)
        yield()
    end
    wait(writer)
    step!(sim; frames = 1)                   # one more frame drains whatever still pends
    @test port(sim, "", :a) === 1000.0       # the last staged level is what stands
    @test port(sim, "", :b) === 0.0          # never staged, never touched
end

@testset "an empty drain is free: the frame top adds no work to a quiet loop (§11.1)" begin
    sim = Simulation(chain3(); h = 1//10)
    init!(sim, fragment(inputs = (u = 0.0,)))
    @test @ballocated(drain!($sim)) == 0
end

@testset "a populated drain is as free as an empty one, whatever the batch touches (§11.4, D-202)" begin
    sim = Simulation(two_root_inputs(); h = 1//10)
    init!(sim, fragment(inputs = (a = 0.0, b = 0.0)))
    stage!(sim, "a" => 1.0); drain!(sim)             # warm the writer's one scatter
    @test @ballocated(drain!($sim), setup = (stage!($sim, "a" => 1.0)), evals = 1) == 0
    # A never-drained sparsity pattern costs the same nothing: the scatter is
    # one specialization per writer, never one per touched-face combination.
    @test @ballocated(drain!($sim), setup = (stage!($sim, "b" => 1.0)), evals = 1) == 0
    @test @ballocated(drain!($sim), setup = (stage!($sim, "a" => 1.0, "b" => 2.0)),
                      evals = 1) == 0
end

# A surface past the 32-wide threshold where Base's tuple `map` leaves its
# inlined small-tuple path: the merge and the scatter must lean on neither.
wide_root_inputs(n) = Group(NamedTuple(Symbol(:s, i) => Sum(sa = 1.0, sb = 1.0) for i in 1:n);
                      inputs = Tuple(vcat(["a$i" => "s$i/a" for i in 1:n],
                                          ["b$i" => "s$i/b" for i in 1:n])))

# Its baseline (§14.6): a generated fixture's full-coverage condition, generated
# beside it. Totality is a precondition of `init!`, and 34 hand-written root-input
# values is exactly the case baselines exist for.
wide_zero(n) = fragment(inputs = NamedTuple(Symbol(p, i) => 0.0
                                           for p in ("a", "b") for i in 1:n))

@testset "a wide surface stages, merges and drains like a narrow one (§11.4, D-202)" begin
    sim = Simulation(wide_root_inputs(17); h = 1//10)      # 34 root faces
    init!(sim, wide_zero(17))
    stage!(sim, "a3" => 1.5)
    stage!(sim, "b9" => -2.0, "a3" => 2.5)           # merge: newest wins, untouched survive
    drain!(sim)
    @test port(sim, "", :a3) === 2.5
    @test port(sim, "", :b9) === -2.0
    @test port(sim, "", :a1) === 0.0                 # never staged, never touched
    dense = () -> stage!(sim, ("a$i" => Float64(i) for i in 1:17)...,
                              ("b$i" => -Float64(i) for i in 1:17)...)
    dense(); drain!(sim)                             # warm the wide scatter
    @test port(sim, "", :b17) === -17.0
    @test @ballocated(drain!($sim), setup = (stage!($sim, "a7" => 1.0)), evals = 1) == 0
    @test @ballocated(drain!($sim), setup = ($dense()), evals = 1) == 0
end
