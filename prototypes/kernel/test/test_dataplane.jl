# --- the data plane's core exchange (§11, increment 9) ---------------------------
# Staging → drain inbound, snapshot publication outbound: §11.1's planes 1–3 with
# the roster empty, so the harness register's derived surface is every root face.

# Two root slots feeding one summing junction: the sparse-batch hazard's shape.
two_slots() = Group((; s = Sum(sa = 1.0, sb = 1.0)), (),
                    ("a" => "children/s/a", "b" => "children/s/b"), ())

# A chain whose published ports are in lockstep at every boundary — g2 computes
# 2·g1 in the same sweep, so any snapshot mixing two boundaries breaks it.
chain3() = Group((; p = Plant(), g1 = Gain(2.0), g2 = Gain(2.0)),
                 ("children/p/y" => "children/g1/e",
                  "children/g1/out" => "children/g2/e"),
                 ("u" => "children/p/u",), ())

@testset "a staged batch lands at its frame top, and nowhere earlier (§11.1, §11.4)" begin
    sim = Simulation(fed(Plant(), "u"); h = 1//10)
    init!(sim)
    run!(sim, 0.3)
    stage!(sim, "in" => 1.0)
    @test port(sim, "", :in) === 0.0         # staging never touches a live slot (§11.1)
    run!(sim, 0.8)
    @test port(sim, "", :in) === 1.0         # the next frame's drain applied it

    # The frame's outcome is a pure function of the drained batch: the same
    # write applied directly at the same stopped point is the same trajectory,
    # bitwise.
    ref = Simulation(fed(Plant(), "u"); h = 1//10)
    init!(ref)
    run!(ref, 0.3)
    set_slot!(ref, "in", 1.0)
    run!(ref, 0.8)
    @test port(sim, "children/c", :y) === port(ref, "children/c", :y)
    @test port(sim, "children/c", :power) === port(ref, "children/c", :power)
end

@testset "a batch staged while stopped waits for the first frame top, not boundary zero (§11.3)" begin
    # The contrast with the direct write: set_slot! before init! makes the
    # predicate hold in the authored state and fire at t₀ (test_events); a
    # *staged* batch is pending, not applied, so boundary zero runs on the
    # un-drained slot and the edge arrives with frame 1's drain.
    sim = Simulation(fed(Trigger(0.5), "sig"); h = 1//10)
    stage!(sim, "in" => 1.0)
    init!(sim)
    @test modes(sim, "children/c").count == 0
    run!(sim, 0.1)
    @test modes(sim, "children/c").count == 1
end

@testset "coalescing: merge only — newest wins per face, untouched faces survive (§11.4)" begin
    sim = Simulation(two_slots(); h = 1//10)
    init!(sim)
    stage!(sim, "a" => 1.0)
    stage!(sim, "b" => 2.0)                  # sparse: must not clobber the pending `a`
    stage!(sim, "a" => 3.0)                  # re-staged: the newest level, the per-face ZOH
    run!(sim, 0.1)
    @test port(sim, "", :a) === 3.0
    @test port(sim, "", :b) === 2.0
    @test port(sim, "children/s", :e) === 5.0
end

@testset "every check runs at staging, on the writer's side; the drain is pure (§11.4)" begin
    sim = Simulation(two_slots(); h = 1//10)
    init!(sim)
    # Each rejection is written into the harness register's diagnostic cell on
    # the staging task (§11.8) — a stopped-sim staging waits there exactly as
    # its surviving entries wait in the staging cell, both drained at the next
    # run's first frame top.
    stage!(sim, "flaps" => 1.0)
    stage!(sim, "a" => "high")
    # A rejected entry is discarded; the rest of its batch stands.
    stage!(sim, "b" => 4.0, "gear" => 1.0)
    run!(sim, 0.1)
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

@testset "the shim converts to the activation's slot types (§11.4)" begin
    sim = Simulation(fed(Plant(), "u"), D8; h = 1//10)
    init!(sim)
    stage!(sim, "in" => 1.0)                 # convert to the slot's declared type: D8
    run!(sim, 0.1)
    @test port(sim, "", :in) === D8(1.0)
end

@testset "publication: one immutable value per frame-top boundary (§11.2)" begin
    sim = Simulation(chain3(); h = 1//10)
    @test latest(sim) === nothing            # nothing is published before init!
    init!(sim)
    snap0 = latest(sim)
    @test snap0.t == 0.0 && snap0.frame == 0 # the boundary-zero snapshot (§14.5)

    set_slot!(sim, "u", 1.0)
    run!(sim, 0.5)
    snap = latest(sim)
    @test snap.frame == 5 && snap.t == sim.clock.t
    # Boundary-consistent and whole-table: every port bitwise the live table's,
    # the root slots riding along as the source cells they are (§11.2).
    for (path, name) in (("children/p", :y), ("children/p", :power),
                         ("children/g1", :out), ("children/g2", :out), ("", :u))
        @test port(snap, path, name) === port(sim, path, name)
    end

    # The binding rule: the run moves on, the published snapshot holds.
    y5 = port(snap, "children/p", :y)
    run!(sim, 1.0)
    @test port(snap, "children/p", :y) === y5
    @test latest(sim).frame == 10

    # Every frame top publishes, the off-tick boundary included.
    simo = Simulation(chain3(); h = 1//20, n = 2)
    init!(simo)
    run!(simo, 0.05)                         # one frame, not a base tick
    @test latest(simo).frame == 1
end

@testset "the exchange is wait-free and coherent: no reader ever sees a torn world (§11.2)" begin
    sim = Simulation(chain3(); h = 1//1000)
    set_slot!(sim, "u", 1.0)
    init!(sim)
    stop = Threads.Atomic{Bool}(false)
    reader = Threads.@spawn begin
        seen, bad, tprev, mono = 0, 0, -1.0, true
        while true
            s = latest(sim)
            if s !== nothing
                seen += 1
                # In-lockstep within one snapshot: g2 was computed as 2·g1 in
                # the same sweep, so any mix of two boundaries breaks it.
                port(s, "children/g2", :out) === 2.0 * port(s, "children/g1", :out) || (bad += 1)
                s.t ≥ tprev || (mono = false)
                tprev = s.t
            end
            stop[] && break                  # one final sample past the stop, so seen ≥ 1
            yield()
        end
        (seen, bad, mono)
    end
    run!(sim, 2.0)
    stop[] = true
    (seen, bad, mono) = fetch(reader)
    @test seen > 0 && bad == 0 && mono
    @test latest(sim).t == sim.clock.t
end

@testset "staging from another task: the CAS merge loses nothing it shouldn't (§11.4)" begin
    sim = Simulation(two_slots(); h = 1//10)
    init!(sim)
    writer = Threads.@spawn for i in 1:1000
        stage!(sim, "a" => Float64(i))
        yield()
    end
    t = 0.0
    while !istaskdone(writer)                # drains race the stager, mid-run
        t += 0.1
        run!(sim, t)
        yield()
    end
    wait(writer)
    run!(sim, t + 0.1)                       # one more frame drains whatever still pends
    @test port(sim, "", :a) === 1000.0       # the last staged level is what stands
    @test port(sim, "", :b) === 0.0          # never staged, never touched
end

@testset "an empty drain is free: the frame top adds no work to a quiet loop (§11.1)" begin
    sim = Simulation(chain3(); h = 1//10)
    init!(sim)
    @test @ballocated(drain!($sim)) == 0
end
