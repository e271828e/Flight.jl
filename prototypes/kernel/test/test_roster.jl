# --- the roster and claims (§11.3, §11.4, §11.6's attach-point slice) -----------
# Increment 10: attach/detach admission, both claim sources, per-device staging
# cells, and the harness register as the derived remainder.
#
# The malformed bindings live at top level for the README's local-scope reason:
# a trait method defined inside a @testset binds a new local function, and the
# conformance check would see a binding declaring nothing at all.

struct NoSides <: AbstractBinding end                # neither side declared

struct NoEnum <: AbstractBinding end                 # input declared, enumeration missing
is_input(::NoEnum) = true

struct GreedyPlus <: AbstractBinding end             # both claim sources at once
is_input(::GreedyPlus) = true
is_greedy(::GreedyPlus) = true
claims(::GreedyPlus) = ("a",)

struct Sourceless <: AbstractBinding end             # a source without its side
is_greedy(::Sourceless) = true

struct Drifted <: AbstractBinding end                # `claims` under a false `is_input`
is_output(::Drifted) = true
claims(::Drifted) = ("a",)

struct Reader <: AbstractBinding end                 # output side only: absent here
is_output(::Reader) = true

@testset "the binding conformance check names every drift at the attach point (§11.6)" begin
    sim = Simulation(two_slots(); h = 1//10)
    d = Pad("d")
    for (b, frag) in ((NoSides(), "neither side"),
                      (NoEnum(), "no `claims` enumeration"),
                      (GreedyPlus(), "alternatives, not layers"),
                      (Sourceless(), "without `is_input`"),
                      (Drifted(), "while `is_input` reads false"))
        err = failure(() -> attach!(sim, d, b))
        @test err isa BuildError && occursin("BindingContractMismatch", err.msg)
        @test occursin(frag, err.msg)
    end
    # The output side is an absence, not a conformance drift, and it is named
    # *after* the conformance clauses — which is why Drifted above reported its
    # drift rather than falling through to this.
    err = failure(() -> attach!(sim, d, Reader()))
    @test err isa BuildError && occursin("output side", err.msg)
    @test isempty(sim.plane.roster)                  # none of the six was rostered
end

@testset "admission is three checks in spec order (§11.3)" begin
    sim = Simulation(two_slots(); h = 1//10)
    d1 = Pad("d1")
    @test attach!(sim, d1, Enumerated("a")).id == 1
    # Identity before claims: the same instance re-attached — even under an
    # overlapping claim — is AlreadyAttached, never a self-ClaimConflict.
    err = failure(() -> attach!(sim, d1, Enumerated("a")))
    @test err isa BuildError && occursin("AlreadyAttached", err.msg)
    # Claims: face exclusivity, always two *distinct* devices named.
    err = failure(() -> attach!(sim, Pad("d2"), Enumerated("b", "a")))
    @test err isa BuildError && occursin("ClaimConflict", err.msg)
    @test occursin("device 1", err.msg)
    # Affinity: the calling task is a single-slot resource.
    @test attach!(sim, Panel("p1"), Enumerated("b")).id == 2
    err = failure(() -> attach!(sim, Panel("p2"), Enumerated()))
    @test err isa BuildError && occursin("CallerTaskConflict", err.msg)
    # An enumeration drifted onto a nonexistent face is a diagnosable anomaly.
    err = failure(() -> attach!(sim, Pad("d3"), Enumerated("flaps")))
    @test err isa BuildError && occursin("AttachUnknownFace", err.msg)
    # Detaching what was never rostered is an error, not a silent no-op.
    err = failure(() -> detach!(sim, Pad("ghost")))
    @test err isa BuildError && occursin("not rostered", err.msg)
end

@testset "a device writes inside its claim, every check at its own staging (§11.3, §11.4)" begin
    sim = Simulation(two_slots(); h = 1//10)
    da, db = Pad("da"), Pad("db")
    ha = attach!(sim, da, Enumerated("a"))           # the handle is the write capability (§11.6)
    hb = attach!(sim, db, Enumerated("b"))
    init!(sim)
    stage!(ha, "a" => 1.0)
    stage!(hb, "b" => 2)                             # the shim converts to the slot's Float64
    @test port(sim, "", :a) === 0.0                  # staged is pending, never applied (§11.1)
    run!(sim, 0.1)
    @test port(sim, "", :a) === 1.0
    @test port(sim, "", :b) === 2.0
    # Out-of-claim is always OutOfClaimEntry for a device — naming the incumbent
    # when the face is claimed elsewhere — and the rest of the batch stands.
    @test_logs (:warn, r"OutOfClaimEntry: `b` is outside device 1 \(Pad\)'s claim.*claimed by device 2") #=
        =# stage!(ha, "b" => 9.0, "a" => 3.0)
    @test_logs (:warn, r"OutOfClaimEntry: `flaps`") stage!(ha, "flaps" => 1.0)
    run!(sim, 0.2)
    @test port(sim, "", :a) === 3.0
    @test port(sim, "", :b) === 2.0
    # The empty enumeration: an honest may-write-nothing degenerate (§11.6).
    hc = attach!(sim, Pad("dc"), Enumerated())
    @test_logs (:warn, r"OutOfClaimEntry") stage!(hc, "a" => 9.0)
end

@testset "the computed claim is the complement at the attach instant (§11.3, §11.6)" begin
    sim = Simulation(two_slots(); h = 1//10)
    d1, g = Pad("d1"), Pad("gui")
    attach!(sim, d1, Enumerated("a"))
    hg = attach!(sim, g, Greedy())                   # greedy last: exactly what is left
    @test sim.plane.roster[2].writer.faces == [:b]
    init!(sim)
    stage!(hg, "b" => 5.0)
    run!(sim, 0.1)
    @test port(sim, "", :b) === 5.0
    # Past the attach point nothing downstream tells the sources apart.
    @test_logs (:warn, r"OutOfClaimEntry") stage!(hg, "a" => 9.0)

    # A rostered greedy claimant empties the harness surface: every harness
    # stage! in such a session is rejected by name (D-192).
    @test isempty(sim.plane.harness.faces)
    @test_logs (:warn, r"ClaimedFaceEntry: `b` is claimed by device 2") stage!(sim, "b" => 9.0)
    @test_logs (:warn, r"ClaimedFaceEntry: `a` is claimed by device 1") stage!(sim, "a" => 9.0)

    # A second greedy stakes the empty remainder: legal, useless, said out loud.
    g2 = Pad("gui2")
    @test_logs (:warn, r"EmptyGreedyClaim") attach!(sim, g2, Greedy())
    @test isempty(sim.plane.roster[3].writer.faces)
end

@testset "the harness surface is the unclaimed complement, recomputed at roster changes (§11.3)" begin
    sim = Simulation(two_slots(); h = 1//10)
    @test sim.plane.harness.faces == [:a, :b]        # the empty roster's complement
    d = Pad("d")
    attach!(sim, d, Enumerated("a"))
    @test sim.plane.harness.faces == [:b]
    init!(sim)
    @test_logs (:warn, r"ClaimedFaceEntry: `a` is claimed by device 1 \(Pad\)") stage!(sim, "a" => 1.0)
    stage!(sim, "b" => 2.0)
    run!(sim, 0.1)
    @test port(sim, "", :a) === 0.0
    @test port(sim, "", :b) === 2.0
    # Detach releases the claims: the surface regains the face from the next run.
    detach!(sim, d)
    @test sim.plane.harness.faces == [:a, :b]
    stage!(sim, "a" => 3.0)
    run!(sim, 0.2)
    @test port(sim, "", :a) === 3.0
end

@testset "the recompilation seam: a pending harness batch is renormalized at attach (§11.4)" begin
    sim = Simulation(two_slots(); h = 1//10)
    stage!(sim, "a" => 1.0, "b" => 2.0)              # staged while stopped, roster still empty
    # The attach reshapes the pending batch through the new schema, discarding
    # the newly claimed face with the incumbent and the site named.
    @test_logs (:warn, r"ClaimedFaceEntry: `a` is claimed by device 1.*renormalization") #=
        =# attach!(sim, Pad("d"), Enumerated("a"))
    init!(sim)
    run!(sim, 0.1)
    @test port(sim, "", :a) === 0.0                  # discarded at the attach, never drained
    @test port(sim, "", :b) === 2.0                  # reshaped, re-staged, drained
    # At detach the surface only broadens: every pending entry survives.
    sim2 = Simulation(two_slots(); h = 1//10)
    d2 = Pad("d2")
    attach!(sim2, d2, Enumerated("a"))
    stage!(sim2, "b" => 4.0)
    detach!(sim2, d2)
    init!(sim2)
    run!(sim2, 0.1)
    @test port(sim2, "", :b) === 4.0
end

@testset "device ids are monotonic per Simulation and never reused (§11.3)" begin
    sim = Simulation(two_slots(); h = 1//10)
    d1 = Pad("d1")
    @test attach!(sim, d1, Enumerated("a")).id == 1
    @test attach!(sim, Pad("d2"), Enumerated("b")).id == 2
    detach!(sim, d1)
    err = failure(() -> attach!(sim, Pad("dx"), Enumerated("b")))     # rejected: ClaimConflict
    @test err isa BuildError
    @test attach!(sim, Pad("d3"), Enumerated("a")).id == 3            # not 1, and no id burned
end

@testset "the roster is frozen per run: attach and detach are stopped-sim operations (§11.3)" begin
    sim = Simulation(chain3(); h = 1//100000)
    init!(sim)
    d = Pad("d")
    @test attach!(sim, d, Enumerated("u")).id == 1   # also warms both compile paths, so
    detach!(sim, d)                                  # the mid-run checks below race no JIT
    t = Threads.@spawn run!(sim, 1.0)                # 100k frames: alive throughout the checks
    while !(@atomic sim.plane.running)
        yield()
    end
    # Inline try/catch rather than `failure`: a fresh closure would JIT-compile
    # mid-run, and the run could end inside that pause.
    err_a = try attach!(sim, d, Enumerated("u")) catch e; e end
    err_d = try detach!(sim, d) catch e; e end
    wait(t)
    @test err_a isa BuildError && occursin("ServiceLifecycle", err_a.msg)
    @test err_d isa BuildError && occursin("ServiceLifecycle", err_d.msg)
    # The freeze lifts with the run: the same operations are legal again.
    @test attach!(sim, d, Enumerated("u")).id == 2
    detach!(sim, d)
end

@testset "the frame's outcome is a pure function of the drained batches, whoever staged them (§11.4)" begin
    sim = Simulation(two_slots(); h = 1//10)
    da, db = Pad("da"), Pad("db")
    ha = attach!(sim, da, Enumerated("a"))
    hb = attach!(sim, db, Enumerated("b"))
    init!(sim)
    run!(sim, 0.3)
    stage!(ha, "a" => 0.7)
    stage!(hb, "b" => -1.3)
    run!(sim, 0.8)
    ref = Simulation(two_slots(); h = 1//10)
    init!(ref)
    run!(ref, 0.3)
    set_slot!(ref, "a", 0.7)
    set_slot!(ref, "b", -1.3)
    run!(ref, 0.8)
    @test port(sim, "children/s", :e) === port(ref, "children/s", :e)
end

@testset "an empty drain stays free with a populated roster (§11.1, §11.4)" begin
    sim = Simulation(two_slots(); h = 1//10)
    attach!(sim, Pad("da"), Enumerated("a"))
    attach!(sim, Pad("gui"), Greedy())
    init!(sim)
    @test @ballocated(drain!($sim)) == 0
end
