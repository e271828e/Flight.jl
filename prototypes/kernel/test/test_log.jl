# --- the log (§11.2, increment 11) ----------------------------------------------
# Publication now follows *every* boundary — t* boundaries included, retiring
# the increment-9 stand-in — and the log is retention policy over exactly those
# publications: the switch, the stride, the bound with its progressive
# re-decimation, and the two unconditional endpoints outside it.

@testset "every boundary publishes: t* included, boundary-consistent (§11.2, §10.6)" begin
    m = Group((; src = Sawtooth(1.0), s = Stamper(0.315)),
              ("children/src/q" => "children/s/sig",), (), ())
    sim = Simulation(m; h = 1//10)
    init!(sim)
    run!(sim; t_end = 0.5)
    snaps = logged(sim)
    t★ = modes(sim, "children/s").t_fired
    i = findfirst(s -> s.t == t★, snaps)   # bitwise: a snapshot published at t* itself
    @test i !== nothing
    @test length(snaps) == 7               # boundary zero + 5 frame tops + one t*
    @test [s.frame for s in snaps] == [0, 1, 2, 3, 4, 4, 5]   # t* shares its frame's ordinal
    # The t* snapshot is the settled boundary's: the re-sweep after the firing
    # is what it captures, so `armed` has already dropped — while the boundary
    # before it still shows the armed value. Boundary-consistency at t*.
    @test port(snaps[i], "children/s", :armed) === false
    @test port(snaps[i-1], "children/s", :armed) === true
    ts = [s.t for s in snaps]
    @test issorted(ts) && allunique(ts)    # t never decreases, t* strictly inside
end

@testset "the log is the publications themselves: full density, zero copies (§11.2)" begin
    sim = Simulation(fed(Plant(), "u"); h = 1//10)
    set_slot!(sim, "in", 1.0)
    init!(sim)
    step!(sim; t_plus = 1.0)
    snaps = logged(sim)
    @test [s.frame for s in snaps] == collect(0:10)      # boundary zero + every frame top
    @test snaps[end] === latest(sim)                     # the same object publication handed out
    ys = [port(s, "children/c", :y) for s in snaps]
    @test issorted(ys) && allunique(ys)                  # the step response, one value per boundary
    @test ys[end] === port(sim, "children/c", :y)        # the terminal endpoint is the live boundary
    step!(sim; t_plus = 1.0)
    @test length(logged(sim)) == 21                      # the session accumulates: one trajectory, one log
end

@testset "log_every thins retention, never publication (§11.2)" begin
    sim = Simulation(fed(Plant(), "u"); h = 1//10, log_every = 3)
    init!(sim)
    step!(sim; frames = 1)
    @test latest(sim).frame == 1                         # published to live readers…
    @test [s.frame for s in logged(sim)] == [0, 1]       # …not retained: `last` alone holds it
    step!(sim; t_plus = 0.9)
    @test [s.frame for s in logged(sim)] == [0, 3, 6, 9, 10]
end

@testset "the endpoints are unconditional and outside the bound (§11.2)" begin
    sim = Simulation(fed(Plant(), "u"); h = 1//10, log_every = 4, log_max = 2)
    init!(sim)
    run!(sim; t_end = 4.0)
    snaps = logged(sim)
    @test [s.frame for s in snaps] == [0, 16, 32, 40]    # two generations in
    @test sim.log.live == 2                              # the bound counts the middle alone
end

@testset "re-decimation: stride doubles, coverage stays global, the bound holds continuously (§11.2, D-137)" begin
    sim = Simulation(fed(Plant(), "u"); h = 1//10, log_max = 8)
    init!(sim)
    ok_bound = ok_ends = ok_sorted = true
    for k in 1:128                                       # one frame at a time: every
        step!(sim; frames = 1)                           # intermediate state is checked
        ok_bound &= sim.log.live ≤ 8
        snaps = logged(sim)
        ok_ends &= snaps[1].frame == 0 && snaps[end].frame == k
        ts = [s.t for s in snaps]
        ok_sorted &= issorted(ts) && allunique(ts)
    end
    @test ok_bound && ok_ends && ok_sorted
    # 128 = 16·8 boundaries, four generations in: the middle sits at exactly
    # stride·(1..max) — coverage global at the effective stride, gap-free —
    # and the retained final boundary dedups against the terminal endpoint.
    @test sim.log.stride == 16
    @test [s.frame for s in sim.log.snaps] == collect(16:16:128)
    @test length(logged(sim)) == 9

    # The effective stride composes with the authored one: log_every · 2^k.
    sim2 = Simulation(fed(Plant(), "u"); h = 1//10, log_every = 2, log_max = 4)
    init!(sim2)
    run!(sim2; t_end = 4.0)
    @test [s.frame for s in logged(sim2)] == [0, 8, 16, 24, 32, 40]
    @test sim2.log.stride == 16                          # 2 · 2³
end

@testset "log = false retains nothing; publication is upstream of the switch (§11.2)" begin
    sim = Simulation(fed(Plant(), "u"); h = 1//10, log = false)
    init!(sim)
    run!(sim; t_end = 0.5)
    @test isempty(logged(sim))
    @test latest(sim).frame == 5
end

@testset "view policies, never trajectory-determining (§11.2)" begin
    # A localized-event trajectory, so the t* machinery is in the loop too:
    # retention differing in every keyword, the trajectory bitwise the same.
    mk(; kw...) = Simulation(single(Bouncer(1.0, 0.315)); h = 1//10, kw...)
    a, b, c = mk(), mk(log = false), mk(log_every = 7, log_max = 3)
    foreach(init!, (a, b, c))
    foreach(s -> run!(s; t_end = 2.0), (a, b, c))
    qa = state(a, "children/c").q
    @test qa === state(b, "children/c").q
    @test qa === state(c, "children/c").q
    ca = modes(a, "children/c").count
    @test ca === modes(b, "children/c").count
    @test ca === modes(c, "children/c").count
end

@testset "a warm restart is a new trajectory: the log starts over (§11.2)" begin
    sim = Simulation(fed(Plant(), "u"); h = 1//10)
    init!(sim)
    run!(sim; t_end = 1.0)
    @test length(logged(sim)) == 11
    init!(sim)
    @test [s.frame for s in logged(sim)] == [0]          # the new boundary zero, alone
end

@testset "logged is a stopped-sim read behind the §11.3 gate" begin
    sim = Simulation(fed(Plant(), "u"); h = 1//100000, log_max = 16)
    init!(sim)
    logged(sim)                                          # warms the compile path, so the
                                                         # mid-run check below races no JIT
    t = Threads.@spawn run!(sim; t_end = 1.0)
    while lifecycle(sim) !== :running
        yield()
    end
    err = try logged(sim) catch e; e end
    wait(t)
    @test err isa BuildError && occursin("ServiceLifecycle", err.msg)
    @test length(logged(sim)) ≥ 2                        # the gate lifts with the run
end

@testset "the retention keywords are validated with their siblings (§11.2)" begin
    m = fed(Plant(), "u")
    @test occursin("log must", failure(() -> Simulation(m; h = 1//10, log = 1)).msg)
    @test occursin("log_every", failure(() -> Simulation(m; h = 1//10, log_every = 0)).msg)
    @test occursin("log_max", failure(() -> Simulation(m; h = 1//10, log_max = 0)).msg)
    @test occursin("log_max", failure(() -> Simulation(m; h = 1//10, log_max = 1.5)).msg)
    sim = Simulation(m; h = 1//10, log_max = Inf)        # the explicit opt-out
    init!(sim)
    run!(sim; t_end = 1.0)
    @test length(logged(sim)) == 11
end
