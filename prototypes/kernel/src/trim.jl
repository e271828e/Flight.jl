# The trim service (§14.7, §14.8): the problem value the aircraft author ships,
# the solver seam it is driven across, the in-house Levenberg–Marquardt behind
# that seam, and the report it comes back as.
#
# Everything the service needs already exists one file up. The write side is the
# condition algebra's two application registers — the dynamic walk for the one
# setup application, the specialized register for the iterations (§14.4, D-066)
# — and the read side is the compiled reader (readers.jl). What trim adds is the
# loop between them, the vectorization at its edges, and the discipline the
# whole section is about: **the simulation's stores have exactly one writer, the
# commit through boundary zero**. Every evaluation happens on scratch executors
# the invocation instantiates and discards (§9.2's one-owner rule, §14.8's
# scratch-store paragraph), so a solve that does not converge leaves the
# simulation bit-for-bit untouched — "never initialized" included.

# --- the problem (§14.7) --------------------------------------------------------

"""
    TrimProblem(; guess, lower, upper, condition, reads, residuals, tolerances)

What the aircraft author ships: what the solver may vary, what those decisions
make of the model, what to read back after each evaluation, and which equations
the readings must satisfy. The field set is **normative and closed** (§14.7),
and every field is required — there is no defaulting a half-specified problem
into a solve.

- `guess`, `lower`, `upper` — same-*named* all-`Float64` NamedTuples. The names
  are the pairing and order carries no semantics (§9.5): the bounds are
  canonicalized to `guess`'s field order at setup, so a permuted spelling is a
  non-event rather than one variable's bound silently applied to another.
- `condition` — the condition-*valued* function over the decisions (§14.9): a
  problem is an implicitly specified condition, and solving makes it explicit.
  Everything fixed per problem is closed over; `d` is the one value that cannot
  be.
- `reads` — the declared read set (`reads(…)`, §14.4), gathered after each
  evaluation.
- `residuals` — `residuals(reads::NamedTuple, d::NamedTuple) → NamedTuple`, the
  residual *system* as named equations, never a scalar cost.
- `tolerances` — an all-`Float64` NamedTuple, same-named as the residual return.
  It rides *in the problem* because a relocated problem carries its own
  convergence test (§14.7).

The value is inert: nothing is validated here, because most of what §14.7
requires is only checkable against a build — a selector resolves against a
model, and the residual key set is observed at the setup guess evaluation.
`trim!` runs the whole list in one collecting pass (`TrimProblemInvalid`).
"""
struct TrimProblem{G,L,U,C,R,F,T}
    guess::G
    lower::L
    upper::U
    condition::C
    reads::R
    residuals::F
    tolerances::T
end

TrimProblem(; guess, lower, upper, condition, reads, residuals, tolerances) =
    TrimProblem(guess, lower, upper, condition, reads, residuals, tolerances)

"""
The service's own seeding tag (§9.4): an activation is keyed by a *concrete*
scalar, and `Dual{TrimTag,Float64,N}` is the one trim seeds its decisions
through. Owning the tag is what keeps a trim activation distinguishable from a
linearization's or a user's own.
"""
struct TrimTag end

# --- the report (§14.8) ---------------------------------------------------------

"""
What `trim!` returns: a value to read, never an exception to catch (§14.8).
Non-convergence is an expected *outcome* — hitting the infeasible edge of an
envelope is information — so it is reported, not thrown; a malformed problem is
the different case, and that one throws at setup.

- `converged` — **the service's** per-residual box test at the point the backend
  returned, evaluated in the residuals' own physical units. That verdict, and
  nothing else, gates the commit (§14.8).
- `solution` — guess-shaped, hence warm-startable by `merge`.
- `residuals`, `tolerances` — the solved-point numbers the verdict was read off,
  beside the thresholds they were read against.
- `committed_residuals` — the same residuals re-gathered from the boundary-zero
  world after the commit, or `nothing` when there was no commit. There is no
  `committed` flag: a converged solve is always committable (§14.8), so the
  absence of the commit *is* the absence of these numbers.
- `status`, `nevals`, `niters` — the backend's, verbatim and authoritative over
  nothing (D-158). `:bypassed` is the zero-decision problem's, where there was
  no backend call at all.
- `saturated` — the decisions sitting at a bound at the returned point, as
  `(name, :lower | :upper)`: the classic CG-limit diagnostic (D-070).
- `fired_events` — the events boundary zero fired at the commit, as
  `(path, name)`; non-empty means the committed stores sit at the post-handler
  point rather than the reported solution, and it raises `TrimCommitEvents`
  beside this field (§14.5, §14.8).
"""
struct TrimReport
    converged::Bool
    solution::NamedTuple
    residuals::NamedTuple
    tolerances::NamedTuple
    committed_residuals::Union{Nothing,NamedTuple}
    status::Symbol
    nevals::Int
    niters::Int
    saturated::Vector{Tuple{Symbol,Symbol}}
    fired_events::Vector{Tuple{String,Symbol}}
end

# --- the backend seam (§14.8) ---------------------------------------------------

"""
    solve(backend, eval!, d0, lower, upper, tol) → (; d, status, nevals, niters)

The backend seam: a **pinned signature, value-passed**, one required method per
backend (§14.8). The backend sees vectors and never names — the declared side is
canonical on both ends and the service has canonicalized before packing — so the
solution it returns unpacks by the order it was given.

`eval!(r, J, d)` is in-place, always fills `r` — the residual vector in
`tolerances`' field order — and fills `J` **iff `J !== nothing`**. The Jacobian
is requested by argument, so a Jacobian-free backend simply always passes
`nothing` and there is still exactly one evaluation method to implement.

`d0`, `lower` and `upper` are packed in `guess`'s field order, `±Inf` meaning
unbounded; a backend that ignores bounds ignores two vectors rather than a
missing argument. `tol` is in `tolerances`' order and is data the backend *may*
stop on, decisive of nothing: the convergence verdict is the service's, re-read
at the returned point (§14.8).

`status` comes from a deliberately **open** set and is recorded verbatim in the
report.
"""
function solve end

"""
    LevenbergMarquardt(; maxiter = 100, λ₀ = 1e-3)

The default backend (§14.8): in-house, dense, ~100 lines — the §10.2 stepper
precedent applied to the solver (a tiny needed core against a heavy dependency),
sharpened by the fact that the per-residual physical tolerances are a
convergence test no external package spells natively.

A damping loop on `(JᵀJ + λ·diag(JᵀJ)) δ = −Jᵀr`, accept/reject on the
**normalized** residual norm `‖r ./ tol‖` with `λ` scaled down by 10 on an
accepted step and up by 10 on a rejected one, each step **projected onto the
box** (§14.8's bound treatment), stopping when `all(abs.(r) .≤ tol)` — the
per-residual test in the service's own units — or at `maxiter`, or when a step
stalls below the `eps` scale of the decisions it would move.

Descent test and stopping rule are measured in the same units on purpose. In
this register "the tolerances *are* the stopping criterion" (§14.8), so LM's
damping loop tests exactly what the service will re-test: a raw `norm(r)` sums
forces against moments, and on a problem whose residuals differ by orders of
magnitude in physical scale it can reject every projected step short of a point
the box test would accept.
"""
struct LevenbergMarquardt
    maxiter::Int
    λ₀::Float64
end

LevenbergMarquardt(; maxiter::Integer = 100, λ₀::Real = 1e-3) =
    LevenbergMarquardt(Int(maxiter), Float64(λ₀))

# The damping bracket. Below `λMIN` the damped normal equations are the
# undamped ones to rounding; above `λMAX` the step is `−Jᵀr/λ·diag` scaled into
# the noise, and a rejection there is a stall rather than a slower descent.
const LM_λMIN = 1e-12
const LM_λMAX = 1e12

_within(r, tol) = all(i -> abs(r[i]) ≤ tol[i], eachindex(r))

# The residual norm in the tolerances' own units — the descent test, measured
# where the box test is (§14.8). Written as a fold rather than `norm(r ./ tol)`
# so it costs no temporary.
_scaled_norm(r, tol) = sqrt(sum(i -> (r[i] / tol[i])^2, eachindex(r)))

function solve(bk::LevenbergMarquardt, eval!, d0::Vector{Float64},
               lower::Vector{Float64}, upper::Vector{Float64}, tol::Vector{Float64})
    n, m = length(d0), length(tol)
    d, r, J = copy(d0), zeros(m), zeros(m, n)
    dt, rt, Jt = similar(d), similar(r), similar(J)
    eval!(r, J, d)
    nevals, λ = 1, bk.λ₀

    for iter in 1:bk.maxiter
        _within(r, tol) && return (; d, status = :converged, nevals, niters = iter - 1)
        # The normal equations of the linearized step, with Marquardt's own
        # scaling: the damping rides the curvature of each column rather than
        # the identity, so a decision the residuals barely respond to is not
        # frozen by the same λ that steadies a stiff one. A zero column — a
        # decision whose Jacobian entries all vanished, the saturated-actuator
        # case §14.7 names — takes the floor and keeps the system definite.
        A, g = J' * J, J' * r
        damp = Diagonal(max.(diag(A), LM_λMIN))
        nrm = _scaled_norm(r, tol)
        accepted = false
        while true
            δ = (A + λ * damp) \ (-g)
            all(isfinite, δ) || break
            @. dt = clamp(d + δ, lower, upper)          # step projection onto the box
            # A step the decisions cannot represent is the end of the descent,
            # whatever the residuals still are: reported as its own status, not
            # dressed up as convergence or as an exhausted iteration count.
            maximum(abs, dt .- d) ≤ eps(maximum(abs, d) + 1.0) && break
            eval!(rt, Jt, dt)
            nevals += 1
            if _scaled_norm(rt, tol) < nrm
                d .= dt; r .= rt; J .= Jt
                λ = max(λ / 10, LM_λMIN)
                accepted = true
                break
            end
            λ *= 10
            λ > LM_λMAX && break
        end
        accepted || return (; d, status = :stalled, nevals, niters = iter)
    end
    (; d, status = _within(r, tol) ? :converged : :maxiter, nevals, niters = bk.maxiter)
end

# --- setup validation (§14.7, §13.1) --------------------------------------------
# One collecting pass over what a build can answer, then the guess evaluation,
# which is where the residual return is observed. Both throws are §13.1's
# register: `TrimProblemInvalid` values for the problem's own fields, with the
# read set's `TapResolution` values spliced in beside them — the read keeps its
# own kind and the *problem* is what the setup is refusing, which is what
# `_resolve_reads` was factored apart for.

function _report_trim!(viol::Vector{Diagnostic})
    isempty(viol) && return nothing
    throw(BuildError(viol))
end

_tviol(field::Symbol, reason::Symbol; kw...) =
    TrimProblemInvalid(; field = field, reason = reason, kw...)

# `guess`, `lower` and `upper`: NamedTuples, one key set between them, all
# fields `Float64`. The key-set comparison is by *set*, order being no
# mismatch — the canonicalization below pairs a permuted spelling by name.
function _check_decisions!(viol::Vector{Diagnostic}, p::TrimProblem)
    named = true
    for (name, v) in ((:guess, p.guess), (:lower, p.lower), (:upper, p.upper))
        v isa NamedTuple && continue
        push!(viol, _tviol(name, :not_a_namedtuple; observed = typeof(v)))
        named = false
    end
    named || return nothing
    for (name, v) in ((:lower, p.lower), (:upper, p.upper))
        Set(keys(v)) == Set(keys(p.guess)) ||
            push!(viol, _tviol(name, :key_set; names = collect(keys(v)),
                               expected = collect(keys(p.guess))))
    end
    for (name, v) in ((:guess, p.guess), (:lower, p.lower), (:upper, p.upper))
        _check_floats!(viol, name, v)
    end
    # And the box has to admit a point at all. Checked per decision, over the
    # pairs that survived the two checks above, because an inverted pair is a
    # box no projection can honor — `clamp` would answer the upper bound and
    # the report would name a saturation nobody authored.
    for k in keys(p.guess)
        (haskey(p.lower, k) && haskey(p.upper, k) &&
         p.lower[k] isa Float64 && p.upper[k] isa Float64) || continue
        p.lower[k] ≤ p.upper[k] ||
            push!(viol, _tviol(:lower, :inverted_box; key = k, value = p.lower[k],
                               bound = p.upper[k]))
    end
    nothing
end

# `tolerances`: a NamedTuple of `Float64`s, each finite and strictly positive.
# Positivity is load-bearing rather than cosmetic. A tolerance is the half-width
# of the box its residual has to sit in, so zero and negative name no box at
# all — and the acceptance test measures `‖r ./ tol`‖ (§14.8), so a non-positive
# one sends the descent test to `Inf`/`NaN`, rejects every trial step and
# returns `:stalled` at the guess. That is a malformed problem, named here.
function _check_tolerances!(viol::Vector{Diagnostic}, p::TrimProblem)
    if !(p.tolerances isa NamedTuple)
        push!(viol, _tviol(:tolerances, :not_a_namedtuple; observed = typeof(p.tolerances)))
        return nothing
    end
    _check_floats!(viol, :tolerances, p.tolerances)
    for k in keys(p.tolerances)
        v = p.tolerances[k]
        v isa Float64 || continue           # the type violation is already named above
        (isfinite(v) && v > 0) ||
            push!(viol, _tviol(:tolerances, :nonpositive_tolerance; key = k, value = v))
    end
    nothing
end

function _check_floats!(viol::Vector{Diagnostic}, name::Symbol, v::NamedTuple)
    bad = Pair{Symbol,Any}[k => typeof(v[k]) for k in keys(v) if !(v[k] isa Float64)]
    isempty(bad) || push!(viol, _tviol(name, :field_types; bad = bad))
    nothing
end

# The read set: a `reads(…)` value whose selectors resolve against this build.
# Its violations are §14.4's own `TapResolution` values, spliced into this list
# as they are — the read is named as authored, with the list in hand — because
# here the *problem* is what is malformed (§14.8) and one throw reports both
# halves of it.
function _check_reads!(viol::Vector{Diagnostic}, p::TrimProblem, b::Build)
    if !(p.reads isa Reads)
        push!(viol, _tviol(:reads, :not_a_read_set; observed = typeof(p.reads)))
        return nothing
    end
    (reader, rviol) = _resolve_reads(p.reads, b, Float64)
    append!(viol, rviol)
    reader
end

# The residual return, observed at the setup guess evaluation (§14.7): its key
# set is `tolerances`' — order free, the return being reordered to it before
# packing — and every field is a real scalar, the residual system being named
# *equations*.
function _check_residuals(r, tolerances::NamedTuple)
    viol = Diagnostic[]
    if !(r isa NamedTuple)
        push!(viol, _tviol(:residuals, :not_a_namedtuple; observed = typeof(r)))
    else
        Set(keys(r)) == Set(keys(tolerances)) ||
            push!(viol, _tviol(:residuals, :key_set; names = collect(keys(r)),
                               expected = collect(keys(tolerances))))
        bad = Pair{Symbol,Any}[k => typeof(r[k]) for k in keys(r) if !(r[k] isa Real)]
        isempty(bad) || push!(viol, _tviol(:residuals, :field_types; bad = bad))
    end
    _report_trim!(viol)
end

# --- the service (§14.8) ---------------------------------------------------------

"""
    trim!(sim, problem; baseline, t₀ = 0.0, backend = LevenbergMarquardt()) → TrimReport

Solve `problem` over `baseline` and, on convergence, commit the solution through
boundary zero. Returns a [`TrimReport`](@ref); non-convergence never throws.

**The scratch world, in two halves (D-213).** Every invocation instantiates its
own executors from the build's cached layouts and discards them with the call.
Setup applies `override(baseline, condition(guess))` to a *nominal* set by the
dynamic walk, checks §14.6's root-input totality against the build's faces
before any evaluation (`UninitializedInputs`, all-or-nothing), and runs one
**establishment round** — boundary zero's sweep with every discrete output
stage admitted, due or not (D-205), with no projection, no guards and no `g`.
The seeded set is then written by the specialized register, and its **frozen
cells are copied from the nominal set** as zero-partial constants: at the
seeded activation the discrete tier never runs (§9.4), so nothing there can
derive a discrete output cell from the authored `s`, and without the copy those
cells would hold the build probe's synthesized values — a fabricated zero being
a fine probe input and a terrible flight condition (§14.6's barrier, reaching
the scratch world). The iterations themselves are untouched: raw write → sweep
→ read cycles at the seeded activation, the continuous chain and `f` alone, no
boundaries and no events (§14.5).

**The verdict is the service's, uniformly.** After the backend returns, the
service evaluates once more at the returned point and reads `converged` off the
per-residual box test in the residuals' own units. That, and never a backend's
`status`, gates the commit (§14.8). No convergence means no commit, and no
commit means the simulation is bit-for-bit untouched, lifecycle included: a
`built` simulation stays `built`.

**The commit is literally an `init!`** — `init!(sim, override(baseline,
condition(solution)); t₀)`, the same composite over the same baseline, so its
totality is setup's and the check is structurally unfailable through this path.
Boundary zero's `project` and any guard already holding then move the committed
stores off the solved point; both movers are surfaced rather than left silent,
as the report's fired-event list with `TrimCommitEvents` beside it, and as the
committed-state residuals with `TrimCommitResiduals` when they leave the box.

**The zero-decision problem is legal** (§14.8): `guess = (;)` bypasses the
solver outright — nothing to pack, no seeded activation, no backend call. The
nominal half's establishment round is the one evaluation, the ordinary box test
decides, and the commit runs as usual. That degenerate problem is the "is this
operating point an equilibrium?" probe, useful in its own right and free.

`trim!` is a stopped-sim service: legal in `built`, `initialized` and
`stopped`, refused while `running` and on an `errored` simulation, exactly as
`init!` is. The prototype spells the anchor `t₀` where the spec writes `t0`.
"""
function trim!(sim::Simulation{Float64}, problem::TrimProblem; baseline,
               t₀::Real = 0.0, backend = LevenbergMarquardt())
    lc = lifecycle(sim)
    lc === :running && throw(BuildError(ServiceLifecycle(op = "trim!", status = :running)))
    lc === :errored && throw(BuildError(ServiceLifecycle(op = "trim!", status = :errored)))

    b = sim.build
    viol = Diagnostic[]
    _check_decisions!(viol, problem)
    _check_tolerances!(viol, problem)
    reader = _check_reads!(viol, problem, b)
    _report_trim!(viol)                       # one collected throw, before any evaluation

    guess, tolerances = problem.guess, problem.tolerances
    K, RK = keys(guess), keys(tolerances)
    N = length(K)
    tol = Float64[tolerances[k] for k in RK]

    # --- the nominal half (D-213) ------------------------------------------------
    ex_nom = _scratch(sim, Float64)
    plan = resolve_condition(override(baseline, problem.condition(guess)), b, Float64)
    assert_total(plan, b.flat, "trim!")       # (§14.6): pre-evaluation, all-or-nothing
    apply!(ex_nom, plan)
    _round!(ex_nom, ESTABLISH)                # every discrete output stage, due or not
    ex_nom.bodies.rhs()
    r0 = problem.residuals(gather(reader, ex_nom), guess)
    _check_residuals(r0, tolerances)          # the return, observed where §14.7 says

    if N == 0
        # The solver is bypassed outright: the establishment round above is the
        # one evaluation, and the ordinary box test decides (§14.8).
        r = Float64[NamedTuple{RK}(r0)[k] for k in RK]
        return _verdict!(sim, problem, baseline, guess, r, tol, :bypassed, 1, 0,
                         Tuple{Symbol,Symbol}[], reader, t₀)
    end

    # --- the seeded half ---------------------------------------------------------
    TD = ForwardDiff.Dual{TrimTag,Float64,N}
    act = activation(b, TD)                   # the cached Stratum-C re-run (§9.4)
    ex = _scratch(sim, TD, act)
    _establish_frozen!(ex, act, ex_nom, sim.build)
    d_dual = _seeded(K, guess, TD)
    plan_d = compile_plan(override(baseline, problem.condition(d_dual)), b, TD)
    reader_d = _compile_reads(problem.reads, b, TD)

    # One evaluation: write the decisions through the shape-compiled plan, sweep,
    # gather, and take `r` off the values and `J` off the partials of one and the
    # same seeded pass (§14.7).
    #
    # §14.7's return check runs a second time on the *first* seeded evaluation.
    # The nominal guess evaluation observes the return at `Float64`; a residual
    # lambda that branches on the scalar it is handed — on `eltype`, on a
    # method that has no `Dual` arm — can answer a different key set here, and
    # without the re-check the reorder below raises a bare `ErrorException`
    # from `NamedTuple{RK}` instead of the collected `TrimProblemInvalid` that
    # names what is wrong. Once only: the per-iteration path takes one `Bool`
    # load and no work at all.
    checked = Ref(false)
    function eval!(r::Vector{Float64}, J, d::Vector{Float64})
        d_nt = _seeded(K, d, TD)
        apply!(ex, plan_d, override(baseline, problem.condition(d_nt)))
        evaluate!(ex)
        raw = problem.residuals(gather(reader_d, ex), d_nt)
        if !checked[]
            checked[] = true
            _check_residuals(raw, tolerances)
        end
        res = NamedTuple{RK}(raw)
        for i in eachindex(r)
            v = res[i]
            r[i] = ForwardDiff.value(v)
            J === nothing && continue
            for j in 1:N
                J[i, j] = ForwardDiff.partials(v, j)
            end
        end
        nothing
    end

    lower = Float64[problem.lower[k] for k in K]
    upper = Float64[problem.upper[k] for k in K]
    # The guess enters the box here, once. §14.8 honors bounds by step
    # projection, and a backend's returns that never *step* — the
    # already-within-tolerance return, a stall at the first iteration — would
    # otherwise hand back the unprojected guess: an out-of-box point committed
    # as converged, with `saturated` empty because it sits at no bound. Every
    # point the backend sees and returns lies in the box because the one it
    # starts from does.
    d0 = clamp.(Float64[guess[k] for k in K], lower, upper)
    out = solve(backend, eval!, d0, lower, upper, tol)

    # The verdict, at the backend's returned point and in the service's own
    # units: one residual evaluation, noise against the solve that produced it.
    r = zeros(Float64, length(tol))
    eval!(r, nothing, out.d)
    _verdict!(sim, problem, baseline, NamedTuple{K}(Tuple(out.d)), r, tol,
              out.status, out.nevals, out.niters, _saturated(K, out.d, lower, upper),
              reader, t₀)
end

# A non-nominal deployment is refused rather than served: the commit runs
# through boundary zero on the simulation's own stores, and those are the
# nominal world's. The seeded activation trim needs is the service's scratch,
# never the deployment's.
trim!(sim::Simulation, ::TrimProblem; kw...) = throw(BuildError(
    ArgumentInvalid(call = :trim!, reason = :non_nominal, value = string(typeof(sim)))))

trim!(::Simulation, other; kw...) = throw(BuildError(
    TrimProblemInvalid(field = :problem, reason = :not_a_problem, observed = typeof(other))))

# --- the pieces the service is built out of --------------------------------------

# One scratch executor: the same buffer set the `Simulation` owns, at whatever
# scalar, from the same cached layouts and the same bound entry data — and it
# dies with the call (§9.2, §14.8, glossary `scratch`).
_scratch(sim::Simulation, ::Type{T}) where {T} = _scratch(sim, T, activation(sim.build, T))
_scratch(sim::Simulation, ::Type{T}, act::Activation{T}) where {T} =
    compile(sim.build, act, sim.D, sim.Φ, sim.Δt; chunk_size = sim.chunk_size)

# D-213's copy: a frozen component's stages are outside the seeded activation's
# executable set (§9.4), so its output cells can only come from the nominal half
# — where the authored `s` has just been published by the establishment round.
# A discrete producer's cells are pinned `Float64` at every activation, so this
# is a value copy; the zero-partial embedding happens where a continuous
# consumer reads them (§14.3).
function _establish_frozen!(ex::Executor, act::Activation{T}, nom::Executor,
                            b::Build) where {T}
    for ci in eachindex(b.flat.comps)
        _frozen(b.tiers, ci, T) || continue
        path = b.flat.paths[ci]
        for name in keys(act.decls[ci].outs)
            scatter!(ex.store, act.layout.addr[(path, name)],
                     gather(nom.store, nom.act.layout.addr[(path, name)]))
        end
    end
    nothing
end

# The decisions at the seeded activation: field `i` carrying the unit partial in
# slot `i`, which is what makes one sweep yield `r` and `J` together (§14.7).
# `v` is indexed positionally, so the packed vector and the guess NamedTuple
# seed through the same code.
_seeded(K::Tuple, v, ::Type{ForwardDiff.Dual{TG,Float64,N}}) where {TG,N} =
    NamedTuple{K}(ntuple(i -> ForwardDiff.Dual{TG}(Float64(v[i]),
                                                   ntuple(j -> Float64(i == j), Val(N))...),
                         Val(N)))

# The decisions sitting at a bound at the returned point (§14.8): the comparison
# is exact because the projection assigns the bound itself, and an infinite
# bound is never met by a finite decision.
function _saturated(K::Tuple, d::Vector{Float64}, lower::Vector{Float64},
                    upper::Vector{Float64})
    out = Tuple{Symbol,Symbol}[]
    for i in eachindex(d)
        d[i] == lower[i] ? push!(out, (K[i], :lower)) :
        d[i] == upper[i] && push!(out, (K[i], :upper))
    end
    out
end

# The verdict, the commit and the report, shared by both registers — the solved
# problem and the bypassed zero-decision one, which differ in how they got their
# residual vector and in nothing after it (§14.8).
function _verdict!(sim::Simulation, p::TrimProblem, baseline, solution::NamedTuple,
                  r::Vector{Float64}, tol::Vector{Float64}, status::Symbol,
                  nevals::Int, niters::Int, saturated::Vector{Tuple{Symbol,Symbol}},
                  reader, t₀)
    RK = keys(p.tolerances)
    residuals = NamedTuple{RK}(Tuple(r))
    converged = _within(r, tol)
    converged || return TrimReport(false, solution, residuals, p.tolerances, nothing,
                                   status, nevals, niters, saturated,
                                   Tuple{String,Symbol}[])

    init!(sim, override(baseline, p.condition(solution)); t₀ = Float64(t₀))

    # The commit's fired events, read off the per-boundary counts right after
    # `init!` returns — they are reset at the next boundary, and there is none
    # (§10.6, §14.5).
    es = sim.exec.events
    fired = Tuple{String,Symbol}[es.names[i] for i in eachindex(es.count) if es.count[i] > 0]
    isempty(fired) || @warn logline(TrimCommitEvents(events = fired))

    # The committed-state residuals, nearly free: that boundary's sweep has just
    # run, so the declared reads need only gather from it — with one `rhs` for
    # the derivative reads, `ẋbuf` being integrator scratch and this a service
    # evaluation (§7.5, §14.8).
    sim.exec.bodies.rhs()
    committed = NamedTuple{RK}(p.residuals(gather(reader, sim.exec), solution))
    off = Tuple{Symbol,Float64,Float64}[(k, Float64(committed[k]), tol[i])
                                        for (i, k) in enumerate(RK)
                                        if !(abs(committed[k]) ≤ tol[i])]
    isempty(off) || @warn logline(TrimCommitResiduals(residuals = off))

    TrimReport(true, solution, residuals, p.tolerances, committed, status, nevals, niters,
               saturated, fired)
end
