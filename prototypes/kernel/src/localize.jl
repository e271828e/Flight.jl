# Localization mechanics (§10.4): the frame loop that carries each grid step
# through `integrate → arrival sweep → trigger → θ = 0 validation → bracket →
# root-find → t* → remainder step`, iterated under `localization_budget`. The
# runtime consultation of `Build.policies` lives here — a sign-form guard's
# crossing is bracketed by trial evaluations over the seam's dense output
# (§10.2, `dense!`) and fired at a `t*` boundary, strictly inside the frame.
#
# Frame ≠ boundary (§10.4): the frame `[tₙ, tₙ₊₁]` is the unit of scheduling —
# tick eligibility keys to it, and `t*` is never a tick — while a boundary is a
# published consistency point. Every `t*` produced here is a boundary that is
# not a frame top.

"""Frame top `k`, computed from the index and `t₀` — never accumulated (§10.4)."""
_grid_time(sim::Simulation, k::Int) = sim.exec.clock.t₀ + k * sim.h

"""
    frame!(sim, k)

Advance through the frame `[tₖ₋₁, tₖ]`, leaving the clock at the indexed frame
top with the state and table at their arrival values — the frame-top boundary
itself is the caller's, exactly as before. A model with no localized events
(every non-nominal activation included, its events having compiled out) takes
the bare step: no arrival machinery, no extra sweep, today's exact path.
The one exception is a §13.5 stop observed at a `t*` publication: the frame's
remainder was abandoned, so the clock stays at `t*` — where the stores are —
and the frame top is never stamped.
"""
function frame!(sim::Simulation{T}, k::Int) where {T}
    t_to = _grid_time(sim, k)
    sim.has_localized ? _localized_frame!(sim, t_to) : step!(sim, T(sim.h))
    sim.policy.hit === nothing && (sim.exec.clock.t = t_to)
    nothing
end

# The localization loop (§10.4). Each turn integrates one segment — the whole
# frame first, then remainder steps from each `t*` — and either returns at the
# frame top or fires one `t*` boundary and goes around. The budget counts
# localizations, one per `t*` boundary produced: per segment, root-finding runs
# are already structurally bounded (at most one per declared event), while the
# segment count is the quantity chattering inflates without bound.
function _localized_frame!(sim::Simulation{T}, t_to) where {T}
    es = sim.exec.events
    n = length(es.prior)
    (x₀, _) = startpoint(sim.stepper)         # the seam's retained pair (§10.2):
    count = 0                                 # x₀ = x(t_seg) after each step!
    fill!(es.loc_warned, false)
    while true
        t_seg = sim.exec.clock.t
        h′ = t_to - t_seg
        step!(sim, h′)

        # The arrival sweep at the segment's end — interior, on the raw
        # unprojected state, before any discrete cell refreshes (§10.4): the
        # sweep that closes the integration step is what raises the trigger.
        sim.exec.bodies.sweep_1(); sim.exec.bodies.sweep_2()
        _guards!(es)
        copyto!(es.σ1, es.σ)

        # The trigger (§10.4): localized policy, prior not-holding at the last
        # boundary's quiescence, holding at arrival — §2.1's directional edge.
        any_trig = false
        for i in 1:n
            es.trig[i] = es.localized[i] && !es.prior[i] && es.σ1[i] ≥ 0
            any_trig |= es.trig[i]
        end
        any_trig || return nothing

        # Budget exhaustion degrades; it does not throw (§10.4): the remainder
        # step has already completed, and this crossing fires in the coming
        # boundary's ordinary iteration — boundary granularity for this frame.
        if count ≥ sim.localization_budget
            for i in 1:n
                (es.trig[i] && !es.loc_warned[i]) || continue
                es.loc_warned[i] = true   # at most one report per event per frame
                (path, name) = es.names[i]
                _report!(sim.loop_diag,   # the loop's own cell (§11.8): folded at the next frame top
                         ChatteringBudget(path, name, Float64(t_to),
                                          sim.localization_budget, count))
            end
            return nothing
        end

        copyto!(sim.xnext, sim.exec.xbuf)          # retain xₙ₊₁ before the trials clobber it

        # The θ = 0 validation (§10.4): x(t_seg) back into the buffer, one
        # interior sweep, no interpolant — x̂(0) = xₙ identically. σ₀ is the left
        # bracket value, and it discriminates the edge's cause: σ₀ holding means
        # the frame-top drain flipped the guard (`u` is the only thing that can
        # differ from the prior's evaluation context), so no in-frame crossing
        # exists — not localizing is the action, and it consumes no budget.
        copyto!(sim.exec.xbuf, x₀)
        sim.exec.clock.t = t_seg
        sim.exec.bodies.sweep_1(); sim.exec.bodies.sweep_2()
        _guards!(es)
        copyto!(es.σ0, es.σ)
        remaining = false
        for i in 1:n
            es.trig[i] = es.trig[i] && es.σ0[i] < 0
            remaining |= es.trig[i]
        end
        if !remaining
            copyto!(sim.exec.xbuf, sim.xnext)   # epoch-caused only: fall through to the frame top
            return nothing
        end

        # ẋₙ₊₁, paid only past a validated trigger (§10.4): one sweep and the
        # RHS block at the arrival state completes the interpolant's data.
        copyto!(sim.exec.xbuf, sim.xnext)
        sim.exec.clock.t = t_seg + h′
        evaluate!(sim)
        copyto!(sim.ẋnext, sim.exec.ẋbuf)

        # Root-find each validated event; the boundary fires at the earliest t*
        # (§10.4). Ties need no decision — every standing edge at θ★ fires in
        # the boundary's own iteration — and a later crossing simply does not
        # hold at θ★, so it re-triggers on the remainder and re-localizes there.
        θ★ = 1.0
        for i in 1:n
            es.trig[i] || continue
            θ★ = min(θ★, _crossing(sim, i, es.σ0[i], es.σ1[i], t_seg, h′))
        end
        if θ★ ≥ 1.0
            # Degenerate: the crossing is the frame top itself (§10.4). The
            # localization is discarded and the event fires inside the frame
            # top's ordinary iteration — one boundary, no zero-length remainder.
            copyto!(sim.exec.xbuf, sim.xnext)
            return nothing
        end

        # The t* boundary (§10.4): interpolated state in, then the §10.6
        # sequence with ticks never due — t* is off the harmonic grid by
        # construction, so this is exactly the off-tick boundary: projection
        # (authority rests here, not with the raw trials), the full event
        # phase iterated to quiescence under a fresh firing budget, priors
        # updated from the settled samples — and the settled boundary publishes
        # before integration resumes (§11.2): every boundary is a published
        # consistency point. The interpolant is invalidated by falling out of
        # scope — the handlers made it a lie for t > t*.
        dense!(sim.stepper, sim.exec.xbuf, sim.xnext, sim.ẋnext, θ★, h′)
        sim.exec.clock.t = t_seg + θ★ * h′
        offtick_boundary!(sim)
        publish!(sim)

        # Every publication is a stop-face sampling point (§13.5): a face
        # holding in the t* snapshot makes it the final one — the frame's
        # remainder is abandoned, and the hit reaches the loop through the
        # policy seam.
        face = _stop_hit(sim, sim.policy)
        if face !== nothing
            sim.policy.hit = face
            return nothing
        end
        count += 1
    end
end

"""
One trial evaluation (§10.4): write x̂(θ) into the state buffer, set the
mid-step clock, run the interior sweep, evaluate the guards. Discrete cells
ZOH-hold through it — the interior sweep has no discrete entries — and the
state is raw: projection's reach is the boundary, not the trial.
"""
function _trial!(sim::Simulation, θ::Float64, t_seg, h′)
    dense!(sim.stepper, sim.exec.xbuf, sim.xnext, sim.ẋnext, θ, h′)
    sim.exec.clock.t = t_seg + θ * h′
    sim.exec.bodies.sweep_1(); sim.exec.bodies.sweep_2()
    _guards!(sim.exec.events)
    nothing
end

"""
Bracketed, derivative-free root-finding over trial evaluations (§10.4, D-018):
ITP (Oliveira & Takahashi 2020) on the θ ∈ [0, 1] bracket, entered only with
σ₀ measured not-holding and σ₁ holding — the observed bracket *is* the
convergence certificate. Stops once `hi - lo ≤ localization_tol` (relative: the
bracket in θ against the segment, D-133) and returns the **holding endpoint of
the final bracket** — every `hi` is an actual observation, so the returned
point is strictly later than the segment start and the guard observably holds
there. A return of exactly 1.0 is the degenerate crossing-at-the-frame-top,
discarded by the caller.
"""
function _crossing(sim::Simulation, i::Int, σ₀::Float64, σ₁::Float64, t_seg, h′)
    es, tol = sim.exec.events, sim.localization_tol
    lo, hi = 0.0, 1.0
    σlo, σhi = σ₀, σ₁
    # ITP constants over the unit bracket: κ₁ = 0.2, κ₂ = 2, n₀ = 1; ε is the
    # half-width target, and past n_max the projection radius hits zero, which
    # degrades each remaining iteration to plain bisection.
    ε = tol / 2
    nmax = max(ceil(Int, -log2(tol)), 0) + 1
    j = 0
    while hi - lo > tol
        xh = 0.5 * (lo + hi)
        (lo < xh < hi) || break                       # the bracket is at float resolution
        xf = (lo * σhi - hi * σlo) / (σhi - σlo)      # regula falsi, well-defined: σlo < 0 ≤ σhi
        s = sign(xh - xf)
        δ = 0.2 * (hi - lo)^2
        xt = δ ≤ abs(xh - xf) ? xf + s * δ : xh       # truncate toward the midpoint
        r = max(ε * exp2(nmax - j) - 0.5 * (hi - lo), 0.0)
        θ = abs(xt - xh) ≤ r ? xt : xh - s * r        # project into the minmax radius
        (lo < θ < hi) || (θ = xh)
        _trial!(sim, θ, t_seg, h′)
        σθ = es.σ[i]
        σθ ≥ 0 ? (hi = θ; σhi = σθ) : (lo = θ; σlo = σθ)
        j += 1
    end
    hi
end
