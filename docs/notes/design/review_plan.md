# Coherence-review edit plan (2026-07 full-document pass)

Working checklist for the fixes agreed after the two independent consistency
reviews of `framework_design.md` (v0.17). Sources: session review + external
review, merged and adjudicated 2026-07-28. Delete this file when the last box
is ticked.

Tags: **[M]** mechanical (apply, then review diff) · **[D]** decision
(discuss fork first, then edit). Line numbers refer to v0.17
(commit 824bc8a7) and will drift as edits land.

---

## WP1 — Editorial batch → v0.18 [M] — DONE (v0.18)

- [x] §13.2 (~1734): restore the lost `**events(::C)**:` bullet header
      (content grafted onto the `locals` bullet by the v0.8 edit).
- [x] Stale §19 pointers → §16: `:1248`, `:1301` (§12.3 slot init / trace
      header → §16.6), `:1596` (§13 intro — services settled), `:2088`
      (§14.3 enforcement → §16.6 `UninitializedSlots`), `:3199`, `:3245`
      (§17.4).
- [x] §14.6 body refresh (`:2173`, `:2186`): BOBYQA-as-register and
      gradient-trim-as-"open option" superseded by rows 69–70 (NLS + AD
      default, BOBYQA an extension fallback); "cost normalizations" →
      per-residual scalings.
- [x] Axis numbering: "axis-7 services" → axis 8 at `:1574` (§12.10) and
      row 31 (`:3530`).
- [x] §16 title (`:2501`): drop "part 1: the condition substrate".
- [x] §17.1 (`:3044`): `@project` → `project`.
- [x] §2 (`:101`): "under rules to be fixed in the execution/periphery
      axes" → settled pointer to §12.
- [x] §5.2 (`:349`): "axis 5 should restate it" → points at §11.3/§11.4.
- [x] Row 24 (`:3523`): add missing "Amended in v0.6 → row 44" marker.
- [x] §12.5 title (`:1344`): "active widgets" → "staging contract".
- [x] §17.5 (`:3440`): init-consistency parenthetical → §16.5, noting
      boundary-zero's due `h` discharges most of it (only the `t₀` sample
      depends on authored `z`).

Excluded on purpose: §16.5's `stop_on`-in-guard-list mislabel (`:2698`) —
fixed once, in WP4.

## WP2 — Signature cluster: what do stage functions receive? [D]

One joint decision; three gaps that must land coherently:

- [ ] Time in stages: outputs cannot depend on `t` (capability gap, not a
      contradiction). Lean: add `t` to `g_s1`/`g_s2` — time is already
      ambient in `f`/`h`/guards/handlers; the Clock-component alternative
      makes stages treat time as a signal while everything else doesn't.
      Either way, qualify §12.10's freshness claim (the anticipatory-`h`
      idiom is exact for timetable scripts but undiscoverable).
- [ ] `Δt` mechanism: instance-carried `Δt` impossible in principle
      (§13.6's `===`-identical-siblings argument — two fields, same
      immutable value, different `rates` keys). Row 19 rejected
      `h`-argument-*only*; forced move is passing it to discrete stages
      and `h`; spelling open.
- [ ] Workspace hand-off: §9.3 promises it is "handed to the update
      function"; no signature shows it. Rides the same decision.
- [ ] Fork to settle: three separate arguments vs. one small context
      argument carrying `t`/`Δt`/workspace.

Edit sites: §5.2, §9.3, §11.5, §13.2 example, §17.5 sampler, row 19
amendment + new row, `sketch_decoder.jl`, `sketch_io.jl`.

## WP3 — Input-face typing: abstract faces as constraints [D]

- [ ] Re-open row 33's subtype rejection under row 54's
      producer-determines insight; consumer abstract `inputs` entries as
      constraints (`producer_face <: declared`), exact equality the rule
      for concrete declarations. Demonstrated consumer: §7 field handles
      (strut ↔ terrain substitutability).
- [ ] Root-slot carve-out: faces surfacing as root slots take their type
      from the consumer declaration (no producer) — staging cells,
      `probe_value`, trace typing need them concrete; abstract entries
      legal only on component-fed inputs.
- [ ] Resolve §17.5 sketch faces (`q = RQuat` bare UnionAll): legalized
      under the new rule, or corrected to `RQuat{Float64}` under status
      quo.

Edit sites: §13.2 inputs bullet, §7 cross-ref, row 33 amendment + new
row, §17.5, `sketch_decoder.jl` if affected.

## WP4 — Boundary cluster: Tier-2, baselines, boundary zero [D]

One discussion; one new paragraph (§11.4 or §11.6) plus corrections:

- [ ] §2.1 fossil (`:114`), the direct contradiction: "degrades
      gracefully to Tier 1 in real-time mode" vs. §11.7's bit-identical
      invariant. Lean: Tier-2 always runs; the cost is pacer debt.
- [ ] `t*` boundary status: state the frame-vs-boundary distinction —
      frames are grid steps (drain, pacing, tick eligibility); boundaries
      are published consistency points, of which `t*` is one. At `t*`:
      snapshot yes, event iteration + once-per-event accounting yes,
      ticks never due, drain no (input timing must not depend on
      localization). To decide: frame counter / §12.8 waiters at `t*`;
      trace replay-pointer indexing.
- [ ] Guard baseline convention: previous sample lives in loop state;
      boundary zero establishes the baseline (fires on true/positive at
      `t₀`, made normative incl. sign convention for continuous guards);
      "newly-fired" semantics within §11.6 rounds; re-baseline on warm
      restart (guard memory correctly not in `z`).
- [ ] `stop_on` at boundary zero: `run!` checks the boundary-zero
      snapshot before the first step; fix §16.5's parenthetical
      mislabeling `stop_on` faces as guard firings (`:2698`).

Edit sites: §2.1, §11.4, §11.6, §12.9, §15.5, §16.5, new rows.

## WP5 — Output-device read addressing [D]

- [ ] Resolve the three-way inconsistency (§12.3 `:1228` faces-only /
      §17.4 `:3252` snapshot paths / §15.5 `:2424` citing a "face-binding
      precedent" §17.4 doesn't contain). Proposal: writes speak the root
      contract; reads see the whole public table (output devices =
      diagnostic observation; GUI deep-reading panels as precedent).
- [ ] Rewrite §15.5's observation-by-path rejection to derive purely from
      the load-bearing/diagnostic distinction (drop the citation); scope
      §12.3's "periphery speaks face names only" to the write surface.
      Explicit sign-off required — this rewrites settled argumentation.

## WP6 — Spackle

- [ ] [M] `Δt_base`/`n` in the `Simulation` ctor spelling (§13.7, §14.1,
      §17.4; default `n = 1`).
- [ ] [M] `capture` explicitly includes root slots (§16.1, §16.10) —
      required by §16.6 totality for capture → apply.
- [ ] [M] Handler asymmetry rationale: `x⁺` complete (flat-buffer
      writeback), `m⁺` partial (per-field cells) — one sentence at
      §5.2/§14.5.
- [ ] [D-small] Unconnected-output warning: retire (lean — fires on every
      observation-oriented port under mandatory `outputs`, poisons the
      sole warning stream) or narrow. §8, §15.2, decision log.
- [ ] [M] §4.2: "non-breaking change" → "non-breaking for consumers;
      visible to the scheduler".
- [ ] [M] `m` is continuous-only: state flatly in §3.2; fix §16.1's
      "modes and discrete state (`m`, `z`)" phrasing; discrete FSM enums
      live in `z`.
- [ ] [D] Collections of children: `NTuple{N, C}` fields as children
      (path segments `"aircraft/1"`; `rates`/`connections` addressing to
      settle) vs. one-field-per-child + programmatic generation. Lean:
      admit homogeneous tuples (§16.9 swarm worlds are a named use case).
      Alternative: park as a migration-outline entry.

## WP7 — Additions (last, in this order)

- [ ] Taught-contracts index: one short section gathering §5.2's stage
      doctrine, boundary-sampling (§17.5), interval alignment (§16.5),
      stop-face sampling (§15.5), levels-never-deltas (§12.3),
      one-home-per-datum, liveness (§12.5), plus WP4's guard baselines.
- [ ] Component test-rig idiom in §15.7's library commitments: one-child
      assembly exporting the child's input faces — `design_world`'s
      little sibling.
- [ ] Log/trace persistence: explicit recorded-deferral in §19 (HDF5
      export scope, field-handle summarization, trace file format —
      settled at migration).
- [ ] API synopsis [D-light]: single-page entry-point listing (`build`,
      `Simulation`, `init!`, `trim!`, `linearize`, `capture`, `attach!`,
      `run!`, control plane, device handle, condition algebra,
      `faces`/`resolve`). Written last — closing coherence check over
      WP2–WP6 and bridge into the migration outline. Placement (new
      section vs. appendix) to decide.

---

Process: each [D] package gets full prose discussion of the fork before
any edit, then per-change approval, then one local commit (one doc
version bump per WP) for review before push.
