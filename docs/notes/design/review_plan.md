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

## WP2 — Signature cluster [D] — DONE (v0.19)

Settled as the named-bundle hand-off (C′) after discussion, plus renames:

- [x] Every function `fn(comp, args)`, destructured by name; bundle law
      (field iff store/tier fact exists; undeclared absent, never
      `nothing`); `t` everywhere, `Δt` discrete-only, `ws` by declaration;
      `project` stays positional. Rejected: positional+clock view, kwargs
      (both `kwarg_decl` and slurp variants), full context object, Clock
      component. (Rows 74.)
- [x] Stage/letter renaming: `f` flow / `g` jump / `h_x`/`h_xu`/`h_z`/
      `h_zu` stages, `y_*` products; bare `h` = step size only;
      wrong-tier stage name = build error. (Row 75.)
- [x] `inputs`/`outputs`/`locals` → `input_types`/`output_types`/
      `local_types` (three-register inventory). (Row 76.)
- [x] Workspace: `init_workspace` → `workspace(::C[, ::Type{T}])`
      allocator, both tiers, `undef` idiom, `T`-generic scratch. (Row 77.)
- [x] Consequentials: §11.5 `comp.Δt` → bundle field (row 19 amended);
      §12.10 claim strengthened; §14.3 probe placeholder period; §14.4
      workspace re-typing; rows 33/35 amendment markers; both sketches
      updated.

## WP3 — Input-face typing: abstract faces as constraints [D] — DONE (v0.20)

Scope grew in discussion far beyond the original three boxes; settled as
the **port-typing cluster** (rows 78–79 + amendments to 33/53/54/55/76):

- [x] Subtype rule `producer_face <: entry` at nominal faces; entries =
      face constraints, not cell types; abstract = structural
      substitutability (§7 handles), never needed for eltype genericity.
- [x] Root-slot carve-out (tight bound determines the type; abstract at
      root = build error) + fan-out sub-rule (unique concrete declaration
      among consumers). `probe_value` never meets an abstract type —
      structural, no new obligation.
- [x] **`T`-signature retired** (user re-opened): `output_types(::C)` /
      `local_types(::C)` concrete nominal, both tiers; activation leaf
      walk (continuous parametrizes `Float64` leaves/`Real` type params,
      `Int`/`Bool`/enum/reference fields pin, discrete pins wholesale);
      §14.5 conformance split — exact at nominal, `{T, Float64}` +
      zero-partial embedding at parametrized leaves (embedding guarantee:
      promotion airtight, no lossy `Dual → Float64` cast).
- [x] Root slots re-typed by the same leaf walk (the §14.1 gap); seeding
      register stated in §16.10; phantom-producer pedagogy in §12.3;
      declarative non-participation recorded as §16.10 door.
- [x] §17.5 tightened to `RQuat{Float64}`; both sketches de-`T`'d;
      consequentials in §6, §9.2, §13.3, §13.5, §13.8, §14.1, §14.3,
      §14.4; header v0.20.

## WP4 — Boundary cluster: Tier-2, baselines, boundary zero [D] — DONE (v0.21)

Settled as rows 80–82 + amendments to rows 1/20/28:

- [x] §2.1 fossil corrected: Tier-2 pace-independent, localization cost =
      §11.7 pacer debt (row 80).
- [x] `t*` = boundary, not frame (§11.4 new paragraph): full §11.6
      iteration w/ per-boundary once-per-event, snapshot + §12.8 counter
      (re-spelled as *boundary* counter) + `stop_on` check; ticks never
      due, no drain, no separate pacing; replay pointers = boundary
      counter, trace frame-indexed (row 81).
- [x] Guard conditions normative (§2.1 + §11.6): positive = holding
      (`g ≥ 0`), not-holding → holding edges vs per-event baselines in
      loop state; boundary-zero baseline = nothing-holds (§16.5 derived);
      warm-restart re-baseline; negated-guard second event for the
      opposite direction (row 82).
- [x] Endpoint policy (user question): holding-endpoint bracket return —
      `t* = tₙ` structurally impossible, `t* = tₙ₊₁` degenerates to the
      grid boundary; earliest-`t*` for multiple guards; grid times
      indexed, never accumulated (remainder targets the grid point)
      (row 82).
- [x] `stop_on` checked at every published boundary from `t₀`; §16.5
      parenthetical fixed (`stop_on` ≠ event firing).

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
