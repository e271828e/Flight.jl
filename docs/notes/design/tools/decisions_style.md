# Decision log style guide

The rules `framework_decisions.md` adheres to. Distilled from the Pass A
conversion brief (2026-08-14) once that pass closed; consulted whenever an
entry is added, amended or rewritten. The counterpart of `spec_style.md`, which
governs the spec.

The log's purpose sets its register: each entry records a **settled decision**
and — the reason the file exists — the alternatives considered and rejected,
with the reasoning that rejected them. Entries state the *current* position;
history lives in git.

## The template

```markdown
### D-nnn — Short imperative title

**Status.** ratified

**Position.** The current decision, stated in one or two sentences, citations
bare (§5.2, D-037).

**Spec.** §a.b, §c.d

**Rationale.** Why — the recorded reasoning.

**Rejected.**
- *Alternative:* why not, as recorded.
- *Superseded position — description:* why it fell, as recorded.

**Divergence.** (rare, optional) Only when the entry explicitly contrasts an
external tool's behavior (Simulink/MATLAB/Modelica does X; we do Y; because Z).
```

## Field rules

- **ID.** `D-nnn`, zero-padded to three digits. Identifiers are stable and
  never reused; they preserve the row numbers of the log's retired table form.
  The heading is exactly `### D-nnn — Title` (spaced em dash) — `linkify.jl`
  and `check_rows.jl` both parse that shape.
- **Title.** Short (up to ~15 words — prefer fewer, but never at the cost of
  clarity), imperative where the decision is an action ("Reject algebraic loops
  at build"), a plain noun phrase where it isn't ("Aggregation mechanism").
  Distinctive enough to identify the entry in an index of ~200. No terminal
  period. For superseded entries, title the *topic*, not the dead position.
- **Status.** `ratified`, or `superseded → D-0NN` when a later entry replaced
  it. Nothing else; no dates.
- **Position.** The headline ruling, one or two sentences. Everything else goes
  to Rationale. Where an entry genuinely bundles parallel rulings with no
  distinct headline, content preservation wins and Position holds the full
  chain — see the Pass C queue.
- **Spec.** The sections where the decided mechanism lives, sorted and
  deduplicated. Omit the field entirely if there are none.
- **Rationale.** The recorded reasoning and any adjudication narrative. Where
  an entry records no rationale beyond its rejections (common in early
  entries): `Recorded only through the rejections below.` For superseded
  entries with nothing else recorded: `See D-0NN for the surviving decision.`
- **Rejected.** One list item per rejected alternative, in the original order:
  italicized name or short description, colon, the recorded reasons. Positions
  that fell *within* a surviving decision are items here, marked
  `*Superseded position — …:*`.
- **Divergence.** Only where the text explicitly contrasts an external tool's
  behavior with ours. External precedent cited *in support* of a decision is
  not divergence — it stays where it is. Ruled 2026-08-14: FlightCore, the
  predecessor in-house framework, does **not** count as external — contrast
  with it is lineage (§7.4 territory) and belongs in Rationale.

## Standing rules

1. **Historical entries are annotated, never rewritten.** An entry records what
   was settled when it was settled. Amend by adding a later entry and marking
   the earlier one superseded; do not retrofit a changed position into old
   prose.
2. **Entries ≤ D-121 keep pre-D-122 vocabulary.** The D-121/D-122 renames apply
   forward, not backward — rewording an old entry into current vocabulary
   falsifies the record of what was decided in what terms.
3. **Citations are written bare, then linked mechanically.** Write `§5.2`,
   `Appendix B`, `D-037`, ranges as `D-166–D-168`; run
   `julia docs/notes/design/tools/linkify.jl` and it converts them to reference
   links and regenerates the definitions block. Never hand-write anchors.
4. **The `## Index` is generated.** `linkify.jl` rebuilds it from the entry
   headings and their `**Status.**` lines, so it cannot drift. Do not hand-edit
   it; add the entry and re-run.
5. **Never resolve doubt by deleting.** Material that looks like it belongs in
   the spec rather than the log gets flagged, not dropped.
6. **The battery is the acceptance test.** `check_refs.jl`, `check_rows.jl`,
   `check_glossary.jl --strict`, and `linkify.jl` as a no-op on re-run.

## Pass C queue — per-entry prose rewrite

Pass A converted the table to entries structurally; Pass C is the optional,
demand-driven prose pass. Nothing here is scheduled — it is the standing list
of known-rough spots, from the Pass A findings adjudication (45/45 closed,
2026-08-14).

- **No-headline bundle** — D-045, D-048, D-050, D-056, D-058, D-059, D-062,
  D-067–D-072, and D-113's five-kinds-in-one-entry structure: Positions holding
  full semicolon-chains that want a real headline.
- **Provenance prefixes** in Positions ("Round-4 gaps finding 7:", …).
- **`D.3`/`D.6`/`D.8` glossary-token spelling** normalization (D-131).
- **Terse early entries** whose Rationale reads "Recorded only through the
  rejections below." — expansion candidates.

## Calibration examples

Four hand-converted entries spanning the difficulty range, kept as the
register's reference points. (Their originals are the retired table's rows;
the row-citation forms they show were swept to `D-nnn` in Pass B.)

### Example 1 — terse early entry (expansion into fields)

Original row 10:

> Structured immutable state over framework-owned flat `Vector{T}` | Mutable
> views (aliasing, silent missing-ẋ); fully structured/no flat vector (same
> machinery needed anyway, loses standard integrator interface)

Converted:

```markdown
### D-010 — Structure continuous state as immutable values over a flat backing

**Status.** ratified

**Position.** Structured immutable state over a framework-owned flat
`Vector{T}`.

**Spec.** §7.1

**Rationale.** Recorded only through the rejections below.

**Rejected.**
- *Mutable views:* aliasing; silent missing-ẋ.
- *Fully structured, no flat vector:* the same machinery is needed anyway;
  loses the standard integrator interface.
```

### Example 2 — superseded entry

Original row 7:

> Superseded by row 37: aggregation is by explicit summing junctions | Earlier
> design of record: reduce-ports with a canonical fold — reversed because it
> was the declaration vocabulary's last wrapper, a three-site census (all
> Newton–Euler, one library file) could not justify canonical-fold,
> multi-connection legality and identity-element machinery, and the aggregate
> wasn't even observable; Σ-junctions as default (arity/positional ceremony —
> objection dissolved by §11's loud declarations); contribution buses
> (invisible dataflow — verdict unchanged); contribution buses also import
> scoping rules and admit accidental contributions — a component contributes
> by being in scope rather than by being wired

Converted:

```markdown
### D-007 — Aggregation mechanism

**Status.** superseded → D-037

**Position.** Aggregation is by explicit summing junctions (D-037).

**Rationale.** See D-037 for the surviving decision.

**Rejected.**
- *Superseded position — reduce-ports with a canonical fold (the earlier
  design of record):* reversed because it was the declaration vocabulary's
  last wrapper, a three-site census (all Newton–Euler, one library file) could
  not justify canonical-fold, multi-connection legality and identity-element
  machinery, and the aggregate wasn't even observable.
- *Σ-junctions as default:* arity/positional ceremony — objection dissolved by
  §11's loud declarations.
- *Contribution buses:* invisible dataflow — verdict unchanged; they also
  import scoping rules and admit accidental contributions — a component
  contributes by being in scope rather than by being wired.
```

### Example 3 — finding entry (narrative stays in Rationale)

Original row 164 (abridged):

> Increment-2 finding (§11.1 extended): **declarations written in a local
> scope silently do not exist** — inside a `let`, a function body or a
> `@testset`, `h_x(::MyComp, (; x)) = …` binds a new *local* function rather
> than adding a method to the global `h_x` […] Mitigation adopted at the other
> end: **a component that declares nothing and defines no stage is a build
> error** […] | Documenting the caveat only (the failure is silent […]); a
> world-age-style check comparing `which(...).primary_world` against the
> caller's world (measures a mechanism that turned out not to be the cause)

Converted (abridged to the same degree):

```markdown
### D-164 — Reject components that declare nothing and define no stage

**Status.** ratified

**Position.** A component that declares nothing and defines no stage is a
build error; the authoring rule is that declarations live at module top level.
An increment-2 finding, extending §11.1.

**Spec.** §11.1

**Rationale.** **Declarations written in a local scope silently do not
exist** — inside a `let`, a function body or a `@testset`,
`h_x(::MyComp, (; x)) = …` binds a new *local* function rather than adding a
method to the global `h_x` […] Mitigation adopted at the other end: an inert
component is unwritable on purpose, which costs a line and catches the
misspelled-declaration family too […]

**Rejected.**
- *Documenting the caveat only:* the failure is silent and its symptom — an
  inert component — is exactly what the check names.
- *A world-age-style check comparing `which(...).primary_world` against the
  caller's world:* measures a mechanism that turned out not to be the cause.
```

Note the sorting: the adopted mitigation is the Position; the trap narrative,
the discovery story and the diagnosis dead-end are Rationale, wording intact.

### Example 4 — monster entry (multi-part ruling)

Original row 173 (abridged):

> **The discrete state letter `z` fuses into `x`** (§3.2, §5.2, §7.3, §11.2,
> §12.3, §12.5, §13.7, §14.1, §15.2, §15.5, Appendix C and all companions
> swept; row 56 untouched and load-bearing […]): `h_z`/`h_zu` → `h_x`/`h_xu`
> […] Motivation: `z` collided with the shift operator […] | Renaming `z` to
> another letter (fixes the collision only […]); `is_discrete(::C)::Bool`
> trait as stateless-leaf tier decider (derivable hence redundant […]); fusing
> `z` into `m` instead (semantically tempting […])

Converted (abridged to the same degree):

```markdown
### D-173 — Fuse the discrete state letter `z` into `x`

**Status.** ratified

**Position.** The discrete state letter `z` fuses into `x`: `h_z`/`h_zu` →
`h_x`/`h_xu`, `init_z` → `init_x`, `y_z` → `y_x`, bundle field `z` → `x`,
`z⁺ = g(z,u,t)` → `x⁺ = g(x,u,t)`; `f`/`g` stay distinct (flow and jump maps).

**Spec.** §3.2, §5.2, §7.3, §11.2, §12.3, §12.5, §13.7, §14.1, §15.2, §15.5,
Appendix C (all companions swept)

**Rationale.** `z` collided with the shift operator (`z⁻¹`, five spec sites —
§13.7's `UnitDelay` bullet had the letter meaning two things in three lines)
and was nonstandard besides […] D-056 untouched and load-bearing — it is what
makes the fusion safe […] Tier doctrine restated: […] (walking at `T` vs
pinned, D-166–D-167) […]

**Rejected.**
- *Renaming `z` to another letter:* fixes the collision only, keeps the
  doubled API; every candidate letter is taken […]
- *An `is_discrete(::C)::Bool` trait as stateless-leaf tier decider:*
  derivable hence redundant […] which D-039 refused for kind […]
- *Fusing `z` into `m` instead:* semantically tempting — discrete state is
  latched memory — but it kills the method fusion, breaks the `x[k]`/`x⁺`
  convention, and §3.2 already ruled the discrete tier's state IS its memory.
```

Note: the sweep list in the original's opening parenthetical becomes the
`**Spec.**` line; the multi-part tier-doctrine ruling stays in Rationale rather
than bloating the Position.
