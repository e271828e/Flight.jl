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
  to Rationale — a finding's narrative, its discovery story and its dead ends
  are Rationale, not Position. Where an entry genuinely bundles parallel rulings
  with no distinct headline, content preservation wins and Position holds the
  full chain — see the Pass C queue.
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

Pass A converted the table to entries structurally; Pass C is the optional
prose pass. From the Pass A findings adjudication (45/45 closed, 2026-08-14),
sized against the log 2026-08-15.

**Sequencing: after Phase 4 of the spec rewrite, not before.** Phase 4's
renumber sweep rewrites the `§` citations inside these same entries, so doing
prose first means reviewing them twice.

**Prerequisite: rebuild the content-preservation audit.** Pass A's
`audit_fragments.jl` — a per-entry token-multiset comparison catching a dropped
citation, code span or named diagnostic — was deleted with the scratch
directory. Any batch claiming to preserve content needs it back first (~60
lines, about half the work of `sweep_rows.jl`).

Items 1–3 are ~21 distinct entries, one working session including the audit
rebuild. Items 1 and 2 are conveniently clustered — D-067–D-072 and
D-135–D-140 are each contiguous runs.

1. **No-headline bundle** — 14 entries: D-045, D-048, D-050, D-056, D-058,
   D-059, D-062, D-067–D-072, and D-113's five-kinds-in-one-entry structure.
   1,208 Position words between them, mean 86 against a corpus mean of 61; the
   worst are D-113 (143 words, 8 semicolons) and D-071 (132, 9). The work is to
   extract a one-sentence headline and move the semicolon chain into Rationale,
   every citation and code span preserved. **Note the tension:** the Position
   field rule above currently blesses these chains. Moving the chain to
   Rationale honours it; if this item lands, revisit that clause so the guide
   and the log agree.
2. **Provenance prefixes** in Positions — 7 entries: D-131 and D-135–D-140,
   each opening "Round-4 gaps finding 5, amending D-052: …". **Ruled
   2026-08-15: delete the prefix outright**, not relocate it to Rationale. The
   Position keeps only the ruling; where the prefix also carries a substantive
   clause ("amending D-052", "closing D-133's §13.2 deferral"), that clause
   survives — it states what the decision does, not where it came from. The
   round provenance stays recoverable from this file's git history.
3. **Bare `D.n` glossary citations** — 9 tokens across 4 entries (D-131 and
   neighbours), spelled `D.3`/`D.6`/`D.8` rather than `§D.3`. A real defect,
   not cosmetics: `check_refs.jl` only recognises the `§D.n` form, so these
   nine are unvalidated and unlinked. Small enough to fold into Phase 4's
   mechanical sweeps instead of running it here.

**Removed from the queue (2026-08-15): "terse early entries as expansion
candidates."** Measured at 71 entries — 37% of the log — whose Rationale reads
`Recorded only through the rejections below.` They read that way because the
original table rows recorded nothing else; the reasoning survives only as the
rejections. Expanding them means sourcing rationale that was never written, out
of the cited spec sections (which state mechanism, not reasoning) or the git
history of `framework_design.md` — historical reconstruction, not prose work,
and in tension with standing rule 1, since writing a rationale into an old
entry today is retrofitting. If ever revisited, first survey which entries have
a real documented source; do not expand wholesale.
