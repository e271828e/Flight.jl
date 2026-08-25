# Decision log style guide

The rules `framework_decisions.md` adheres to. Distilled from the brief that
converted the log's retired table form into entries; consulted whenever an
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
bare (§5.2, D-037). Where the entry bundles parallel rulings, a headline
sentence and then one bullet per ruling.

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
  are Rationale, not Position. Where an entry bundles parallel rulings, the
  headline sentence is followed by a bulleted list, one ruling per bullet: the
  rulings stay in the ruling field, and the reader gets a headline instead of a
  semicolon chain. No provenance — which review round or finding raised a
  decision belongs to this file's git history, not to the Position.
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
2. **Entries keep the vocabulary of their day.** The D-121/D-122 renames and
   every later one (D-172, D-185–D-187, D-194–D-196, the 2026-08-15
   "generic-holding" / "trial evaluation" renames) apply forward, not
   backward — a rename sweeps the spec and its companions, never the log — and
   rewording an old entry into current vocabulary falsifies the record of what
   was decided in what terms; rule 7 licenses prose rewrites, not vocabulary
   modernization. The current name is always one spec read away, or in the
   superseding entry.
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
7. **A prose rewrite is audited, not trusted.** `audit_fragments.jl [rev]`
   compares every entry's multiset of code spans, inline math, citations and
   named identifiers against a git revision — losing one is fatal. Rewording an
   entry is licensed by that check, not by care.
8. **Do not expand the terse early entries.** 71 entries — 37% of the log —
   have `Recorded only through the rejections below.` for a Rationale, because
   the original table rows recorded nothing else and the reasoning survives
   only as the rejections. Writing one today means sourcing rationale that was
   never written, out of the cited spec sections (which state mechanism, not
   reasoning) or out of `framework_design.md`'s git history: historical
   reconstruction, not prose work, and retrofitting under rule 1. If ever
   revisited, first survey which entries have a real documented source; do not
   expand wholesale.
