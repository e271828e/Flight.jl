# framework_spec.md restructure — normative plan

Executed 2026-08-03 in three phases. This file is the single source of truth for
the target outline, the section-number map, and the conventions every phase
follows. It is committed for audit and can be deleted once the restructure is
accepted.

## Files

Working directory: `/Users/miguel/.julia/dev/Flight.jl-core-redesign-2/docs/notes/design/`

- **Restructured:** `framework_spec.md`.
- **Reference-swept only** (their `§N.M` citations rewritten to the new
  numbering, nothing else): `framework_decisions.md`,
  `event_visibility_walkthrough.md`, `inbound_periphery_walkthrough.md`,
  `review3_2026-07-30_gaps.md`, `review3_2026-07-30_improvements.md`,
  `review3_2026-07-30_inconsistencies.md`.
- **Never touched, never `git add`-ed:** `framework_design.md` (historical,
  self-referential), `framework_spec_pre_restructure.md` (comparison snapshot),
  `framework_spec.pdf`. The three `review3_*` files are deliberately untracked:
  rewrite them in place but do not `git add` them.

## Target outline (after phase 2)

Parts are new `#`-level headings (the document title remains the first `#`).
Chapters stay `##`, sections `###`. "← old X" gives the source; everything not
listed as merged/split moves verbatim.

```
# Part I — Foundations
## 1. Purpose and method            ← old preamble + §1 (deduplicated) + old §18's pointer
## 2. Formalism                     (2.1, 2.2 unchanged)
## 3. Component taxonomy            (3.1–3.3 unchanged)
## 4. Ports and signals             (4.1–4.3 unchanged)
### 4.4 Function-valued signals: environment access   ← old §7
## 5. Evaluation order and feedthrough                (5.1–5.4 unchanged)
### 5.5 Diagnostics: feedthrough tracing              ← old §10
## 6. Composition: connections, aggregation and hierarchy   (new lede, 1–3 sentences)
### 6.1 Connections and hierarchy                     ← old §8
### 6.2 Aggregation: explicit summing junctions       ← old §6
## 7. State and data representation ← old §9 (7.1–7.5 ← 9.1–9.5)

# Part II — Execution
## 8. Time and execution            ← old §11 (8.1–8.7 ← 11.1–11.7)
## 9. Runtime periphery: the data plane        ← old §12 lede + 12.1–12.5 → 9.1–9.5
## 10. Runtime periphery: lifecycle and orchestration  ← old 12.6–12.12 → 10.1–10.7
                                                (new lede, 1–3 sentences)

# Part III — Authoring and build
## 11. The declaration layer: components and assemblies   ← old §13 (11.1–11.8)
## 12. The build pipeline           ← old §14 (12.1–12.7)

# Part IV — Failure and services
## 13. Error discipline             ← old §15 (13.1–13.7)
## 14. Stopped-sim services         ← old §16 (14.1–14.10)

# Part V — Grounding
## 15. Case studies                 ← old §17 (15.1–15.5)
## 16. Open axes                    ← old §19

## Appendix A–D                     (unchanged, incl. D.1–D.10)
```

Old §18 (Decision log) is dissolved: its content is a pointer to
`framework_decisions.md` already present in the introduction; make sure the
merged §1 retains exactly one copy of that pointer and of the formerly
duplicated preamble/§1 opening paragraph. Chapter titles keep their old
wording except where noted above; subsection titles are unchanged in phase 2.

## Reference map (phase 2)

Applied to every `§`-citation in `framework_spec.md` **and** the six swept
companion files, including citations inside fenced code blocks and inside
headings (e.g. old "17.3 … the §12 staging shapes" → "15.3 … the §9 staging
shapes").

| old | new | | old | new |
|---|---|---|---|---|
| §1–§5 (all) | unchanged | | §12.6–§12.12 | §10.1–§10.7 |
| §6 | §6.2 | | §13.x | §11.x |
| §7 | §4.4 | | §14.x | §12.x |
| §8 | §6.1 | | §15.x | §13.x |
| §9.x, §9 | §7.x, §7 | | §16.x | §14.x |
| §10 | §5.5 | | §17.x | §15.x |
| §11.x, §11 | §8.x, §8 | | §18 | see below |
| §12.1–§12.5 | §9.1–§9.5 | | §19 | §16 |

Appendix citations unchanged. Bare `§12` (5 occurrences in the spec, 1 in
`review3_2026-07-30_inconsistencies.md`) is adjudicated per occurrence: prefer
the specific new chapter (§9 or §10) the context means; only where the whole
periphery is genuinely meant, cite both ("§9, §10"). The single `§18` citation
is resolved in context (reword to name `framework_decisions.md` or point at
§1). Record every adjudication in the final report.

**Mechanics are non-negotiable:** the rewrite is a single table-driven pass —
match each token `§N(.M)?` once and map it through the table (old §13 → new
§11 while old §11 → new §8: sequential search-and-replace would double-map).
Verify per file: total `§`-token count unchanged (modulo the documented §12/§18
adjudications), and every cited section exists in the new outline.

## Phase 3 — interior splits (renumbering on top of phase 2)

- **§5.2** (194 lines) → split at the natural seam into `### 5.2` (the
  two-stage output contract) and a new `### 5.3` (structural feedthrough
  mechanics: the view-set survey with its code blocks). Old 5.3→5.4, 5.4→5.5,
  5.5→5.6. Re-adjudicate existing §5.2 citations: re-point those that clearly
  mean the view-set material; leave ambiguous ones at §5.2.
- **§9.3** (old §12.3, 318 lines) → split into ~3 sections along its internal
  structure (root input slots / per-device staging / the input trace — exact
  seams from the text). Subsequent sections shift (devices, GUI write path).
  Re-adjudicate the ~90 citations the same way (default: first section).
- **Unnumbered `####` headings** (no numbering impact, no citation impact)
  inside: the devices section (old §12.4), §11.2 (declaration inventory — one
  heading per declaration group), and case studies §15.4, §15.5.
- Same sweep + verification discipline as phase 2, across all seven files.

## Phase 4 — navigation

- Convert `§N(.M)?` and `Appendix X` citations in the spec's **prose** (never
  inside fenced code blocks or inline code spans) to markdown links against
  GitHub-style auto-slugs of the actual headings (lowercase, punctuation
  stripped, spaces→hyphens). Companion files link cross-file:
  `framework_spec.md#slug` — except `framework_decisions.md`, deliberately
  left plain: it is one large table whose rows are edited constantly, and
  links would clutter every cell's source (exemption honored by
  `linkify.jl`). `row N` citations stay plain text everywhere (the decisions
  table has no per-row anchors; row numbers are stable by contract).
- Insert a generated, linked table of contents (`## Contents`, unnumbered)
  between the document's front matter and the `# Part I` heading, listing
  parts, chapters and sections.
- Scripts live in `docs/notes/design/tools/` and are committed: the one-shot
  renumber script(s), `linkify.jl`/TOC generator, and `check_refs.jl` — the
  durable checker that extracts headings and verifies (a) every `§`-citation in
  every design file resolves to an existing heading, (b) every markdown anchor
  link resolves. Julia, plain stdlib.

## Conventions (all phases)

- Moves are **verbatim**: body text relocates unchanged except the citation
  sweep, the §1 deduplication, and explicitly allowed stitching (chapter ledes
  for new §6 and §10, ≤3 sentences each; minimal connective edits at merge
  seams). No rewording, no "improvements".
- Verify after each phase: extracted outline matches this plan; per-file
  citation-count conservation; checker passes; `git diff --stat` reviewed.
- One commit per logical step (phase 2 atomic: moves + renumber + sweep).
  Commit messages: a single concise sentence, subject line only, **no body, no
  Co-Authored-By, no session links, no attribution of any kind**.
