# `framework_spec.md` rewrite plan

*Drafted 2026-08-13. Supersedes the abandoned `rewrite_test/` attempt and the
deleted `tools/restructure_plan.md`.*

## Goal and non-goals

**Goal.** Make the spec readable for its actual use case — a designer returning
cold to one section — without losing a single normative claim. Five defects to
fix, plus supporting improvements:

1. Overloaded sentences (rule + qualification + rejected alternative + citations
   in one breath); two named idioms to eliminate: the *possessive citation*
   ("the arity distinction … is §12.7's") and the *aphoristic lead-in* ("the
   absolute register anchors, and anchoring severs").
2. Rule-first ordering where background-first would serve; Part II (Execution)
   before Part III (Authoring and build) despite depending on it.
3. Adversarial rationale (rejections) inlined instead of living in the decision
   log.
4. Definite-article coinages ("the drain", "the epoch rule") used far from
   their definitions with no local gloss.
5. No links into the glossary (and no anchors in the glossary to link to).

**Non-goals.** No design changes of any kind. The decision log changes only to
*receive* migrated rationale and to track renumbered section references.
Companion documents change only where they cite spec section numbers.

## Ground rules

- **Content preservation is verified, not assumed.** Every style-pass section
  gets a normative-claim inventory extracted *before* rewriting and checked
  *after* (Phase 3 procedure). Cheap global invariants are scripted: the set of
  `row N` citations, the set of link labels, glossary coverage.
- **The spec stays ground truth at every commit.** No half-migrated states
  across commits; each commit is a complete, self-consistent pass over a stated
  scope (a section, a chapter, or one mechanical transformation).
- **Per-batch approval.** Claude proposes each batch as a diff; Miguel approves;
  then commit. One-sentence commit subjects, no attribution trailers.
  *Amended for P3.2 (2026-08-13):* per-section approval is replaced by
  end-of-pass review — Miguel reads the finished spec side by side with the
  untracked pre-P3.2 baseline copy (`framework_spec_p32_baseline.md`, the spec
  at `007484ea`; deleted at Phase 3 close). In-flight content preservation
  moves to an independent verifier agent per section (P3.2 step 5); decisions
  that would have surfaced at approval time queue into chapter-end digests.
- **Sequencing: mechanical → migration → style → structure.** The reorder comes
  last because it renumbers sections and touches every citing document; doing
  it over stable content makes it one mechanical move instead of a moving
  target.

## Phase 0 — Preparation

**P0.1 Style guide.** Write `tools/spec_style.md` (short, ~1 page): the rules
the style pass enforces, so the pass stays consistent across sessions and
context windows. Contents:

- One burden per sentence: a sentence states a rule, *or* qualifies one, *or*
  cites — not all three. Target register: rigorous but pedagogical, economical
  but not terse.
- Topic sentences are findable claims, not aphorisms.
- Ban the possessive-citation idiom; spell out what the cited section provides.
- Coined terms: gloss (5–10 words) + glossary link on first use *per section*;
  bare term afterward within the section.
- Section template: **context** (what problem this section solves) →
  **rule** (normative, visually marked) → **mechanism / example** (code sketch
  where prose runs long) → **consequences and pointers**.
- Labeling convention for the template roles (proposed: sparing bold
  **Rule.** / **Why.** / **Example.** markers — see open decision 5).
- Constructive vs. adversarial rationale: constructive (why the rule is shaped
  as it is) stays in the spec under *Why*; adversarial (why alternative X
  loses) lives in the decision log, cited as "(row N)".

**P0.2 Coinage inventory.** Build the list of coined terms: start from the
glossary's bold terms, then sweep the spec for definite-article terms not in
the glossary. Terms found missing get glossary entries added (flagged for
approval — these are the one place Phase 1 adds content).

**P0.3 Audit and adapt the existing tooling.** `tools/` already has
`check_refs.jl`, `linkify.jl`, `sweep_refs.jl`, `slugs.jl` from the abandoned
attempt. Audit them; keep what fits, rewrite what doesn't. Required checks
(run before every commit of every phase):

- *Link check*: every reference-style label is defined exactly once; every
  defined label is used; every anchor target exists in the file.
- *Row-citation set*: the multiset of `row N` citations in the spec, diffed
  against the pre-phase baseline; any change must be explained by the batch.
- *Glossary coverage*: every inventory term has an anchor; every glossary
  anchor is linked from at least one body section.

**Done when:** style guide approved, inventory list approved, checks run green
on the untouched spec (baseline established).

## Phase 1 — Mechanical passes (content-neutral)

**P1.1 Reference-style links.** Convert every inline
`[§8.6](#86-event-iteration-…)` to `[§8.6][s8-6]` with a link-definition table
at file bottom. Scriptable end-to-end (adapt `linkify.jl`); one commit for the
whole file. This is what makes Phase 4's renumbering tractable — cross-ref
updates become edits to one table.

**P1.2 Glossary anchors.** Give every glossary entry an anchor. Proposed
scheme: `<a id="g-drain"></a>**drain** — …` (explicit anchors keep entries as
bold terms rather than promoting them to headings and bloating the TOC).
Scriptable from the glossary's own bold-term structure; one commit.

**P1.3 Glossary links and first-use glosses.** Semi-manual, batched per Part:
at each coined term's first use in a section, add the gloss and the
`[drain][g-drain]` link. Overlap warning: Phase 3 will reword some of these
glosses — accepted, they're cheap, and the readability win is immediate for
ongoing consultation of the spec.

**Done when:** checks green; spec renders identically in content, differs only
in link plumbing, anchors, glosses.

## Phase 2 — Rationale migration

Scope: the ~76 inline "rejected" passages. Procedure per site:

1. Classify: *adversarial* (litigates an alternative) or *constructive*
   (explains the rule's shape). Constructive stays.
2. For adversarial: find the row it cites (or should cite). Compare depth.
   - Row already carries the argument (the common case — spot-checks show rows
     17–24 are already dense): reduce the spec text to "(row N)".
   - Spec's inline version is richer: enrich the row first (approval batch on
     `framework_decisions.md`), then reduce the spec text.
   - No row exists: propose a new row (numbers append-only), then reduce.
3. Batch per chapter; run the row-citation check — the citation *set* must not
   shrink (every deleted argument leaves its "(row N)" behind).

**Done when:** no adversarial litigation remains inline; decision log is the
single home of "what we rejected and why"; checks green.

## Phase 3 — Style pass (the bulk)

**P3.1 Triage.** Skim all ~70 sections and rate each: *ok* (leave alone),
*light* (sentence surgery only), *full* (rewrite to template). Produces the
work list and a realistic effort picture before committing to it. Expected
worst offenders from sampling: §8.5–8.6, §9.3–9.5, parts of §14.

**P3.2 Per-section procedure** (the content-preservation core):

1. Extract the section's normative claims as a bullet inventory (untracked
   `p32_inventories/` beside the spec, one file per section; deleted at Phase 3
   close — kept on disk so the verifier and the end-of-pass review can consult
   it).
2. Rewrite to the template: context first, rule marked, mechanism with a code
   sketch where the section is prose-heavy (candidate sketches: a
   `sample_times` declaration next to its compiled `(D, Φ)` pairs for §8.5; a
   staging/drain timeline for §9.4).
3. Re-verify every inventory claim against the rewrite; diff the section's row
   citations and glossary links.
4. Commit (one commit per full section or per chapter light-cluster), checks
   green. No new section numbers during P3.2: proposed splits (§9.2, §10.4)
   are realized as `####` sub-headings; real splits, if still wanted, wait for
   Phase 4.
5. An independent verifier agent (fresh context) compares the baseline section
   against the rewrite using the claim inventory; findings are fixed (follow-up
   commit) before the next section launches. Decisions needing Miguel — R1
   borderline keep/reduce calls, unclear-meaning flags, structural proposals —
   queue into a chapter-end digest rather than blocking per section.

**P3.3 Work order.** Recommended: chapter-linear within a part, starting with
Part II (worst offender, and the part most consulted during kernel
increments), then Part III, I, IV, V, appendices. See open decision 4.

**Done when:** every *light*/*full* section rewritten and verified; glossary
entries re-checked against their (possibly reworded) owning sections — the
owning-section-wins precedence makes stale glossary compression a real risk
here.

## Phase 4 — Structural reorder

*P4.1 and P4.2 applied 2026-08-15 (`7217b632`). Chapters were permuted*
**8→10, 9→11, 10→12, 11→8, 12→9**, *subsection numbers preserved (old §8.6 =
new §10.6), and the parts renumbered so that Authoring and build is Part II
and Execution is Part III. This is the old→new map open decision 2 asks to be
kept for archaeology; the one-shot scripts that applied it were deleted once
spent. **Every § number in this plan is pre-swap — apply the map to read it.**
In particular, P4.2's "§8–10 ↔ §11–12" below describes the move itself, not a
citation to be updated. The corpus (spec plus the eight companions) and
`tools/gloss_table.md` are fully swept.*

*P4.3 applied 2026-08-15 (`8ddf5643`..`493f602d`). Each Part now opens with a
roadmap. **No cross-section moves survived**: a fresh survey of all 81 section
openings found Phase 3 had already absorbed them, and the one candidate — §9.6,
which announces itself as an aside inside the build pipeline — was left in place
rather than pay another renumbering sweep for a grounding aside. The Phase 3
triage files were untracked and are gone, so that survey, not the original flag
list, is the record. Two things were fixed that the phase did not ask for: the
appendices were lifted out of Part V into their own back matter, since an index,
an API synopsis, a kind set and a glossary are not grounding; and the trim-commit
account was deduplicated between §14.5 and §14.8.*

**P4.1 The part swap.** Move Part III (Authoring and build) before Part II
(Execution). Stopped-sim services *stay* after Execution: §14.5 defines
boundary zero in terms of §8.6's boundary machinery, and the services'
semantics are parameterized by run semantics — chronology loses to dependency
direction here. Recorded as the rationale so the question stays settled.

**P4.2 Renumbering fallout.** Sections renumber (§8–10 ↔ §11–12). Update, in
order: the spec's own headings and TOC; the link-definition table (Phase 1's
payoff — body text barely changes); the decision log's ~140 § citations;
`sample_time_proposal.md` (~32); the walkthrough docs
(`localization_validation_walkthrough.md`, `event_visibility_walkthrough.md`,
`inbound_periphery_walkthrough.md`, `frozen_discrete_walkthrough.md`,
`trim_environment_walkthrough.md`) and `framework_extensions.md`. Scriptable
as an old→new section-number map applied everywhere (adapt `sweep_refs.jl`);
one commit, plus a manual skim for prose that *describes* order ("as Part II
established…").

**P4.3 Within-part ordering and roadmaps.** Fix section-internal
rule-before-background orderings that Phase 3 flagged but couldn't fix locally
(most background-first fixes land in Phase 3; this catches cross-section
moves). Add a short roadmap paragraph at the top of each Part: what it covers,
in what order, what it assumes from earlier parts.

**Done when:** every § citation in every document resolves to the section it
meant before the move, and each Part opens with a roadmap naming what it
assumes from earlier parts.

Forward-reference *count* is deliberately not the criterion. Measured at the
swap: cross-part citation traffic between the two swapped parts runs 30 one
way against 31 the other, so the move trades one set of forward references for
a set the same size. The swap stands on conceptual dependency instead — a
reader meets how a model is authored before how it executes — and the P4.3
roadmaps are what carry the forward references that remain.

## Phase 5 — Final verification

*Applied 2026-08-15 (`5100f5a9`) — the rewrite is closed. The full linear
read found no structural, ordering or vocabulary defects; what survived was
citation-language staleness from Pass B (§1's and the glossary's descriptions
of the decision log still said "row N"; one line-break-split "rows 106–107"
citation had escaped both the sweep and `check_rows.jl`'s line-by-line scan,
now hardened) plus three wording nits. The battery ran green throughout; the
walkthroughs quote no spec text verbatim, so the quote check was vacuous;
`rewrite_test/` was already deleted 2026-08-13; `glossarize.jl` was the one
one-shot tool retired (`slugs.jl` is a live include, `sweep_rows.jl` is
check_rows' named fixer); the row-citation baseline was regrown 113 → 121.*

- Full linear read-through of the spec (fresh session, reader's stance).
- All Phase 0 checks green; row-citation set identical to the Phase 2 exit
  baseline; glossary bidirectional check (every entry linked from body, every
  body coinage in glossary).
- Confirm companion walkthroughs' quoted spec passages (if any quote verbatim)
  still match.
- Delete `rewrite_test/` (see open decision 3) and retire any `tools/` scripts
  that were one-shot.

## Open decisions (Miguel)

1. **Single file vs. split per part.** Recommend keeping the single file:
   whole-document search is the dominant workflow, and the link plumbing
   assumes one anchor namespace. Revisit only if the file becomes hostile to
   editors.
2. **Renumbering policy.** Recommend plain renumbering in Phase 4 (numbers
   stay meaningful and ordered) over stable section IDs (permanent confusion
   between order and identity). The old→new map is kept in the Phase 4 commit
   message-adjacent notes for archaeology.
3. **Fate of `rewrite_test/`.** Recommend deleting it in Phase 5 once this
   plan has delivered; it was to be ignored, not mined.
4. **Phase 3 work order.** Recommend Part II first (worst text, highest
   consultation value during kernel work) rather than strictly linear.
5. **Template markers.** Recommend literal sparing bold markers (**Rule.**,
   **Why.**, **Example.**) over an unmarked convention — skimmability is the
   point, and unmarked conventions drift.

## Effort shape (rough)

- Phase 0–1: one to two sessions (mostly scripted).
- Phase 2: one to two sessions (~76 sites, batched per chapter).
- Phase 3: the bulk — several sessions; triage (P3.1) gives the real number.
- Phase 4: one session (scripted sweep + manual skim).
- Phase 5: one session.

Kernel increment 3 stays blocked behind nothing here in principle — the spec
is normatively unchanged throughout — but Phases 2–4 churn the text badly
enough that interleaving spec-citing design work mid-phase is not worth the
merge pain. Finish a phase, then interleave if needed.
