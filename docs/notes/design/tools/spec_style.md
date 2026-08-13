# Spec style guide

The rules the `framework_spec.md` rewrite enforces (plan: `spec_rewrite_plan.md`,
Phase 3). Consulted at the top of every rewrite batch so the pass stays
consistent across sessions. The target register: **rigorous but pedagogical,
economical but not terse**. Growth in word count is acceptable; a norm buried
in pedagogy is not.

## Sentences

- **One burden per sentence.** A sentence states a rule, *or* qualifies it,
  *or* explains it, *or* cites. Prefer two short sentences over one loaded one.
- At most one parenthetical or dash-delimited clause per sentence; never
  nested. Em-dashes mark genuine appositions, not a stack of pending clauses.
- **Topic sentences are findable claims, not aphorisms.** A returning reader
  skims bold lead-ins; each must state its paragraph's content.
  - Bad: "The absolute register anchors, and anchoring severs."
  - Good: "An `Absolute` entry detaches its child from the enclosing scope's
    grid."
- **No possessive citations.** Spell out what the cited section provides when
  the dependency is load-bearing; a bare "([§12.7][s12-7])" is fine as a
  pointer.
  - Bad: "The arity distinction that carries it into the phase-body surface —
    zero-arg interior, tick-indexed boundary — is §12.7's."
  - Good: "The two sweep variants surface in the phase-body signatures:
    interior bodies take no arguments, boundary bodies take the tick index
    ([§12.7][s12-7])."

## Terms

- **First use per section: gloss + glossary link; bare term afterward.** The
  gloss is 5–10 words, parenthesized: "the [drain][g-drain] (the frame-top
  swap that publishes staged device inputs into the root slots)". Every
  section must be locally re-readable by a reader arriving cold.
- Coinages are governed by `tools/coinage_inventory.md`: a term used in three
  or more places, or far from its definition, gets a glossary entry; a local
  one- or two-use coinage gets an inline gloss only. New coinages introduced
  during the rewrite join the inventory.

## Section template

Content order within a section, mirroring how a reader reasons into the
problem:

1. **Context** — what problem this section exists to solve; the background a
   cold reader needs. Comes first, before any rule.
2. **Rule** — the normative content, visually marked.
3. **Mechanism / example** — how it plays out; a small code sketch wherever
   prose runs long (sketches may be simplified or incomplete).
4. **Consequences and pointers** — implications, interactions, where the
   machinery lives.

Markers: sparing bold **Rule.** / **Why.** / **Example.** labels. Mark every
rule a skimmer must be able to find; do not label every paragraph — ubiquity
would destroy the signal. *Why* holds constructive rationale only.

## Rationale

- **Constructive rationale stays** (why the rule is shaped as it is — what
  makes it believable and memorable). It lives under *Why*.
- **Adversarial rationale goes** (why alternative X loses). It lives in
  `framework_decisions.md`, cited as "(row N)". Before deleting an inline
  argument, verify the row carries it; enrich the row first if the inline
  version is richer (Phase 2 procedure).

## Content preservation

Style edits never change what the spec norms. Every Phase 3 section rewrite
follows the claim-inventory procedure (plan, P3.2): extract normative claims
before, verify each survives after, diff row citations and glossary links.
When a sentence resists rewriting because its meaning is unclear, that is a
finding, not an obstacle — flag it for discussion instead of guessing.
