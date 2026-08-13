# Coinage inventory

The working list of coined terms for the spec rewrite (plan:
`spec_rewrite_plan.md`, P0.2). It drives Phase 1.3 (first-use gloss + glossary
link per section) and feeds `tools/check_glossary.jl`. A working document: the
initial list comes from the Appendix D extraction plus a heuristic
head-noun sweep of the body ("the X rule/sweep/law/…"), and it grows as the
Phase 2/3 reading uncovers coinages the sweep missed.

**Triage rule.** A coinage used in three or more places, or used far from its
definition, gets a glossary entry (bucket C until approved). A one- or two-use
local coinage gets an inline gloss at the use site only (bucket D). Aliases of
an existing entry link to that entry (bucket B).

## A. Canonical terms

The 158 bold entries of Appendix D. Not duplicated here —
`check_glossary.jl` extracts them live from the spec, so the glossary remains
the single home. Multi-alias entries ("stage function / two-stage outputs")
split on "/" for matching.

## B. Aliases — link to the owning entry

Coined spellings that name an existing glossary concept. Phase 1.3 links them
to the owning entry's anchor; no new entries needed.

| Coinage in body | Owning entry |
|---|---|
| frame-top drain, un-pause drain | drain |
| frozen roster | roster |
| stepper seam, backend seam | seam |
| honest priors | prior |
| epoch rule | input epoch (verify the rule text lives there when linking) |
| levels doctrine | the levels-never-deltas doctrine — verify owning entry; row 24 carries it |
| interior sweep, boundary sweep | sweep (both variants are bold-defined inside that one entry; Phase 1.2 may give such in-entry sub-terms their own anchors) |
| bundle law | bundle (the entry states the law: a name is present iff the store exists) |

Compositional uses of the D.10 meta-vocabulary ("didactic register",
"inspection register", "what-if register", "device contract", "face
contract", "authoring contract", "mid-run mutation doctrine") are ordinary
compositions of *register* / *contract* / *doctrine* with a qualifier, not
separate coinages: link the meta-term or the owning section, whichever the
context wants.

## C. Proposed new glossary entries (approval pending; written in Phase 1)

| Term | Uses | Notes |
|---|---|---|
| one-hop law | 5 | §5.2; never defined in glossary prose |
| handler return law | 3 | §5.2; never defined in glossary prose |

Both are §5.2 hand-off laws, siblings of the bundle law (bucket B). Whether
they join the `bundle` / `stage function` entries as bold sub-terms with their
own anchors, or get short entries of their own, is decided when the entries
are drafted in Phase 1.

## D. Low-frequency coinages — inline gloss only

One or two uses each; promote to bucket C only if usage grows during the
rewrite: arrival sweep, post-transition sweep, sampling seam, getfield walk,
branch-shape rule, locality law, active-widget contract, by-allocation
register.
