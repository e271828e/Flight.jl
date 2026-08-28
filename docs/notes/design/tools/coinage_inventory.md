# Coinage inventory

The working list of the spec's coined terms. It governs the first-use gloss +
glossary link per section and feeds `tools/check_glossary.jl`. A working
document: the initial list came from the Appendix D extraction plus a heuristic
head-noun sweep of the body ("the X rule/sweep/law/…"), and it grows as later
reading uncovers coinages the sweep missed.

**Triage rule.** A coinage used in three or more places, or used far from its
definition, gets a glossary entry (bucket C until approved). A one- or two-use
local coinage gets an inline gloss at the use site only (bucket D). Aliases of
an existing entry link to that entry (bucket B).

## A. Canonical terms

The bold entries of Appendix D. Not duplicated here —
`check_glossary.jl` extracts them live from the spec, so the glossary remains
the single home. Multi-alias entries ("stage function / two-stage outputs")
split on "/" for matching.

## B. Aliases — link to the owning entry

Coined spellings that name an existing glossary concept. They link to the
owning entry's anchor; no new entries needed.

| Coinage in body | Owning entry |
|---|---|
| frame-top drain, un-pause drain | drain |
| frozen roster | roster |
| stepper seam, backend seam | seam |
| honest priors | prior |
| epoch rule | input epoch (verify the rule text lives there when linking) |
| levels doctrine | the levels-never-deltas doctrine — verify owning entry; D-024 carries it |
| interior sweep, boundary sweep | sweep (both variants are bold-defined inside that one entry; such in-entry sub-terms could be given their own anchors) |
| bundle law | bundle (the entry states the law: a name is present iff the store exists) |
| name-transparent container | container children (the entry states the naming opt-out `transparent_container` declares) |

Compositional uses of the §D.10 meta-vocabulary ("didactic register",
"inspection register", "what-if register", "device contract", "face
contract", "authoring contract", "mid-run mutation doctrine") are ordinary
compositions of *register* / *contract* / *doctrine* with a qualifier, not
separate coinages: link the meta-term or the owning section, whichever the
context wants.

## C. Promoted coinages — settled

| Term | Uses | Notes |
|---|---|---|
| handler return law | 3 | §5.2; never defined in glossary prose |
| termination record | 13 | defined in §13.5; also §11.8, §12.1, §12.4, §12.6, §13.2, §13.6 and Appendix C |
| scratch | 45 | promoted to an Appendix D entry (`g-scratch`) with D-213; bucket A from here on |

The handler return law is a §5.2 hand-off law, sibling of the bundle law
(bucket B). Settled the third way: it got its own `####` sub-heading in §5.2,
and the body cites the section rather than a glossary anchor.

The termination record went the ordinary way: promoted by D-203 to an
Appendix D entry (`g-termination-record`), so it is bucket A from here on and
this row is the record of the promotion.

## D. Low-frequency coinages — inline gloss only

One or two uses each; promote to bucket C only if usage grows: arrival sweep, post-transition sweep, sampling seam, getfield walk,
branch-shape rule, locality law, active-widget contract, by-allocation
register.
