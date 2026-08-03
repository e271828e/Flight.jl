#!/usr/bin/env julia
#
# Phase-2 restructure of framework_spec.md (see tools/restructure_plan.md).
#
# One-shot. Run AFTER tools/sweep_refs.jl, which rewrites the citations; this
# script only moves blocks, renumbers headings and inserts the five `# Part`
# headings, the two new chapter ledes and the merged §1.
#
# Body text is relocated verbatim. The only prose written here is the merged §1
# (assembled from the old preamble, old §1 and old §18's pointer, with the
# duplicated opening sentence dropped) and the two new ledes.
#
# Run from anywhere:  julia docs/notes/design/tools/phase2_restructure.jl

const SPEC = normpath(joinpath(@__DIR__, "..", "framework_spec.md"))

lines = readlines(SPEC)

"Index of the single line starting with `prefix`."
function anchor(prefix)
    idx = findall(startswith(prefix), lines)
    length(idx) == 1 || error("$(length(idx)) lines start with $(repr(prefix))")
    return only(idx)
end

"Lines `a:b` with leading/trailing blanks and horizontal rules removed."
function block(a, b)
    out = lines[a:b]
    while !isempty(out) && (isempty(strip(out[end])) || strip(out[end]) == "---")
        pop!(out)
    end
    while !isempty(out) && isempty(strip(out[1]))
        popfirst!(out)
    end
    return out
end

"The chapter block starting at heading `prefix` and ending before `nextprefix`."
chapter(prefix, nextprefix) = block(anchor(prefix), anchor(nextprefix) - 1)

"Replace a block's own heading line with `heading`."
retitle(blk, heading) = [heading; blk[2:end]]

"Renumber `### old.N` subsection headings to `### new.(N + shift)`."
function renumber(blk, old, new; shift = 0)
    map(blk) do line
        m = match(Regex("^### $old\\.(\\d+) (.*)\$"), line)
        m === nothing ? line : "### $new.$(parse(Int, m[1]) + shift) $(m[2])"
    end
end

# ---------------------------------------------------------------- old blocks

title    = lines[1]
old1     = chapter("## 1. ", "## 2. ")
old2     = chapter("## 2. ", "## 3. ")
old3     = chapter("## 3. ", "## 4. ")
old4     = chapter("## 4. ", "## 5. ")
old5     = chapter("## 5. ", "## 6. ")
old6     = chapter("## 6. ", "## 7. ")
old7     = chapter("## 7. ", "## 8. ")
old8     = chapter("## 8. ", "## 9. ")
old9     = chapter("## 9. ", "## 10. ")
old10    = chapter("## 10. ", "## 11. ")
old11    = chapter("## 11. ", "## 12. ")
old12ld  = block(anchor("## 12. "), anchor("### 12.1 ") - 1)
old12a   = block(anchor("### 12.1 "), anchor("### 12.6 ") - 1)   # 12.1–12.5
old12b   = block(anchor("### 12.6 "), anchor("## 13. ") - 1)     # 12.6–12.12
old13    = chapter("## 13. ", "## 14. ")
old14    = chapter("## 14. ", "## 15. ")
old15    = chapter("## 15. ", "## 16. ")
old16    = chapter("## 16. ", "## 17. ")
old17    = chapter("## 17. ", "## 18. ")
old19    = chapter("## 19. ", "## Appendix A.")
appendix = block(anchor("## Appendix A."), length(lines))

# ------------------------------------------------------------- merged §1
#
# Paragraph 1 = the old preamble's opening paragraph verbatim, followed by the
# second sentence of old §1's opening paragraph (their shared first sentence is
# kept once). Then old §1's ground rules verbatim, the preamble's status
# paragraph verbatim, and the preamble's decision-log pointer verbatim with old
# §18's row-stability sentence appended — old §18 says nothing else.

new1_head = ["## 1. Purpose and method", "",
"""
This document specifies a modeling and simulation framework intended to
replace `FlightCore` as the substrate for `FlightPhysics` and `FlightApps`.
It is the normative statement of the design: what the framework *is*, in
present tense. It is derived from `framework_design.md`, which remains the
historical original and carries the full record of how the design was
reached. The new framework must match or surpass `FlightCore` in
functionality, performance and flexibility, while being more rigorous and
explicit — reducing the learning curve and the number of latent footguns for
model authors.""" |> strip |> String, ""]

# old §1's ground rules: everything from "Ground rules adopted" to the end.
ground_rules = old1[findfirst(startswith("Ground rules adopted"), old1):end]

new1_tail = ["",
"""
All design axes are settled — the formalism, the component taxonomy, the
signal and scheduling model, time and execution, the runtime periphery, the
declaration layer, the build pipeline, error discipline and the stopped-sim
services. Only §16's items — the migration outline, the GUI panel authoring
API and the log/trace persistence deferral — remain open.""" |> strip |> String,
"",
"""
Decision rationale, including the alternatives considered and the reasons
they were rejected, lives in `framework_decisions.md`, cited throughout as
"row N": one row per settled decision. Row numbers are stable, so a citation
here always names the same row there.""" |> strip |> String]

new1 = [new1_head; ground_rules; new1_tail]

# ------------------------------------------------------------- new ledes

lede6 = """
Components become a system through wiring: connections that route signals
across the assembly hierarchy, and ordinary junction components wherever
several signals must combine into one. §6.1 gives the connection and
hierarchy rules; §6.2 gives the aggregation idiom they force."""

lede10 = """
Where §9 fixes how data crosses the loop boundary, this chapter covers the
machinery that drives the loop itself: the control plane and the scheduling
primitives, the shutdown protocol, and the run lifecycle from `init!` through
replay."""

# ------------------------------------------------------- assembled chapters

new4 = [old4; ""; retitle(old7, "### 4.4 Function-valued signals: environment access")]
new5 = [old5; ""; retitle(old10, "### 5.5 Diagnostics: feedthrough tracing")]
new6 = ["## 6. Composition: connections, aggregation and hierarchy", "",
        split(strip(lede6), "\n")...,
        "",
        retitle(old8, "### 6.1 Connections and hierarchy")...,
        "",
        retitle(old6, "### 6.2 Aggregation: explicit summing junctions")...]
new7 = renumber(retitle(old9, "## 7. State and data representation"), 9, 7)
new8 = renumber(retitle(old11, "## 8. Time and execution"), 11, 8)

# The old §12 lede becomes chapter 9's; its "§12.1–§12.6" range straddled the
# chapter split, so the dash left by the sweep becomes a comma.
lede9 = replace.(retitle(old12ld, "## 9. Runtime periphery: the data plane"),
                 "(§9–§10.1)" => "(§9, §10.1)")
new9  = [lede9; ""; renumber(old12a, 12, 9)]
new10 = ["## 10. Runtime periphery: lifecycle and orchestration", "",
         split(strip(lede10), "\n")...,
         "",
         renumber(old12b, 12, 10; shift = -5)...]

new11 = renumber(retitle(old13, "## 11. The declaration layer"), 13, 11)
new12 = renumber(retitle(old14, "## 12. The build pipeline"), 14, 12)
new13 = renumber(retitle(old15, "## 13. Error discipline"), 15, 13)
new14 = renumber(retitle(old16, "## 14. Stopped-sim services"), 16, 14)
new15 = renumber(retitle(old17, "## 15. Case studies"), 17, 15)
new16 = retitle(old19, "## 16. Open axes")

# ------------------------------------------------------------- assembly

part(heading, blk) = ["# $heading"; ""; blk]

units = [[title],
         part("Part I — Foundations", new1),
         old2, old3, new4, new5, new6, new7,
         part("Part II — Execution", new8), new9, new10,
         part("Part III — Authoring and build", new11), new12,
         part("Part IV — Failure and services", new13), new14,
         part("Part V — Grounding", new15), new16,
         appendix]

doc = join((join(u, "\n") for u in units), "\n\n---\n\n") * "\n"
write(SPEC, doc)
println("framework_spec.md restructured: ", count("\n", doc) + 1, " lines")
