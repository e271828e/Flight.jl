#!/usr/bin/env julia
#
# Phase-4 navigation for the design docs (see tools/restructure_plan.md).
#
# Two jobs, both idempotent — re-run this after any renumbering or renaming and
# the links and the table of contents regenerate against the new headings:
#
#   1. Turn `§N`, `§N.M`, `§X.N` (appendix) and `Appendix X` citations into
#      markdown links against GitHub-style heading anchors. In-document links
#      in framework_spec.md, cross-file links (`framework_spec.md#slug`) in the
#      walkthroughs and the review files.
#   2. Regenerate the spec's `## Contents` between the front matter and
#      `# Part I`.
#
# Never linkified: fenced code blocks, inline code spans, heading lines (so the
# computed slugs stay clean), and text already inside a markdown link — which is
# what makes re-running a no-op. `row N` citations are left alone: the decision
# table has no per-row anchors.
#
# framework_decisions.md is deliberately excluded: it is one giant table whose
# rows are edited constantly, and links would clutter every cell's source.
#
# Run from anywhere:  julia docs/notes/design/tools/linkify.jl

include(joinpath(@__DIR__, "slugs.jl"))

const DESIGN = normpath(joinpath(@__DIR__, ".."))
const SPEC = "framework_spec.md"

const COMPANIONS = ["event_visibility_walkthrough.md",
                    "inbound_periphery_walkthrough.md"]

# Scanned left to right. The first two alternatives are the exclusion zones —
# they are matched so they can be skipped whole; the rest are the citations.
const SCAN = r"`+[^`]*`+|\[[^\]]*\]\([^)]*\)|§\d+(?:\.\d+)?|§[A-D]\.\d+|Appendix [A-D](?![\w–—-])"

iscitation(tok) = startswith(tok, "§") || startswith(tok, "Appendix ")

"Citation key for a matched token: `§9.4` → `9.4`, `Appendix C` → `C`."
key(tok) = startswith(tok, "§") ? tok[nextind(tok, 1):end] : tok[end:end]

"""
    linkify!(file, targets, prefix) -> nlinks

Rewrite every citation in `file` into `[citation](prefix#slug)`, skipping code
blocks, code spans, headings and existing links.
"""
function linkify!(file, targets, prefix)
    path = joinpath(DESIGN, file)
    lines = readlines(path; keep = true)
    infence = false
    n = 0
    for (i, line) in pairs(lines)
        startswith(line, "```") && (infence = !infence; continue)
        (infence || startswith(line, "#")) && continue
        lines[i] = replace(line, SCAN => function (tok)
            iscitation(tok) || return tok                  # code span or link
            slug = get(targets, key(tok), nothing)
            slug === nothing && return tok                 # dangling: check_refs reports it
            n += 1
            return "[$tok]($prefix#$slug)"
        end)
    end
    write(path, join(lines))
    return n
end

"""
    toc(hs) -> Vector{String}

The `## Contents` block. Parts and appendices sit at the outer level, chapters
and the `D.x` glossary groups one level in, sections two. `####` headings are
not listed.
"""
function toc(hs)
    out = ["## Contents", ""]
    for (i, h) in enumerate(hs)
        (i == 1 || h.level >= 4 || h.text == "Contents") && continue
        isappendix = h.number !== nothing && occursin(r"^[A-D]", h.number)
        indent = h.level == 1 ? "" :               # parts
                 h.level == 2 ? (isappendix ? "" : "  ") :
                 isappendix ? "  " : "    "        # D.x groups vs. N.M sections
        push!(out, "$indent- [$(h.text)](#$(h.slug))")
    end
    push!(out, "")
    return out
end

"Regenerate the spec's `## Contents` between the front matter and `# Part I`."
function retoc!(hs)
    path = joinpath(DESIGN, SPEC)
    lines = readlines(path)
    part1 = findfirst(startswith("# Part I"), lines)
    part1 === nothing && error("`# Part I` not found")
    old = findfirst(isequal("## Contents"), lines)
    if old !== nothing
        old < part1 || error("`## Contents` is not in the front matter")
        deleteat!(lines, old:part1-1)
        part1 = findfirst(startswith("# Part I"), lines)
    end
    block = [toc(hs); "---"; ""]
    splice!(lines, part1:part1-1, block)
    write(path, join(lines, "\n") * "\n")
    return count(l -> startswith(strip(l), "- ["), block)
end

hs = headings(joinpath(DESIGN, SPEC))
tg = targets(hs)
dup = collisions(hs)
isempty(dup) || println("slug collisions resolved by suffix: ", dup)

println("linkify")
println("  ", rpad(SPEC, 42), lpad(linkify!(SPEC, tg, ""), 5), " links")
for f in COMPANIONS
    println("  ", rpad(f, 42), lpad(linkify!(f, tg, SPEC), 5), " links")
end
println("  contents: ", retoc!(hs), " entries")
