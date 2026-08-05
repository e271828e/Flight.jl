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
# The walkthroughs cite their own numbered sections, so a bare `§N` in one of
# them is ambiguous: it may mean the walkthrough's own §N or the spec's. Such
# tokens are left alone and warned about — the author resolves them by writing
# the link explicitly, local (`#slug`) or cross-file, and the explicit link then
# survives every later run untouched.
#
# framework_decisions.md is deliberately excluded: it is one giant table whose
# rows are edited constantly, and links would clutter every cell's source.
#
# Run from anywhere:  julia docs/notes/design/tools/linkify.jl

include(joinpath(@__DIR__, "slugs.jl"))

const DESIGN = normpath(joinpath(@__DIR__, ".."))
const SPEC = "framework_spec.md"

const COMPANIONS = ["event_visibility_walkthrough.md",
                    "inbound_periphery_walkthrough.md",
                    "trim_environment_walkthrough.md"]

# Companions that cite their own numbered sections (the walkthroughs; the review
# files' numbered findings are never cited with `§`). Only here is a bare `§N`
# ambiguous between the file and the spec.
const SELF_CITING = ["event_visibility_walkthrough.md",
                     "inbound_periphery_walkthrough.md"]

# Scanned left to right. The first two alternatives are the exclusion zones —
# they are matched so they can be skipped whole; the rest are the citations.
const SCAN = r"`+[^`]*`+|\[[^\]]*\]\([^)]*\)|§\d+(?:\.\d+)?|§[A-D]\.\d+|Appendix [A-D](?![\w–—-])"

iscitation(tok) = startswith(tok, "§") || startswith(tok, "Appendix ")

"Citation key for a matched token: `§9.4` → `9.4`, `Appendix C` → `C`."
key(tok) = startswith(tok, "§") ? tok[nextind(tok, 1):end] : tok[end:end]

"The file's own top-level section numbers — `§N` tokens it may be citing itself with."
function ownsections(file)
    hs = headings(joinpath(DESIGN, file))
    return Set(h.number for h in hs if h.number !== nothing && !occursin('.', h.number))
end

"""
    linkify!(file, targets, prefix; own) -> (nlinks, ambiguous)

Rewrite every citation in `file` into `[citation](prefix#slug)`, skipping code
blocks, code spans, headings and existing links. Bare tokens naming one of
`own`'s section numbers are skipped and returned as `(lineno, token)` pairs.
"""
function linkify!(file, targets, prefix; own = Set{String}())
    path = joinpath(DESIGN, file)
    lines = readlines(path; keep = true)
    infence = false
    n = 0
    ambiguous = Tuple{Int,String}[]
    for (i, line) in pairs(lines)
        startswith(line, "```") && (infence = !infence; continue)
        (infence || startswith(line, "#")) && continue
        lines[i] = replace(line, SCAN => function (tok)
            iscitation(tok) || return tok                  # code span or link
            k = key(tok)
            if k in own                                    # self vs. spec: author decides
                push!(ambiguous, (i, tok))
                return tok
            end
            slug = get(targets, k, nothing)
            slug === nothing && return tok                 # dangling: check_refs reports it
            n += 1
            return "[$tok]($prefix#$slug)"
        end)
    end
    write(path, join(lines))
    return n, ambiguous
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
println("  ", rpad(SPEC, 42), lpad(first(linkify!(SPEC, tg, "")), 5), " links")
for f in COMPANIONS
    n, ambiguous = linkify!(f, tg, SPEC;
                            own = f in SELF_CITING ? ownsections(f) : Set{String}())
    println("  ", rpad(f, 42), lpad(n, 5), " links")
    for (lineno, tok) in ambiguous
        println("    warn: $f:$lineno bare $tok also names this file's own section",
                " — write the link explicitly")
    end
end
println("  contents: ", retoc!(hs), " entries")
