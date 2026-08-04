#!/usr/bin/env julia
#
# Cross-reference checker for the framework spec and its companion files.
#
# Two checks, both over framework_spec.md and the six companions:
#
#   1. Citations — every `§N` / `§N.M` / `§X.N` / `Appendix X` names a heading
#      that exists in framework_spec.md. Citations inside fenced code blocks,
#      code spans and headings count too: the phase-2/3 sweeps rewrote them, and
#      framework_decisions.md keeps all of its citations plain by design.
#   2. Anchors — every markdown link target (`#slug` in-file,
#      `framework_spec.md#slug` cross-file) resolves to a real heading anchor.
#
# `row N` citations are out of scope: they name rows of framework_decisions.md,
# not sections.
#
# Usage:  julia docs/notes/design/tools/check_refs.jl
# Exits nonzero if anything dangles.

include(joinpath(@__DIR__, "slugs.jl"))

const DESIGN = normpath(joinpath(@__DIR__, ".."))
const SPEC = "framework_spec.md"

const COMPANIONS = ["framework_decisions.md",
                    "event_visibility_walkthrough.md",
                    "inbound_periphery_walkthrough.md",
                    "review4_2026-08-03_dryrun.md",
                    "review4_2026-08-03_gaps.md",
                    "review4_2026-08-03_inconsistencies.md"]

const CITATION = r"§([A-D]|\d+)(?:\.(\d+))?|Appendix ([A-D])(?![\w–—-])"
const ANCHOR = r"\]\(([^)#]*)#([^)]+)\)"

function main()
    hs = headings(joinpath(DESIGN, SPEC))
    numbers = keys(targets(hs))
    slugs = Set(h.slug for h in hs)
    dup = collisions(hs)

    println("outline of $SPEC: ", length(numbers), " citable headings (",
            count(h -> h.level == 1, hs) - 1, " parts, ",
            count(n -> !occursin('.', n), numbers), " chapters/appendices, ",
            count(n -> occursin('.', n), numbers), " sections)")
    isempty(dup) || println("  slug suffixes in use: ", dup)

    bad = Tuple{String,Int,String,String}[]
    tc = ta = 0
    for file in [SPEC; COMPANIONS]
        path = joinpath(DESIGN, file)
        if !isfile(path)
            println("  skipped (absent): ", file)
            continue
        end
        own = Set(h.slug for h in headings(path))
        nc = na = 0
        for (lineno, line) in enumerate(eachline(path))
            for m in eachmatch(CITATION, line)
                nc += 1
                id = m[3] !== nothing ? m[3] :
                     m[2] === nothing ? m[1] : "$(m[1]).$(m[2])"
                id in numbers || push!(bad, (file, lineno, m.match, "unknown section"))
            end
            for m in eachmatch(ANCHOR, line)
                na += 1
                dest, slug = m[1], m[2]
                pool = dest == "" ? own : dest == SPEC ? slugs : nothing
                if pool === nothing
                    push!(bad, (file, lineno, m.match, "unknown link destination"))
                elseif !(slug in pool)
                    push!(bad, (file, lineno, m.match, "unknown anchor"))
                end
            end
        end
        tc += nc; ta += na
        println("  ", rpad(file, 42), lpad(nc, 5), " citations ", lpad(na, 5), " anchors")
    end
    println("total: ", tc, " citations, ", ta, " anchors")

    if isempty(bad)
        println("OK — every citation and every anchor resolves.")
        return 0
    end
    println("\nDANGLING (", length(bad), "):")
    for (file, lineno, tok, why) in bad
        println("  $file:$lineno: $tok — $why")
    end
    return 1
end

exit(main())
