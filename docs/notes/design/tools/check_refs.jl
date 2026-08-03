#!/usr/bin/env julia
#
# Cross-reference checker for the framework spec and its companion files.
#
# Extracts the chapter/section outline of framework_spec.md and verifies that
# every `§N` / `§N.M` citation — in the spec and in the companion files — names
# a heading that exists. Appendix citations (`§D.4`) are checked against the
# appendix headings the same way. `row N` citations are out of scope: they name
# rows of framework_decisions.md, not sections.
#
# Citations inside fenced code blocks are checked too (the phase-2 sweep
# rewrote them as well). Headings are read outside fenced blocks only.
#
# Usage:  julia docs/notes/design/tools/check_refs.jl
# Exits nonzero if any citation dangles.

const DESIGN = normpath(joinpath(@__DIR__, ".."))
const SPEC = "framework_spec.md"

const COMPANIONS = ["framework_decisions.md",
                    "event_visibility_walkthrough.md",
                    "inbound_periphery_walkthrough.md",
                    "review3_2026-07-30_gaps.md",
                    "review3_2026-07-30_improvements.md",
                    "review3_2026-07-30_inconsistencies.md"]

const CITATION = r"§([A-D]|\d+)(?:\.(\d+))?"

"Outline of `framework_spec.md`: the set of citable heading numbers."
function outline(path)
    ids = String[]
    infence = false
    for line in eachline(path)
        startswith(line, "```") && (infence = !infence; continue)
        infence && continue
        for (re, id) in ((r"^## (\d+)\. ", m -> m[1]),
                         (r"^### (\d+)\.(\d+) ", m -> "$(m[1]).$(m[2])"),
                         (r"^## Appendix ([A-D])[.:] ", m -> m[1]),
                         (r"^### ([A-D])\.(\d+) ", m -> "$(m[1]).$(m[2])"))
            m = match(re, line)
            m === nothing || push!(ids, id(m))
        end
    end
    return ids
end

function main()
    ids = outline(joinpath(DESIGN, SPEC))
    valid = Set(ids)
    length(valid) == length(ids) || error("duplicate heading numbers in $SPEC")

    chapters = sort([parse(Int, i) for i in ids if !occursin('.', i) && all(isdigit, i)])
    println("outline of $SPEC: ", length(ids), " citable headings (",
            length(chapters), " chapters, ",
            count(i -> occursin('.', i), ids), " sections)")

    dangling = Tuple{String,Int,String,String}[]
    total = 0
    for file in [SPEC; COMPANIONS]
        path = joinpath(DESIGN, file)
        isfile(path) || (println("  skipped (absent): ", file); continue)
        n = 0
        for (lineno, line) in enumerate(eachline(path))
            for m in eachmatch(CITATION, line)
                n += 1
                id = m[2] === nothing ? m[1] : "$(m[1]).$(m[2])"
                id in valid || push!(dangling, (file, lineno, m.match, strip(line)))
            end
        end
        total += n
        println("  ", rpad(file, 42), lpad(n, 5), " citations")
    end
    println("total: ", total, " citations")

    if isempty(dangling)
        println("OK — every citation resolves.")
        return 0
    end
    println("\nDANGLING (", length(dangling), "):")
    for (file, lineno, tok, ctx) in dangling
        println("  $file:$lineno: $tok — ", first(ctx, 90))
    end
    return 1
end

exit(main())
