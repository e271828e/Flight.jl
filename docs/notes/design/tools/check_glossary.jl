#!/usr/bin/env julia
#
# Glossary coverage checker for framework_spec.md (rewrite plan P0.3).
#
# Extracts Appendix D's bold entry terms and, in the default report mode,
# prints the advisories that are expected to be nonempty before Phase 1 lands:
#
#   - entries with no `<a id="g-…">` anchor (Phase 1.2 adds them);
#   - entries whose term text never appears in the body — drift or pruning
#     candidates, advisory only (some meta-entries legitimately never match).
#
# In strict mode (--strict, the pre-commit check from Phase 1 on), errors on:
#
#   - an entry without an anchor, or duplicate anchor ids;
#   - an anchor no body link points at (`#g-…` inline or `[g-…]` reference
#     label) — every glossary entry must be reachable from the text.
#
# Usage:  julia docs/notes/design/tools/check_glossary.jl [--strict]

const DESIGN = normpath(joinpath(@__DIR__, ".."))
const SPEC = joinpath(DESIGN, "framework_spec.md")

"Split the file at Appendix D. Returns (body, glossary) line vectors."
function split_at_glossary()
    lines = readlines(SPEC)
    i = findfirst(startswith("## Appendix D"), lines)
    i === nothing && error("Appendix D heading not found")
    return lines[1:i-1], lines[i:end]
end

"Matchable alias strings for an entry term: strip code ticks, split on \" / \",
drop a trailing parenthetical."
function aliases(term)
    t = replace(term, "`" => "")
    parts = split(t, " / ")
    return [strip(replace(p, r"\s*\([^)]*\)$" => "")) for p in parts]
end

function main(strict::Bool)
    body, glossary = split_at_glossary()
    bodytext = lowercase(join(body, "\n"))

    # entries: optional anchor prefix, then the bold term opening the line
    entries = Tuple{String,Union{String,Nothing}}[]   # (term, anchor id)
    for line in glossary
        m = match(r"^(?:<a id=\"(g-[a-z0-9_-]+)\"></a>\s*)?\*\*(.+?)\*\*(?= —| –|—)", line)
        m === nothing && continue
        push!(entries, (m[2], m[1]))
    end
    println("glossary: ", length(entries), " entries")

    noanchor = [t for (t, a) in entries if a === nothing]
    ids = [a for (_, a) in entries if a !== nothing]
    dup = [id for id in unique(ids) if count(==(id), ids) > 1]

    unused = String[]
    for (t, _) in entries
        any(occursin(lowercase(al), bodytext) for al in aliases(t)) || push!(unused, t)
    end

    # anchors nothing links to: `#g-…` inline or `[g-…]:`/`[…][g-…]` reference use
    unlinked = [id for id in ids
                if !occursin("#$id", bodytext) && !occursin("[$id]", bodytext)]

    isempty(dup) || println("duplicate anchor ids: ", dup)
    if !isempty(noanchor)
        println(length(noanchor), " entries without anchors",
                strict ? ":" : " (expected before Phase 1.2):")
        strict && foreach(t -> println("  ", t), noanchor)
    end
    if !isempty(unused)
        println(length(unused), " entry terms never matched in the body (advisory):")
        foreach(t -> println("  ", t), unused)
    end
    if strict && !isempty(unlinked)
        println(length(unlinked), " anchors with no body link:")
        foreach(id -> println("  ", id), unlinked)
    end

    fail = !isempty(dup) || (strict && (!isempty(noanchor) || !isempty(unlinked)))
    println(fail ? "FAIL" : "OK")
    return fail ? 1 : 0
end

exit(main("--strict" in ARGS))
