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
#     label) — every glossary entry must be reachable from the text, unless
#     whitelisted below (strict-flip adjudication, 2026-08-13). A whitelisted
#     anchor that *acquires* a body link is reported as a stale whitelist entry.
#
# Usage:  julia docs/notes/design/tools/check_glossary.jl [--strict]

const DESIGN = normpath(joinpath(@__DIR__, ".."))
const SPEC = joinpath(DESIGN, "framework_spec.md")

# Anchors adjudicated as deliberately unlinked (strict skips them):
const WHITELIST = Dict(
    # heading-defined: the term is verbatim a section title; the TOC links the
    # section and the entry's own text cites it — no prose use to link
    "g-periodic-discrete-component" => "heading-defined (§3.2)",
    "g-immutable-value-semantics" => "heading-defined (§4.1)",
    "g-next-snapshot-wait" => "heading-defined (§10.3)",
    "g-declaration-inventory" => "heading-defined (§11.2)",
    "g-always-on-conformance-check" => "heading-defined (§12.5)",
    "g-collect-the-checks-fail-the-evaluations-fast" => "heading-defined (§13.1)",
    "g-carrier-exception" => "heading-defined (§13.2)",
    # meta / never used bare in body prose
    "g-entry" => "disambiguation entry; never used bare by its own text",
    "g-row" => "'row N' citations are mechanical; linking them was rejected",
    "g-the-letters" => "only body occurrences are §5's own definition passage",
    "g-schema-vs-layout" => "glossary-coined recall key; body says schema/layout separately",
    "g-fragment-tree" => "glossary-coined recall key; body says fragment functions",
    "g-service-lifecycle" => "glossary-coined recall key; body says run lifecycle",
    # coined by the polysemy pass (2026-08-15) to separate the assembly/Build's
    # derived checkable surface from the declared-interface sense of "contract";
    # the seven body sites read as ordinary prose and are deliberately unlinked
    "g-derived-contract" => "polysemy coinage; disambiguates prose, no body link",
    # orphaned by the class-B glossary sense audit (2026-08-15): the anchor's
    # sole body link was a subject-verb mis-parse ("the trace records...",
    # verb "to record"), not the noun compound; delinked, leaving no body use
    "g-trace-record" => "sense audit: sole link was a mis-parsed verb use, not the noun",
)

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
        # hyphen/space variants both count as a body match
        variants(al) = (lowercase(al), lowercase(replace(al, " " => "-")))
        any(occursin(v, bodytext) for al in aliases(t) for v in variants(al)) ||
            push!(unused, t)
    end

    # anchors nothing links to: `#g-…` inline or `[g-…]:`/`[…][g-…]` reference use
    linked(id) = occursin("#$id", bodytext) || occursin("[$id]", bodytext)
    unlinked = [id for id in ids if !linked(id) && !haskey(WHITELIST, id)]
    stale = [id for id in keys(WHITELIST) if linked(id)]

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
    if strict && !isempty(stale)
        println(length(stale), " whitelisted anchors now linked (trim WHITELIST):")
        foreach(id -> println("  ", id), sort!(collect(stale)))
    end

    fail = !isempty(dup) ||
           (strict && (!isempty(noanchor) || !isempty(unlinked) || !isempty(stale)))
    println(fail ? "FAIL" : "OK")
    return fail ? 1 : 0
end

exit(main("--strict" in ARGS))
