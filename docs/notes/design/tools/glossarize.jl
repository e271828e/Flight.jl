#!/usr/bin/env julia
#
# Assisted glossary linking for rewrite-plan Phase 1.3.
#
# Turns the first use of each glossary term in each section into an inline link
# to its Appendix D anchor — `the [drain](#g-drain) …` — within one Part of the
# spec at a time. The output is a *proposal*: every insertion is printed for
# review, and the reviewing human fixes wrong-sense links and adds glosses; the
# ambiguous vocabulary below is never auto-linked at all.
#
# Policy, encoding tools/coinage_inventory.md:
#
#   - (id, term) pairs are read live from the glossary's own anchor lines;
#     aliases split on " / ", trailing parentheticals stripped. The EXTRA table
#     adds the inventory's bucket-B aliases (frame-top drain → g-drain …).
#   - An alias that is code-styled in the glossary (`Build`, `at`, `w`) only
#     ever matches its backticked spelling — bare prose "build" or "at" never
#     links.
#   - Prose aliases in AMBIGUOUS — words that double as ordinary English or
#     have competing technical senses (entry, log, row, anchor, phase, class,
#     condition, register, …) — are skipped entirely; the reviewer links the
#     technical uses by hand.
#   - First use per section (### / ## resets), tracked per alias, so "interior
#     sweep" and "boundary sweep" each get their first link even though both
#     resolve to g-sweep.
#   - Never inside: fenced code, inline code spans (except as the whole match
#     of a code-styled alias), existing links, headings.
#   - Plurals: alias + "s" and y→ies match; the link text keeps what matched.
#
# Usage:  julia docs/notes/design/tools/glossarize.jl '# Part I' '# Part II'
# (start marker inclusive, end marker exclusive; both matched by startswith)

const SPEC = normpath(joinpath(@__DIR__, "..", "framework_spec.md"))

const EXTRA = ["frame-top drain" => "g-drain",
               "interior sweep"  => "g-sweep",
               "boundary sweep"  => "g-sweep",
               "bundle law"      => "g-bundle",
               "stepper seam"    => "g-seam",
               "frozen roster"   => "g-roster",
               "honest priors"   => "g-prior",
               "epoch rule"      => "g-input-epoch"]

const AMBIGUOUS = Set(["entry", "kind", "register", "capture", "merge",
                       "override", "log", "row", "anchor", "phase", "prior",
                       "condition", "class", "store", "view", "frame", "debt",
                       "due", "blessed", "worked", "nominal", "baseline",
                       "binding", "batch", "iff", "the letters", "normative"])

"(pattern alias, id, iscode) triples, longest alias first."
function aliases()
    out = Tuple{String,String,Bool}[]
    for line in eachline(SPEC)
        m = match(r"^<a id=\"(g-[a-z0-9_-]+)\"></a>\*\*(.+?)\*\*(?= —| –|—)", line)
        m === nothing && continue
        for raw in split(m[2], " / ")
            a = strip(replace(raw, r"\s*\([^)]*\)$" => ""))
            iscode = startswith(a, "`")
            a = replace(a, "`" => "")
            (!iscode && lowercase(a) in AMBIGUOUS) && continue
            push!(out, (a, m[1], iscode))
        end
    end
    append!(out, [(a, id, false) for (a, id) in EXTRA])
    return sort!(out; by = t -> -length(t[1]))
end

"Byte ranges of code spans and links in a line — no insertions inside these."
function protected(line)
    out = UnitRange{Int}[]
    for m in eachmatch(r"`+[^`]*`+|\[[^\]]*\]\([^)]*\)|\[[^\]]*\]\[[^\]]*\]|<a id=\"[^\"]*\"></a>", line)
        push!(out, m.offset:(m.offset + ncodeunits(m.match) - 1))
    end
    return out
end

"Case-insensitive whole-word pattern for a prose alias, plural-tolerant."
function pattern(alias, iscode)
    e = replace(alias, r"([\\.^$|?*+()\[\]{}])" => s"\\\1")
    iscode && return Regex("`$e`")
    p = endswith(alias, "y") ? "(?:$(e)|$(chop(e))ies)" : "(?:$(e)s?)"
    return Regex("(?i)\\b$p\\b")
end

function main(from, to)
    al = aliases()
    pats = [(a, id, pattern(a, iscode), iscode) for (a, id, iscode) in al]
    lines = readlines(SPEC)
    lo = findfirst(l -> startswith(l, from), lines)
    hi = findfirst(l -> startswith(l, to), lines)
    (lo === nothing || hi === nothing) && error("part markers not found")

    seen = Set{String}()          # aliases already linked in the current section
    infence = false
    n = 0
    for i in lo:hi-1
        line = lines[i]
        startswith(line, "```") && (infence = !infence; continue)
        infence && continue
        startswith(line, "#") && (empty!(seen); continue)
        for (a, id, pat, iscode) in pats
            a in seen && continue
            # an existing link to this id in the section counts as its first use
            occursin("(#$id)", line) && any(occursin(pat, x)
                for x in [m[1] for m in eachmatch(r"\[([^\]]+)\]\(#[^)]*\)", line)]) &&
                (push!(seen, a); continue)
            m = match(pat, line)
            m === nothing && continue
            lohi = m.offset:(m.offset + ncodeunits(m.match) - 1)
            any(!isempty(intersect(lohi, r)) for r in protected(line)) && continue
            head = m.offset == 1 ? "" : line[1:prevind(line, m.offset)]
            tail = line[nextind(line, last(lohi)):end]
            line = head * "[" * m.match * "](#$id)" * tail
            push!(seen, a)
            n += 1
            println("  $(lpad(i, 5)): $(m.match) → $id")
        end
        lines[i] = line
    end
    write(SPEC, join(lines, "\n") * "\n")
    println("$n links inserted in $from … $to")
end

main(ARGS[1], ARGS[2])
