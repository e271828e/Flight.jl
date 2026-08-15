#!/usr/bin/env julia
#
# Row-citation sweep for the framework spec and its companions (decisions-log
# rewrite, Pass B step 1). Converts the prose `row N` family to the entry form
# the log itself already speaks — `row 37` → `D-037`, `rows 106–107` →
# `D-106–D-107`, `rows 44, 104 and 106` → `D-044, D-104 and D-106` — dropping
# the now-redundant `row`/`rows` word, as Pass A did inside the log.
#
# framework_decisions.md is not swept: Pass A already converted it (708 D-refs,
# no row form left). It is scanned anyway, so a regression there would show up.
#
# Only the citation form changes. Every other character on a rewritten line is
# preserved, and that is *checked*, not merely intended: each line is reduced to
# a canonical form in which every citation group becomes `⟪n,m,…⟫`, and the
# rewrite is rejected unless the before- and after-canonicalizations are byte
# equal. Separators inside a group (en-dash ranges, `, `, ` and `, `/`) are
# carried over verbatim, so the canonical forms line up group for group.
#
# Two further guards, both fatal:
#   * every number swept must name an entry that exists in framework_decisions.md
#     (so a `row 5` that meant something else would stop the run, not be mangled);
#   * the exotic forms the runbook warned about (`row-79`, `rows ≤145`) are not
#     silently passed through — they abort, because their target spelling is a
#     judgement call, not a substitution.
#
# Fenced code blocks are skipped. As of the sweep there are no citations inside
# them, nor inside code spans; the skip keeps that true for later runs.
#
# Usage:  julia docs/notes/design/tools/sweep_rows.jl [--apply]
# Dry run by default: prints per-file counts and every distinct rewrite.
# Exits nonzero if any guard trips.

const DESIGN = normpath(joinpath(@__DIR__, ".."))
const DECISIONS = "framework_decisions.md"

const FILES = ["framework_spec.md",
               "framework_extensions.md",
               "sample_time_proposal.md",
               "event_visibility_walkthrough.md",
               "inbound_periphery_walkthrough.md",
               "trim_environment_walkthrough.md",
               "frozen_discrete_walkthrough.md",
               "localization_validation_walkthrough.md",
               DECISIONS]

# A citation group in either spelling. The body is numbers joined by range
# dashes and list separators; both regexes capture the whole group so the
# canonicalizer can collapse it to its number set.
const ROWGROUP = r"\brows?[\s-]+(?:[≤≥]\s*)?(\d+(?:\s*[–-]\s*\d+)?(?:(?:,\s*|\s+and\s+|/)\d+(?:\s*[–-]\s*\d+)?)*)"
const DGROUP = r"\bD-\d+(?:\s*[–-]\s*D-\d+)?(?:(?:,\s*|\s+and\s+|/)D-\d+(?:\s*[–-]\s*D-\d+)?)*"

# Link targets are generated from heading text, never authored, so the sweep
# leaves them alone; a heading whose own citation is swept gets its slug (and
# the links to it) regenerated afterwards, with check_refs.jl as the guard.
const TARGET = r"\]\([^)]*\)"

# A row citation, unless it is inside a link target — scanned left to right so
# the target alternative wins.
const SCAN = Regex(TARGET.pattern * "|" * ROWGROUP.pattern)

# Forms whose conversion is not a plain substitution — see the header.
const EXOTIC = r"\brows?-\d|\brows?\s+[≤≥]"

"A line with its link targets blanked, so slugs read as neither citation nor prose."
blank(line) = replace(line, TARGET => "](…)")

"Entry numbers defined by the log's `### D-nnn — Title` headings."
function defined()
    ns = Set{Int}()
    for line in eachline(joinpath(DESIGN, DECISIONS))
        m = match(r"^### D-(\d+) — ", line)
        m === nothing || push!(ns, parse(Int, m[1]))
    end
    isempty(ns) && error("no entry headings found in $DECISIONS")
    return ns
end

dref(n) = "D-" * lpad(n, 3, '0')

"Numbers named by a group body, ranges expanded. Entry prefixes are dropped
first, so a `D-185–D-187` range reads as the range it is rather than as a pair."
function numbers(body::AbstractString)
    out = Int[]
    for m in eachmatch(r"(\d+)(?:\s*[–-]\s*(\d+))?", replace(body, "D-" => ""))
        lo = parse(Int, m[1])
        hi = m[2] === nothing ? lo : parse(Int, m[2])
        lo <= hi || error("descending range \"$(m.match)\"")
        append!(out, lo:hi)
    end
    return out
end

"A line with every citation group of `re` collapsed to its number set."
canonical(line, re) =
    replace(line, re => m -> "⟪" * join(sort(numbers(m)), ",") * "⟫")

"""
    sweepline(line, ok) -> (newline, nrewrites, badnums)

Rewrite every row-form citation on one line. `badnums` collects numbers that
name no entry; when it is non-empty the caller aborts rather than write.
"""
function sweepline(line, ok)
    n = 0
    bad = Int[]
    out = replace(line, SCAN => function (m)
        startswith(m, "](") && return m              # link target
        body = match(ROWGROUP, m)[1]
        ns = numbers(body)
        append!(bad, filter(x -> !(x in ok), ns))
        n += 1
        return replace(body, r"\d+" => d -> dref(parse(Int, d)))
    end)
    return out, n, bad
end

function main(apply::Bool)
    ok = defined()
    println("decision log: ", length(ok), " entries (max ", maximum(ok), ")")

    fatal = String[]
    rewrites = Dict{String,Int}()          # "before → after" => count
    total = 0
    for file in FILES
        path = joinpath(DESIGN, file)
        lines = readlines(path)
        n = 0
        infence = false
        for (i, line) in pairs(lines)
            startswith(line, "```") && (infence = !infence; continue)
            infence && continue
            e = match(EXOTIC, blank(line))
            e === nothing ||
                push!(fatal, "$file:$i: exotic form \"$(e.match)\" — convert by hand")
            new, k, bad = sweepline(line, ok)
            k == 0 && continue
            isempty(bad) ||
                push!(fatal, "$file:$i: cites undefined entries $(sort(unique(bad)))")
            if canonical(blank(line), ROWGROUP) != canonical(blank(new), DGROUP)
                push!(fatal, "$file:$i: rewrite altered more than the citation form")
                continue
            end
            for m in eachmatch(ROWGROUP, blank(line))
                after = replace(m.match, ROWGROUP =>
                                x -> replace(match(ROWGROUP, x)[1],
                                             r"\d+" => d -> dref(parse(Int, d))))
                rewrites["$(m.match) → $after"] = get(rewrites, "$(m.match) → $after", 0) + 1
            end
            lines[i] = new
            n += k
        end
        total += n
        println("  ", rpad(file, 42), lpad(n, 4), " groups rewritten")
        apply && isempty(fatal) && n > 0 && write(path, join(lines, "\n") * "\n")
    end
    println("total: ", total, " citation groups, ", length(rewrites), " distinct forms")

    if !isempty(fatal)
        println("\nABORTED (", length(fatal), "):")
        foreach(f -> println("  ", f), fatal)
        return 1
    end

    println("\ndistinct rewrites:")
    for k in sort!(collect(keys(rewrites)))
        println("  ", lpad(rewrites[k], 3), "×  ", k)
    end
    println(apply ? "\napplied." : "\ndry run — pass --apply to write.")
    return 0
end

exit(main("--apply" in ARGS))
