#!/usr/bin/env julia
#
# Row-citation checker for the framework spec and its companions (rewrite plan
# P0.3). Complements check_refs.jl, which owns `§` citations and anchors and
# declares `row N` out of scope.
#
# Two checks:
#
#   1. Validity — every `row N` / `rows N–M` / `rows N and M` citation, and
#      every `D-nnn` / `D-nnn–D-mmm` entry citation, in every scanned file,
#      names an entry that exists in framework_decisions.md (`### D-nnn`
#      headings; D-nnn preserves the retired table's row N).
#   2. Coverage (spec only) — the set of rows framework_spec.md cites must be a
#      superset of the committed baseline tools/row_baseline.txt. Phase 2 may
#      consolidate duplicate citations of a row, but a row dropping out of the
#      spec entirely means an argument lost its pointer. Rows newly cited are
#      reported, not errors. After an approved change to the covered set, run
#      with --rebaseline to regenerate the file.
#
# Usage:  julia docs/notes/design/tools/check_rows.jl [--rebaseline]
# Exits nonzero if a citation dangles or the spec's coverage shrank.

const DESIGN = normpath(joinpath(@__DIR__, ".."))
const SPEC = "framework_spec.md"
const DECISIONS = "framework_decisions.md"
const BASELINE = joinpath(@__DIR__, "row_baseline.txt")

const FILES = [SPEC,
               DECISIONS,   # rows cite each other; those pointers dangle too
               "framework_extensions.md",
               "sample_time_proposal.md",
               "event_visibility_walkthrough.md",
               "inbound_periphery_walkthrough.md",
               "trim_environment_walkthrough.md",
               "frozen_discrete_walkthrough.md",
               "localization_validation_walkthrough.md"]

# One citation group: "row 80", "rows 185–187", "rows 7 and 37", "rows 33, 34
# and 55". The inner token is a number or an (en- or hyphen-dashed) range.
const GROUP = r"\brows?\s+(\d+(?:\s*[–-]\s*\d+)?(?:(?:,\s*|,?\s+and\s+)\d+(?:\s*[–-]\s*\d+)?)*)"
const TOKEN = r"(\d+)(?:\s*[–-]\s*(\d+))?"

# Entry citations in the post-table format: "D-037", "D-166–D-168".
const DGROUP = r"\bD-(\d+)(?:\s*[–-]\s*D-(\d+))?"

"All row numbers cited by one group match, ranges expanded."
function expand(group::AbstractString)
    out = Int[]
    for m in eachmatch(TOKEN, group)
        lo = parse(Int, m[1])
        hi = m[2] === nothing ? lo : parse(Int, m[2])
        lo <= hi || error("descending range \"$(m.match)\"")
        append!(out, lo:hi)
    end
    return out
end

"Row numbers defined by the decision log's entry headings."
function definedrows()
    rows = Set{Int}()
    for line in eachline(joinpath(DESIGN, DECISIONS))
        m = match(r"^### D-(\d+) — ", line)
        m === nothing || push!(rows, parse(Int, m[1]))
    end
    isempty(rows) && error("no entry headings found in $DECISIONS")
    return rows
end

function main(rebaseline::Bool)
    defined = definedrows()
    println("decision log: ", length(defined), " rows (max ", maximum(defined), ")")

    bad = Tuple{String,Int,String}[]
    speccited = Set{Int}()
    for file in FILES
        path = joinpath(DESIGN, file)
        cited = Set{Int}()
        n = 0
        for (lineno, line) in enumerate(eachline(path))
            for g in eachmatch(GROUP, line)
                for r in expand(g[1])
                    n += 1
                    push!(cited, r)
                    r in defined || push!(bad, (file, lineno, "row $r"))
                end
            end
            for g in eachmatch(DGROUP, line)
                lo = parse(Int, g[1])
                hi = g[2] === nothing ? lo : parse(Int, g[2])
                lo <= hi || error("descending range \"$(g.match)\"")
                for r in lo:hi
                    n += 1
                    push!(cited, r)
                    r in defined || push!(bad, (file, lineno, "D-$r"))
                end
            end
        end
        file == SPEC && (speccited = cited)
        println("  ", rpad(file, 42), lpad(n, 4), " citations, ",
                lpad(length(cited), 3), " distinct rows")
    end

    if rebaseline
        write(BASELINE, join(sort(collect(speccited)), "\n") * "\n")
        println("baseline regenerated: ", length(speccited), " rows")
    elseif isfile(BASELINE)
        baseline = Set(parse(Int, l) for l in eachline(BASELINE) if !isempty(l))
        lost = sort(collect(setdiff(baseline, speccited)))
        gained = sort(collect(setdiff(speccited, baseline)))
        isempty(gained) || println("newly cited by the spec (fine): ", gained)
        if !isempty(lost)
            println("COVERAGE SHRANK — rows no longer cited by the spec: ", lost)
            return 1
        end
    else
        println("no baseline (", BASELINE, ") — run with --rebaseline to create it")
    end

    if isempty(bad)
        println("OK — every row citation names an existing row.")
        return 0
    end
    println("\nDANGLING (", length(bad), "):")
    for (file, lineno, tok) in bad
        println("  $file:$lineno: $tok")
    end
    return 1
end

exit(main("--rebaseline" in ARGS))
