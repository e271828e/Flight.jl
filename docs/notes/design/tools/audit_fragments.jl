#!/usr/bin/env julia
#
# Content-preservation audit for the decision log (Pass C gate). Rebuilt from
# the description of Pass A's deleted original: a per-entry token-multiset
# comparison that catches a dropped citation, code span or named diagnostic
# when an entry's prose is rewritten.
#
# A prose rewrite is allowed to change every word; it is not allowed to lose a
# fact. So the audit ignores prose entirely and compares, per `D-nnn` entry,
# the multiset of four token classes:
#
#   code  — code spans, `like_this`
#   math  — inline math, $like this$
#   cite  — §N, §N.M, §D.n, Appendix X, D-nnn
#   name  — CamelCase identifiers with an internal hump (`TrimProblemInvalid`,
#           `NamedTuple`), the shape a named diagnostic takes; deliberately not
#           all-caps words (GUI, LM), which are prose.
#
# Multisets, not sets: dropping one of an entry's three `§14.6` citations is a
# loss and shows up as one.
#
# Reference links are unwrapped first (`[§9.4][s9-4]` → `§9.4`), so a linkify
# run between the two snapshots is invisible here — the citation is the same
# citation whether or not it is linked. The generated definitions block at EOF
# is dropped for the same reason, and the `## Index` never enters, being ahead
# of the first entry heading.
#
# Losses are fatal; gains are printed and tolerated, since a rewrite may
# legitimately add a citation and may never legitimately drop one. An entry
# that disappears between snapshots is a loss of the whole entry.
#
# Usage:  julia docs/notes/design/tools/audit_fragments.jl [rev]
# Compares the working-tree file against `rev` (default HEAD), so the natural
# call before a rewrite batch names the commit the batch started from.
# Exits nonzero if anything was lost.

const DESIGN = normpath(joinpath(@__DIR__, ".."))
const DECISIONS = "framework_decisions.md"
const RELPATH = "docs/notes/design/" * DECISIONS

const CLASSES = ["code" => r"`[^`]+`",
                 "math" => r"\$[^\$\n]+\$",
                 "cite" => r"§(?:[A-D]|\d+)(?:\.\d+)?|Appendix [A-D]\b|\bD-\d+\b",
                 "name" => r"\b[A-Z][A-Za-z0-9_]*[a-z][A-Z][A-Za-z0-9_]*\b"]

const REFLINK = r"\[([^\]]+)\]\[[^\]]+\]"      # [§9.4][s9-4]
const INLINK = r"\[([^\]]+)\]\([^)]*\)"        # [text](target)
const REFDEF = r"^\[[sd][^\]]*\]:\s"           # generated definitions block

"Entry bodies keyed by `D-nnn`, headings included, generated blocks dropped."
function entries(text::AbstractString)
    out = Dict{String,String}()
    id = nothing
    buf = IOBuffer()
    for line in split(text, '\n')
        occursin(REFDEF, line) && continue
        m = match(r"^### (D-\d+) — ", line)
        if m !== nothing
            id === nothing || (out[id] = String(take!(buf)))
            id = m[1]
        end
        id === nothing || println(buf, line)
    end
    id === nothing || (out[id] = String(take!(buf)))
    isempty(out) && error("no entry headings found")
    return out
end

"Token counts of one entry, keyed by (class, token)."
function tokens(body::AbstractString)
    plain = replace(replace(body, REFLINK => s"\1"), INLINK => s"\1")
    counts = Dict{Tuple{String,String},Int}()
    for (class, re) in CLASSES, m in eachmatch(re, plain)
        key = (class, m.match)
        counts[key] = get(counts, key, 0) + 1
    end
    return counts
end

function main(rev::AbstractString)
    old = entries(read(`git show $rev:$RELPATH`, String))
    new = entries(read(joinpath(DESIGN, DECISIONS), String))
    println("audit: $DECISIONS — $rev (", length(old), " entries) vs working tree (",
            length(new), ")")

    lost = String[]
    gained = String[]
    for id in sort!(collect(keys(old)))
        if !haskey(new, id)
            push!(lost, "$id — entry absent from the working tree")
            continue
        end
        before, after = tokens(old[id]), tokens(new[id])
        for key in union(keys(before), keys(after))
            d = get(after, key, 0) - get(before, key, 0)
            line = "$id  $(key[1])  $(key[2])" * (abs(d) == 1 ? "" : "  ×$(abs(d))")
            d < 0 && push!(lost, line)
            d > 0 && push!(gained, line)
        end
    end
    for id in sort!(collect(setdiff(keys(new), keys(old))))
        push!(gained, "$id — new entry")
    end

    if !isempty(gained)
        println("\nadded (", length(gained), "):")
        foreach(l -> println("  ", l), sort!(gained))
    end
    if !isempty(lost)
        println("\nLOST (", length(lost), "):")
        foreach(l -> println("  ", l), sort!(lost))
        return 1
    end
    println("\nOK — every entry keeps every code span, citation and named identifier.")
    return 0
end

exit(main(get(ARGS, 1, "HEAD")))
