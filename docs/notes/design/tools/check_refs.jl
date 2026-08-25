#!/usr/bin/env julia
#
# Cross-reference checker for the framework spec and its companion files.
#
# Two checks, both over framework_spec.md and the companions:
#
#   1. Citations — every `§N` / `§N.M` / `§X.N` / `Appendix X` names a heading
#      that exists in framework_spec.md. Citations inside fenced code blocks,
#      code spans and headings count too: the phase-2/3 sweeps rewrote them, and
#      framework_decisions.md keeps all of its citations plain by design.
#   2. Anchors — every markdown link target (`#slug` in-file,
#      `framework_spec.md#slug` or `framework_decisions.md#slug` cross-file)
#      resolves to a real heading anchor. In-file anchors are checked against the
#      file's own headings, which is what keeps the companions' self-references
#      (a walkthrough citing its own `§N`) honest.
#   3. Reference links (Phase 1.1 form, `[§9.4][s9-4]` + a generated definitions
#      block) — every label used outside code has a definition in the same file,
#      every definition's target resolves like an inline anchor, and every
#      definition is used (an unused one means the block is stale: re-run
#      linkify.jl). Labels are `s`- or `d`-prefixed by construction, and fenced
#      code is skipped, so array indexing in sketches never reads as a reference
#      link.
#
# Plus one advisory, not an error: in the walkthroughs — the companions that cite
# their own numbered sections — a link labelled `§N` that points at the spec while
# the walkthrough itself has a section `N`. Often legitimate (the spec is usually
# what is meant), but it is the shape a mis-resolved self-reference takes, so the
# set is printed and any growth in it deserves a look. Current known set, all
# hand-verified as correct spec citations: event_visibility_walkthrough.md l.20
# (§7). Anything beyond that line is new and needs the same hand check.
#
# Decision citations are ordinary links after the Pass B sweep — `D-037` is a
# `[D-037][d-037]` reference resolving into framework_decisions.md — so checks 2
# and 3 cover them with no special case. check_rows.jl remains the guard on
# citation *existence* for both the `D-nnn` and the retired `row N` spellings.
#
# Usage:  julia docs/notes/design/tools/check_refs.jl
# Exits nonzero if anything dangles.

include(joinpath(@__DIR__, "slugs.jl"))

const DESIGN = normpath(joinpath(@__DIR__, ".."))
const SPEC = "framework_spec.md"
const DECISIONS = "framework_decisions.md"

const COMPANIONS = ["framework_decisions.md",
                    "framework_extensions.md",
                    "event_visibility_walkthrough.md",
                    "inbound_periphery_walkthrough.md",
                    "trim_environment_walkthrough.md",
                    "frozen_discrete_walkthrough.md",
                    "localization_validation_walkthrough.md",
                    "sample_time_proposal.md"]

# The companions that cite their own numbered sections (see the advisory above).
const SELF_CITING = ["event_visibility_walkthrough.md",
                     "inbound_periphery_walkthrough.md"]

const CITATION = r"§([A-D]|\d+)(?:\.(\d+))?|Appendix ([A-D])(?![\w–—-])"
const ANCHOR = r"\]\(([^)#]*)#([^)]+)\)"
const LABELLED = r"\[§(\d+)\](?:\(([^)#]*)#|\[)"
const REFUSE = r"\]\[(s[A-D0-9][A-Za-z0-9-]*|d-\d+)\]"
const REFDEF = r"^\[(s[A-D0-9][A-Za-z0-9-]*|d-\d+)\]:\s*(\S*)#(\S+)\s*$"

"Explicit HTML anchor ids in a file (`<a id=\"…\">` — the glossary's g- anchors),
which are link targets exactly like heading slugs."
htmlids(path) = Set(m[1] for line in eachline(path)
                    for m in eachmatch(r"<a id=\"([^\"]+)\">", line))

function main()
    hs = headings(joinpath(DESIGN, SPEC))
    numbers = keys(targets(hs))
    slugs = union(Set(h.slug for h in hs), htmlids(joinpath(DESIGN, SPEC)))
    logslugs = union(Set(h.slug for h in headings(joinpath(DESIGN, DECISIONS))),
                     htmlids(joinpath(DESIGN, DECISIONS)))
    dup = collisions(hs)

    "Anchors a link destination may name, or `nothing` if it names no known file."
    anchors(dest, own) = dest == "" ? own :
                         dest == SPEC ? slugs :
                         dest == DECISIONS ? logslugs : nothing

    println("outline of $SPEC: ", length(numbers), " citable headings (",
            count(h -> startswith(h.text, "Part "), hs), " parts, ",
            count(n -> !occursin('.', n), numbers), " chapters/appendices, ",
            count(n -> occursin('.', n), numbers), " sections)")
    isempty(dup) || println("  slug suffixes in use: ", dup)

    bad = Tuple{String,Int,String,String}[]
    ambiguous = Tuple{String,Int,String}[]
    tc = ta = tr = 0
    for file in [SPEC; COMPANIONS]
        path = joinpath(DESIGN, file)
        if !isfile(path)
            println("  skipped (absent): ", file)
            continue
        end
        ownhs = headings(path)
        own = union(Set(h.slug for h in ownhs), htmlids(path))
        ownnums = file in SELF_CITING ?
                  Set(h.number for h in ownhs
                      if h.number !== nothing && !occursin('.', h.number)) :
                  Set{String}()
        nc = na = nr = 0
        defs = Dict{String,Tuple{Int,String,String}}()   # label => (lineno, dest, slug)
        uses = Dict{String,Int}()                        # label => first-use line
        infence = false
        for (lineno, line) in enumerate(eachline(path))
            startswith(line, "```") && (infence = !infence)
            for m in eachmatch(CITATION, line)
                nc += 1
                id = m[3] !== nothing ? m[3] :
                     m[2] === nothing ? m[1] : "$(m[1]).$(m[2])"
                id in numbers || push!(bad, (file, lineno, m.match, "unknown section"))
            end
            for m in eachmatch(ANCHOR, line)
                na += 1
                dest, slug = m[1], m[2]
                pool = anchors(dest, own)
                if pool === nothing
                    push!(bad, (file, lineno, m.match, "unknown link destination"))
                elseif !(slug in pool)
                    push!(bad, (file, lineno, m.match, "unknown anchor"))
                end
            end
            for m in eachmatch(LABELLED, line)
                (m[2] === nothing || m[2] == SPEC) && m[1] in ownnums &&
                    push!(ambiguous, (file, lineno, "§$(m[1])"))
            end
            infence && continue                          # ref links: prose only
            d = match(REFDEF, line)
            d !== nothing && (defs[d[1]] = (lineno, d[2], d[3]))
            for m in eachmatch(REFUSE, line)
                nr += 1
                haskey(uses, m[1]) || (uses[m[1]] = lineno)
            end
        end
        for (l, lineno) in uses
            haskey(defs, l) ||
                push!(bad, (file, lineno, "[$l]", "undefined reference label"))
        end
        for (l, (lineno, dest, slug)) in defs
            haskey(uses, l) ||
                push!(bad, (file, lineno, "[$l]", "unused definition — re-run linkify.jl"))
            pool = anchors(dest, own)
            if pool === nothing
                push!(bad, (file, lineno, "[$l]: $dest#$slug", "unknown link destination"))
            elseif !(slug in pool)
                push!(bad, (file, lineno, "[$l]: #$slug", "unknown anchor"))
            end
        end
        tc += nc; ta += na; tr += nr
        println("  ", rpad(file, 42), lpad(nc, 5), " citations ", lpad(na, 5),
                " anchors ", lpad(nr, 5), " refs")
    end
    println("total: ", tc, " citations, ", ta, " anchors, ", tr, " reference links")

    if !isempty(ambiguous)
        println("\nself-vs-spec advisory (", length(ambiguous),
                " walkthrough citations resolved to the spec while naming a section",
                " number the walkthrough also has):")
        for (file, lineno, tok) in ambiguous
            println("  $file:$lineno: $tok")
        end
    end

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
