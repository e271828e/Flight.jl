#!/usr/bin/env julia
#
# Phase-4 reorder for framework_spec.md and its corpus.
#
# Three transformations, run in order:
#   1. Move the "Part III — Authoring and build" block so it sits before
#      "Part II — Execution", then swap the two part numerals (II <-> III)
#      so the parts stay ordered I..V.
#   2. Renumber the spec's own chapter/section headings for chapters 8-12
#      through a fixed permutation (subsection numbers are preserved).
#   3. Sweep every visible `§N` / `§N.M` citation token across the 9-file
#      corpus through the same permutation.
#
# Mechanics for steps 2-3: every chapter number is mapped ONCE through a
# pure Dict lookup (get(CHMAP, ch, ch)), not by sequential search-and-replace
# -- sequential replace would double-map (old §8 -> new §10, then that §10
# would be re-mapped again to §12 on a later pass over the same text).
#
# Run from anywhere:  julia docs/notes/design/tools/phase4_reorder.jl

const DESIGN = normpath(joinpath(@__DIR__, ".."))
const SPEC = joinpath(DESIGN, "framework_spec.md")

const CORPUS = ["framework_spec.md",
                "framework_decisions.md",
                "framework_extensions.md",
                "event_visibility_walkthrough.md",
                "inbound_periphery_walkthrough.md",
                "trim_environment_walkthrough.md",
                "frozen_discrete_walkthrough.md",
                "localization_validation_walkthrough.md",
                "sample_time_proposal.md"]

# Chapter permutation. Everything outside this map (1-7, 13-16) is identity.
const CHMAP = Dict(8 => 10, 9 => 11, 10 => 12, 11 => 8, 12 => 9)
remap_chapter(ch::Int) = get(CHMAP, ch, ch)

const TOKEN = r"§(\d+)(?:\.(\d+))?"

# Per-occurrence adjudications for the citation sweep, keyed by
# (file, 1-based pre-sweep line number, matched token) -- see sweep_refs.jl
# for the architecture this mirrors. The chapter permutation here is total
# and unambiguous, so this table is intentionally empty; if a citation ever
# seems to need special-casing, that's a signal to stop and report it rather
# than adjudicate it here.
const OVERRIDES = Dict{Tuple{String,Int,String},String}()

# ---------------------------------------------------------------------------
# Step 1: move the Part III block before Part II, then swap the two numerals.
# ---------------------------------------------------------------------------
function move_part_block!()
    lines = readlines(SPEC; keep = true)
    before_multiset = sort(copy(lines))

    idx(text) = findfirst(l -> strip(l) == text, lines)
    i2 = idx("# Part II — Execution")
    i3 = idx("# Part III — Authoring and build")
    i4 = idx("# Part IV — Failure and services")
    i2 === nothing && error("Part II heading not found")
    i3 === nothing && error("Part III heading not found")
    i4 === nothing && error("Part IV heading not found")
    i2 < i3 < i4 || error("unexpected part heading order: $i2, $i3, $i4")

    block = lines[i3:i4-1]                       # Part III heading .. line before Part IV
    remainder = vcat(lines[1:i3-1], lines[i4:end])  # i2 unaffected: i2 < i3
    newlines = vcat(remainder[1:i2-1], block, remainder[i2:end])

    if sort(newlines) != before_multiset
        error("line-multiset assertion FAILED: block move changed file contents")
    end
    println("block move: line-multiset assertion PASSED (", length(lines), " lines before and after)")

    # Numeral swap so parts stay ordered I..V.
    for (i, l) in pairs(newlines)
        s = strip(l)
        if s == "# Part III — Authoring and build"
            newlines[i] = replace(l, "Part III — Authoring and build" => "Part II — Authoring and build")
        elseif s == "# Part II — Execution"
            newlines[i] = replace(l, "Part II — Execution" => "Part III — Execution")
        end
    end

    write(SPEC, join(newlines))
end

# ---------------------------------------------------------------------------
# Step 2: renumber the spec's own `## N. Title` / `### N.M Title` headings.
# ---------------------------------------------------------------------------
function rewrite_headings!()
    lines = readlines(SPEC; keep = true)
    chapter_re = r"^## (\d+)\. (.*)$"
    section_re = r"^### (\d+)\.(\d+) (.*)$"
    changes = Tuple{String,String}[]

    for (i, l) in pairs(lines)
        nl = endswith(l, "\n") ? "\n" : ""
        s = rstrip(l, ['\n'])

        m = match(chapter_re, s)
        if m !== nothing
            old_ch = parse(Int, m[1])
            new_ch = remap_chapter(old_ch)
            if new_ch != old_ch
                new_s = "## $(new_ch). $(m[2])"
                push!(changes, (s, new_s))
                lines[i] = new_s * nl
            end
            continue
        end

        m2 = match(section_re, s)
        if m2 !== nothing
            old_ch = parse(Int, m2[1])
            new_ch = remap_chapter(old_ch)
            if new_ch != old_ch
                new_s = "### $(new_ch).$(m2[2]) $(m2[3])"
                push!(changes, (s, new_s))
                lines[i] = new_s * nl
            end
        end
    end

    write(SPEC, join(lines))
    return changes
end

# ---------------------------------------------------------------------------
# Step 3: sweep every visible §N / §N.M citation token across the corpus.
# ---------------------------------------------------------------------------
function remap_token(file, i, tok)
    key = (file, i, tok)
    haskey(OVERRIDES, key) && return OVERRIDES[key]
    m = match(TOKEN, tok)
    ch = parse(Int, m[1])
    ch > 16 && error("unmapped citation $tok in $file:$i (chapter > 16)")
    sub = m[2]
    new_ch = remap_chapter(ch)
    return sub === nothing ? "§$new_ch" : "§$new_ch.$sub"
end

function sweep(file)
    path = joinpath(DESIGN, file)
    lines = readlines(path; keep = true)
    before = 0
    used = Set{Tuple{String,Int,String}}()
    for (i, line) in pairs(lines)
        lines[i] = replace(line, TOKEN => function (tok)
            before += 1
            key = (file, i, tok)
            haskey(OVERRIDES, key) && push!(used, key)
            return remap_token(file, i, tok)
        end)
    end
    unused = setdiff(Set(k for k in keys(OVERRIDES) if k[1] == file), used)
    isempty(unused) || error("stale overrides for $file: $(sort(collect(unused)))")
    write(path, join(lines))
    after = sum(length(collect(eachmatch(TOKEN, l))) for l in lines)
    println(rpad(file, 42), " tokens: ", before, " → ", after)
    return before, after
end

# ---------------------------------------------------------------------------

function chapter_distribution()
    counts = Dict{Int,Int}()
    for f in CORPUS
        for line in readlines(joinpath(DESIGN, f))
            for m in eachmatch(TOKEN, line)
                ch = parse(Int, m[1])
                counts[ch] = get(counts, ch, 0) + 1
            end
        end
    end
    return counts
end

function print_distribution_table(before, after)
    chapters = sort(collect(union(keys(before), keys(after))))
    println(rpad("chapter", 10), rpad("before", 10), "after")
    for ch in chapters
        b = get(before, ch, 0)
        a = get(after, ch, 0)
        println(rpad(string(ch), 10), rpad(string(b), 10), a)
    end
end

function main()
    println("phase-4 reorder: block move + heading renumber + citation sweep\n")

    before_dist = chapter_distribution()   # captured before any mutation

    println("--- step 1: block move + numeral swap ---")
    move_part_block!()

    println("\n--- step 2: heading renumber ---")
    changes = rewrite_headings!()
    for (old, new) in changes
        println("  ", old, "  →  ", new)
    end
    println(length(changes), " heading lines rewritten")

    println("\n--- step 3: citation sweep ---")
    tb = ta = 0
    for f in CORPUS
        b, a = sweep(f)
        tb += b
        ta += a
    end
    println(rpad("TOTAL", 42), " tokens: ", tb, " → ", ta)

    println("\n--- chapter distribution (corpus-wide §N tokens) ---")
    after_dist = chapter_distribution()
    print_distribution_table(before_dist, after_dist)
end

main()
