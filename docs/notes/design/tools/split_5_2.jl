#!/usr/bin/env julia
#
# Phase-3 interior split of §5.2 (see tools/restructure_plan.md).
#
# One-shot. Splits "5.2 Two-stage outputs and structural feedthrough" at the
# seam before the "**The letters**" paragraph:
#
#   ### 5.2 Two-stage outputs: signatures, bundles and the hand-off laws
#         the two output stages, the signature block, the named-bundle hand-off,
#         the bundle law, the handler return law, the views
#   ### 5.3 Structural feedthrough: stage roles, schedule and step boundaries
#         the letters and their dependence classes, the per-stage roles, the
#         schedule, step-boundary semantics, the orthodox departure, shared
#         expensive computations
#
# Old 5.3 → 5.4, 5.4 → 5.5, 5.5 → 5.6.
#
# Citations are rewritten in ONE table-driven pass over all seven files (the
# shifts) plus a per-occurrence override table (the §5.2 re-adjudications: a
# citation moves to §5.3 only where it clearly targets material that moved;
# ambiguous ones stay at §5.2, which remains valid).
#
# Run from anywhere:  julia docs/notes/design/tools/split_5_2.jl

const DESIGN = normpath(joinpath(@__DIR__, ".."))
const SPEC = "framework_spec.md"

const FILES = [SPEC,
               "framework_decisions.md",
               "event_visibility_walkthrough.md",
               "inbound_periphery_walkthrough.md",
               "review3_2026-07-30_gaps.md",
               "review3_2026-07-30_improvements.md",
               "review3_2026-07-30_inconsistencies.md"]

# Uniform shifts, applied to every §-token exactly once.
const SHIFT = Dict("§5.3" => "§5.4", "§5.4" => "§5.5", "§5.5" => "§5.6")

# Per-occurrence re-adjudication of §5.2, keyed by (file, pre-edit line number).
# Every entry points at material that moved into the new §5.3; each line holds
# exactly one §5.2 token (asserted below).
const REPOINT = Dict{Tuple{String,Int},String}(
    (SPEC, 130)  => "guards/handlers read the fresh boundary table",
    (SPEC, 1059) => "the boundary sequence",
    (SPEC, 1125) => "the boundary sweep in the sequence",
    (SPEC, 1308) => "the event-iteration question deferred to §8.6",
    (SPEC, 1330) => "the cross-component same-snapshot claim",
    (SPEC, 2765) => "auto-publication at stage-1 position",
    (SPEC, 2775) => "event declaration order",
    (SPEC, 2920) => "event declaration order",
    (SPEC, 2986) => "auto-publication of unproduced state-matching outputs",
    (SPEC, 4889) => "the law runs once in a stage, `f`/`g` copy",
    (SPEC, 4909) => "the derivative/output-overlap claim (orthodox departure)",
    (SPEC, 5311) => "auto-published state at stage-1 position",
    (SPEC, 5474) => "handler-phase visibility",
    (SPEC, 5896) => "bare `h` is the step size (the letters)",
    (SPEC, 5951) => "stage-2 outputs presumed input-dependent",
    (SPEC, 6013) => "no incremental `f`-only re-evaluation (the schedule)",
    (SPEC, 6033) => "iteration to quiescence",
    ("framework_decisions.md", 116) => "row 100: within-round event visibility",
    ("event_visibility_walkthrough.md", 4)   => "step-boundary semantics",
    ("event_visibility_walkthrough.md", 12)  => "the same-boundary-snapshot answer",
    ("event_visibility_walkthrough.md", 182) => "\"their own component's refreshed ports\"",
    ("review3_2026-07-30_gaps.md", 244) => "event declaration order",
    ("review3_2026-07-30_gaps.md", 257) => "the same-boundary-snapshot claim",
    ("review3_2026-07-30_gaps.md", 285) => "step-boundary semantics",
    ("review3_2026-07-30_inconsistencies.md", 418) => "\"Departure from the orthodox formalism\"",
    ("review3_2026-07-30_inconsistencies.md", 423) => "the `h_x`-may-be-absent bullet",
)

const TOKEN = r"§(\d+)(?:\.(\d+))?"

function sweep!()
    repointed = Dict(f => 0 for f in FILES)
    kept = Dict(f => 0 for f in FILES)
    used = Set{Tuple{String,Int}}()
    for file in FILES
        path = joinpath(DESIGN, file)
        lines = readlines(path; keep = true)
        for (i, line) in pairs(lines)
            if haskey(REPOINT, (file, i))
                n = count(m -> m.match == "§5.2", collect(eachmatch(TOKEN, line)))
                n == 1 || error("$file:$i carries $n §5.2 tokens, expected 1")
            end
            lines[i] = replace(line, TOKEN => function (tok)
                if tok == "§5.2"
                    if haskey(REPOINT, (file, i))
                        push!(used, (file, i))
                        repointed[file] += 1
                        return "§5.3"
                    end
                    kept[file] += 1
                    return tok
                end
                return get(SHIFT, tok, tok)
            end)
        end
        write(path, join(lines))
    end
    stale = setdiff(Set(keys(REPOINT)), used)
    isempty(stale) || error("stale re-pointings: $(sort(collect(stale)))")
    for file in FILES
        (repointed[file] + kept[file]) == 0 && continue
        println("  ", rpad(file, 42), " §5.2: ", rpad(repointed[file], 3),
                "→ §5.3, ", kept[file], " kept")
    end
    println("  total re-pointed: ", sum(values(repointed)),
            ", kept: ", sum(values(kept)))
end

function split_heading!()
    path = joinpath(DESIGN, SPEC)
    lines = readlines(path)
    old = findfirst(startswith("### 5.2 "), lines)
    old === nothing && error("§5.2 heading not found")
    seam = findfirst(l -> startswith(l, "**The letters**"), lines)
    seam === nothing && error("seam paragraph not found")
    old < seam || error("seam is not inside §5.2")

    lines[old] = "### 5.2 Two-stage outputs: signatures, bundles and the hand-off laws"
    for (from, to) in ("### 5.5 " => "### 5.6 ", "### 5.4 " => "### 5.5 ",
                       "### 5.3 " => "### 5.4 ")
        i = findfirst(startswith(from), lines)
        i === nothing && error("heading $from not found")
        lines[i] = to * lines[i][length(from)+1:end]
    end
    splice!(lines, seam:seam-1,
            ["### 5.3 Structural feedthrough: stage roles, schedule and step boundaries", ""])
    write(path, join(lines, "\n") * "\n")
    println("  split at line ", seam, " (\"**The letters**\")")
end

println("phase 3 — §5.2 split")
sweep!()
split_heading!()
