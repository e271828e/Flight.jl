#!/usr/bin/env julia
#
# Phase-3 interior split of §9.3 (see tools/restructure_plan.md).
#
# One-shot. Splits "9.3 Inbound: root input slots and per-device staging" into
# three sections along its own internal structure:
#
#   ### 9.3 Inbound: root input slots, claims and the frozen roster
#         write surface, slot exclusivity, claims, the two surface classes,
#         slot initial values, the roster freeze, device identity and
#         admission, device death
#   ### 9.4 Inbound: per-device staging, representation and the drain
#         the staging cell and CAS merge, the attach-compiled representation,
#         staging-site diagnostics, the levels doctrine, the drain, mappings
#   ### 9.5 Inbound: the input trace
#         the trace, record density, the header, on-by-default, and the
#         closing note on rejected inbound shapes
#
# Old 9.4 (devices) → 9.6, old 9.5 (GUI write path) → 9.7.
#
# Citations are rewritten in ONE table-driven pass over all seven files (the
# shifts) plus per-occurrence re-adjudication of §9.3: a citation moves only
# where it clearly targets material that moved; ambiguous ones stay at §9.3,
# which remains valid.
#
# Run from anywhere:  julia docs/notes/design/tools/split_9_3.jl

const DESIGN = normpath(joinpath(@__DIR__, ".."))
const SPEC = "framework_spec.md"
const DEC = "framework_decisions.md"
const GAPS = "review3_2026-07-30_gaps.md"
const INC = "review3_2026-07-30_inconsistencies.md"

const FILES = [SPEC, DEC,
               "event_visibility_walkthrough.md",
               "inbound_periphery_walkthrough.md",
               GAPS,
               "review3_2026-07-30_improvements.md",
               INC]

# Uniform shifts (old 9.4/9.5 move down by two).
const SHIFT = Dict("§9.4" => "§9.6", "§9.5" => "§9.7")

# §9.3 citations whose target moved into the new §9.4 (staging, representation,
# staging-site diagnostics, levels, drain, mappings), keyed by (file, line).
const TO_9_4 = Set([
    (SPEC, 186),   # glossary-style aside: the per-device *staging cell*
    (SPEC, 1192),  # "input drain at frame top"
    (SPEC, 1534),  # "devices submit pending input writes … never touching live slots"
    (SPEC, 2119),  # "`map_input` by the staging checks"
    (SPEC, 2144),  # "§9.3's shared pure conditioning helper" (mappings)
    (SPEC, 2293),  # "§9.3's drain source"
    (SPEC, 3909),  # "mappings run on the device task"
    (SPEC, 4960),  # "the exercise that selected per-device cells"
    (SPEC, 5034),  # "conditioning upstream as mapping data"
    (SPEC, 5487),  # taught contract: "Levels, never deltas"
    (SPEC, 5595),  # "shape and normalization shim compiled"
    (SPEC, 5820),  # `EntryTypeMismatch` — staging-site only
    (SPEC, 5989),  # glossary "staging cell"
    (SPEC, 6220),  # glossary "batch"
    (SPEC, 6227),  # glossary "binding" — `TableBinding`'s generic `map_input`
    (SPEC, 6246),  # glossary "coalescing" — the CAS merge
    (SPEC, 6266),  # glossary "drain"
    (DEC, 119),    # row 103: the conditioning helper
    (DEC, 120),    # row 104: staged representation and coalescing
    (GAPS, 38), (GAPS, 106), (GAPS, 125), (GAPS, 136), (GAPS, 224),
    (GAPS, 364), (GAPS, 389),
    (INC, 76), (INC, 246), (INC, 255),
])

# §9.3 citations whose target moved into the new §9.5 (the input trace, record
# density, the header, the on-by-default rule).
const TO_9_5 = Set([
    (SPEC, 1045),  # "the trace holds staged inputs"
    (SPEC, 1587),  # "staging fed from a recording"
    (SPEC, 1617),  # trace/bit-identical replay
    (SPEC, 1636),  # "the plain on/off switch the trace has"
    (SPEC, 2571),  # "the entry point the §9.3 trace exists for"
    (SPEC, 2589),  # "§9.3's capture placement" (trace header)
    (SPEC, 2635),  # "the register's sparse records"
    (SPEC, 4400),  # trace-header capture placement
    (SPEC, 4423),  # "records the resolved pre-sequence stores and slots"
    (SPEC, 5051),  # "log/trace (§9.2–§9.3)" — see DASH_FIX below
    (SPEC, 5416),  # "header of initial stores and slot values"
    (SPEC, 5589),  # "never on the trace"
    (SPEC, 5708),  # glossary "trace" as reproduction tool
    (SPEC, 6359),  # "the log is recomputable from the trace"
    (SPEC, 6364),  # "before the boundary-zero sequence runs"
    (SPEC, 6370),  # "information rather than surface width"
    (DEC, 123),    # row 107: trace records match batch density
    (GAPS, 22), (GAPS, 35), (GAPS, 47), (GAPS, 68), (GAPS, 478),
])

# The one range citation whose endpoints now straddle the split: "log/trace
# (§9.2–§9.3)" becomes "log/trace (§9.2, §9.5)". Token count is preserved.
const DASH_FIX = "§9.2–§9.5" => "§9.2, §9.5"

const TOKEN = r"§(\d+)(?:\.(\d+))?"

function sweep!()
    to4 = Dict(f => 0 for f in FILES)
    to5 = Dict(f => 0 for f in FILES)
    kept = Dict(f => 0 for f in FILES)
    used = Set{Tuple{String,Int}}()
    for file in FILES
        path = joinpath(DESIGN, file)
        lines = readlines(path; keep = true)
        for (i, line) in pairs(lines)
            key = (file, i)
            if key in TO_9_4 || key in TO_9_5
                n = count(m -> m.match == "§9.3", collect(eachmatch(TOKEN, line)))
                n == 1 || error("$file:$i carries $n §9.3 tokens, expected 1")
            end
            lines[i] = replace(line, TOKEN => function (tok)
                if tok == "§9.3"
                    if key in TO_9_4
                        push!(used, key); to4[file] += 1; return "§9.4"
                    elseif key in TO_9_5
                        push!(used, key); to5[file] += 1; return "§9.5"
                    end
                    kept[file] += 1
                    return tok
                end
                return get(SHIFT, tok, tok)
            end)
        end
        file == SPEC && (lines[5051] = replace(lines[5051], DASH_FIX))
        write(path, join(lines))
    end
    stale = setdiff(union(TO_9_4, TO_9_5), used)
    isempty(stale) || error("stale re-pointings: $(sort(collect(stale)))")
    for file in FILES
        (to4[file] + to5[file] + kept[file]) == 0 && continue
        println("  ", rpad(file, 42), " §9.3: ", rpad(to4[file], 3), "→ §9.4, ",
                rpad(to5[file], 3), "→ §9.5, ", kept[file], " kept")
    end
    println("  total: ", sum(values(to4)), " → §9.4, ", sum(values(to5)),
            " → §9.5, ", sum(values(kept)), " kept")
end

function split_headings!()
    path = joinpath(DESIGN, SPEC)
    lines = readlines(path)
    head = findfirst(startswith("### 9.3 "), lines)
    seam_b = findfirst(startswith("**Staging: one atomic cell per attached device"), lines)
    seam_c = findfirst(startswith("**The input trace** is the sequence"), lines)
    (head === nothing || seam_b === nothing || seam_c === nothing) &&
        error("anchors not found")
    head < seam_b < seam_c || error("seams out of order")

    for (from, to) in ("### 9.5 " => "### 9.7 ", "### 9.4 " => "### 9.6 ")
        i = findfirst(startswith(from), lines)
        i === nothing && error("heading $from not found")
        lines[i] = to * lines[i][length(from)+1:end]
    end
    lines[head] = "### 9.3 Inbound: root input slots, claims and the frozen roster"
    splice!(lines, seam_c:seam_c-1, ["### 9.5 Inbound: the input trace", ""])
    splice!(lines, seam_b:seam_b-1,
            ["### 9.4 Inbound: per-device staging, representation and the drain", ""])
    write(path, join(lines, "\n") * "\n")
    println("  seams at lines ", seam_b, " and ", seam_c)
end

println("phase 3 — §9.3 split")
sweep!()
split_headings!()
