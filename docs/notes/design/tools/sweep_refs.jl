#!/usr/bin/env julia
#
# Phase-2 citation sweep for the framework_spec.md restructure.
#
# One-shot, table-driven rewrite of every `§N` / `§N.M` citation in the spec and
# its six companion files, from the pre-restructure numbering to the phase-2
# numbering defined in tools/restructure_plan.md.
#
# Mechanics: every token is matched ONCE by a single regex pass and mapped
# through the table below. Sequential search-and-replace would double-map
# (old §13 → new §11 while old §11 → new §8).
#
# Run from anywhere:  julia docs/notes/design/tools/sweep_refs.jl
#
# Appendix citations (§D.4 …) are not matched and stay as they are.

const DESIGN = normpath(joinpath(@__DIR__, ".."))

const FILES = ["framework_spec.md",
               "framework_decisions.md",
               "event_visibility_walkthrough.md",
               "inbound_periphery_walkthrough.md",
               "review3_2026-07-30_gaps.md",
               "review3_2026-07-30_improvements.md",
               "review3_2026-07-30_inconsistencies.md"]

# Per-occurrence adjudications, keyed by (file, 1-based line number in the
# pre-sweep file, matched token). Bare §12 is ambiguous after old §12 split
# into new §9 (data plane) and §10 (lifecycle and orchestration); bare §18 is
# dissolved into §1.
const OVERRIDES = Dict{Tuple{String,Int,String},String}(
    # "the staging rules settled in §12" — staging is the data plane.
    ("framework_spec.md", 62, "§12") => "§9",
    # "The rule, binding for the periphery (§12): external readers observe the
    # signal table only at step boundaries." — a read/write-path rule.
    ("framework_spec.md", 1123, "§12") => "§9",
    # §12.1–§12.6 range in the old §12 lede: the endpoints now straddle the
    # chapter split, so the range is respelled as chapter 9 plus §10.1 (the
    # trailing endpoint is mapped by the table; see the lede stitch in
    # phase2_restructure.jl, which turns the dash into a comma).
    ("framework_spec.md", 1500, "§12.1") => "§9",
    # "the staging rules — and the concurrency model generally — belong to §12."
    # Staging is §9, the concurrency model is §10.
    ("framework_spec.md", 1493, "§12") => "§9, §10",
    # heading "Torture test for the §12 staging shapes" — staging, so §9.
    ("framework_spec.md", 4938, "§12") => "§9",
    # "the §12 periphery decisions were forced by this cast" — the cast spans
    # devices, control plane, shutdown and lifecycle: the whole periphery.
    ("framework_spec.md", 5033, "§12") => "§9, §10",
    # glossary "periphery": staged writes, snapshot reads AND the control
    # surface — the whole periphery.
    ("framework_spec.md", 6288, "§12") => "§9, §10",
    # glossary "torture test": "filter/joystick/GUI against §12's staging
    # shapes" — staging, so §9.
    ("framework_spec.md", 6544, "§12") => "§9",
    # glossary "row": the decision-log pointer, now carried by §1.
    ("framework_spec.md", 6535, "§18") => "§1",
    # "The defects cluster in the periphery (§12)" — the whole periphery.
    ("review3_2026-07-30_inconsistencies.md", 499, "§12") => "§9, §10",
)

"Map one old citation to its phase-2 spelling. `sub === nothing` means a bare
chapter citation."
function remap(ch::Int, sub::Union{Int,Nothing})
    ch <= 5 && return sub === nothing ? "§$ch" : "§$ch.$sub"
    if ch in (6, 7, 8, 10)
        sub === nothing || error("old §$ch has no subsections, got §$ch.$sub")
        return ch == 6 ? "§6.2" : ch == 7 ? "§4.4" : ch == 8 ? "§6.1" : "§5.5"
    end
    ch == 9  && return sub === nothing ? "§7" : "§7.$sub"
    ch == 11 && return sub === nothing ? "§8" : "§8.$sub"
    if ch == 12
        sub === nothing && error("bare §12 needs an adjudication")
        return sub <= 5 ? "§9.$sub" : "§10.$(sub - 5)"
    end
    ch == 18 && error("§18 needs an adjudication")
    ch == 19 && return "§16"
    if 13 <= ch <= 17                       # 13→11, 14→12, 15→13, 16→14, 17→15
        new = ch - 2
        return sub === nothing ? "§$new" : "§$new.$sub"
    end
    error("unmapped citation §$ch")
end

const TOKEN = r"§(\d+)(?:\.(\d+))?"

function sweep(file)
    path = joinpath(DESIGN, file)
    lines = readlines(path; keep = true)
    before = 0
    used = Set{Tuple{String,Int,String}}()
    for (i, line) in pairs(lines)
        lines[i] = replace(line, TOKEN => function (tok)
            before += 1
            key = (file, i, tok)
            if haskey(OVERRIDES, key)
                push!(used, key)
                return OVERRIDES[key]
            end
            m = match(TOKEN, tok)
            ch = parse(Int, m[1])
            sub = m[2] === nothing ? nothing : parse(Int, m[2])
            return remap(ch, sub)
        end)
    end
    unused = setdiff(Set(k for k in keys(OVERRIDES) if k[1] == file), used)
    isempty(unused) || error("stale overrides for $file: $(sort(collect(unused)))")
    write(path, join(lines))
    after = sum(length(collect(eachmatch(TOKEN, l))) for l in lines)
    println(rpad(file, 42), " tokens: ", before, " → ", after)
    return before, after
end

function main()
    println("phase-2 citation sweep")
    tb = ta = 0
    for f in FILES
        b, a = sweep(f)
        tb += b; ta += a
    end
    println(rpad("TOTAL", 42), " tokens: ", tb, " → ", ta)
end

main()
