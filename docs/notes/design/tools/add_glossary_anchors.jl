#!/usr/bin/env julia
#
# One-shot for rewrite-plan Phase 1.2: give every Appendix D glossary entry an
# HTML anchor, same-line, so body text can link terms:
#
#     <a id="g-drain"></a>**drain** — the frame-top …
#
# The id is derived from the entry's primary alias (the part before " / "),
# with backticks and a trailing parenthetical stripped and Δ transliterated to
# `d`, then slugified exactly as the heading slugger does. Entry anchors only:
# bold sub-terms defined inside an entry (interior/boundary sweep, the bundle
# law) are reached through their owning entry, per tools/coinage_inventory.md
# bucket B.
#
# Refuses to run if any `g-` anchor already exists (it is a one-shot; the
# anchors live in the source afterwards and check_glossary.jl guards them).
# Errors on id collisions before writing anything.
#
# Run from anywhere:  julia docs/notes/design/tools/add_glossary_anchors.jl

include(joinpath(@__DIR__, "slugs.jl"))

const SPEC = normpath(joinpath(@__DIR__, "..", "framework_spec.md"))

"Anchor id for an entry term."
function anchorid(term)
    t = split(term, " / ")[1]
    t = replace(t, "`" => "", "Δ" => "d")
    t = strip(replace(t, r"\s*\([^)]*\)$" => ""))
    return "g-" * slugify(t)
end

lines = readlines(SPEC)
start = findfirst(startswith("## Appendix D"), lines)
start === nothing && error("Appendix D heading not found")

any(occursin("<a id=\"g-", l) for l in lines) &&
    error("g- anchors already present — one-shot already applied")

seen = Dict{String,String}()
n = 0
for i in start:length(lines)
    m = match(r"^\*\*(.+?)\*\*(?= —| –|—)", lines[i])
    m === nothing && continue
    id = anchorid(m[1])
    haskey(seen, id) && error("id collision: \"$(m[1])\" and \"$(seen[id])\" both → $id")
    seen[id] = m[1]
    lines[i] = "<a id=\"$id\"></a>" * lines[i]
    global n += 1
end

write(SPEC, join(lines, "\n") * "\n")
println("anchored $n entries")
