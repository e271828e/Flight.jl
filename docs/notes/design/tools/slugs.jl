# Shared heading/slug machinery for the design-doc tooling.
#
# `included` by tools/linkify.jl and tools/check_refs.jl so both compute the
# same anchors. Plain stdlib.
#
# The slugger reproduces GitHub's: strip markdown formatting, lowercase, drop
# punctuation and symbols (keeping letters, digits, spaces, `-` and `_`), then
# spaces → hyphens. Duplicates get `-1`, `-2`, … in document order, exactly as
# GitHub's slugger does.

"GitHub-style anchor slug for a heading's text."
function slugify(text::AbstractString)
    s = replace(text, r"`([^`]*)`" => s"\1")            # inline code
    s = replace(s, r"\*\*([^*]+)\*\*" => s"\1")         # bold
    s = replace(s, r"\*([^*]+)\*" => s"\1")             # italics
    s = replace(s, r"\[([^\]]*)\]\([^)]*\)" => s"\1")   # links → their text
    s = lowercase(s)
    s = replace(s, r"[^\p{L}\p{Nd}\p{M} _-]" => "")
    return replace(s, ' ' => '-')
end

"""
    headings(path) -> Vector{NamedTuple}

Every ATX heading outside fenced code blocks, in document order, as
`(level, text, number, slug)`. `number` is the citable number — `"9"`, `"9.4"`,
`"A"` (an appendix), `"D.7"` — or `nothing` for unnumbered headings (the title,
the parts, `## Contents`, every `####`).
"""
function headings(path)
    out = NamedTuple[]
    seen = Dict{String,Int}()
    infence = false
    for line in eachline(path)
        startswith(line, "```") && (infence = !infence; continue)
        infence && continue
        m = match(r"^(#{1,6}) +(.*?)\s*$", line)
        m === nothing && continue
        level, text = length(m[1]), String(m[2])
        number = nothing
        for (re, f) in ((r"^(\d+)\.(\d+) ", g -> "$(g[1]).$(g[2])"),
                        (r"^(\d+)\. ", g -> String(g[1])),
                        (r"^Appendix ([A-D])[.:] ", g -> String(g[1])),
                        (r"^([A-D])\.(\d+) ", g -> "$(g[1]).$(g[2])"))
            g = match(re, text)
            g === nothing || (number = f(g); break)
        end
        base = slugify(text)
        n = get(seen, base, 0)
        seen[base] = n + 1
        push!(out, (; level, text, number, slug = n == 0 ? base : "$base-$n"))
    end
    return out
end

"Map from citation key (`\"9.4\"`, `\"A\"`, `\"D.7\"`) to anchor slug."
function targets(hs)
    d = Dict{String,String}()
    for h in hs
        h.number === nothing && continue
        haskey(d, h.number) && error("two headings claim number $(h.number)")
        d[h.number] = h.slug
    end
    return d
end

"Headings whose slug needed a `-N` suffix (empty is the good case)."
collisions(hs) = [(h.text, h.slug) for h in hs if h.slug != slugify(h.text)]
