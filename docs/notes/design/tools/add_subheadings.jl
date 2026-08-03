#!/usr/bin/env julia
#
# Phase-3 unnumbered `####` headings (see tools/restructure_plan.md).
#
# One-shot. Inserts navigation headings inside the three long sections the plan
# names: the devices section (§9.6), the declaration inventory (§11.2) and the
# two big case studies (§15.4, §15.5). Every heading is derived from the text's
# own structure — the bold lead-in or the bullet's declaration name — and no
# body text is touched. `####` headings are unnumbered, so they carry no
# citation impact and do not enter the §N.M outline.
#
# Run from anywhere:  julia docs/notes/design/tools/add_subheadings.jl

const SPEC = normpath(joinpath(@__DIR__, "..", "framework_spec.md"))

# anchor (a unique line prefix) => heading text inserted just above it
const HEADINGS = [
    # §9.6 Devices
    "**Every attached device receives the same handle**" =>
        "Every attached device receives the same handle",
    "**The authoring contract: four functions" =>
        "The authoring contract: four functions, one optional, one trait",
    "**The author owns the loop body" =>
        "The author owns the loop body; the framework owns the bracket",
    "**The binding: framework-legible" =>
        "The binding: framework-legible by enumeration, opaque in its mappings",
    "**One shipped binding type" =>
        "One shipped binding type; conditioning has an owner",
    "**Bad datum versus bug" =>
        "Bad datum versus bug: two classes, two fates",
    # §11.2 The declaration inventory — one heading per declaration group
    "- **State, modes, discrete state**" => "State, modes, discrete state",
    "- **`input_types(::C)`**"           => "`input_types(::C)`",
    "- **`output_types(::C)`**"          => "`output_types(::C)`",
    "- **`local_types(::C)`**"           => "`local_types(::C)`",
    "- **`events(::C)`**"                => "`events(::C)`",
    "- **No stage tags anywhere.**"      => "No stage tags anywhere",
    "- **Custom structs as port types**" => "Custom structs as port types",
    "- **Completeness of the declaration set**" =>
        "Completeness of the declaration set",
    # §15.4 The interactive C172X demo
    "**Architectures examined here and rejected**" =>
        "Architectures examined here and rejected",
    "**Surface walkthrough.**" => "Surface walkthrough",
    "**Frame anatomies.**"    => "Frame anatomies",
    # §15.5 The strapdown IMU
    "**Why the merge buys nothing.**" => "Why the merge buys nothing",
    "**The counterexample**"          => "The counterexample",
    "**The idiom: integrate-and-difference.**" =>
        "The idiom: integrate-and-difference",
    "**Exactness condition, stated once**" => "Exactness condition, stated once",
    "**Why `u.V` is fresh" => "Why `u.V` is fresh — the line that would silently zero",
    "**Author-knowledge note**" => "Author-knowledge note",
]

lines = readlines(SPEC)

sites = map(HEADINGS) do (anchor, heading)
    idx = findall(startswith(anchor), lines)
    length(idx) == 1 || error("$(length(idx)) lines start with $(repr(anchor))")
    (only(idx), heading)
end

for (i, heading) in sort(sites; by = first, rev = true)   # bottom-up: indices hold
    pad = isempty(strip(lines[i-1])) ? String[] : [""]
    splice!(lines, i:i-1, [pad; "#### $heading"; ""])
end

write(SPEC, join(lines, "\n") * "\n")
println("inserted ", length(sites), " `####` headings")
