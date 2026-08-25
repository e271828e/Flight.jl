# --- helpers shared across the test files ---------------------------------------

const D8 = ForwardDiff.Dual{Nothing,Float64,8}

# The entries a phase body walks, per variant (§10.5): the boundary variant is
# the full list — its discrete entries wearing their `Gated` wrapper, unwrapped
# here — the interior one the continuous entries alone.
walked(body, variant = :boundary) =
    [e isa Gated ? e.e : e for c in getfield(body, variant) for e in c.entries]

# The gated entries of a body's boundary variant.
gated(body) = count(e isa Gated for c in body.boundary for e in c.entries)

# A model of one component, wrapped in the assembly every model's root must be.
single(c) = Group((; c = c), (), (), ())

# The same, with the component's one input face handed up to a root slot `in`.
fed(c, face) = Group((; c = c), (), ("in" => "children/c/$face",), ())

# The error a build raises, for the tests that read the message.
failure(f) =
    try
        f()
        nothing
    catch e
        e
    end
