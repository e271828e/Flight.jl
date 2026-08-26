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

# One writer's record in a snapshot's framework status (§11.8), by name: the
# devices as "device 1 (Pad)", the harness register as "harness", the loop as
# "loop".
writer_status(snap, who::String) =
    only(w for w in snap.status.writers if w.who == who)

# A device-task diagnostic's single account, wherever timing put it (§11.8,
# §12.4): the reporting task may or may not beat the run's last frame top, so
# the record is either in the terminal status's totals or presented by the
# run's-end sweep — exactly one of the two, never both, never neither.
accounted(sim, logs, who::String, field::Symbol, kind::String) =
    (getfield(writer_status(latest(sim), who).totals, field) ≥ 1) ⊻
    any(occursin("$kind from $who, past the final", string(l.message)) for l in logs)

crash_accounted(sim, logs, who::String) = accounted(sim, logs, who, :crash, "DeviceCrash")
