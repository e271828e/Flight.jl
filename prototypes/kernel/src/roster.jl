# The roster and claims (§11.3), with the attach-point slice of §11.6's device
# contract: the periphery's two roots, the declared sides, the enumeration
# contract with its bidirectional conformance check, the three-part admission,
# both claim sources, and the harness register as the derived remainder. The
# task side of the contract — the handle, the wrapper, the run's bracket and
# tail — lives in devices.jl; what stands further out in the spec —
# the trace — is deliberately absent (README); the binding conventions and
# the compiled gather live in bindings.jl.
#
# This file holds the types and the pure mechanics; `attach!`/`detach!` and
# the staging entry points live in sim.jl, beside the loop whose frames they
# configure.

"""
The periphery's two roots (§11.6): a device is a user type subtyping
`AbstractDevice`, a binding a value subtyping `AbstractBinding` — one
mandatory word each, buying `attach!`'s dispatch gate. Device identity is the
*instance*: the same object (`===`) may occupy at most one roster entry, and
two instances of one type — two joysticks — are two devices. (Which is why
the library's stub devices are mutable structs: an immutable stub would be
egal to its field-for-field twin, and `===` would stop meaning identity.)
The binding stays an `attach!` argument, never a device field: the same
device binds differently per deployment, and narrowing the binding narrows
the claim (§11.3).
"""
abstract type AbstractDevice end
abstract type AbstractBinding end

"""
The declared sides and the claim source (§11.6): sides are *declared*, never
inferred — the root carries the false defaults, so silence = absent — and
`is_greedy` selects the computed claim source *within* the input side (§11.3).
`needs_calling_task` is the affinity trait: the calling task is a single-slot
resource, so the roster admits at most one holder (§11.1, §11.3).
"""
is_input(::AbstractBinding) = false
is_output(::AbstractBinding) = false
is_greedy(::AbstractBinding) = false
needs_calling_task(::AbstractDevice) = false

"""
The enumeration contract (§11.6): an input-side binding's `claims(b)` is
called once at attach, and what it names is staked. The root's fallback is
error-throwing, so a declared side whose method was never written fails
loudly at the attach point rather than degrading into silence. `claims` never
returns a sentinel meaning "compute it for me": the enumeration contract has
one meaning and the `is_greedy` trait carries the other, so the maximal
surface stays reachable only through an explicit declaration.
"""
claims(b::AbstractBinding) = throw(BuildError(
    "BindingContractMismatch: $(typeof(b)) declares an input side and defines no " *
    "`claims` enumeration — the enumeration is the interface (§11.6)"))

"""
The output side's enumeration (§11.6, §14.4): an output-side binding's
`reads(b)` is called once at attach and compiled to one gather — the labeled
NamedTuple of selectors that fixes what `map_output` receives. The root's
fallback is error-throwing for the same reason `claims`'s is: a declared side
whose method was never written fails loudly at the attach point rather than
degrading into silence.
"""
reads(b::AbstractBinding) = throw(BuildError(
    "BindingContractMismatch: $(typeof(b)) declares an output side and defines no " *
    "`reads` enumeration — the enumeration is the interface (§11.6)"))

"""
§11.6's bidirectional conformance check over each (trait, method) pair, run
at the attach point: the same fact stated twice, in two registers, with the
framework paid to compare them. The declared-but-missing direction is the
`claims`/`reads` fallback firing when the attach calls it; the
defined-but-undeclared direction is one `which` against that fallback — the
reflection class, run at a stopped-sim service point, never inside a frame.
Greediness stays orthogonal to `reads`: a greedy front end driving a compiled
output gather is legal and currently uninstantiated (§11.6).
"""
function check_binding(b::AbstractBinding)
    T = typeof(b)
    isin, isout, greedy = is_input(b), is_output(b), is_greedy(b)
    drifted = which(claims, Tuple{T}) !== which(claims, Tuple{AbstractBinding})
    rdrifted = which(reads, Tuple{T}) !== which(reads, Tuple{AbstractBinding})
    greedy && !isin && throw(BuildError(
        "BindingContractMismatch: $T declares `is_greedy` without `is_input` — " *
        "greediness is a claim source within the input side, and a source without " *
        "its side is meaningless (§11.6)"))
    isin || isout || throw(BuildError(
        "BindingContractMismatch: $T declares neither side — a binding that touches " *
        "nothing is a configuration mistake, not a degenerate (§11.6)"))
    isin && greedy && drifted && throw(BuildError(
        "BindingContractMismatch: $T declares `is_greedy` and defines its own `claims` — " *
        "the two claim sources are alternatives, not layers (§11.6)"))
    isin || !drifted || throw(BuildError(
        "BindingContractMismatch: $T defines `claims` while `is_input` reads false — a " *
        "method written and never reached is exactly the drift this check catches (§11.6)"))
    isout || !rdrifted || throw(BuildError(
        "BindingContractMismatch: $T defines `reads` while `is_output` reads false — a " *
        "method written and never reached is exactly the drift this check catches (§11.6)"))
    nothing
end

"""
One roster entry (§11.3): the device instance, its binding, the stable id
assigned at `attach!` — monotonic per `Simulation`, never reused, living
exactly as long as the entry — the compiled writer over its claim set, the
entry's drain thunk, the per-attachment `should_abort` policy (§11.6: never
a device property — the same joystick is advisory in one deployment and
load-bearing in another), and the handle `attach!` constructed and returned.
Past the attach point nothing distinguishes a computed claim from a returned
one: the source is exhausted there, and validation, storage and the drain
treat the two identically.
"""
struct RosterEntry
    dev::AbstractDevice
    binding::AbstractBinding
    id::Int
    writer::Writer
    drain::Function                 # the compiled (cell, scatter) pair as a callable, see _drain_thunk
    should_abort::Bool              # the per-attachment failure policy (§11.6, §12.4)
    diag::DiagCell                  # the device's diagnostic cell (§11.8), shared with the handle
    acct::WriterAccount             # the loop's private account behind the cell (§11.8)
    handle::Any                     # the DeviceHandle; `Any` for include order only, read off the frame path
end

_who(e::RosterEntry) = "device $(e.id) ($(typeof(e.dev)))"

# The stopped-sim compile of one writer's drain (§11.4): a zero-argument thunk
# capturing the store and the writer *concretely* — this dynamic dispatch is
# what specializes the closure — so the frame top's call boxes no argument and
# the drain, empty or populated, stays allocation-free (D-202). The per-entry
# dispatch this leaves is
# the iteration §11.4 licenses in place of the tuple specialization the freeze
# would also permit.
_drain_thunk(store, w::Writer) = () -> _drain!(store, w)

"""
The data plane's mutable holder: the roster in attachment order — which is the
drain's application order — the harness register's writer and drain thunk, the
exclusivity index behind the `ClaimedFaceEntry` payload, the store the thunks
compile against, and the id counter; the §11.3 freeze itself is the
lifecycle's `:running` state (devices.jl). Mutable and
abstractly typed deliberately: the harness writer's *type* changes at every
roster change (its schema is recompiled), so it cannot be a `Simulation` type
parameter, and everything here is stopped-sim configuration read behind
function barriers.

The harness register is a writer, so it owns a diagnostic cell and its
account like any other (§11.8, D-200: the harness register is a
diagnostic writer, its staging diagnostics needing a single-writer home).
Unlike the writer, the cell
survives roster changes: its diagnostics are the harness's, whatever surface
it speaks for. `run_tasks` is the run's device-id → `Task` registry, filled
at spawn and emptied at run end, which is what lets `publish!` read
`task_state` off the handles the loop owns (D-193); while stopped it is
empty, and every device reads `:none` — device tasks are run-scoped
observables (§12.4).
"""
mutable struct DataPlane
    roster::Vector{RosterEntry}     # attachment order (§11.3): the drain applies in it
    harness::Writer                 # the derived-surface writer, recompiled at roster changes
    harness_drain::Function         # its drain thunk, recompiled with it
    harness_diag::DiagCell          # the harness register's diagnostic cell (§11.8)
    harness_acct::WriterAccount     # and the loop's account behind it
    run_tasks::Dict{Int,Task}       # the run's device tasks, by device id (§12.2, D-193)
    claimedby::Dict{Symbol,String}  # face → incumbent: the exclusivity index
    store::Any                      # the model's store bundle, captured into the drain thunks
    next_id::Int
end

function DataPlane(layout::Layout, store)
    w = Writer(layout, Symbol[f for (f, _) in layout.root_inputs])
    DataPlane(RosterEntry[], w, _drain_thunk(store, w), DiagCell(EMPTY_DIAG),
              WriterAccount(), Dict{Int,Task}(), Dict{Symbol,String}(), store, 1)
end

"""
The claim, from its source (§11.3): computed — `is_greedy`, the unclaimed
complement at this instant, disjoint from every incumbent by construction and
never recomputed, which is what makes attachment order load-bearing — or
returned, the enumeration called once, each face validated against the root
face set (`AttachUnknownFace`: a mapping drifted onto a nonexistent face is a
diagnosable anomaly, never a silent write). Duplicates within one enumeration
collapse; the empty enumeration is the honest may-write-nothing degenerate.
"""
function _claim(plane::DataPlane, layout::Layout, b::AbstractBinding)
    faceset = Symbol[f for (f, _) in layout.root_inputs]
    is_greedy(b) && return Symbol[f for f in faceset if !haskey(plane.claimedby, f)]
    claim = Symbol[]
    for f in claims(b)
        s = Symbol(f)
        s in faceset || throw(BuildError(
            "AttachUnknownFace: $(typeof(b)) claims `$s`, which names no root input " *
            "face — the root faces are $(_facelist(faceset)) (§11.3)"))
        s in claim || push!(claim, s)
    end
    claim
end

"""
Recompute the exclusivity index and recompile the harness register's writer
after a roster change (§11.3, §11.4): the harness surface is the unclaimed
complement — the faces no rostered device speaks for — so it moves exactly
when the roster does, and is as fixed within a run as any claim set. Being
always present, the harness cell gets the compilation unasked.

The recompilation has one seam (§11.4): a pending harness batch staged before
a stopped-sim `attach!` may hold the old shape, or may name a face the new
claim covers. The batch is taken, reshaped through the new schema and
re-staged — newly claimed faces discarded with `ClaimedFaceEntry` — so the run
always starts with cells matching the run's schemas. At `detach!` the surface
only broadens and every pending entry survives the reshape.
"""
function reclaim!(plane::DataPlane, layout::Layout)
    empty!(plane.claimedby)
    for e in plane.roster, f in e.writer.faces
        plane.claimedby[f] = _who(e)
    end
    old = plane.harness
    pending = @atomicswap old.cell.pending = nothing
    w = Writer(layout, Symbol[f for (f, _) in layout.root_inputs if !haskey(plane.claimedby, f)])
    plane.harness = w
    plane.harness_drain = _drain_thunk(plane.store, w)
    if pending !== nothing
        batch = pending[]
        entries = [old.faces[i] => batch.vals[i] for i in 1:length(old.faces) if batch.mask[i]]
        renorm = _normalize(plane.harness, entries, plane.claimedby, plane.harness_diag;
                            site = :renormalization)
        renorm === nothing || _stage!(plane.harness, renorm)
    end
    nothing
end
