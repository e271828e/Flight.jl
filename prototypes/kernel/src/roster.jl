# The roster and claims (§11.3), with the attach-point slice of §11.6's device
# contract: the periphery's two roots, the declared sides, the enumeration
# contract with its bidirectional conformance check, the three-part admission,
# both claim sources, and the harness register as the derived remainder. The
# task side of the contract — the handle, the wrapper, the run's bracket and
# tail — lives in devices.jl, the trace whose schema list the roster's writers
# grow in trace.jl; the binding conventions and the compiled gather live in
# bindings.jl.
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
    BindingContractMismatch(binding = string(typeof(b)), reason = :claims_missing)))

"""
The output side's enumeration (§11.6, §14.4): an output-side binding's
`reads(b)` is called once at attach and compiled to one gather — the labeled
NamedTuple of selectors that fixes what `map_output` receives. The root's
fallback is error-throwing for the same reason `claims`'s is: a declared side
whose method was never written fails loudly at the attach point rather than
degrading into silence.
"""
reads(b::AbstractBinding) = throw(BuildError(
    BindingContractMismatch(binding = string(typeof(b)), reason = :reads_missing)))

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
        BindingContractMismatch(binding = string(T), reason = :greedy_without_input)))
    isin || isout || throw(BuildError(
        BindingContractMismatch(binding = string(T), reason = :neither_side)))
    isin && greedy && drifted && throw(BuildError(
        BindingContractMismatch(binding = string(T), reason = :greedy_with_claims)))
    isin || !drifted || throw(BuildError(
        BindingContractMismatch(binding = string(T), reason = :claims_without_input)))
    isout || !rdrifted || throw(BuildError(
        BindingContractMismatch(binding = string(T), reason = :reads_without_output)))
    nothing
end

"""
§11.6's device twin of `check_binding`, run at the same attach point: a
device with no `loop` method of its own is a build-time configuration
mistake — `loop`'s fallback exists only as the comparison target this `which`
check needs, and refusing here (before the device is ever rostered, let alone
spawned) is what makes `DeviceContractMismatch` `service, fail-fast` rather
than a task-side `DeviceCrash` (Appendix C).
"""
function check_device(dev::AbstractDevice)
    T = typeof(dev)
    # against the handle type the wrapper calls with: a `loop(::T, ::DeviceHandle)`
    # is a method, and `Tuple{T,Any}` would not see it
    which(loop, Tuple{T,DeviceHandle}) === which(loop, Tuple{AbstractDevice,DeviceHandle}) &&
        throw(BuildError(DeviceContractMismatch(device = string(T), reason = :no_loop)))
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
# would also permit. The trace register and this writer's index into its schema
# list ride along (§11.5), which is why every thunk is compiled at one site
# alone — `_install_writers!` (trace.jl) — and recompiled whenever the indices
# move.
_drain_thunk(store, w::Writer, reg, widx::Int) = () -> _drain!(store, w, reg, widx)

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

# The register is an argument because the drain thunks close over it (§11.5):
# with the roster empty the harness register is the sole writer, index 1 of the
# first set the header captures, and `_install_writers!` re-fixes that at the
# capture and at every roster change.
function DataPlane(layout::Layout, store, reg)
    w = Writer(layout, Symbol[f for (f, _) in layout.root_inputs])
    DataPlane(RosterEntry[], w, _drain_thunk(store, w, reg, 1), DiagCell(EMPTY_DIAG),
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
            AttachUnknownFace(binding = string(typeof(b)), face = s, candidates = faceset)))
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

The roster change is also where the trace's schema list grows (§11.5): the
current writer set is appended and every drain thunk recompiled against its
new index, which is what `_install_writers!` does at the tail here.
"""
function reclaim!(plane::DataPlane, layout::Layout, reg)
    empty!(plane.claimedby)
    for e in plane.roster, f in e.writer.faces
        plane.claimedby[f] = _who(e)
    end
    old = plane.harness
    pending = @atomicswap old.cell.pending = nothing
    w = Writer(layout, Symbol[f for (f, _) in layout.root_inputs if !haskey(plane.claimedby, f)])
    plane.harness = w
    _install_writers!(reg, plane)          # §11.5: the schema list grows, the thunks follow
    if pending !== nothing
        batch = pending[]
        entries = [old.faces[i] => batch.vals[i] for i in 1:length(old.faces) if batch.mask[i]]
        renorm = _normalize(plane.harness, entries, plane.claimedby, plane.harness_diag;
                            site = :renormalization)
        renorm === nothing || _stage!(plane.harness, renorm)
    end
    nothing
end
