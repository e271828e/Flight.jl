# Cell-store representation bench

Increment 1 of the kernel prototype. Settles the one representation question
the spec deliberately left open: **how the signal table's cells are stored and
addressed** by the compiled executor (§12.7). Everything else on the path is
already decided and is held fixed here — flat state backing with compile-time
offsets (§7.1), scatter/gather as the whole protocol (§4.3), a concretely-typed
entry tuple traversed by a chunked unrolled walk (§12.7).

## The two candidates

Both are pure; the hybrid was rejected.

- **C1 — heterogeneous store, type-domain indices.** One store enumerating every
  cell in its own type; an entry carries each cell's position as a *type*
  parameter, so access is `getfield` at a constant. Risk: the store type grows
  with the model and appears in every entry's signature, so N identical
  components may still compile N distinct bodies.
- **C2 — per-eltype homogeneous stores, field offsets.** Cells are flattened by
  the §7.1 leaf walk into one contiguous buffer per element type; an entry
  carries plain integer offsets in *fields*, with the port type as the address
  token's parameter. Hoped-for win: two instances of one component type differ
  only in offset values, so one compiled body serves all N. Risk (the reason
  this is measured, not argued): inlining and const-prop may respecialize on
  the offsets anyway and silently restore N bodies.

The bench isolates exactly that axis. `src/protocol.jl` holds the entry shape,
the gather/scatter/stage-call body, the chunked walk and the sweep builder —
identical for both. A candidate supplies only `build_store`, `cell_addr`,
`gather`, `scatter!` and `snapshot`.

## Gates, in order

1. **`@ballocated(sweep()) == 0`** — mandatory. `sweep()` is the zero-arg
   interior variant of `sweep_hxu` (§12.7). A candidate that cannot pass is out
   whatever else it does.
2. **Compile cost vs. N identical instances** — the discriminating measurement.
   A flat curve for C2 means the code sharing is real; a slope in the class of
   C1's means it is not, and C2's whole argument goes with it.
3. **Tiebreakers** — sweep runtime, and snapshot-capture cost and allocation
   (§9.2 publishes by copying out of the store, so the store's shape prices it).

The decision row is written after the numbers, not before. Permission to
discard the whole prototype is reserved.

## Result: C2, on every gate (2026-08-08)

Apple Silicon, Julia 1.12.6, chunk 16, one cold process per point. Raw data in
`results/` (untracked); regenerate with `run_bench.jl`.

| cand | scalar | N | `run!` bodies | compile s | sweep ns | ns/entry | alloc | snap ns |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| C1 | Float64 | 25 | 28 | 1.87 | 334 | 11.9 | 0 | 52 |
| C1 | Float64 | 100 | 103 | 10.16 | 2616 | 25.4 | 0 | 235 |
| C1 | Float64 | 400 | 403 | 56.18 | 19709 | 48.9 | 0 | 1138 |
| C2 | Float64 | 25 | 4 | 0.89 | 150 | 5.4 | 0 | 56 |
| C2 | Float64 | 100 | 4 | 1.00 | 722 | 7.0 | 0 | 278 |
| C2 | Float64 | 400 | 4 | 1.14 | 2847 | 7.1 | 0 | 1129 |
| C1 | Dual8 | 50 | 53 | 9.74 | 3271 | 61.7 | 0 | 1483 |
| C1 | Dual8 | 200 | 203 | 35.11 | 17250 | 85.0 | 0 | 8222 |
| C2 | Dual8 | 50 | 4 | 8.58 | 2310 | 43.6 | 0 | 1654 |
| C2 | Dual8 | 200 | 4 | 10.04 | 10000 | 49.3 | 0 | 8445 |
| C2 | Dual8 | 400 | 4 | 8.26 | 19375 | 48.1 | 0 | 16958 |

- **Gate 1** ties: zero allocation for both candidates, sweep and snapshot, at
  every N and both activations.
- **Gate 2** decides: C1's nominal compile is superlinear (0.60 → 56.2 s over
  N = 1 → 400); C2's is flat (0.64 → 1.14 s). The sharing is structural, not a
  timing artifact — compiled MethodInstances at N = 100 are 4 `run!` / 3
  `gather` / 3 `scatter!` for C2 against 103 / 105 / 306 for C1. That count is
  also what closes the const-prop risk: offsets survive as runtime data across
  `invoke_sweep`'s barrier, so no per-instance respecialization happens.
- On `Dual8` the *ratio* narrows (3.5× at N = 200) but the shape is the point:
  C2 saturates near 9 s from N ≈ 50 — bounded by chunk-type count, not by model
  size — while C1 keeps climbing. C2's nominal and Dual activations for a
  C172X-scale model both land in single-digit seconds, which makes §12.7's
  mitigation ladder optional rather than load-bearing.
- **Gate 3** follows: C2's per-entry cost is flat (~7.1 ns nominal, ~48 ns
  Dual); C1's degrades with model size (5.0 → 48.9 ns nominal) as 400+ distinct
  bodies and a pointer chase per gather bite. Snapshot runtime ties.

Recorded as decision row 162; §12.7 amended.

## The model

`chain_spec(N)`: `Source` → `Filter` and `Junction` → a chain of N `Workhorse`
instances, each consuming the previous one's `pose`. The chain forces a real
topological order and keeps the sweep from folding away; the fixed head keeps
the model heterogeneous at every N. Ports cover the three shapes a store has to
carry — a scalar, an `SVector{3}`, and a custom isbits struct (`Pose3`) — and
`Workhorse`'s body is ~20 ops, the size class of the 2026-07 compile-cost
anchors.

Deliberate simplification: declarations take the activation scalar `T`
explicitly instead of being walked out of nominal `Float64` declarations
(§11.2). The walk is not what this bench measures.

## Files

| file | role |
| --- | --- |
| `src/leaves.jl` | flatten/reconstruct over the closed leaf vocabulary — shared (C1 needs it for state, C2 also for cells) |
| `src/model.jl` | synthetic components, stage functions, spec generator, state layout |
| `src/protocol.jl` | candidate interface, entry, chunked walk, `build_sweep` |
