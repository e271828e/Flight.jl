# One (candidate, N) measurement, in a fresh process. Compile cost is only
# meaningful cold, so the driver spawns one of these per point rather than
# looping in-session.
#
#   julia --project=. bench_point.jl C1 100 16 Float64
#   julia --project=. bench_point.jl C1 25 16 Dual8

using StaticArrays, BenchmarkTools, ForwardDiff

include("src/leaves.jl")
include("src/model.jl")
include("src/protocol.jl")
include("src/c1_hetero.jl")
include("src/c2_flat.jl")

const CANDIDATES = Dict("C1" => C1, "C2" => C2)

# The two activations of interest (§7.2): the nominal one, and the 8-partial
# `Dual` that linearization and trim run through — the instruction multiplier
# §12.7's anchors put at ~20×.
const SCALARS = Dict("Float64" => Float64,
                     "Dual8" => ForwardDiff.Dual{Nothing,Float64,8})

const name = ARGS[1]
const cand = CANDIDATES[name]
const N = parse(Int, ARGS[2])
const chunk = parse(Int, length(ARGS) >= 3 ? ARGS[3] : "16")
const scalar = length(ARGS) >= 4 ? ARGS[4] : "Float64"
const T = SCALARS[scalar]

spec = chain_spec(N)

# --- gate 2: cold compile of the sweep, then of the snapshotter --------------
# Measured separately: the snapshot's flat NamedTuple is a large type in its own
# right, and letting it into the sweep number would price the wrong thing.

Base.cumulative_compile_timing(true)

c0 = Base.cumulative_compile_time_ns()[1]
t0 = time_ns()
sweep, store, layout, xbuf = build_sweep(cand, spec, T; chunk_size = chunk)
invoke_sweep(sweep)
t1 = time_ns()
c1 = Base.cumulative_compile_time_ns()[1]

s0 = Base.cumulative_compile_time_ns()[1]
ts0 = time_ns()
snap = build_snapshotter(cand, spec, layout, store, T)
invoke_snapshot(snap)
ts1 = time_ns()
s1 = Base.cumulative_compile_time_ns()[1]

Base.cumulative_compile_timing(false)

# --- gates 1 and 3 ------------------------------------------------------------

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 0.5

sweep_alloc = @ballocated invoke_sweep($sweep)
sweep_s = @belapsed invoke_sweep($sweep)
snap_alloc = @ballocated invoke_snapshot($snap)
snap_s = @belapsed invoke_snapshot($snap)

result = (candidate = name, scalar = scalar, N = N, chunk = chunk,
          entries = length(spec), cells = length(cell_inventory(spec, T)),
          entry_types = entry_type_count(cand, spec, T),
          chunk_types = length(unique(map(typeof, sweep.chunks))),
          sweep_compile_s = (c1 - c0) / 1e9, sweep_wall_s = (t1 - t0) / 1e9,
          snap_compile_s = (s1 - s0) / 1e9, snap_wall_s = (ts1 - ts0) / 1e9,
          sweep_alloc = sweep_alloc, sweep_ns = sweep_s * 1e9,
          snap_alloc = snap_alloc, snap_ns = snap_s * 1e9)

println("RESULT ", result)
