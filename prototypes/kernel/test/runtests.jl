# Acceptance tests for the kernel prototype: the continuous-tier walking
# skeleton (increment 2), the discrete tier at one rate (increment 3), the
# hierarchy (increment 4), the multi-rate grid (increment 5), events at
# boundaries (increment 6), localization (increment 7), the stepper seam
# with its second backend (increment 8) and the data plane's core exchange
# (increment 9).
#
#   julia --project=. test/runtests.jl

using Test, StaticArrays, LinearAlgebra, BenchmarkTools, ForwardDiff

include("../src/leaves.jl")
include("../src/declare.jl")
include("../src/assembly.jl")
include("../src/store.jl")
include("../src/executor.jl")
include("../src/build.jl")
include("../src/dataplane.jl")
include("../src/stepper.jl")
include("../src/sim.jl")
include("../src/localize.jl")
include("../src/library.jl")

include("utils.jl")
include("test_continuous.jl")
include("test_structure.jl")
include("test_discrete.jl")
include("test_hierarchy.jl")
include("test_multirate.jl")
include("test_events.jl")
include("test_localization.jl")
include("test_stepper.jl")
include("test_dataplane.jl")
