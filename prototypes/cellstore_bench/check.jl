# Correctness and gate 1, for both candidates. Run before trusting any number
# out of `run_bench.jl`: a fast wrong sweep would win every gate.
#
#   julia --project=. check.jl

using Test, StaticArrays, BenchmarkTools

include("src/leaves.jl")
include("src/model.jl")
include("src/protocol.jl")
include("src/c1_hetero.jl")
include("src/c2_flat.jl")

const T = Float64

"""
Independent evaluation of the model: plain values threaded by name, no store,
no schedule compilation. Slow and type-unstable on purpose — it shares nothing
with the machinery under test except the stage functions themselves, so
agreement with it is evidence the scatter/gather protocol is wired correctly
and not merely self-consistent.
"""
function reference_snapshot(spec::ModelSpec, ::Type{T}) where {T}
    produced = Dict{Tuple{Int,Symbol},Any}()
    names, vals = Symbol[], Any[]
    for (ci, cs) in enumerate(spec.comps)
        u = NamedTuple{keys(in_types(cs, T))}(tuple((produced[src] for (_, src) in cs.inputs)...))
        y = h_xu(cs.comp, x_value(cs, T), u)
        for port in keys(out_types(cs, T))
            produced[(ci, port)] = getfield(y, port)
            push!(names, Symbol(cs.path, :_, port))
            push!(vals, getfield(y, port))
        end
    end
    NamedTuple{tuple(names...)}(tuple(vals...))
end

"""
Leaf-wise approximate equality. Exact equality is the wrong test here: the
compiled sweeps contract `q*k - cmd*e` into an FMA and the reference (running
type-unstable, in a different body) does not, so results differ in the last
ulp and accumulate along the chain. Same class as the muladd/FMA test failure
already settled in this repo — never `==` on separately compiled float products.
"""
function ≊(a::NamedTuple, b::NamedTuple)
    keys(a) == keys(b) &&
        isapprox(collect(Float64, _leaf_values(a)), collect(Float64, _leaf_values(b)); rtol = 1e-12)
end

function build(::Type{C}, N::Int; chunk_size = 16) where {C<:Candidate}
    spec = chain_spec(N)
    sweep, store, layout, _ = build_sweep(C, spec, T; chunk_size)
    snap = build_snapshotter(C, spec, layout, store, T)
    (; spec, sweep, snap)
end

@testset "candidates compute the model" begin
    for N in (1, 3, 17), C in (C1, C2)
        b = build(C, N)
        invoke_sweep(b.sweep)
        @test invoke_snapshot(b.snap) ≊ reference_snapshot(b.spec, T)
    end
end

@testset "candidates agree with each other" begin
    for N in (1, 3, 17)
        b1, b2 = build(C1, N), build(C2, N)
        invoke_sweep(b1.sweep)
        invoke_sweep(b2.sweep)
        @test invoke_snapshot(b1.snap) ≊ invoke_snapshot(b2.snap)
    end
end

@testset "chunk size does not change results" begin
    for cs in (1, 4, 64), C in (C1, C2)
        spec = chain_spec(17)
        sweep, store, layout, _ = build_sweep(C, spec, T; chunk_size = cs)
        snap = build_snapshotter(C, spec, layout, store, T)
        invoke_sweep(sweep)
        @test invoke_snapshot(snap) ≊ reference_snapshot(spec, T)
    end
end

@testset "gate 1: the interior sweep does not allocate" begin
    for N in (1, 17), C in (C1, C2)
        b = build(C, N)
        invoke_sweep(b.sweep)
        invoke_snapshot(b.snap)
        @test @ballocated(invoke_sweep($(b.sweep))) == 0
        @test @ballocated(invoke_snapshot($(b.snap))) == 0
    end
end

@testset "gate 2, structurally: distinct entry bodies per candidate" begin
    # C1 addresses in the type domain, so every instance is its own type;
    # C2 addresses in fields, so instances of one component type share a body.
    for N in (10, 40)
        @test entry_type_count(C1, chain_spec(N), T) == N + 3
        @test entry_type_count(C2, chain_spec(N), T) == 4
    end
end
