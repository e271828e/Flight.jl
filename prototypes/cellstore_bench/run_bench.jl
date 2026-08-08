# Driver: one child process per (candidate, N), collected into a table.
#
#   julia --project=. run_bench.jl                          # default sweep
#   julia --project=. run_bench.jl 1,10,100                 # explicit Ns
#   julia --project=. run_bench.jl 1,10,50 16 Dual8 C1 240  # Ns chunk scalar cands budget
#
# A candidate that blows its per-point budget stops scaling — that is a result,
# not a failure, and the table records how far it got. The budget matters on the
# `Dual8` activation, where C1's compile cost is the thing being measured and
# can run to minutes per point.

using Printf, Dates

const HERE = @__DIR__
arg(i, default) = length(ARGS) >= i ? ARGS[i] : default

const NS = parse.(Int, split(arg(1, "1,2,5,10,25,50,100,200,400"), ','))
const CHUNK = parse(Int, arg(2, "16"))
const SCALAR = arg(3, "Float64")
const CANDS = split(arg(4, "C1,C2"), ',')
const BUDGET_S = parse(Float64, arg(5, "900"))

function run_point(cand, N::Int, chunk::Int, scalar::String)
    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$HERE $(joinpath(HERE, "bench_point.jl")) $cand $N $chunk $scalar`
    t = time()
    out = read(cmd, String)
    wall = time() - t
    line = only(filter(l -> startswith(l, "RESULT "), split(out, '\n')))
    merge(Main.eval(Meta.parse(line[length("RESULT ")+1:end])), (; process_wall_s = wall))
end

rows = []
for cand in CANDS
    for N in NS
        @printf("  %s/%s N=%-4d ... ", cand, SCALAR, N)
        flush(stdout)
        r = run_point(cand, N, CHUNK, SCALAR)
        push!(rows, r)
        @printf("compile %7.2fs  sweep %9.1fns  alloc %d  bodies %d\n",
                r.sweep_compile_s, r.sweep_ns, r.sweep_alloc, r.entry_types)
        if r.sweep_compile_s > BUDGET_S
            println("  $cand exceeded the per-point budget at N=$N; stopping its scaling")
            break
        end
    end
end

mkpath(joinpath(HERE, "results"))
stamp = Dates.format(now(), "yyyy-mm-dd_HHMM")
csv = joinpath(HERE, "results", "cellstore_$(SCALAR)_$(stamp).csv")
open(csv, "w") do io
    ks = keys(rows[1])
    println(io, join(ks, ","))
    for r in rows
        println(io, join((getfield(r, k) for k in ks), ","))
    end
end

println("\n| cand | scalar |    N | bodies | compile s | snap compile s | sweep ns | ns/entry | alloc | snap ns |")
println("| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
for r in rows
    @printf("| %s | %s | %d | %d | %.2f | %.2f | %.0f | %.1f | %d | %.0f |\n",
            r.candidate, r.scalar, r.N, r.entry_types, r.sweep_compile_s, r.snap_compile_s,
            r.sweep_ns, r.sweep_ns / r.entries, r.sweep_alloc, r.snap_ns)
end
println("\nwrote $csv")
