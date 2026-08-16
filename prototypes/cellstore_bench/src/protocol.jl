# The candidate-agnostic half of the executor: entry shape, chunked unrolled
# walk, sweep construction (§9.7). A candidate supplies only its store and its
# cell *addressing*; everything else here is shared, so the measurement isolates
# the one axis under test.
#
# A candidate `C` implements:
#
#   build_store(C, spec, T)   -> (store, layout)   # layout is candidate-private
#   cell_addr(C, layout, ci, port) -> address token
#   gather(store, addr)       -> port value
#   scatter!(store, addr, v)  -> nothing
#
# Four functions. Snapshot capture (gate 3) is shared machinery on top of
# `gather`, so it is not part of the interface.
#
# The address token is where the candidates part company: C1 puts the cell
# index in the *type* domain (a singleton), C2 puts a byte-free integer offset
# in a *field* and carries the port type as the parameter. Everything the
# entry does with it — the gather, the stage call, the scatter — is this file's
# code in both cases.

abstract type Candidate end

# --- entries ------------------------------------------------------------------

"""
One schedule entry: one component's interior `h_xu` stage. Carries what selects
code (component type, port names, address token types) in type parameters and
what is plain data (the state offset, and whatever the address tokens hold) in
fields — §9.7's rule.
"""
struct Entry{Comp,XT,IN,ON,IA<:Tuple,OA<:Tuple}
    comp::Comp
    in_addrs::IA
    out_addrs::OA
    x_off::Int
end

@inline function run!(e::Entry{Comp,XT,IN,ON}, store, xbuf) where {Comp,XT,IN,ON}
    x = reconstruct(XT, xbuf, e.x_off)
    u = gather_inputs(Val(IN), e.in_addrs, store)
    y = h_xu(e.comp, x, u)
    scatter_outputs!(store, e.out_addrs, Val(ON), y)
end

@generated function gather_inputs(::Val{Ns}, addrs::Tuple, store) where {Ns}
    args = [:(gather(store, addrs[$i])) for i in 1:length(Ns)]
    quote
        $(Expr(:meta, :inline))
        NamedTuple{$Ns}(($(args...),))
    end
end

@generated function scatter_outputs!(store, addrs::Tuple, ::Val{Ns}, y::NamedTuple) where {Ns}
    stmts = [:(scatter!(store, addrs[$i], getfield(y, $(QuoteNode(Ns[i]))))) for i in 1:length(Ns)]
    quote
        $(Expr(:meta, :inline))
        $(stmts...)
        nothing
    end
end

function make_entry(::Type{C}, spec::ModelSpec, layout, xl, ci::Int, ::Type{T}) where {C<:Candidate,T}
    cs = spec.comps[ci]
    XT = typeof(x_value(cs, T))
    IN = keys(in_types(cs, T))
    ON = keys(out_types(cs, T))
    @assert Tuple(first(p) for p in cs.inputs) == IN "input order must match the declaration"
    in_addrs = Tuple(cell_addr(C, layout, src[1], src[2]) for (_, src) in cs.inputs)
    out_addrs = Tuple(cell_addr(C, layout, ci, name) for name in ON)
    Entry{typeof(cs.comp),XT,IN,ON,typeof(in_addrs),typeof(out_addrs)}(
        cs.comp, in_addrs, out_addrs, xl.offsets[ci])
end

# --- the walk -----------------------------------------------------------------

"""
A chunk: a statically typed tuple of entries behind a non-inlined function
barrier, bound over the store and state buffer. Chunk size is the executor's
only representation freedom (§9.7) and is held equal across candidates.
"""
struct Chunk{E<:Tuple,S,X}
    entries::E
    store::S
    xbuf::X
end

@noinline (c::Chunk)() = _walk(c.entries, c.store, c.xbuf)

@inline _walk(::Tuple{}, store, xbuf) = nothing
@inline function _walk(t::Tuple, store, xbuf)
    run!(t[1], store, xbuf)
    _walk(Base.tail(t), store, xbuf)
end

"""The zero-arg interior variant of `sweep_hxu` (§9.7) — what gate 1 measures."""
struct Sweep{C<:Tuple}
    chunks::C
end

@inline (s::Sweep)() = _walkchunks(s.chunks)

@inline _walkchunks(::Tuple{}) = nothing
@inline function _walkchunks(t::Tuple)
    t[1]()
    _walkchunks(Base.tail(t))
end

# --- construction -------------------------------------------------------------
# Type-opaque: entries are built into an untyped buffer and splatted once per
# chunk. The compiled tuple's only consumer is the unrolled walk (§9.7).

"""
    build_sweep(C, spec, T; chunk_size = 16)

Returns `(sweep, store, layout, xbuf)`. `sweep()` is the interior `sweep_hxu`.
"""
function build_sweep(::Type{C}, spec::ModelSpec, ::Type{T};
                     chunk_size::Int = 16) where {C<:Candidate,T}
    store, layout = build_store(C, spec, T)
    xl = state_layout(spec, T)
    xbuf = xl.buf
    entries = Any[make_entry(C, spec, layout, xl, ci, T) for ci in eachindex(spec.comps)]
    chunks = Any[]
    for lo in 1:chunk_size:length(entries)
        hi = min(lo + chunk_size - 1, length(entries))
        es = tuple(entries[lo:hi]...)
        push!(chunks, Chunk(es, store, xbuf))
    end
    Sweep(tuple(chunks...)), store, layout, xbuf
end

"""
    entry_type_count(C, spec, T)

Number of distinct `Entry` types in the schedule — the number of `run!` bodies
the executor must compile. The structural form of gate 2: it states the sharing
claim as a count, independently of any timing.
"""
function entry_type_count(::Type{C}, spec::ModelSpec, ::Type{T}) where {C<:Candidate,T}
    _, layout = build_store(C, spec, T)
    xl = state_layout(spec, T)
    length(unique(typeof(make_entry(C, spec, layout, xl, ci, T))
                  for ci in eachindex(spec.comps)))
end

"""
    invoke_sweep(sweep)

The call every measurement goes through. Constant propagation is disabled and
inlining blocked *on purpose*: benchmark interpolation (`\$sweep`) hands the
compiler a constant object whose entries are immutable and whose offsets are
therefore foldable, which would hand C2 a specialization the real loop — where
the simulation is an ordinary runtime value — never gets. This barrier is the
difference between measuring the representation and measuring BenchmarkTools.
"""
Base.@constprop :none @noinline invoke_sweep(s::Sweep) = s()

# --- snapshot capture (gate 3) ------------------------------------------------
# §11.2 publishes by copying the table out; the store's shape prices that copy.
# Flat naming (`path_port`) rather than the real nested-by-path shape — enough
# to price the gathers, which is all this gate asks.

struct Snapshotter{Ns,A<:Tuple,S}
    addrs::A
    store::S
end

@generated function (s::Snapshotter{Ns})() where {Ns}
    args = [:(gather(s.store, s.addrs[$i])) for i in 1:length(Ns)]
    quote
        $(Expr(:meta, :inline))
        NamedTuple{$Ns}(($(args...),))
    end
end

Base.@constprop :none @noinline invoke_snapshot(s::Snapshotter) = s()

function build_snapshotter(::Type{C}, spec::ModelSpec, layout, store,
                           ::Type{T}) where {C<:Candidate,T}
    names = Symbol[]
    addrs = Any[]
    for (ci, cs) in enumerate(spec.comps)
        for port in keys(out_types(cs, T))
            push!(names, Symbol(cs.path, :_, port))
            push!(addrs, cell_addr(C, layout, ci, port))
        end
    end
    a = tuple(addrs...)
    Snapshotter{tuple(names...),typeof(a),typeof(store)}(a, store)
end
