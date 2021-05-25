using BlockArrays

struct LBVector{K, L} <: AbstractVector{Float64}
    #K: number of blocks
    #L: labels
    data::AbstractArray{Float64}
    labels::NamedTuple{L, NTuple{K, Int}}
    # LBVector{K, L}(blocks, labels) where {K, L} = new(mortar())
    #here, we enforce the number of blocks and length of L to be equal
end

Base.size(x::LBVector) = size(x.data)
Base.length(x::LBVector) = length(x.data)
Base.firstindex(x::LBVector) = 1
Base.lastindex(x::LBVector) = lastindex(x.data)
Base.getindex(x::LBVector, i) = getindex(x.data, i)
Base.setindex!(x::LBVector, v, i) = setindex(x.data, v, i)

x = LBVector{2, (:a, :b)}(ones(2), (a = 1, b = 2))
# println(x)
