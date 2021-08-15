using ComponentArrays

const WrenchAxes = getaxes(ComponentVector(F = zeros(3), M = zeros(3)))

struct Wrench{D <: AbstractVector{Float64}} <: AbstractVector{Float64}
    data::ComponentVector{Float64, D, typeof(WrenchAxes)}
end

#avoid ComponentVector(F = F, M = M), which forwards to ComponentArray{T}
#(::NamedTuple), which needs to create the axes from scratch and is very slow.
#much faster to provide the axes directly

function Wrench(input::AbstractVector{<:Real})
    if length(input) != 6
        throw(ArgumentError("Got input length $(length(input)), expected 6"))
    end
    data = ComponentVector{Float64}(undef, WrenchAxes)
    data .= input
    Wrench(data)
end

function Wrench(; F = SVector(0,0,0), M = SVector(0,0,0))
    data = ComponentVector{Float64}(undef, WrenchAxes)
    data.F = F; data.M = M
    Wrench(data)
end
ComponentArrays.ComponentVector(wr::Wrench) = copy(getfield(wr,:data))

Base.size(::Wrench) = (6,)
Base.length(::Wrench) = 6

#this is type unstable! it returns ::Any (?!!!)
# Base.getproperty(wr::Wrench, i::Symbol) = getproperty(wr, Val(i))
# Base.getproperty(wr::Wrench, ::Val{:data}) = getfield(wr, :data)
# Base.getproperty(wr::Wrench, i::Val{S} where {S}) = getproperty(wr.data, i)

Base.getproperty(wr::Wrench, i::Symbol) = getproperty(getfield(wr,:data), i)
Base.setproperty!(wr::Wrench, i::Symbol, v) = setproperty!(wr, Val(i), v)
Base.setproperty!(wr::Wrench, ::Val{:F}, v) = (wr.F .= v) #gets the F block, then broadcasts
Base.setproperty!(wr::Wrench, ::Val{:M}, v) = (wr.M .= v) #gets the M block, then broadcasts

Base.getindex(wr::Wrench, i) = getindex(getfield(wr,:data), i)
Base.setindex!(wr::Wrench, v, i) = setindex!(getfield(wr,:data), v, i)

Base.eltype(::Wrench) = Float64 #helps with allocation efficiency
Base.similar(::Wrench) = Wrench(ComponentVector{Float64}(undef, WrenchAxes))
Base.similar(::Type{<:Wrench}) = Wrench(ComponentVector{Float64}(undef, WrenchAxes))

Base.show(io::IO, wr::Wrench) = print(io, "Wrench(F = $(wr.F), M = $(wr.M))")
Base.show(io::IO, ::MIME"text/plain", wr::Wrench) = print(io, "Wrench(F = $(wr.F), M = $(wr.M))")

#since Wrench <: AbstractVector, broadcasting works out of the box (and with
#it the non-broadcast operators +, -, etc), but it falls back to the default
#broadcast implementation in Base, which returns a generic Vector. we want
#broadcasted operations to return a Wrench. so we do this:

struct WrenchStyle{D} <: Broadcast.AbstractArrayStyle{1} end

WrenchStyle{D}(::Val{1}) where {D} = WrenchStyle{D}()
Base.BroadcastStyle(::Type{Wrench{D}}) where {D} = WrenchStyle{D}()
function Base.similar(::Broadcast.Broadcasted{WrenchStyle{D}}, ::Type{ElType}) where {D,ElType}
    similar(Wrench{D})
end
function Base.BroadcastStyle(::WrenchStyle{D1}, ::WrenchStyle{D2}) where {D1,D2}
    WrenchStyle{promote_type(D1, D2)}()
end