using ComponentArrays, StaticArrays, LinearAlgebra
using UnPack

using Flight

#considerar hacer x0, y0, etc ComponentVectors basados en SVectors. de esa
#manera, no tengo que guardar sus longitudes como type parameters adicionales, y
#puedo echar mano de ellas con length cuando lo necesite en el codigo, porque
#van a estar estaticamente determinadas.

#lo mismo con MVectors. puedo hacer que x, ẋ y u sean
#ComponentVector(MVector)s. eso hace que sus tamanos esten determinados
#estaticamente.

#pero la putada es que tanto SVector como MVector hacen loop unrolling y todas
#estas historias, y para matrices de 16x16 y similares que voy a tener yo, la
#cosa empieza a salirse de madre. lo que necesito es un SizedArray, cuya
#longitud es conocida en tiempo de compilacion, pero que no anade optimizaciones
#que para matrices grandes pueden ser contraproducentes

#FieldArray no es una alternativa, porque luego pierdo la capacidad de h

function StaticArrays.SVector{L}(x::ComponentVector) where {L}
    ComponentVector(SVector{L}(getdata(x)), getaxes(x))
end

const tV = AbstractVector{<:Float64}
const tM = AbstractMatrix{<:Float64}

Base.@kwdef struct LinearSystemDescriptor3{
                    LX, LU, LY,
                    tX <: tV, tU <: tV, tY <: tV,
                    tA <: tM, tB <: tM, tC <: tM, tD <: tM} <: SystemDescriptor

    ẋ0::tX; x0::tX; u0::tU; y0::tY;
    A::tA; B::tB; C::tC; D::tD;
    x_cache::tX; y_cache::tY; y_out::tY;
    Δx_cache::tX; Δu_cache::tU

    function LinearSystemDescriptor3(ẋ0, x0, u0, y0, A, B, C, D)

        # ẋ0, x0, u0, y0= map(x->SizedVector{length(x)}(x), (ẋ0, x0, u0, y0))
        # A, B, C, D= map(M->SizedMatrix{size(M)...}(M), (A, B, C, D))

        lengths = map(length, (x0, u0, y0))
        types = map(typeof, (x0, u0, y0, A, B, C, D))

        println("OK")
        new{(lengths..., types...)...}(
            ẋ0, x0, u0, y0, A, B, C, D, copy(x0), copy(y0), copy(y0), copy(x0), copy(u0))

    end

end

Systems.init(::SystemX, desc::LinearSystemDescriptor3) = copy(desc.x0)
Systems.init(::SystemU, desc::LinearSystemDescriptor3) = copy(desc.u0)
Systems.init(::SystemY, desc::LinearSystemDescriptor3) = SVector{length(desc.y0)}(desc.y0)


function test()

    x0 = ComponentVector(a = 1.0, b = 0.5, c = 0.3)
    ẋ0 = similar(x0)
    u0 = ComponentVector(e = 0.1, a = 0.2)
    y0 = ComponentVector(p = 0.3, q = 0.8, r = 2.0, h = 3.0)

    ẋ_buffer = similar(x0) #preallocated buffers to avoid allocations in f_ode!
    y_buffer1 = copy(y0)#preallocated buffers to avoid allocations in f_ode!
    y_buffer2 = copy(y0)#preallocated buffers to avoid allocations in f_ode!

    A = (x0).^0.9 * x0'
    B = (x0).^0.9 * u0'
    C = (y0).^0.9 * x0'
    D = (y0).^0.9 * u0'

    ẋ = rand(3)
    x = rand(3)
    u = rand(2)
    y = rand(4)

    @show ẋ = x0 + A * (x - x0) + B * (u - u0)
    @show y = y0 + C * (x - x0) + D * (u - u0)

    lsd = LinearSystemDescriptor3(ẋ0, x0, u0, y0, A, B, C, D)
    sys = System(lsd)

    return sys

end

#########################
#
function Systems.f_ode!(sys::System{<:LinearSystemDescriptor3{LX, LU, LY}}) where {LX, LU, LY}

    @unpack ẋ, x, u, y, params = sys
    @unpack ẋ0, x0, u0, y0, A, B, C, D, x_cache, y_cache, y_out, Δx_cache, Δu_cache = params

    #to avoid allocations, we cannot simply do:
    #ẋ = ẋ0 + A * (x - x0) + B * (u - u0)
    #y = y0 + C * (x - x0) + D * (u - u0)

    @. Δx_cache = x - x0
    @. Δu_cache = u - u0

    ẋ .= ẋ0
    mul!(x_cache, A, Δx_cache)
    ẋ .+= x_cache
    mul!(x_cache, B, Δu_cache)
    ẋ .+= x_cache

    y_out .= y0
    mul!(y_cache, C, Δx_cache)
    y_out .+= y_cache
    mul!(y_cache, D, Δu_cache)
    y_out .+= y_cache

    sys.y = SVector{LY}(y_out)

    return nothing

end

# function StaticArrays.SizedVector{L}(x::ComponentVector) where {L}
#     ComponentVector(SizedVector{L}(getdata(x)), getaxes(x))
# end

# function StaticArrays.SizedMatrix{M,N}(x::ComponentMatrix) where {M,N}
#     ComponentMatrix(SizedMatrix{M,N}(getdata(x)), getaxes(x))
# end