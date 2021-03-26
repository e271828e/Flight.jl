module TestQuaternions

include("../src/quaternions.jl") #from Revise

using Test
using LinearAlgebra

using Main.Quaternions: Quat, UnitQuat, AbstractQuat
using Main.Quaternions: Quat16, Quat32, Quat64
using Main.Quaternions: UnitQuat16, UnitQuat32, UnitQuat64

export test_quaternions

function test_quaternions()
    test_Quat()
    test_UnitQuat()
end

function test_Quat()

    @testset "Quat constructors ($T)" for T in [Float16, Float32, Float64]
        constructors(Quat, T)
    end

    #type parameter promotion goes here

    @testset "Quat basics" begin basics(Quat) end

    @testset "Quat functions" begin
        v = [2.0,1,0,-1]
        q = Quat(v)
        q_copy = copy(q)
        q_2 = q[2]
        q[2] = 0
        @test q[2] == 0.0
        @test q_copy[2] == q_2
    end

    @testset "Quat operators" begin
        #see Python tests
    end

    # inferred goes here

end

function test_UnitQuat()

    @testset "UnitQuat constructors ($T)" for T in [Float16, Float32, Float64]
        constructors(UnitQuat, T)
    end

    @testset "UnitQuat basics" begin basics(UnitQuat) end

    #type parameter promotion

    #Quat promotion

    @testset "UnitQuat norm" begin
        #check constructor input normalization
        v = [1.0, 2, 3, 4]
        q = UnitQuat(v)
        @test q[:] == normalize(v)
        #break norm constraint by bypassing inner constructor input normalization
        q = UnitQuat{Float64}(Quat(v), enforce_norm = false)
        @test q[:] != normalize(v)
        #see that it is restored by normalize!
        normalize!(q)
        @test q[:] == normalize(v)
    end

    @testset "UnitQuat functions" begin
        #include also UnitQuat - Quat operators
    end

    @testset "UnitQuat operators" begin
        #include also UnitQuat - Quat operators
    end
end


function constructors(::Type{Q}, ::Type{T}) where {Q, T}

    #explicit
    @test isa(Q{T}(), Q{T})
    @test isa(Q{T}([1, 2, 3, 4]), Q{T})
    @test isa(Q{T}(1), Q{T})
    @test isa(Q{T}(Q{Float32}()), Q{T}) #type should be Q{T}, not Q32
    @test isa(Q{T}(real = Float64(6), imag = Vector{Float16}([1., 2, 3])), Q{T})
    @test isa(Q{T}(real = Float16(6)), Q{T})
    @test isa(Q{T}(imag = [2, 1, -3]), Q{T})

    #implicit
    @test isa(Q(), Q{Float64}) #should default to Float64
    @test isa(Q(Vector{T}([1, 2, 3, 4])), Q{T})
    @test isa(Q(T(1)), Q{T})
    @test isa(Q(Q{T}()), Q{T}) #should infer parameter type from input Quat
    @test isa(Q(real = Float32(6), imag = Vector{T}([1, 2, 3])), Q{Float64})
    @test isa(Q(real = Float16(6)), Q{Float64})
    @test isa(Q(imag = Vector{T}([1, 2, 3])), Q{Float64})

end

function basics(::Type{T}) where T

    v = [1.0, 2, 3, 4]
    q = T(v)
    @test q.real == q[1]
    @test q.imag == q[2:end]
    @test norm(q) == norm(q[:])
    normalize!(q)
    @test norm(q) ≈ 1

end

#=
function upcasting(::Type{T}) where T
    println("\nTesting upcasting between parameter types for type $T")
    q64 = T(rand(Float64, 4))
    q32 = T(rand(Float32, 4))
    q16 = T(rand(Float16, 4))
    println(typeof(promote(q16, q32)))
    println(typeof(promote(q16, q64)))
    println(typeof(promote(q32, q64)))
end

function inner_operators(::Type{T}) where T
    println("\nTesting inner operators for type $T")
    q64 = T(rand(Float64, 4))
    q32 = T(rand(Float32, 4))
    q16 = T(rand(Float16, 4))

    println(normalize!(q))
    println(q16 * q64)
    println(q64 * q64)
    #Quat, Real
end

function type_parameter_promotion(::Type{Q}, ::Type{T1}, ::Type{T2}) where {Q, T1, T2}
    println("\nTesting promotion for types $T1 and $T2")
    q64 = T1{Float64}(1)
    q32 = T1{Float32}(1)
    q16 = T1{Float16}(1)
    u64 = T2{Float64}(1)
    u32 = T2{Float32}(1)
    u16 = T2{Float16}(1)
    println(promote(q16, u16))
    println(promote(q32, u32))
    println(promote(q64, u64))
    println(promote(q32, u64))
    println(promote(q64, u32))
    println(promote(q32, u16))
    println(promote(q16, u32))
end

=#

end #module