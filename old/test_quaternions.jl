module TestQuaternions

include("../src/quaternions.jl")

using Test
using LinearAlgebra
using .Quaternions

export test_quaternions

function test_quaternions()
    @testset "Quat" begin test_Quat() end
    @testset "UnitQuat" begin test_UnitQuat() end
end

function test_Quat()

    @testset "Constructors ($T)" for T in [Float16, Float32, Float64]
        constructors(Quat, T)
    end

    @testset "Promotion and conversion" begin
        q64 = Quat(rand(Float64, 4))
        q32 = Quat(rand(Float32, 4))
        q16 = Quat(rand(Float16, 4))
        q16_q32 = promote(q16, q32)
        q16_q64 = promote(q16, q64)
        q32_q64 = promote(q32, q64)
        q32_f64 = promote(q32, 3.)

        @test isa(q16_q32, Tuple{Quat32, Quat32}) && q16_q32[1] ≈ q16
        @test isa(q16_q64, Tuple{Quat64, Quat64}) && q16_q64[1] ≈ q16
        @test isa(q32_q64, Tuple{Quat64, Quat64}) && q32_q64[1] ≈ q32
        @test isa(q32_f64, Tuple{Quat32, Quat32}) && q32_f64[2] ≈ Quat(3.)
    end

    @testset "Basics" begin
        v = [1.0, 2, 3, 4];
        q = Quat(v)
        @test all(x[1]==x[2] for x in zip(q, v))

        @test q.real == v[1]
        @test q.imag == v[2:end]

        q_copy = copy(q)
        q_2 = q[2]
        q[2] = 0
        @test q[2] == 0.0
        @test q_copy[2] == q_2

        @test q'.imag ≈ -q.imag

        @test norm(q) == norm(q[:])
        normalize!(q)
        @test norm(q) ≈ 1

    end

    @testset "Operators" begin
        v1 = rand(4); v2 = rand(4); v3 = rand(4)
        q1 = Quat(v1); q2 = Quat(v2); q3 = Quat(v3)

        #unary minus, addition and substraction
        @test (-q1)[:] == -v1
        @test (q1 + q2)[:] == v1 + v2
        @test (q1 - q2)[:] == v1 - v2

        #multiplication, inverse and division by quaternion
        @test q1 * Quat(1.0) ≈ q1
        @test (q1 * q2) * q3 ≈ (q1 * (q2 * q3))
        @test (q1 + q2) * q3 ≈ (q1 * q3 + q2 * q3)
        @test q1 * inv(q1) ≈ Quat(1.0)
        @test q1 / q2 ≈ q1 * inv(q2)
        @test q1 \ q2 ≈ inv(q1) * q2
        @test q1 \ q2 != q2 / q1 #inv(q1)*q2 != q2*inv(q1)

        #and by real
        @test 5q1 ≈ Quat64(5*q1[:])
        @test q1/5 ≈ Quat64(q1[:]/5)
        @test 5/q1 ≈ Quat(5.) * inv(q1)
        @test 5\q1 ≈ inv(Quat(5.)) * q1
        #see Python tests
    end

    @testset "Type stability" begin
        #exercise different constructors, type parameters and operators
        @inferred norm(Quat32(real=4) * Quat16([3,2,1,4])' / Quat(rand(4))+1)
    end

end

function test_UnitQuat()

    @testset "Constructors ($T)" for T in [Float16, Float32, Float64]
        constructors(UnitQuat, T)
    end

    @testset "Promotion and conversion" begin
        u64 = UnitQuat(rand(Float64, 4))
        u32 = UnitQuat(rand(Float32, 4))
        q64 = Quat(rand(Float64, 4))
        q32 = Quat(rand(Float32, 4))
        u64_u32 = promote(u64, u32)
        u64_q32 = promote(u64, q32)
        u64_q64 = promote(u64, q64)
        u32_q32 = promote(u32, q32)
        u32_q64 = promote(u32, q64)

        @test isa(u64_u32, Tuple{UnitQuat64, UnitQuat64}) && u64_u32[2] ≈ u32
        @test isa(u64_q32, Tuple{Quat64, Quat64}) && u64_q32[1] ≈ u64 && u64_q32[2] ≈ q32
        @test isa(u64_q64, Tuple{Quat64, Quat64}) && u64_q64[1] ≈ u64 && u64_q64[2] ≈ q64
        @test isa(u32_q32, Tuple{Quat32, Quat32}) && u32_q32[1] ≈ u32 && u32_q32[2] ≈ q32
        @test isa(u32_q64, Tuple{Quat64, Quat64}) && u32_q64[1] ≈ u32 && u32_q64[2] ≈ q64
    end

    @testset "Basics" begin
        #default constructor should return identity quaternion
        @test UnitQuat() == Quat(1.0)

        #check constructor input normalization
        v = [1.0, 2, 3, 4]
        v_norm = normalize(v)
        u = UnitQuat(v)

        @test u[:] == v_norm
        #break norm constraint by bypassing inner constructor input normalization
        u = UnitQuat{Float64}(Quat(v), enforce_norm = false)

        @test u[:] != v_norm
        #see that it is restored by normalize!

        normalize!(u)
        @test u[:] == v_norm

        @test all(x[1]==x[2] for x in zip(u, v_norm))
        @test u.real == v_norm[1]
        @test u.imag == v_norm[2:end]
        @test u'.imag ≈ -u.imag
    end

    @testset "Operators" begin
        u1 = UnitQuat(rand(4)); u2 = UnitQuat(rand(4))
        q = Quat(rand(4))

        #unary minus
        @test (-u1)[:] == -u1

        #multiplication, inverse and division by UnitQuat
        @test u1 * UnitQuat() ≈ u1
        @test u1' == inv(u1)
        @test u1 * u1' ≈ UnitQuat()
        @test u1 / u2 ≈ u1 * u2'
        @test u1 \ u2 ≈ u1' * u2
        @test u1 \ u2 != u2 / u1 #inv(u1)*u2 != u2*inv(u1)

        #compatibility with Quat
        @test u1 * Quat(1.0) ≈ u1
        @test q * u1' ≈ q / Quat(u1)
        @test u1 / q ≈ Quat(u1) / q
        @test u1 \ q ≈ Quat(u1) \ q
    end

    @testset "Type stability" begin
        #exercise different constructors, type parameters and operators
        @inferred norm(UnitQuat32(imag=[1,2,3])' * UnitQuat16())
    end

end

#required constructors for AbstractQuat subtypes
function constructors(::Type{Q}, ::Type{T}) where {Q, T}

    #explicit
    t1 = (
        Q{T}(),
        Q{T}([1, 2, 3, 4]),
        Q{T}(1),
        Q{T}(Q{Float32}()),
        Q{T}(real = Float16(6)),
        Q{T}(imag = [2, 1, -3]),
        Q{T}(real = Float64(6), imag = Vector{Float16}([1., 2, 3]))
        )
    @test all(isa.(t1, Q{T}))

    #implicit, inferred type parameter
    t2 = (
        Q(Vector{T}([1, 2, 3, 4])),
        Q(T(1)),
        Q(Q{T}()), #should infer parameter type from input Quat
    )
    @test all(isa.(t1, Q{T}))

    #implicit, default type parameter
    t3 = (
        Q(),
        Q(real = Float32(6), imag = Vector{T}([1, 2, 3])),
        Q(real = Float16(6)),
        Q(imag = Vector{T}([1, 2, 3])),
    )
    @test all(isa.(t3, Q{Float64}))

end

end #module