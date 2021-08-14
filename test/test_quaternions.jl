module TestQuaternions

using Test
using LinearAlgebra
using StaticArrays
using Flight.Quaternions

export test_quaternions

function test_quaternions()
    @testset verbose = true "Quaternions" begin
        @testset verbose = true "Quat" begin test_Quat() end
        @testset verbose = true "UnitQuat" begin test_UnitQuat() end
    end
end

function test_Quat()

    @testset "Basics" begin

        #use inner constructor only, do not rely on outer constructors
        data = SVector{4,Float64}([-3, 2, 4, 1])
        q = Quat(data)

        #getindex
        @test q[:] == data[:]

        #getproperty
        @test q.real == data[1]
        @test q.imag == data[2:end]

        #copy
        q_copy = copy(q)
        @test isa(q_copy, Quat)

        #normalization
        q_copy = copy(q)
        @test norm(normalize(q)) ≈ 1 #returns a normalized copy of q

    end

    @testset "Constructors" begin

        v = [2, 3, 1, -4]
        q = Quat(v)

        @test q[:] == v[:] #Vector
        @test Quat(q)[:] == q[:] #AbstractQuat
        @test Quat(2.0)[:] == [2.0, 0, 0, 0] #scalar
        @test Quat()[:] == zeros(4) #zero arguments
        @test Quat(real = 2, imag = [4,3,-2])[:] == [2, 4, 3, -2] #keyword
        @test Quat(real = 2)[:] == [2, 0, 0, 0] #keyword
        @test Quat(imag = [3,4,5])[:] == [0, 3, 4, 5] #keyword
        @test_throws DimensionMismatch Quat([1, 3, 2, 5, 0])
        @test_throws DimensionMismatch Quat([1, 5, 0])

    end

    @testset "Operators" begin

        v1 = rand(4); v2 = rand(4); v3 = rand(4)
        q1 = Quat(v1); q2 = Quat(v2); q3 = Quat(v3)

        #equality (inherited from AbstractVector)
        # @test q1 != v1 #not equal to an AbstractVector with the same components
        @test q1 != Quat(-v1)
        @test q1 == q1

        #approximate equality
        @test q1 ≈ Quat(v1 + [1e-9, 0, 0, 0])
        @test !(q1 ≈ q2)

        #promotion and conversion from real scalar and AbstractVector
        @test promote(q1, 3.213)[2] == Quat(3.213)
        @test convert(Quat, [1, 2, 3, 4]) == Quat([1,2,3,4])

        #conjugation
        @test q1'.imag ≈ -q1.imag

        #unary minus
        @test isa(-q1, Quat)
        @test (-q1)[:] == -(q1[:])

        #addition and substraction
        @test isa(q1 + q2, Quat)
        @test (q1 + q2)[:] == v1 + v2
        @test isa(q1 - q2, Quat)
        @test (q1 - q2)[:] == v1 - v2

        #multiplication, inverse and division by quaternion
        @test isa(q1 * q2, Quat)
        @test q1 * Quat(1.0) ≈ q1
        @test (q1 * q2) * q3 ≈ (q1 * (q2 * q3))
        @test (q1 + q2) * q3 ≈ (q1 * q3 + q2 * q3)
        @test q1 * inv(q1) ≈ Quat(1.0)
        @test q1 / q2 ≈ q1 * inv(q2)
        @test q1 \ q2 ≈ inv(q1) * q2
        @test q1 \ q2 != q2 / q1 #inv(q1)*q2 != q2*inv(q1)

        #and by real
        @test 5q1 ≈ Quat(5*q1[:])
        @test q1/5 ≈ Quat(q1[:]/5)
        @test 5/q1 ≈ Quat(5.) * inv(q1)
        @test 5\q1 ≈ inv(Quat(5.)) * q1
    end

end

function test_UnitQuat()

    @testset "Basics" begin

        q = Quat([-3, 2, 4, 1])
        u_badnorm = UnitQuat(q, normalization = false)
        u = UnitQuat(q)

        #built-in normalization
        @test abs(norm(u_badnorm) - 1) > 0.5 #check normalization has been bypassed
        @test norm(u) ≈ 1

        #getindex
        @test u_badnorm[:] == q[:]

        #getproperty
        @test u_badnorm.real == q[1]
        @test u_badnorm.imag == q[2:end]

        #setindex & setproperty disallowed
        @test_throws ErrorException u[4] = -5
        @test_throws ErrorException u.real = 0
        @test_throws ErrorException u.imag = [1, 2, 3]

        #copy & normalization
        u_badnorm_copy = copy(u_badnorm)
        @test isa(u_badnorm_copy, UnitQuat)
        @test abs(norm(u_badnorm_copy) - 1) > 0.5 #check normalization is bypassed on copy
        @test norm(normalize(u_badnorm)) ≈ 1 #returns a normalized copy of u
        # normalize!(u_badnorm)
        # @test norm(u_badnorm) ≈ 1 #returns a normalized copy of u

    end

    @testset "Constructors" begin

        v = [2.0, 3, 1, -4]
        u = UnitQuat(Quat(v)) #inner constructor
        q = Quat(v)

        #outer constructors
        @test u[:] == normalize(v)[:] #Vector
        @test UnitQuat(q)[:] == normalize(q)[:] #AbstractQuat
        @test UnitQuat(2.0)[:] == [1, 0, 0, 0] #scalar
        @test UnitQuat()[:] == [1, 0, 0, 0] #zero arguments
        @test UnitQuat(real = v[1], imag = v[2:4])[:] == normalize(v)[:]

    end

    @testset "Operators" begin

        v1 = normalize(rand(4)); v2 = normalize(rand(4))
        u1 = UnitQuat(v1); u2 = UnitQuat(v2)

        #equality (inherited from AbstractVector)
        # @test u1 != v1 #not equal to an AbstractVector with the same components
        @test u1 != UnitQuat(-v1)
        @test u1 == UnitQuat(v1)
        @test u1 == Quat(u1) #equal to a Quat with the same components

        #approximate equality
        @test u1 ≈ UnitQuat(Quat(v1 + [1e-9, 0, 0, 0]), normalization = false)
        @test !(u1 ≈ u2)

        #promotion and conversion to Quat
        @test promote(u1, Quat())[1] == Quat(u1)

        #conjugation
        @test u1'.imag ≈ -u1.imag

        #unary minus
        @test isa(-u1, UnitQuat)
        @test (-u1)[:] == -(u1[:])

        #multiplication, inverse and division by UnitQuat
        @test u1 * UnitQuat() ≈ u1
        @test u1' == inv(u1)
        @test u1 * u1' ≈ UnitQuat()
        @test u1 / u2 ≈ u1 * u2'
        @test u1 \ u2 ≈ u1' * u2
        @test u1 \ u2 != u2 / u1 #inv(u1)*u2 != u2*inv(u1)

        #compatibility with Quat
        q = Quat(rand(4))
        @test u1 * Quat(1.0) ≈ u1
        @test q * u1' ≈ q / Quat(u1)
        @test u1 / q ≈ Quat(u1) / q
        @test u1 \ q ≈ Quat(u1) \ q

    end

end

end #module