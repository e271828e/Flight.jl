module TestAttitude

using Test
using LinearAlgebra
using Flight.Attitude
using Flight.Quaternions

export test_attitude

function test_attitude()
    @testset verbose = true "Attitude" begin
        @testset verbose = true "RQuat" begin test_RQuat() end
        @testset "AxisAng" begin test_AxisAng() end
    end
end

function test_RQuat()
    #construction from UnitQuat or anything convertible to UnitQuat, other
    #constructors will be provided by other Rotation subtypes
    #construction from no arguments
    v_ab = rand(4)
    v_bc = rand(4)
    q_ab = RQuat(v_ab)
    q_bc = RQuat(v_bc)
    q_ab_neg = RQuat(-v_ab)
    x_c = rand(3)

    #RQuat accepts any AbstractVector of length 4 and normalizes it, no need to
    #pass a UnitQuat instance
    @test q_ab._quat ≈ RQuat(UnitQuat(v_ab))._quat
    @test_throws DimensionMismatch RQuat([1,2,3])

    #no argument constructor returns unit rotation
    @test RQuat()._quat[:] == [1,0,0,0]

    #operators
    #transformation, inversion, composition, equality, renormalization
    #equality and must account for quaternion double cover
    r_bad = RQuat(UnitQuat([1,2,3,4], enforce_norm = false))
    @test abs(norm(r_bad._quat)-1) > 1e-2
    @test norm(normalize(r_bad)) ≈ 1
    normalize!(r_bad)
    @test norm(r_bad) ≈ 1

    #equality
    @test q_ab != q_bc
    @test q_ab == q_ab_neg
    @test q_ab == q_ab

    #approximate equality
    @test q_ab ≈ q_ab_neg
    @test q_ab ≈ q_ab

    #overloads for inversion, transformation and composition
    @test q_ab' == invert(q_ab)
    @test q_ab ∘ q_bc == compose(q_ab, q_bc)
    @test q_bc * x_c == transform(q_bc, x_c)

    #transformation & composition
    @test q_ab * (q_bc * x_c) ≈ (q_ab ∘ q_bc) * x_c

    #inversion
    @test x_c ≈ q_bc' * (q_bc * x_c)

end

function test_AxisAng()

end

end #module