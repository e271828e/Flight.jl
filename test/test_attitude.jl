module TestAttitude

using Test
using LinearAlgebra
using Flight.Attitude
using Flight.Quaternions

export test_attitude

logrange(x1, x2, n) = (10^y for y in range(log10(x1), log10(x2), length=n))

function test_attitude()
    @testset verbose = true "Attitude" begin
        @testset verbose = true "RQuat" begin test_RQuat() end
        @testset verbose = true "AxisAng" begin test_RAxAng() end
    end
end

function test_RQuat()
    #construction from UnitQuat or anything convertible to UnitQuat, other
    #constructors will be provided by other Rotation subtypes
    v_ab = [0.7, -1, 0.1, 1]
    v_bc = [1.7, 2, -1, 1]
    q_ab = RQuat(v_ab)
    q_bc = RQuat(v_bc)
    q_ab_neg = RQuat(-v_ab)
    x_c = [2, 3, -1.1]

    @testset "Constructors" begin
        #RQuat accepts any AbstractVector of length 4 and normalizes it, no need to
        #pass a UnitQuat instance
        @test q_ab._quat ≈ RQuat(UnitQuat(v_ab))._quat
        @test_throws DimensionMismatch RQuat([1,2,3])
        @test RQuat()._quat[:] == [1,0,0,0]
    end

    @testset "Normalization" begin
        r_bad = RQuat(UnitQuat([1,2,3,4], normalization = false))
        @test abs(norm(r_bad._quat)-1) > 1e-2
        @test norm(normalize(r_bad)) ≈ 1
        normalize!(r_bad)
        @test norm(r_bad) ≈ 1
    end

    @testset "Operators" begin
        #operators: transformation, inversion, composition, equality
        #equality and isapprox must account for quaternion double cover

        #equality
        @test q_ab != q_bc
        @test q_ab == q_ab_neg
        @test q_ab == q_ab

        #isapprox
        @test q_ab ≈ q_ab_neg
        @test q_ab ≈ q_ab
        @test !(q_ab ≈ q_bc)

        # transformation & composition
        @test q_ab * (q_bc * x_c) ≈ (q_ab ∘ q_bc) * x_c

        #inversion
        @test x_c ≈ q_bc' * (q_bc * x_c)

        #time derivative
        ω_ab_b = [10, -4, 2]
        @test dt(q_ab, ω_ab_b) == 0.5 * (q_ab._quat * Quat(imag = ω_ab_b))

    end

end

function test_RAxAng()

    u_ab_bad = [5, -2, 1]
    u_ab = normalize(u_ab_bad)
    μ_ab = 1.214
    u_μ_ab = (u_ab, μ_ab)

    r_ab = RAxAng(u_ab, μ_ab)
    r_bc = RAxAng([-2,-3, -1], 0.72)
    x_c = [-1, 2, 3]

    @testset "Constructors" begin

        #inner constructor
        @test RAxAng(u_ab_bad, μ_ab, normalization = false).axis == u_ab_bad
        @test RAxAng(u_ab_bad, μ_ab).axis ≈ normalize(u_ab)

        #tuple constructor
        @test r_ab.axis == RAxAng(u_μ_ab).axis && r_ab.angle == RAxAng(u_μ_ab).angle

        #zero argument constructor
        @test RAxAng().axis == [1, 0, 0] && RAxAng().angle == 0

    end

    @testset "Conversions" begin

        #sanity checks
        @test RAxAng(RQuat(r_ab)) ≈ r_ab
        @test RQuat(RAxAng([1, 2, -1], 0)) ≈ RQuat([1,0,0,0])
        @test RQuat(RAxAng([3, 1, 1], 2π)) ≈ RQuat([1,0,0,0])

        #[-2π, 2π] including small angles
        u = normalize([7, -5, 2])
        μ_gen = logrange(1e-10, 2π, 10)
        μ_array = [-reverse(collect(μ_gen)); collect(μ_gen)]
        r_array = collect(RAxAng(u, μ) for μ in μ_array)
        r_array_test = collect(RAxAng(RQuat(r)) for r in r_array)
        @test all(r_array .≈ r_array_test)

    end


    #these actually build upon quaternion conversion, but that's an
    #implementation detail we should not care about, so we can test these before
    @testset "Homogeneous Operators" begin

        #equality must account for different axis and angle sign combinations, and
        #2kπ increments
        @test RAxAng(u_ab, μ_ab) == RAxAng(u_ab, μ_ab)
        @test RAxAng(u_ab, μ_ab) ≈ RAxAng(-u_ab, -μ_ab)
        @test RAxAng(-u_ab, μ_ab) ≈ RAxAng(u_ab, -μ_ab)
        @test RAxAng(u_ab, μ_ab) ≈ RAxAng(u_ab, μ_ab + 4π)
        @test RAxAng(u_ab, 0) ≈ RAxAng([1, 2, 5], 0)
        @test RAxAng(u_ab, μ_ab) != RAxAng(u_ab, -μ_ab)
        @test !(RAxAng(u_ab, μ_ab) ≈ RAxAng(u_ab, -μ_ab))

        #approximate equality
        @test RAxAng(u_ab, μ_ab) ≈ RAxAng(u_ab, μ_ab + 1e-8)
        @test !(RAxAng(u_ab, μ_ab) ≈ RAxAng(u_ab, -μ_ab))


        #transformation & composition
        @test r_ab * (r_bc * x_c) ≈ (r_ab ∘ r_bc) * x_c

        #inversion
        @test x_c ≈ r_bc' * (r_bc * x_c)

    end

    @testset "Heterogeneous Operators" begin

        #equality must account for different axis and angle sign combinations, and
        #2kπ increments
        q_ab = RQuat(r_ab); q_bc = RQuat(r_bc)

        #approximate equality
        @test r_ab ≈ q_ab #some loss of precision should be expected
        @test !(r_ab ≈ q_bc) #some loss of precision should be expected

        #transformation & composition
        @test q_ab * (r_bc * x_c) ≈ r_ab * (q_bc * x_c)
        @test (q_ab ∘ r_bc) * x_c ≈ (r_ab ∘ q_bc) * x_c

        #inversion
        @test r_ab' ≈ q_ab'
    end
end

function test_REuler()

    @testset "Constructors" begin

    end

    @testset "Operators" begin

    end

    @testset "Conversions" begin

    end
end


#RAxAng: to and from Quat
#REuler: to and from Quat
#RMatrix: to and from RQuat #to be tested in test_RMatrix

#indirect equality!!
#RQuat == RAxAng(r) == REuler(r) == RMatrix(r)

#indirect conversions, in an ad_hoc test method. create torture test from RAxAng
#then do: RAxAng |> REuler |> RAxAng |> RMatrix |> REuler |> RMatrix |> RAxAng
#and test equality
#RAxAng <-> REuler OK
#RAxAng <-> RMatrix OK
#RMatrix <-> REuler OK



end #module