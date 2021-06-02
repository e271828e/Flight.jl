module TestAttitude

using Test
using LinearAlgebra
using Flight.Attitude
using Flight.Quaternions

export test_attitude

function test_attitude()
    @testset verbose = true "Attitude" begin
        @testset verbose = true "RQuat" begin test_RQuat() end
        @testset verbose = true "AxisAng" begin test_AxAng() end
    end
end

function test_RQuat()
    #construction from UnitQuat or anything convertible to UnitQuat, other
    #constructors will be provided by other Rotation subtypes
    v_ab = rand(4)
    v_bc = rand(4)
    q_ab = RQuat(v_ab)
    q_bc = RQuat(v_bc)
    q_ab_neg = RQuat(-v_ab)
    x_c = rand(3)

    @testset "Constructors" begin
        #RQuat accepts any AbstractVector of length 4 and normalizes it, no need to
        #pass a UnitQuat instance
        @test q_ab._quat ≈ RQuat(UnitQuat(v_ab))._quat
        @test_throws DimensionMismatch RQuat([1,2,3])
        @test RQuat()._quat[:] == [1,0,0,0]
    end

    @testset "Normalization" begin
        r_bad = RQuat(UnitQuat([1,2,3,4], enforce_norm = false))
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

        #overloads for inversion, transformation and composition
        @test q_ab' == invert(q_ab)
        @test q_ab ∘ q_bc == compose(q_ab, q_bc)
        @test q_bc * x_c == transform(q_bc, x_c)

        #transformation & composition
        @test q_ab * (q_bc * x_c) ≈ (q_ab ∘ q_bc) * x_c

        #inversion
        @test x_c ≈ q_bc' * (q_bc * x_c)

        #time derivative
        ω_ab_b = rand(3)
        @test dt(q_ab, ω_ab_b) == 0.5 * (q_ab._quat * Quat(imag = ω_ab_b))

    end


end

function test_AxAng()

    @testset "Constructors" begin

        u_ab = rand(3)
        μ_ab = rand()
        u_μ_ab = (u_ab, μ_ab)

        #inner constructor
        r_ab = AxAng(u_ab, μ_ab, enforce_norm = false)
        @test r_ab.axis == u_ab
        r_ab = AxAng(u_ab, μ_ab)
        @test r_ab.axis == normalize(u_ab)

        #tuple constructor
        @test r_ab.axis == AxAng(u_μ_ab).axis && r_ab.angle == AxAng(u_μ_ab).angle

        #zero argument constructor
        @test AxAng().axis == [1, 0, 0] && AxAng().angle == 0

    end

    #these actually build upon quaternion conversion, but that's an
    #implementation detail we should not care about, so we can test these before
    @testset "Operators" begin

        #equality must account for different axis and angle sign combinations, and
        #2kπ increments
        @test_broken AxAng(u_ab, μ_ab) == AxAng(-u_ab, -μ_ab)
        @test_broken AxAng(-u_ab, μ_ab) == AxAng(u_ab, -μ_ab)
        @test_broken AxAng(u_ab, μ_ab) == AxAng(u_ab, μ_ab + 4π)
        @test_broken AxAng(u_ab, 0) == AxAng(rand(3), 0)
        @test_broken AxAng(u_ab, μ_ab) != AxAng(u_ab, -μ_ab)

        #approximate equality
        @test_broken AxAng(u_ab, μ_ab) ≈ AxAng(u_ab, μ_ab + 1e-8)
        @test_broken !(AxAng(u_ab, μ_ab) ≈ AxAng(u_ab, -μ_ab))

        r_ab = AxAng(rand(3), rand())
        r_bc = AxAng(rand(3), rand())

        #overloads for inversion, transformation and composition
        @test_broken r_ab' == invert(r_ab) #just switch angle sign!
        @test_broken r_ab ∘ r_bc == compose(r_ab, r_bc)
        @test_broken r_bc * x_c == transform(r_bc, x_c)

        #transformation & composition
        @test_broken r_ab * (r_bc * x_c) ≈ (r_ab ∘ r_bc) * x_c

        #inversion
        @test_broken x_c ≈ r_bc' * (r_bc * x_c)

    end

    @testset "RQuat Conversion" begin

        #sanity checks
        @test_broken AxAng(RQuat(r_ab)) ≈ RQuat([cos(r_ab.angle/2), r_ab.axis * sin(r_ab.angle/2)])
        @test_broken RQuat(AxAng(rand(3), 0)) == RQuat([1,0,0,0])
        @test_broken RQuat(AxAng(rand(3), 2π)) == RQuat([1,0,0,0])

        #[-2π, 2π] including small angles
        u = normalize(rand(3))
        logrange(x1, x2, n) = (10^y for y in range(log10(x1), log10(x2), length=n))
        μ_gen = logrange(1e-10, 2π, 10)
        μ_array = [-reverse(collect(μ_gen)); collect(μ_gen)]
        r_array = [AxAng(u, μ) for μ in μ_array]
        @test_broken r_array ≈ [AxAng(RQuat(r)) for r in r_array]

    end

end

end #module