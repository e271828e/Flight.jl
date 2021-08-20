module TestAttitude

using Test
using LinearAlgebra
using Flight.Attitude
using Flight.Quaternions

export test_attitude

logrange(x1, x2, n) = (10^y for y in range(log10(x1), log10(x2), length=n))

function lograngesym(x1, x2, n)
    log_semi = collect(logrange(x1, x2, n))
    log_sym = [-reverse(log_semi); log_semi]
    return (i for i in log_sym)
end


function test_attitude()
    @testset verbose = true "Attitude" begin
        @testset verbose = true "RQuat" begin test_RQuat() end
        @testset verbose = true "RMatrix" begin test_RMatrix() end
        @testset verbose = true "RAxAng" begin test_RAxAng() end
        @testset verbose = true "REuler" begin test_REuler() end
    end
end

function test_RQuat()

        v_ab = [0.7, -1, 0.1, 1]
        v_bc = [1.7, 2, -1, 1]
        q_ab = RQuat(v_ab)
        q_bc = RQuat(v_bc)
        q_ab_neg = RQuat(-v_ab)
        x_c = [2, 3, -1.1]

    @testset "Constructors" begin
        #RQuat accepts any AbstractVector of length 4 and normalizes it, no need to
        #pass a UnitQuat instance
        @test q_ab._u ≈ RQuat(UnitQuat(v_ab))._u
        @test_throws DimensionMismatch RQuat([1,2,3])
        @test RQuat()._u[:] == [1,0,0,0]
    end

    @testset "Normalization" begin
        r_bad = RQuat(UnitQuat([1,2,3,4], normalization = false))
        @test abs(norm(r_bad._u)-1) > 1e-2
        @test norm(normalize(r_bad)) ≈ 1
        # normalize!(r_bad)
        # @test norm(r_bad) ≈ 1
    end

    @testset "Operators" begin
        #equality
        @test q_ab != q_bc
        @test q_ab == q_ab_neg
        @test q_ab == q_ab

        #isapprox
        @test q_ab ≈ q_ab_neg
        @test q_ab ≈ q_ab
        @test !(q_ab ≈ q_bc)

        #transformation, inversion & composition
        @test q_ab * (q_bc * x_c) ≈ (q_ab ∘ q_bc) * x_c
        @test x_c ≈ q_bc' * (q_bc * x_c)

        #time derivative
        ω_ab_b = [10, -4, 2]
        @test dt(q_ab, ω_ab_b) == 0.5 * (q_ab._u * Quat(imag = ω_ab_b))

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

    #these actually build upon quaternion conversion, but that's an
    #implementation detail we should not care about, so we can test these before
    @testset "Operators" begin

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

        #transformation, composition & inversion
        @test r_ab * (r_bc * x_c) ≈ (r_ab ∘ r_bc) * x_c
        @test x_c ≈ r_bc' * (r_bc * x_c)

    end

    @testset "Conversions" begin

        #conversions to and from RQuat (REuler tested in the REuler testset)

        #null rotation
        @test RQuat(RAxAng([1, 2, -1], 0)) ≈ RQuat([1,0,0,0])
        @test RQuat(RAxAng([3, 1, 1], 2π)) ≈ RQuat([1,0,0,0])

        #complete [-2π, 2π] angle range focusing on  small angles
        u = normalize([7, -5, 2])
        μ_array = lograngesym(Attitude.ε_null, 2π, 10)
        r_array = collect(RAxAng(u, μ) for μ in μ_array)
        r_array_test = collect(r |> RQuat |> RAxAng for r in r_array)
        @test all(r_array .≈ r_array_test)

    end

    @testset "Mixed Operators" begin

        q_ab = RQuat(r_ab); q_bc = RQuat(r_bc)

        #approximate equality
        @test r_ab ≈ q_ab #some loss of precision should be expected
        @test !(r_ab ≈ q_bc) #some loss of precision should be expected

        #transformation, composition & inversion
        @test q_ab * (r_bc * x_c) ≈ r_ab * (q_bc * x_c)
        @test (q_ab ∘ r_bc) * x_c ≈ (r_ab ∘ q_bc) * x_c
        @test r_ab' ≈ q_ab'
    end
end

function test_REuler()

    @testset "Constructors" begin

        ψ = 1.459; θ = -1.122; φ = 0.454
        r1 = REuler(ψ, θ, φ) #individual arguments
        r2 = REuler((ψ, θ, φ)) #tuple
        r3 = REuler(θ = θ) #tuple
        r4 = REuler() #zero argument
        @test r1.ψ == r2.ψ && r1.θ == r2.θ && r1.φ == r2.φ
        @test r3.θ == θ && r3.ψ == r3.φ == 0
        @test r4.ψ == r4.θ == r4.φ == 0

        @test_throws AssertionError REuler(θ = 1.7)
        @test_throws AssertionError REuler(ψ = 1.65π)
        @test_throws AssertionError REuler(φ = 1.65π)

    end

    @testset "Operators" begin

        r_ab = REuler(-1.234, √2, √3)
        r_bc = REuler(1.34, -√2, √2)
        x_c = [-1, 2, 3]

        #equality
        @test r_ab == r_ab
        @test r_ab != r_bc

        #approximate equality
        @test r_bc ≈ REuler(1.34000001, -√2, √2)
        @test REuler(ψ = π) ≈ REuler(ψ = -π)
        @test REuler(φ = π) ≈ REuler(φ = -π)
        @test !(r_ab ≈ r_bc)

        #transformation, composition & inversion
        @test r_ab * (r_bc * x_c) ≈ (r_ab ∘ r_bc) * x_c
        @test x_c ≈ r_bc' * (r_bc * x_c)

    end

    @testset "Conversions" begin

        #conversions to and from RQuat and RAxAng

        #null rotation
        @test RQuat(REuler()) ≈ RQuat([1,0,0,0])

        #cover Euler angle ranges including interval bounds but avoiding gimbal lock
        ψ_range = range(-π, π, length = 10)
        θ_range = range(-(π/2 - 0.001), π/2 - 0.001, length = 10)
        φ_range = range(-π, π, length = 10)
        r_array = vec([REuler(i) for i in Iterators.product(ψ_range, θ_range, φ_range)])

        r_array_test = [r |> RQuat |> REuler |> RAxAng |> REuler for r in r_array]
        # r_failed = r_array[r_array .≈ r_array_test]
        # r_test_failed = r_array_test[r_array .≈ r_array_test]
        @test all(r_array .≈ r_array_test)

    end

    @testset "Mixed Operators" begin

        #equality
        r_ab = REuler(0.1, -1, 0.6); r_bc = REuler(-1.8, 0.2, -1)
        q_ab = RQuat(r_ab); q_bc = RQuat(r_bc)
        x_c = [2,3,10]

        #approximate equality
        @test r_ab ≈ q_ab #some loss of precision should be expected
        @test !(r_ab ≈ q_bc) #some loss of precision should be expected

        #transformation, composition & inversion
        @test q_ab * (r_bc * x_c) ≈ r_ab * (q_bc * x_c)
        @test (q_ab ∘ r_bc) * x_c ≈ (r_ab ∘ q_bc) * x_c
        @test r_ab' ≈ q_ab'
    end

end

function test_RMatrix()

    A = [0.590969 0.288519 0.902137;
        0.657604 0.927693 0.218535;
        0.370707  0.481724  0.10268]

    A_orth = qr(A).Q
    r_ab = RMatrix(A_orth)
    r_bc = RMatrix(sin.(A_orth)) #generate a different matrix and let the constructor normalize
    x_c = [-1, 2, 3]

    @testset "Constructors" begin

        A_dist = A_orth + 1e-6*sin.(A)
        r_orth = RMatrix(A_orth, normalization = false)
        r_dist = RMatrix(A_dist, normalization = false)
        @test det(r_orth) ≈ 1
        @test !(det(r_dist) ≈ 1)
        @test det(RMatrix(A_dist)) ≈ 1 #automatic normalization on construction

        r_norm = normalize(r_dist)
        @test det(r_norm) ≈ 1
        # r_norm_copy = RMatrix(A_dist, normalization = false)
        # normalize!(r_norm_copy)
        # @test det(r_norm_copy) ≈ 1

        #error after renormalization
        r_err = r_orth ∘ r_norm'
        @test RAxAng(r_err).angle < 1e-7

    end

    @testset "Operators" begin

        #equality
        @test r_ab == r_ab
        @test r_ab != r_bc

        #approximate equality
        @test r_ab ≈ RMatrix(r_ab._mat .+ 1e-10, normalization = false)
        @test !(r_ab ≈ r_bc)

        #transformation, composition & inversion
        @test r_ab * (r_bc * x_c) ≈ (r_ab ∘ r_bc) * x_c
        @test r_ab * (r_bc * x_c) ≈ (r_ab * r_bc) * x_c
        @test x_c ≈ r_bc' * (r_bc * x_c)

        #time derivative
        ω_ab_b = [10, -4, 2]
        @test dt(r_ab, ω_ab_b) == r_ab._mat * Attitude.skew(ω_ab_b)
    end

    @testset "Conversions" begin

        #cover Euler angle ranges including interval bounds but avoiding gimbal lock
        ψ_range = range(-π, π, length = 10)
        θ_range = range(-(π/2 - 0.001), π/2 - 0.001, length = 10)
        φ_range = range(-π, π, length = 10)
        r_array = vec([REuler(i) for i in Iterators.product(ψ_range, θ_range, φ_range)])

        #this should exercise all cases of RMatrix -> RQuat conversion
        r_array_test = [r |> RMatrix |> RQuat |> RMatrix |> REuler |> RMatrix |> RAxAng |> RMatrix for r in r_array]
        @test all(r_array .≈ r_array_test)

    end

    @testset "Mixed Operators" begin

        #equality
        q_ab = RQuat(r_ab)
        q_bc = RQuat(r_bc)

        #approximate equality
        @test r_ab ≈ q_ab #some loss of precision should be expected
        @test !(r_ab ≈ q_bc) #some loss of precision should be expected

        #transformation, composition & inversion
        @test q_ab * (r_bc * x_c) ≈ r_ab * (q_bc * x_c)
        @test (q_ab ∘ r_bc) * x_c ≈ (r_ab ∘ q_bc) * x_c
        @test r_ab' ≈ q_ab'
    end


end

end #module