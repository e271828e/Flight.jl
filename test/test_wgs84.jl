module TestWGS84

using Test
using LinearAlgebra
using Flight.WGS84

export test_wgs84

function test_wgs84()
    @testset verbose = true "WGS84" begin
        @testset verbose = true "WGS84Pos Essentials" begin test_essentials() end
        @testset verbose = true "WGS84Pos Conversions" begin test_conversions() end
        # @testset verbose = true "WGS84 Models" begin test_models() end #ellipsoid, angular velocity, gravity, gravitation
    end
end

function test_NVectorAlt()

    n_e = [1, 2, -3]
    h = 1000
    p1 = NVectorAlt(n_e, h)

    #inner constructor
    @test p1.n_e ≈ normalize(n_e) #normalizes by default
    @test NVectorAlt(n_e, h, normalization = false).n_e == n_e
    @test p1 === NVectorAlt(p1)

    #equality and isapprox
    @test p1 == p1
    @test p1 != NVectorAlt(-n_e, h)
    @test p1 ≈ NVectorAlt(n_e, h+1e-10)
    @test !(p1 ≈ NVectorAlt(-n_e, h))

    #outer constructors
    @test p1 == NVectorAlt((n_e, h))  #accept also a tuple
    @test NVectorAlt(n_e).n_e ≈ normalize(n_e) #accept also a single unit vector
    @test NVectorAlt(h = 324).n_e == [1, 0, 0] && NVectorAlt(n_e = n_e).h == 0
    @test NVectorAlt() == NVectorAlt(n_e = [1, 0, 0], h = 0)
    @test_throws AssertionError NVectorAlt(h = -500) #ADJUST THE LOWER ALT THRESHOLD

end

function test_Cartesian()

    n_e = [1, 2, -3]
    h = 1000

    #constructors and equality
    p1 = Cartesian(NVectorAlt(n_e, h))
    @test p1 == Cartesian(p1.r)
    @test p1 === Cartesian(p1)

    #AbstractArray interface


    #equality and isapprox
    @test P1 == P1
    @test P1 != NVectorAlt(-n_e, h)
    @test P1 ≈ NVectorAlt(n_e, h+1e-10)
    @test !(P1 ≈ NVectorAlt(-n_e, h))

    #outer constructors
    @test P1 == NVectorAlt((n_e, h))  #accept also a tuple
    @test NVectorAlt(n_e).n_e ≈ normalize(n_e) #accept also a single unit vector
    @test NVectorAlt(h = 324).n_e == [1, 0, 0] && NVectorAlt(n_e = n_e).h == 0
    @test NVectorAlt() == NVectorAlt(n_e = [1, 0, 0], h = 0)
    @test_throws AssertionError NVectorAlt(h = -500) #ADJUST THE LOWER ALT THRESHOLD

end

function test_conversions()

    @show pN = NVectorAlt([1e-16, 1e-16, 1], 10000)
    @show NVectorAlt(Cartesian(pN))
    # @show pS = NVectorAlt([0.1, 0.2, -100], 10000)
    # @show NVectorAlt(Cartesian(pS))
    # @show pEqN = NVectorAlt([1, 2, 0.001], 10000)
    # @show NVectorAlt(Cartesian(pEqN))
    # @show pEqS = NVectorAlt([1, 2, -0.001], 10000)
    # @show NVectorAlt(Cartesian(pEqS))

end

end


#CHECK VALUES AGAINST PYTHON!

    # @test
    # r_nl = RQuat([1,-2,4,3]) #ECEF to LTF (e)
    # LTF(r_nl)._r_nl ≈ LTF(RAxAng(r_nl)) ≈ LTF(REuler(r_nl)) ≈
    # P3 = Geodetic(1, 2, 1000)


    #only allow to construct R from other Location subtypes

    #test convergence for low altitudes