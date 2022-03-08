module TestPropellers

using Test
using Zygote
using BenchmarkTools
using LinearAlgebra
using Interpolations

using Flight

using Flight.Propellers: DefaultAirfoil, cL, cD, cL_α
using Flight.Propellers: Blade, Coefficients, Dataset

export test_propellers

function test_propellers()
    @testset verbose = true "Propellers" begin
        test_default_airfoil()
        test_coefficients()
        test_dataset()
        test_propeller()
    end
end

function test_default_airfoil()

    @testset verbose = true "DefaultAirfoil" begin

        airfoil = DefaultAirfoil()

        α = range(-π/6, π/3, length = 10)
        M = range(0, 1.5, length = 5)

        cL_α_auto_array = Array{Float64}(undef, (length(α), length(M)))
        cL_α_analytic_array = similar(cL_α_auto_array)

        for (j, M) in enumerate(M)
            cL_α_auto = let airfoil = airfoil, M = M
                α -> cL(airfoil, α, M)
            end
            for (i, α) in enumerate(α)
                cL_α_analytic_array[i,j] = cL_α(airfoil, α, M)
                cL_α_auto_array[i,j] = gradient(cL_α_auto, α)[1]
            end
        end

    @test all(cL_α_analytic_array .≈ cL_α_auto_array)

    end #testset

end

function test_coefficients()

    @testset verbose = true "Coefficients" begin

        n = 2
        b = Blade()
        # b = Blade(c̃ = Propellers.ConstantFunction(0.0573))

        coeffs_static = Coefficients(b, n; J=0, M_tip=0, Δβ = 0.0)
        coeffs_moving = Coefficients(b, n; J=0.5, M_tip=0, Δβ = 0.0)

        @test coeffs_static.η_p == 0
        @test coeffs_static.C_Fx > 0
        @test coeffs_static.C_Mx < 0
        @test coeffs_static.C_Fz_α == 0
        @test coeffs_static.C_Mz_α == 0
        @test coeffs_static.C_P < 0

        @test coeffs_moving.η_p > 0
        @test coeffs_moving.C_Fx < coeffs_static.C_Fx
        @test abs(coeffs_moving.C_Mx) < abs(coeffs_static.C_Mx)
        @test coeffs_moving.C_Fz_α < 0
        @test coeffs_moving.C_Mz_α < 0
        @test abs(coeffs_moving.C_P) < abs(coeffs_static.C_P)

    end

end

function test_dataset()

    @testset verbose = true "Dataset" begin

        @testset verbose = true "FixedPitch" begin

            fpd = Dataset(FixedPitch(), Blade(), 2; n_J = 20, n_M_tip = 5)

            J_bounds, M_tip_bounds = bounds(fpd)
            J = range(J_bounds[1], J_bounds[2], length = 100)
            M_tip = range(M_tip_bounds[1], stop = M_tip_bounds[2], step = 0.4)

            data = [fpd(J, M_tip, 0) for (J, M_tip) in Iterators.product(J, M_tip)]

            @test @ballocated($fpd(0.1, 0.3, 0)) == 0
            @test (data isa Array{Coefficients{Float64}})

        end #testset

        @testset verbose = true "VariablePitch" begin

            vpd = Dataset(FixedPitch(), Blade(), 2; n_J = 20, n_M_tip = 5, n_Δβ = 5)

            J_bounds, M_tip_bounds, Δβ_bounds = bounds(vpd)
            J = range(J_bounds[1], J_bounds[2], length = 100)
            M_tip = range(M_tip_bounds[1], stop = M_tip_bounds[2], step = 0.4)
            Δβ = range(Δβ_bounds[1], stop = Δβ_bounds[2], length = 5)

            data = [vpd(J, M_tip, Δβ) for (J, M_tip, Δβ) in Iterators.product(J, M_tip, Δβ)]

            @test @ballocated($vpd(0.1, 0.3, 0.05)) == 0
            @test (data isa Array{Coefficients{Float64}})

        end #testset

    end #testset

end

function test_propeller()

    t_bp = FrameTransform(r = [1.0, 0, 0])
    kin = KinInit(v_eOb_b = [50, 0, 5]) |> KinData #positive α
    atm = AtmosphereDescriptor() |> System
    air = AirData(kin, atm)
    ω = 300

    @testset verbose = true "Propeller" begin

        @testset verbose = true "FixedPitch" begin

            pitch = FixedPitch()
            sense = Propellers.CCW
            fp_sys = Propeller(; pitch, sense, t_bp) |> System

            @test_throws AssertionError f_cont!(fp_sys, kin, air, ω)

            f_cont!(fp_sys, kin, air, -ω)

            wr_p = fp_sys.y.wr_p
            @test wr_p.F[1] > 0
            @test wr_p.F[3] < 0
            @test wr_p.M[1] > 0 #in a CCW propeller should be positive along x
            @test wr_p.M[3] > 0 #in a CCW propeller should be positive along z

            @test @ballocated(f_cont!($fp_sys, $kin, $air, -$ω)) == 0
            @test @ballocated(f_disc!($fp_sys)) == 0


        end #testset

        @testset verbose = true "VariablePitch" begin

            pitch = VariablePitch((-deg2rad(5), deg2rad(10)))
            sense = Propellers.CW
            vp_sys = Propeller(; pitch, sense, t_bp) |> System

            vp_sys.u[] = 0
            f_cont!(vp_sys, kin, air, ω)
            @test vp_sys.y.Δβ ≈ vp_sys.params.pitch.bounds[1]
            Fx_0 = vp_sys.y.wr_p.F[1]

            vp_sys.u[] = 1
            f_cont!(vp_sys, kin, air, ω)
            @test vp_sys.y.Δβ ≈ vp_sys.params.pitch.bounds[2]
            Fx_1 = vp_sys.y.wr_p.F[1]

            @test Fx_1 > Fx_0

            @test @ballocated(f_cont!($vp_sys, $kin, $air, $ω)) == 0
            @test @ballocated(f_disc!($vp_sys)) == 0

        end #testset

    end #testset

end


end