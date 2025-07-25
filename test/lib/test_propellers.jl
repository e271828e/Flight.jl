module TestPropellers

using Test
using BenchmarkTools
using LinearAlgebra
using Interpolations: bounds
using FiniteDiff: finite_difference_derivative

using Flight.FlightCore
using Flight.FlightLib

#non-exported stuff
using Flight.FlightLib.Propellers: DefaultAirfoil, cL, cD, cL_α
using Flight.FlightLib.Propellers: Blade, Coefficients, Lookup

export test_propellers

function test_propellers()
    @testset verbose = true "Propellers" begin
        test_default_airfoil()
        test_coefficients()
        test_lookup()
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
                cL_α_auto_array[i,j] = finite_difference_derivative(cL_α_auto, α)
            end
        end

    @test all(cL_α_analytic_array .≈ cL_α_auto_array)

    end #testset

end

function test_coefficients()

    @testset verbose = true "Coefficients" begin

        n = 2
        b = Blade()
        # b = Blade(c̃ = Propellers.ConstantDistribution(0.0573))

        coeffs_static = Coefficients(n, b, 0, 0, 0.0)
        coeffs_moving = Coefficients(n, b, 0.5, 0, 0.0)

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

function test_lookup()

    @testset verbose = true "Lookup" begin

        J_range = range(0, 1.5, length = 31)
        Mt_range = range(0, 1.5, length = 21)
        Δβ_range = range(0, 0.1, length = 11)
        lookup = Lookup(3, Blade(); J_range, Mt_range, Δβ_range)

        #flat extrapolation
        @test lookup.C_Fx(1.5, 1.5, 0.1) == lookup.C_Fx(3, 2, 2)

        #saving and loading lookup data
        @test size(lookup.data.C_Fx) === (31, 21, 11)
        fname = tempname()
        Propellers.save_lookup(lookup, fname)
        lookup_loaded = Propellers.load_lookup(fname)

        test_point = (J = 0.123, Mt = 0.643, Δβ = 0.04)
        foreach(fieldnames(Coefficients)) do name
            @test getproperty(lookup, name)(test_point...) ==
                    getproperty(lookup_loaded, name)(test_point...)
        end

        #coefficient lookup must not allocate
        @test lookup(test_point...) isa Coefficients{Float64}
        # @test @ballocated($lookup($test_point...)) == 0

        #make sure there are no issues with the singleton Δβ dimension in a
        #fixed pitch lookup
        lookup_fp = Lookup(3, Blade();
                           J_range, Mt_range, Δβ_range = range(0,0,length=1))
        @test lookup_fp(0.1, 0.3, 0) == lookup(0.1, 0.3, 0)

    end #testset

end

function test_propeller()

    t_bp = FrameTransform(r = [1.0, 0, 0])
    kin = KinInit(v_eb_n = [50, 0, 0]) |> KinData
    atm = AtmosphericData()
    air = AirflowData(atm, kin)
    ω = 300

    @testset verbose = true "Propeller" begin

        @testset verbose = true "FixedPitch" begin

            sense = Propellers.CCW
            fp_mdl = Propeller(FixedPitch(); sense, t_bp) |> Model
            @test fp_mdl.u |> isnothing

            f_ode!(fp_mdl, kin, air, -ω)

            wr_p = fp_mdl.y.wr_p
            @test wr_p.F[1] > 0
            @test wr_p.F[3] < 0
            @test wr_p.τ[1] > 0 #in a CCW propeller should be positive along x
            @test wr_p.τ[3] > 0 #in a CCW propeller should be positive along z

            @test @ballocated(f_ode!($fp_mdl, $kin, $air, -$ω)) == 0
            @test @ballocated(f_step!($fp_mdl)) == 0


        end #testset

        @testset verbose = true "VariablePitch" begin

            sense = Propellers.CW
            vp_mdl = Propeller(VariablePitch(-0.1, 0.2); sense, t_bp) |> Model

            vp_mdl.u[] = 0
            f_ode!(vp_mdl, kin, air, ω)
            @test vp_mdl.y.Δβ ≈ bounds(vp_mdl.lookup)[3][1] #lower Δβ bound
            Fx_0 = vp_mdl.y.wr_p.F[1]

            vp_mdl.u[] = 1
            f_ode!(vp_mdl, kin, air, ω)
            @test vp_mdl.y.Δβ ≈ bounds(vp_mdl.lookup)[3][2] #upper Δβ bound
            Fx_1 = vp_mdl.y.wr_p.F[1]

            @test Fx_1 > Fx_0

            @test @ballocated(f_ode!($vp_mdl, $kin, $air, $ω)) == 0
            @test @ballocated(f_step!($vp_mdl)) == 0

        end #testset

    end #testset

end


end