module TestGeneric

using Test
using BenchmarkTools
using UnPack
using ComponentArrays

using Flight

export test_generic

function test_generic()
    @testset verbose = true "Generic" begin
        test_state_space()
        test_pi_compensator()
    end
end

function test_state_space()

    function build_ss(x0, u0, y0)

        ẋ0 = copy(x0)

        A = x0 * x0'
        B = x0 * u0'
        C = y0 * x0'
        D = y0 * u0'

        return StateSpaceModel(ẋ0, x0, u0, y0, A, B, C, D)

    end

    function test_system(ss)

        @unpack ẋ0, x0, u0, y0, A, B, C, D = ss

        sys = System(ss)

        x = 2sys.params.x0
        u = 3sys.params.u0
        sys.x .= x
        sys.u .= u
        f_ode!(sys)

        @test sys.ẋ == ẋ0 + A * (x - x0) + B * (u - u0)
        @test sys.y == y0 + C * (x - x0) + D * (u - u0)

        @test @ballocated(f_ode!($sys)) === 0
        @test @ballocated(f_step!($sys)) === 0
        @test @ballocated(f_disc!($sys, 1)) === 0

    end

    @testset verbose = true "StateSpaceModel" begin

        x0 = ComponentVector(V = 1.0, q = 0.5, θ = 0.3, α = 5.0)
        u0 = ComponentVector(e = 0.1, a = 0.2)
        y0 = ComponentVector(V = 0.3, q = 0.8, θ = 2.0, α = 3.0, f_z = -9.8)

        @testset verbose = true "System" begin
            build_ss(x0, u0, y0) |> test_system #with ComponentVectors
            build_ss(x0[:], u0[:], y0[:]) |> test_system #with plain Vectors
        end

        @testset verbose = true "Filtering" begin

            #with ComponentVectors
            cmp = build_ss(x0, u0, y0)
            cmp = filter(cmp; x = (:V, :q), y = (:V, :q, :f_z))
            @test (length(cmp.x0) == 2 && length(cmp.y0) == 3 && size(cmp.A) == (2,2) &&
                size(cmp.B) == (2, 2) && size(cmp.C) == (3, 2) && size(cmp.D) == (3,2))

            #with plain Vectors
            cmp = build_ss(x0[:], u0[:], y0[:])
            cmp = filter(cmp; x = [1, 3], u = [1], y = [1, 3, 5])
            @test (length(cmp.x0) == 2 && length(cmp.u0) == 1 && length(cmp.y0) == 3 &&
            size(cmp.A) == (2,2) && size(cmp.B) == (2, 1) && size(cmp.C) == (3, 2) && size(cmp.D) == (3,1))

        end

    end

end


function test_pi_compensator(save = false)

    @testset verbose = true "PICompensator" begin

        comp = PICompensator{3}(k_p = 1.0, k_i = 1.0, k_l = 0.0, bounds = (-1, 2));
        sys = System(comp)
        sim = Simulation(sys)

        sys.u.input .= 1.0
        sys.u.sat_enable[2:3] .= false
        step!(sim, 2, true)
        @test sys.y.out[1] == 2.0
        @test sys.y.sat_status[1] == 1
        @test sys.y.out[2] == sys.y.out[3] > sys.y.out[1]

        sys.u.input .= -1.0
        step!(sim, 3, true)
        @test sys.y.out[1] == -1.0
        @test sys.y.sat_status[1] == -1

        sys.u.reset[2] = true
        step!(sim, 2, true)
        @test sys.y.out[2] != 0 #integrator disabled, but we still get proportional output

        sys.u.input[3] = 0
        sys.u.reset[3] = true
        @test f_step!(sys) == true
        @test sys.x[3] == 0 #sys.x changes immediately
        f_ode!(sys)
        @test sys.y.state[3] == 0 #but sys.y needs f_ode! to update
        @test sys.y.out[3] == 0 #but sys.y needs f_ode! to update
        @test f_step!(sys) == false #once reset, no further changes to sys.x[3]

        @test @ballocated($f_ode!($sys)) == 0
        @test @ballocated($f_step!($sys)) == 0

        plots = make_plots(TimeHistory(sim); Plotting.defaults...)
        save ? save_plots(plots, save_folder = joinpath("tmp", "pi_test")) : nothing

    end


end



end #module