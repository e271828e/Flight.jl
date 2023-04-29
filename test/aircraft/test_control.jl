module TestControl

using Test
using BenchmarkTools
using UnPack
using ComponentArrays

using Flight

export test_control

function test_control()
    @testset verbose = true "Control" begin
        test_state_space()
        test_pi_continuous()
        test_pid_discrete()
    end
end

function test_state_space()

    function build_ss(x0, u0, y0)

        ẋ0 = copy(x0)

        A = x0 * x0'
        B = x0 * u0'
        C = y0 * x0'
        D = y0 * u0'

        return LinearStateSpace(ẋ0, x0, u0, y0, A, B, C, D)

    end

    function test_system(ss)

        @unpack ẋ0, x0, u0, y0, A, B, C, D = ss

        sys = System(ss)

        x = 2sys.params.x0
        u = 3sys.params.u0
        sys.x .= x
        sys.u .= u
        f_ode!(sys)

        @test all(sys.ẋ .== ẋ0 + A * (x - x0) + B * (u - u0))
        @test all(sys.y .== y0 + C * (x - x0) + D * (u - u0))

        @test @ballocated(f_ode!($sys)) === 0
        @test @ballocated(f_step!($sys)) === 0
        @test @ballocated(f_disc!($sys, 1)) === 0

    end

    @testset verbose = true "LinearStateSpace" begin

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

function test_pi_continuous(save = false)

    @testset verbose = true "PIContinuous" begin

        comp = PIContinuous{2}(k_p = 1.0, k_i = 1.0, k_l = 0.0);
        sys = System(comp)
        sim = Simulation(sys)

        sys.u.setpoint .= 0.0
        sys.u.feedback .= 1.0
        sys.u.bound_lo[1] = -1
        sys.u.bound_hi[1] = 1
        step!(sim, 2, true)

        @test sys.y.sat_out[1] == -1
        @test sys.y.int_halt[1] #integrator 1 should have been halted
        @test abs(sys.y.y_i[1]) < 0.1
        @test sys.y.out[1] ≈ -1.0

        @test sys.y.sat_out[2] == 0
        @test !sys.y.int_halt[2] #integrator 2 should have not
        @test sys.y.y_i[2] ≈ -2.0 atol = 1e-2
        @test sys.y.out[2] ≈ -3.0 atol = 1e-2
        y_i0 = sys.y.y_i

        sys.u.sat_ext[2] = -sign(sys.y.u_i[2]) #set opposite external saturation
        step!(sim, 1, true)
        @test !sys.y.int_halt[2] #integrator 2 should have not halted
        sys.u.sat_ext[2] = sign(sys.y.u_i[1]) #set same sign saturation
        step!(sim, 1, true)
        @test sys.y.int_halt[2] #integrator 2 should have halted

        sys.u.anti_windup[1] = false #disable anti-windup in the first component
        step!(sim, 2, true)
        @test sys.y.sat_out[1] == -1 #output still saturated
        @test !sys.y.int_halt[1] #but integrator no longer halted
        @test abs(sys.y.y_i[1]) > abs(y_i0[1]) #integrator should have kept accumulating

        sys.u.reset[2] = true
        @test f_step!(sys) == true
        @test sys.x[2] == 0 #sys.x changes immediately
        f_ode!(sys)
        @test sys.y.y_i[2] == 0 #but sys.y needs f_ode! to update
        @test sys.y.out[2] == 0 #idem
        @test f_step!(sys) == false #once reset, no further changes to sys.x[3]

        @test @ballocated($f_ode!($sys)) == 0
        @test @ballocated($f_disc!($sys, 1)) == 0
        @test @ballocated($f_step!($sys)) == 0

        plots = make_plots(TimeHistory(sim); Plotting.defaults...)
        save && save_plots(plots, save_folder = joinpath("tmp", "pi_continuous_test"))

    end #testset

end #function


function test_pid_discrete(save = false)

    @testset verbose = true "PIDDiscrete" begin

        sys = PIDDiscrete{2}(k_p = 1.0, k_i = 1.0, k_d = 0.0) |> System;
        sim = Simulation(sys; Δt = 0.01)

        sys.u.setpoint .= 0.0
        sys.u.feedback .= 1.0

        sys.u.bound_lo[1] = -1
        sys.u.bound_hi[1] = 1
        step!(sim, 2, true)

        @test sys.y.sat_out[1] == -1
        @test sys.y.int_halt[1] #integrator 1 should have been halted
        @test abs(sys.y.y_i[1]) < 0.1
        @test sys.y.out[1] ≈ -1.0

        @test sys.y.sat_out[2] == 0
        @test !sys.y.int_halt[2] #integrator 2 should have not
        @test sys.y.y_i[2] ≈ -2.0
        @test sys.y.out[2] ≈ -3.0

        sys.u.sat_ext[2] = -sign(sys.y.u_i[2]) #set opposite external saturation
        step!(sim, 1, true)
        @test !sys.y.int_halt[2] #integrator 2 should have not halted
        sys.u.sat_ext[2] = sign(sys.y.u_i[1]) #set same sign saturation
        step!(sim, 1, true)
        @test sys.y.int_halt[2] #integrator 2 should have halted

        y_i0 = sys.y.y_i

        sys.u.anti_windup[1] = false #disable anti-windup in the first component
        step!(sim, 2, true)
        @test sys.y.sat_out[1] == -1 #output still saturated
        @test !sys.y.int_halt[1] #but integrator no longer halted
        @test abs(sys.y.y_i[1]) > abs(y_i0[1]) #integrator should have kept accumulating

        sys.u.anti_windup[1] = true
        sys.u.reset .= true
        step!(sim) #let it propagate

        @test all(sys.s.x_i0 .== 0) #state must have been reset
        @test all(sys.y.out .== 0) #output is nulled
        sys.u.reset .= false

        @test @ballocated($f_ode!($sys)) == 0
        @test @ballocated($f_disc!($sys, 1)) == 0
        @test @ballocated($f_step!($sys)) == 0

        plots = make_plots(TimeHistory(sim); Plotting.defaults...)
        save && save_plots(plots, save_folder = joinpath("tmp", "pid_discrete_test"))

        #operate PID as a filtered derivative
        sys = PIDDiscrete{1}(k_p = 0.0, k_i = 0.0, k_d = 1.0, τ_d = 0.2, β_d = 1.0) |> System;
        sim = Simulation(sys; Δt = 0.01)

        step!(sim, 1, true)
        sys.u.setpoint .= 0.0
        sys.u.feedback .= 1.0
        sys.u.bound_lo[1] = -1
        sys.u.bound_hi[1] = 1

        step!(sim, 0.02, true)
        @test sys.y.y_d[1] < 0.0 #feedback is positive, error has decreased, derivative path output must be negative
        step!(sim, 5, true)
        @test sys.y.y_d[1] ≈ 0.0 atol = 1e-6 #must have returned to zero

        # plots = make_plots(TimeHistory(sim); Plotting.defaults...)
        # save && save_plots(plots, save_folder = joinpath("tmp", "pid_discrete_test"))

        return

    end #testset

end #function


end #module