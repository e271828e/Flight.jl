module TestControl

using Test
using BenchmarkTools
using UnPack
using ComponentArrays
using ControlSystems
using StaticArrays

using Flight.FlightCore.Systems
using Flight.FlightCore.Sim
using Flight.FlightCore.Plotting

using Flight.FlightComponents.Control

export test_control

function test_control()
    @testset verbose = true "Control" begin
        test_state_space()
        test_pi_continuous()
        test_pid_discrete()
        test_discrete_integrator()
        test_discrete_lead()
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

        @testset verbose = true "Submodel" begin

            #with ComponentVectors
            cmp = build_ss(x0, u0, y0)
            cmp = submodel(cmp; x = (:V, :q), y = (:V, :q, :f_z))
            @test (length(cmp.x0) == 2 && length(cmp.y0) == 3 && size(cmp.A) == (2,2) &&
                size(cmp.B) == (2, 2) && size(cmp.C) == (3, 2) && size(cmp.D) == (3,2))

            #with plain Vectors
            cmp = build_ss(x0[:], u0[:], y0[:])
            cmp = submodel(cmp; x = [1, 3], u = [1], y = [1, 3, 5])
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

        #test the numerical correctness of the PID discretization

        #first, define an arbitrary PID through its transfer function, convert
        #it to a continuous LinearStateSpace System and simulate it for a unit
        #step input
        k_p = 1
        k_i = 1
        k_d = 0.2
        τ_d = 0.1
        pid_tf = k_p + k_i * tf(1, [1,0]) + k_d * tf([1, 0], [τ_d, 1])

        pid_ss = ss(pid_tf)
        pid_lss = LinearStateSpace(pid_ss) |> System
        pid_lss.u .= 1
        sim = Simulation(pid_lss; dt = 0.0001, t_end = 2)
        Sim.run!(sim)
        th_lss = TimeHistory(sim)
        th_y_lss = (Sim.get_components(th_lss) |> collect)[1]
        y_lss_last = Sim.get_data(th_y_lss)[end]

        #then, define the equivalent discrete PID and simulate it for a unit
        #step input
        pid_disc = PIDDiscrete{1}(; k_p, k_i, k_d, τ_d) |> System
        pid_disc.u.setpoint .= 1
        sim = Simulation(pid_disc; Δt = 0.0001, t_end = 2)
        Sim.run!(sim)
        th_disc = TimeHistory(sim)
        th_y_disc = (Sim.get_components(th_disc.out) |> collect)[1]
        y_disc_last = Sim.get_data(th_y_disc)[end]

        #compare the final values
        @test y_lss_last ≈ y_disc_last atol=1e-4

        end #testset

end #function


function test_discrete_integrator(save = false)

    @testset verbose = true "DiscreteIntegrator" begin

        sys = DiscreteIntegrator() |> System;
        sim = Simulation(sys; Δt = 0.01)

        sys.u.bound_lo = -1
        sys.u.bound_hi = 2

        sys.u.u1 = -1
        step!(sim, 2, true)

        @test sys.s.x0 <= -1
        @test sys.y.y1 ≈ -1.0
        @test sys.y.sat_y1 == -1
        @test sys.y.halted

        sys.u.u1 = 1
        step!(sim, 2, true)
        @test sys.y.sat_y1 == 0
        @test !sys.y.halted

        step!(sim, 2, true)
        @test sys.y.sat_y1 == 1
        @test sys.y.halted

        sys.u.u1 = -1
        step!(sim, 1, true)
        @test sys.y.sat_y1 == 0
        @test !sys.y.halted

        sys.u.sat_ext = -sign(sys.u.u1)
        step!(sim, 1, true)
        @test !sys.y.halted
        sys.u.sat_ext = sign(sys.u.u1)
        step!(sim, 1, true)
        @test sys.y.halted

        Control.reset!(sys)

        @test sys.s.x0 == 0
        @test sys.s.sat_y0 == 0
        @test sys.y.x1 == 0
        @test sys.y.y1 == 0
        @test sys.y.sat_y1 == 0
        @test !sys.y.halted
        @test sys.u.bound_lo != 0
        @test sys.u.bound_hi != 0

        @test @ballocated($f_ode!($sys)) == 0
        @test @ballocated($f_disc!($sys, 1)) == 0
        @test @ballocated($f_step!($sys)) == 0

        end #testset

end #function


function test_discrete_lead(save = false)

    @testset verbose = true "DiscreteLead" begin

        z, p, k = -1, -10, 2.5
        sys = DiscreteLead() |> System;
        @pack! sys.u = z, p, k

        sys_io! = let
            function (sys)
                t = sys.t[]
                sys.u.u1 = sin(t)
            end
        end

        sim = Simulation(sys; Δt = 0.001, t_end = 10, sys_io!)
        Sim.run!(sim)
        th = TimeHistory(sim)

        lead_cont = zpk([z], [p], k)
        step_result = lsim(lead_cont, (x,t)->SVector(sin(t),), 0:0.001:10)
        @test Sim.get_data(th.y1[end])[1] ≈ step_result.y[end] atol = 1e-3

        Control.reset!(sys)
        @test sys.s.u0 == 0
        @test sys.s.x0 == 0
        @test sys.y.u1 == 0
        @test sys.y.y1 == 0
        @test sys.y.p != 0
        @test sys.y.z != 0

        @test @ballocated($f_ode!($sys)) == 0
        @test @ballocated($f_disc!($sys, 1)) == 0
        @test @ballocated($f_step!($sys)) == 0

        end #testset

end #function

end #module