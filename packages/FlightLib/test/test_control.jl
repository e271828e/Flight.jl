module TestControl

using Test
using BenchmarkTools
using ControlSystems
using StaticArrays

using FlightCore
using FlightLib

#import individual components to avoid potential namespace conflicts with ControlSystems
import FlightLib.Linearization: LinearizedSS
import FlightLib.Control.Continuous: PIVector as PIContinuous
import FlightLib.Control.Discrete: PID as PIDDiscrete
import FlightLib.Control.Discrete: PIDVector as PIDDiscreteVector
import FlightLib.Control.Discrete: LeadLag as LeadLagDiscrete
import FlightLib.Control.Discrete: Integrator as IntegratorDiscrete
import FlightLib.Control.Discrete: IntegratorVector as IntegratorDiscreteVector
import FlightLib.Control.Discrete: LQR as DiscreteLQR

export test_control

function test_control()
    @testset verbose = true "Control" begin
        test_continuous_pi()
        test_discrete_integrator()
        test_discrete_integrator_vector()
        test_discrete_leadlag()
        test_discrete_pid()
        test_discrete_pid_vector()
        test_discrete_lqr()
    end
end


function test_continuous_pi(save = false)

    @testset verbose = true "Continuous PIVector" begin

        comp = PIContinuous{2}();
        mdl = Model(comp)
        sim = Simulation(mdl)

        mdl.parameters.k_p .= 1.0
        mdl.parameters.k_i .= 1.0
        mdl.parameters.bound_lo[1] = -1
        mdl.parameters.bound_hi[1] = 1
        mdl.u.input .= -1.0
        step!(sim, 2, true)

        @test mdl.y.sat_out[1] == -1
        @test mdl.y.int_halted[1] #integrator 1 should have been halted
        @test abs(mdl.y.y_i[1]) < 0.1
        @test mdl.y.output[1] ≈ -1.0

        @test mdl.y.sat_out[2] == 0
        @test !mdl.y.int_halted[2] #integrator 2 should have not
        @test mdl.y.y_i[2] ≈ -2.0 atol = 1e-2
        @test mdl.y.output[2] ≈ -3.0 atol = 1e-2

        mdl.u.sat_ext[2] = -sign(mdl.y.u_i[2]) #set opposite external saturation
        step!(sim, 1, true)
        @test !mdl.y.int_halted[2] #integrator 2 should have not halted
        mdl.u.sat_ext[2] = sign(mdl.y.u_i[1]) #set same sign saturation
        step!(sim, 1, true)
        @test mdl.y.int_halted[2] #integrator 2 should have halted

        @test mdl.x[2] != 0
        @test mdl.u.sat_ext[2] != 0
        @test mdl.y.output[2] != 0
        init!(mdl)
        @test mdl.x[2] == 0
        @test mdl.u.sat_ext[2] == 0
        f_ode!(mdl) #but output takes one call to f_ode! to update
        @test mdl.y.output[2] == 0

        @test @ballocated(f_ode!($mdl)) == 0
        @test @ballocated(f_step!($mdl)) == 0
        @test @ballocated(f_periodic!(NoScheduling(), $mdl)) == 0

        save && save_plots(TimeSeries(sim), normpath("tmp/test_control/test_continuous_pi"); Plotting.defaults...)

    end #testset

end #function


function test_discrete_integrator()

    @testset verbose = true "Discrete Integrator" begin

        mdl = IntegratorDiscrete() |> Model;
        sim = Simulation(mdl; Δt = 0.01)

        mdl.parameters.bound_lo[] = -1
        mdl.parameters.bound_hi[] = 2

        mdl.u.input = -1
        step!(sim, 2, true)

        @test mdl.s.x0 <= -1
        @test mdl.y.output ≈ -1.0
        @test mdl.y.sat_out == -1
        @test mdl.y.halted

        mdl.u.input = 1
        step!(sim, 2, true)
        @test mdl.y.sat_out == 0
        @test !mdl.y.halted

        step!(sim, 2, true)
        @test mdl.y.sat_out == 1
        @test mdl.y.halted

        mdl.u.input = -1
        step!(sim, 1, true)
        @test mdl.y.sat_out == 0
        @test !mdl.y.halted

        mdl.u.sat_ext = -sign(mdl.u.input)
        step!(sim, 1, true)
        @test !mdl.y.halted
        mdl.u.sat_ext = sign(mdl.u.input)
        step!(sim, 1, true)
        @test mdl.y.halted

        init!(mdl)

        @test mdl.s.x0 == 0
        @test mdl.s.sat_out_0 == 0
        @test mdl.parameters.bound_lo != -Inf #should be preserved
        @test mdl.parameters.bound_hi != Inf #should be preserved

        f_periodic!(mdl)
        @test mdl.y.x1 == 0
        @test mdl.y.output == 0
        @test mdl.y.sat_out == 0
        @test !mdl.y.halted

        @test @ballocated(f_ode!($mdl)) == 0
        @test @ballocated(f_step!($mdl)) == 0
        @test @ballocated(f_periodic!(NoScheduling(), $mdl)) == 0

        end #testset

end #function



function test_discrete_integrator_vector()

    @testset verbose = true "Discrete Vector Integrator" begin

        mdl = IntegratorDiscreteVector{2}() |> Model;
        sim = Simulation(mdl; Δt = 0.01)

        mdl.parameters.bound_lo .= -1
        mdl.parameters.bound_hi .= 2

        mdl.u.input .= -1
        step!(sim, 2, true)

        @test mdl.s.x0[2] <= -1
        @test mdl.y.output[2] ≈ -1.0
        @test mdl.y.sat_out[2] == -1
        @test mdl.y.halted[2]

        mdl.u.input .= 1
        step!(sim, 2, true)
        @test mdl.y.sat_out[1] == 0
        @test !mdl.y.halted[1]

        step!(sim, 2, true)
        @test mdl.y.sat_out[2] == 1
        @test mdl.y.halted[2]

        mdl.u.input .= -1
        step!(sim, 1, true)
        @test mdl.y.sat_out[1] == 0
        @test !mdl.y.halted[1]

        mdl.u.sat_ext[1] = -sign(mdl.u.input[1])
        step!(sim, 1, true)
        @test !mdl.y.halted[1]
        mdl.u.sat_ext[1] = sign(mdl.u.input[1])
        step!(sim, 1, true)
        @test mdl.y.halted[1]

        init!(mdl)

        @test mdl.s.x0[2] == 0
        @test mdl.s.sat_out_0[2] == 0
        @test mdl.parameters.bound_lo[2] != -Inf #should be preserved
        @test mdl.parameters.bound_hi[2] != Inf #should be preserved

        f_periodic!(mdl)
        @test mdl.y.x1[2] == 0
        @test mdl.y.output[2] == 0
        @test mdl.y.sat_out[2] == 0
        @test !mdl.y.halted[2]

        @test @ballocated(f_ode!($mdl)) == 0
        @test @ballocated(f_step!($mdl)) == 0
        @test @ballocated(f_periodic!(NoScheduling(), $mdl)) == 0

        end #testset

end #function


function test_discrete_leadlag(save = false)

    @testset verbose = true "Discrete LeadLag" begin

        z, p, k = -1, -10, 2.5
        mdl = LeadLagDiscrete() |> Model;
        mdl.parameters.z[] = z
        mdl.parameters.p[] = p
        mdl.parameters.k[] = k

        user_callback! = let
            function (mdl)
                t = mdl.t[]
                mdl.u.u1 = sin(t)
            end
        end

        sim = Simulation(mdl; Δt = 0.001, t_end = 10, user_callback!)
        Sim.run!(sim)
        ts = TimeSeries(sim)

        lead_cont = zpk([z], [p], k)
        step_result = lsim(lead_cont, (x,t)->SVector(sin(t),), 0:0.001:10)
        @test Sim.get_data(ts.y1[end])[1] ≈ step_result.y[end] atol = 1e-3

        init!(mdl)
        @test mdl.s.u0 == 0
        @test mdl.s.x0 == 0

        f_periodic!(mdl)
        @test mdl.y.u1 == 0
        @test mdl.y.y1 == 0
        @test mdl.y.p != 0
        @test mdl.y.z != 0

        @test @ballocated(f_ode!($mdl)) == 0
        @test @ballocated(f_step!($mdl)) == 0
        @test @ballocated(f_periodic!(NoScheduling(), $mdl)) == 0

        end #testset

end #function

function test_discrete_pid(save = false)

    @testset verbose = true "Discrete PID" begin

        mdl = PIDDiscrete() |> Model;
        sim = Simulation(mdl; Δt = 0.01)

        mdl.u.input = 1.0
        mdl.parameters.k_p[] = 1.0
        mdl.parameters.k_i[] = 1.0
        mdl.parameters.k_d[] = 0.1
        step!(sim, 1, true)

        @test mdl.y.y_p == 1.0
        @test mdl.y.y_i ≈ 1.0
        @test mdl.y.out_free ≈ 2.0
        @test mdl.y.output ≈ 2.0
        @test mdl.y.sat_out == 0
        @test !mdl.y.int_halted

        mdl.parameters.bound_lo[] = -1
        mdl.parameters.bound_hi[] = 1
        step!(sim, 1, true)
        @test mdl.y.out_free > 2.0
        @test mdl.y.output ≈ 1
        @test mdl.y.sat_out == 1
        @test mdl.y.int_halted

        mdl.u.input = -1.0
        step!(sim, 2, true)
        @test mdl.y.sat_out == -1
        @test mdl.y.int_halted

        mdl.u.input = 0.1
        step!(sim, 1, true)
        @test mdl.y.sat_out == 0
        @test !mdl.y.int_halted

        mdl.u.sat_ext = -sign(mdl.y.u_i) #set opposite external saturation
        step!(sim)
        @test !mdl.y.int_halted #integrator should have not halted

        mdl.u.sat_ext = sign(mdl.y.u_i) #set same sign saturation
        step!(sim)
        @test mdl.y.int_halted #integrator 2 should have halted

        init!(mdl)

        @test mdl.u.input == 0
        @test mdl.u.sat_ext == 0
        @test mdl.parameters.bound_lo[] != -Inf #should be preserved
        @test mdl.parameters.bound_hi[] != Inf #should be preserved
        @test mdl.s.x_i0 == 0
        @test mdl.s.x_d0 == 0
        @test mdl.s.sat_out_0 == 0

        f_periodic!(mdl)
        @test mdl.y.y_i == 0
        @test mdl.y.y_d == 0
        @test mdl.y.y_p == 0
        @test !mdl.y.int_halted

        @test @ballocated(f_ode!($mdl)) == 0
        @test @ballocated(f_step!($mdl)) == 0
        @test @ballocated(f_periodic!(NoScheduling(), $mdl)) == 0

        save && save_plots(TimeSeries(sim), normpath("tmp/test_control/test_discrete_pid"); Plotting.defaults...)

        #operate PID as a filtered derivative
        init!(mdl)
        mdl.parameters.k_p[] = 0.0
        mdl.parameters.k_i[] = 0.0
        mdl.parameters.k_d[] = 1.0
        mdl.parameters.τ_f[] = 0.2
        mdl.parameters.bound_lo[] = -Inf
        mdl.parameters.bound_hi[] = Inf

        sim = Simulation(mdl; Δt = 0.01)
        mdl.u.input = 1.0
        step!(sim)
        @test mdl.y.y_d > 0.0
        step!(sim, 5, true)
        @test mdl.y.y_d ≈ 0.0 atol = 1e-6 #must have returned to zero

        #test the numerical correctness of the PID discretization

        #define an arbitrary PID through its transfer function, convert it to a
        #continuous LinearizedSS Model and simulate it for a unit step input
        k_p = 1
        k_i = 1
        k_d = 0.2
        τ_f = 0.1
        pid_tf = k_p + k_i * tf(1, [1,0]) + k_d * tf([1, 0], [τ_f, 1])

        pid_ss = ss(pid_tf)
        pid_lss = LinearizedSS(pid_ss) |> Model
        pid_lss.u .= 1
        sim = Simulation(pid_lss; dt = 0.0001, t_end = 2)
        Sim.run!(sim)
        ts_lss = TimeSeries(sim)
        ts_y_lss = (Sim.get_components(ts_lss) |> collect)[1]
        y_lss_last = Sim.get_data(ts_y_lss)[end]

        #define the equivalent discrete PID and simulate it for a unit step input
        init!(mdl)
        mdl.parameters.k_p[] = k_p
        mdl.parameters.k_i[] = k_i
        mdl.parameters.k_d[] = k_d
        mdl.parameters.τ_f[] = τ_f

        mdl.u.input = 1
        sim = Simulation(mdl; Δt = 0.0001, t_end = 2)
        Sim.run!(sim)
        ts_disc = TimeSeries(sim)
        ts_y_disc = ts_disc.output
        y_disc_last = Sim.get_data(ts_y_disc)[end]

        #compare the final values
        @test y_lss_last ≈ y_disc_last atol=1e-3

        end #testset

    end #function

function test_discrete_pid_vector(save = false)

    @testset verbose = true "Discrete Vector PID" begin

        mdl = PIDDiscreteVector{2}() |> Model;
        sim = Simulation(mdl; Δt = 0.01)

        mdl.u.input .= 1.0
        mdl.parameters.k_p .= 1.0
        mdl.parameters.k_i .= 1.0
        mdl.parameters.k_d .= 0.1
        step!(sim, 1, true)

        @test mdl.y.y_p[1] == 1.0
        @test mdl.y.y_i[1] ≈ 1.0
        @test mdl.y.out_free[1] ≈ 2.0
        @test mdl.y.output[1] ≈ 2.0
        @test mdl.y.sat_out[1] == 0
        @test !mdl.y.int_halted[1]

        mdl.parameters.bound_lo .= -1
        mdl.parameters.bound_hi .= 1
        step!(sim, 1, true)
        @test mdl.y.out_free[1] > 2.0
        @test mdl.y.output[1] ≈ 1
        @test mdl.y.sat_out[1] == 1
        @test mdl.y.int_halted[1]

        mdl.u.input .= -1.0
        step!(sim, 2, true)
        @test mdl.y.sat_out[1] == -1
        @test mdl.y.int_halted[1]

        mdl.u.input .= 0.1
        step!(sim, 1, true)
        @test mdl.y.sat_out[1] == 0
        @test !mdl.y.int_halted[1]

        mdl.u.sat_ext[1] = -sign(mdl.y.u_i[1]) #set opposite external saturation
        step!(sim)
        @test !mdl.y.int_halted[1] #integrator should have not halted

        mdl.u.sat_ext[1] = sign(mdl.y.u_i[1]) #set same sign saturation
        step!(sim)
        @test mdl.y.int_halted[1] #integrator 2 should have halted

        init!(mdl)

        @test mdl.u.input[1] == 0
        @test mdl.u.sat_ext[1] == 0
        @test mdl.parameters.bound_lo[1] != 0
        @test mdl.parameters.bound_hi[1] != 0
        @test mdl.s.x_i0[1] == 0
        @test mdl.s.x_d0[1] == 0
        @test mdl.s.sat_out_0[1] == 0

        f_periodic!(mdl)
        @test mdl.y.y_i[1] == 0
        @test mdl.y.y_d[1] == 0
        @test mdl.y.y_p[1] == 0
        @test !mdl.y.int_halted[1]

        @test @ballocated(f_ode!($mdl)) == 0
        @test @ballocated(f_step!($mdl)) == 0
        @test @ballocated(f_periodic!(NoScheduling(), $mdl)) == 0

        save && save_plots(TimeSeries(sim), normpath("tmp/test_control/test_discrete_pid_vector"); Plotting.defaults...)

        #operate PID as a filtered derivative
        init!(mdl)
        mdl.u.input .= 1.0
        mdl.parameters.k_p .= 0.0
        mdl.parameters.k_i .= 0.0
        mdl.parameters.k_d .= 1.0
        mdl.parameters.τ_f .= 0.2
        mdl.parameters.bound_lo .= -Inf
        mdl.parameters.bound_hi .= Inf

        sim = Simulation(mdl; Δt = 0.01)
        step!(sim)
        @test mdl.y.y_d[1] > 0.0
        step!(sim, 5, true)
        @test mdl.y.y_d[1] ≈ 0.0 atol = 1e-6 #must have returned to zero

        #test the numerical correctness of the PID discretization

        #define an arbitrary PID through its transfer function, convert it to a
        #continuous LinearizedSS Model and simulate it for a unit step input
        k_p = 1
        k_i = 1
        k_d = 0.2
        τ_f = 0.1
        pid_tf = k_p + k_i * tf(1, [1,0]) + k_d * tf([1, 0], [τ_f, 1])

        pid_ss = ss(pid_tf)
        pid_lss = LinearizedSS(pid_ss) |> Model
        pid_lss.u .= 1
        sim = Simulation(pid_lss; dt = 0.0001, t_end = 2)
        Sim.run!(sim)
        ts_lss = TimeSeries(sim)
        ts_y_lss = (Sim.get_components(ts_lss) |> collect)[1]
        y_lss_last = Sim.get_data(ts_y_lss)[end]

        #define the equivalent discrete PID and simulate it
        pid_disc = PIDDiscreteVector{1}() |> Model
        pid_disc.u.input .= 1.0
        pid_disc.parameters.k_p .= k_p
        pid_disc.parameters.k_i .= k_i
        pid_disc.parameters.k_d .= k_d
        sim = Simulation(pid_disc; Δt = 0.0001, t_end = 2)
        Sim.run!(sim)
        ts_disc = TimeSeries(sim)
        ts_y_disc = (Sim.get_components(ts_disc.output) |> collect)[1]
        y_disc_last = Sim.get_data(ts_y_disc)[end]

        #compare the final values
        @test y_lss_last ≈ y_disc_last atol=1e-3

        end #testset

end #function

function test_discrete_lqr()
    @testset verbose = true "Discrete LQR" begin

        mdl = DiscreteLQR{10, 2, 2}() |> Model

        @test @ballocated(f_ode!($mdl)) === 0
        @test @ballocated(f_step!($mdl)) === 0
        @test @ballocated(f_periodic!(NoScheduling(), $mdl)) === 0

        #functionality verified via C172Xv1 integration tests
    end
end




end #module