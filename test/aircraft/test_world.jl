module TestWorld

using Test
using UnPack
using BenchmarkTools
using Sockets

using Flight

export test_world

function test_world()
    @testset verbose = true "World" begin

        test_system_methods()

    end
end

function test_system_methods()

    h_trn = HOrth(608.55);

    ac = Cessna172Rv0();
    env = SimpleEnvironment(trn = HorizontalTerrain(altitude = h_trn))
    world = SimpleWorld(ac, env) |> System

    kin_init = KinematicInit(
        v_eOb_n = [0, 0, 0],
        ω_lb_b = [0, 0, 0],
        q_nb = REuler(ψ = 0, θ = 0.0, φ = 0.),
        loc = LatLon(ϕ = deg2rad(40.503205), λ = deg2rad(-3.574673)),
        h = h_trn + 1.8 + 0);

    init_kinematics!(world, kin_init)

    @testset verbose = true "System Methods" begin

        #get a copy of the initial system output
        y0 = world.y
        @test y0 == world.y

        #modify inputs and call f_ode!
        world.u.env.atm.wind.v_ew_n[1] = -5

        #if the right method is called, the modified inputs must propagate to
        #the outputs
        f_ode!(world)
        @test world.y.env.atm.wind.v_ew_n[1] == -5
        @test world.y.ac.air.v_ew_n[1] == -5

        #reset outputs to their initial value
        world.y = y0
        @test y0 == world.y

        #and repeat for f_disc!, which must also propagate inputs to outputs
        f_disc!(world, 0.02)
        @test world.y.env.atm.wind.v_ew_n[1] == -5
        @test world.y.ac.air.v_ew_n[1] == -5

        #mess up a quaternion norm
        world.ac.x.kinematics.pos.q_lb[1] *= 2

        #and make sure the call to f_step! restores it
        x_mod = f_step!(world)
        @test x_mod == true
        @test world.ac.x.kinematics.pos.q_lb[1] ≈ 1

        #make sure we are on the ground to ensure landing gear code coverage
        @test world.ac.y.airframe.ldg.left.strut.wow == true
        @test @ballocated(f_ode!($world)) == 0
        @test @ballocated(f_step!($world)) == 0
        @test @ballocated(f_disc!($world, 0.2)) == 0

        # return nothing
    end

end



end #module