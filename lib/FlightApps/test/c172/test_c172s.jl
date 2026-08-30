module TestC172S

using Test, BenchmarkTools

using FlightCore
using FlightPhysics
using FlightApps

export test_c172s


function test_c172s()
    @testset verbose = true "Cessna 172S" begin

        test_trimming()
        test_linearization()
        test_update_methods()
        test_ground_start()

    end
end

function test_trimming()

    @testset verbose = true "Trimming" begin

        vehicle = Model(C172S.Vehicle())
        atmosphere = Model(SimpleAtmosphere())
        terrain = Model(UniformTerrain())
        params = C172.TrimParameters()
        state = C172.TrimState()

        f_target = C172.get_f_target(vehicle, params, atmosphere, terrain)

        @test @ballocated($f_target($state)) === 0

        success, _ = f_init!(vehicle, params)

        @test success

    end #testset

end #function

function test_linearization()

    @testset verbose = true "Linearization" begin
        @test_nowarn Cessna172Sv0(NED()) |> Model |> linearize
    end #testset

end #function

function test_update_methods()

    @testset verbose = true "Update Methods" begin

        atmosphere = SimpleAtmosphere() |> Model
        terrain = UniformTerrain() |> Model

        location = NVector()
        trn_data = TerrainData(terrain, location)
        vehicle_init = KinInit( h = HOrth(trn_data) + 1.8) |> C172.Init

        aircraft = Model(Cessna172Sv0());

        f_init!(aircraft, vehicle_init, atmosphere, terrain)

        #ensure we are on the ground for full landing gear code coverage
        @test aircraft.y.vehicle.systems.ldg.left.strut.wow == true

        @test @ballocated(f_ode!($aircraft, $atmosphere, $terrain)) == 0
        @test @ballocated(f_step!($aircraft, $atmosphere, $terrain)) == 0
        @test @ballocated(f_periodic!(Unconditional(), $aircraft, $atmosphere, $terrain)) == 0

    end

    return nothing

end

#the aircraft must settle on its landing gear when started on the ground
function test_ground_start()

    @testset verbose = true "Ground Start" begin

        h_trn = HOrth(0.0)
        world = SimpleWorld(Cessna172Sv0(), SimpleAtmosphere(), UniformTerrain(h_trn)) |> Model
        sim = Simulation(world; dt = 0.02, t_end = 10)
        init!(sim, C172.Init(KinInit(h = h_trn + C172.Δh_to_gnd)))
        Sim.run!(sim)

        (; kinematics, systems) = world.aircraft.y.vehicle
        @test all(leg.strut.wow for leg in systems.ldg)
        @test Float64(kinematics.h_o) ≈ C172.Δh_to_gnd atol = 0.05
        @test maximum(abs, kinematics.v_eb_n) < 1e-3

    end

end


end #module
