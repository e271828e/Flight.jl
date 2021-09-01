module TestKinematics

using Test
using LinearAlgebra
using BenchmarkTools

using Flight.Kinematics
using Flight.Attitude
using Flight.WGS84

export test_kinematics

function test_kinematics()
    @testset verbose = true "WGS84 Kinematics" begin
        @testset verbose = true "Initialization" begin test_init() end
        @testset verbose = true "Position Update" begin test_init() end
    end
end

function test_init()

    init = KinInit(
        q_nb = RQuat([1, 2, 3, -2]),
        Ob = NVectorAlt(ϕ = π/3, λ = -π/6, h = 1500),
        ω_lb_b = [0.1, 0.01, -0.4],
        v_eOb_b = [100, 5, -10])

    x = x0(init)
    #now rebuild the initializer from the kinematic state vector and check it
    #against the original
    q_lb = RQuat(x.pos.q_lb, normalization = false)
    q_el = RQuat(x.pos.q_el, normalization = false)
    h = x.pos.h[1]
    ω_eb_b = x.vel.ω_eb_b
    v_eOb_b = x.vel.v_eOb_b

    Ob = NVectorAlt(NVector(q_el), h)
    (R_N, R_E) = radii(Ob)
    q_nb = q_lb
    v_eOb_n = q_nb * v_eOb_b
    ω_el_n = [
        v_eOb_n[2] / (R_E + h),
        -v_eOb_n[1] / (R_N + h),
        0.0]
    ω_el_b = q_nb' * ω_el_n
    ω_lb_b = ω_eb_b - ω_el_b

    init_test = KinInit(q_nb = q_lb, Ob = Ob, ω_lb_b = ω_lb_b, v_eOb_b = v_eOb_b )

    @test init.q_nb ≈ init_test.q_nb
    @test init.Ob ≈ init_test.Ob
    @test init.ω_lb_b ≈ init_test.ω_lb_b
    @test init.v_eOb_b ≈ init_test.v_eOb_b

end


function test_fpos()

    init = KinInit(
        q_nb = RQuat([1, 2, 3, -2]),
        Ob = NVectorAlt(ϕ = π/3, λ = -π/6, h = 1500),
        ω_lb_b = [0.1, 0.01, -0.4],
        v_eOb_b = [100, 5, -10])

    x_kin = x0(init)

    y_kin = Y(Kin())
    ẋ_pos = x0(Pos())

    @btime f_kin!($y_kin, $ẋ_pos, $x_kin)

end

end #module