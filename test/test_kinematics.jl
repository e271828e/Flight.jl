module TestKinematics

using Test
using LinearAlgebra
using Flight.Kinematics
using Flight.Attitude
using Flight.WGS84

export test_kinematics

function test_kinematics()
    @testset verbose = true "WGS84 Kinematics" begin
        @testset verbose = true "Initialization" begin test_init() end
    end
end

function test_init()

    init = KinInit()
    x = XKinWGS84(init)

    init = KinInit(
        n_b = RQuat([1, 2, 3, -2]),
        Ob = WGS84Pos(ϕ = π/3, λ = -π/6, h = 1500),
        ω_lb_b = [0.1, 0.01, -0.4],
        v_eOb_b = [100, 5, -10])
    x = XKinWGS84(init)

    #now rebuild the initializer from the kinematic state vector and check it
    #against the original
    l_b = RQuat(x.att.l_b)
    e_l = RQuat(x.pos.e_l)
    h = x.pos.h[1]
    ω_eb_b = x.vel.ω_eb_b
    v_eOb_b = x.vel.v_eOb_b

    Ob = WGS84Pos(NVector(e_l), h)
    (R_N, R_E) = radii(Ob)
    n_b = l_b
    v_eOb_n = n_b * v_eOb_b
    ω_el_n = [
        v_eOb_n[2] / (R_E + h),
        -v_eOb_n[1] / (R_N + h),
        0.0]
    ω_el_b = n_b' * ω_el_n
    ω_lb_b = ω_eb_b - ω_el_b

    init_test = KinInit(n_b = l_b, Ob = Ob, ω_lb_b = ω_lb_b, v_eOb_b = v_eOb_b )

    @test init.n_b ≈ init_test.n_b
    @test init.Ob ≈ init_test.Ob
    @test init.ω_lb_b ≈ init_test.ω_lb_b
    @test init.v_eOb_b ≈ init_test.v_eOb_b

end



end #module