using Flight
using Sockets

function test_xp12()

    # https://developer.x-plane.com/article/using-the-correct-wing-datarefs/

    t = 0.02

    δe = 0.5
    δa = 0.5
    δr = 0.5
    δf = 0.5

    ac.y.vehicle.components.pwp.propeller.ω
    ω_prop = 100 #rad/s
    ϕ_prop = ω_prop * t

    ac.y.vehicle.components.ldg.nose.strut.ψ_sw
    ψ_sw = 0.5 #rad

    #first create a generic XPlaneOutput, then specialize it for C172 with the
    #appropriate initialization overrides and extract_output

    # TODO: add engine behaviour based on RPM threshold (disc / not disc). when
    # non disc, maybe rotate it according to simulation time

    #TODO: add landing gear override
    #add engine RPM!!

    override_pose = "sim/operation/override/override_planepath[0]"
    override_surf = "sim/operation/override/override_control_surfaces[0]"
    override_prop = "sim/flightmodel2/engines/prop_disc/override[0]"
    override_nws = "sim/operation/override/override_wheel_steer[0]"

    #all of these must be handled by the specific extract_output method
    rele_pos = "sim/flightmodel2/wing/elevator1_deg[8]"
    lele_pos = "sim/flightmodel2/wing/elevator1_deg[9]"
    lflap_pos = "sim/flightmodel2/wing/flap1_deg[0]"
    rflap_pos = "sim/flightmodel2/wing/flap1_deg[1]"
    rud_pos = "sim/flightmodel2/wing/rudder1_deg[10]"
    lail_pos = "sim/flightmodel2/wing/aileron1_deg[0]"
    rail_pos = "sim/flightmodel2/wing/aileron1_deg[1]"

    prop_is_disc = "sim/flightmodel2/engines/prop_is_disc[0]"
    prop_angle = "sim/flightmodel2/engines/prop_rotation_angle_deg[0]"

    nws_angle = "sim/flightmodel2/gear/tire_steer_actual_deg[0]"

    address = IPv4("192.168.1.2")
    port = 49000
    xpc = XPlaneOutput(; address, port)
    IODevices.init!(xpc)

    #first four must go in init! the rest must go in

    # msg_tuple = (
    #     Network.xpmsg_set_dref(override_pose, 1),
    #     Network.xpmsg_set_dref(override_surf, 1),
    #     Network.xpmsg_set_dref(override_prop, 1),
    #     Network.xpmsg_set_dref(override_nws, 1),
    #     Network.xpmsg_set_dref(lele_pos, rad2deg(δe)),
    #     Network.xpmsg_set_dref(rele_pos, rad2deg(δe)),
    #     Network.xpmsg_set_dref(lail_pos, rad2deg(δa)),
    #     Network.xpmsg_set_dref(rail_pos, rad2deg(-δa)),
    #     Network.xpmsg_set_dref(lflap_pos, rad2deg(δf)),
    #     Network.xpmsg_set_dref(rflap_pos, rad2deg(δf)),
    #     Network.xpmsg_set_dref(rud_pos, rad2deg(δr)),
    #     Network.xpmsg_set_dref(prop_is_disc, 0),
    #     Network.xpmsg_set_dref(prop_angle, rad2deg(ϕ_prop)),
    #     Network.xpmsg_set_dref(nws_angle, rad2deg(ψ_sw)),
    # )
    # IODevices.handle_data!(xpc, msg_tuple)

    IODevices.shutdown!(xpc)

end