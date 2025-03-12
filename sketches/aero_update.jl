using LinearAlgebra, StaticArrays, Interpolations

"""
    aero_update()

Recreate the aerodynamic interpolations for the Cessna 172, with focus on the
flap-deflection related coefficients. This serves as a basis for extending the
flap deflection range from 0-30 degrees to 0-40 degrees.

Returns a NamedTuple with the eight relevant interpolations that depend on flap deflection.
"""
function aero_update()
    # Data copied from generate_aero_lookup() function

    ################################ C_D data ##################################

    # Flap deflection contribution to drag
    C_D_δf = (
        δf = deg2rad.(range(0, 40, step = 10)),
        data = [0.0000, 0.0070, 0.0120, 0.0180, 0.0240]
    )

    # Angle of attack and flap deflection combined effect on drag
    C_D_α_δf = (
        α = deg2rad.(range(-5, 20, step = 1)),
        δf = deg2rad.(range(0, 40, step = 10)),
        data = [
            0.0041 0.0013 0.0001 0.0003 0.0020 0.0052 0.0099 0.0162 0.0240 0.0334 0.0442 0.0566 0.0706 0.0860 0.0962 0.1069 0.1180 0.1298 0.1424 0.1565 0.1727 0.1782 0.1716 0.1618 0.1475 0.1097;
            0.0000 0.0004 0.0023 0.0057 0.0105 0.0168 0.0248 0.0342 0.0452 0.0577 0.0718 0.0874 0.1045 0.1232 0.1353 0.1479 0.1610 0.1746 0.1892 0.2054 0.2240 0.2302 0.2227 0.2115 0.1951 0.1512;
            0.0005 0.0025 0.0059 0.0108 0.0172 0.0251 0.0346 0.0457 0.0583 0.0724 0.0881 0.1053 0.1240 0.1442 0.1573 0.1708 0.1849 0.1995 0.2151 0.2323 0.2521 0.2587 0.2507 0.2388 0.2214 0.1744;
            0.0014 0.0041 0.0084 0.0141 0.0212 0.0299 0.0402 0.0521 0.0655 0.0804 0.0968 0.1148 0.1343 0.1554 0.1690 0.1830 0.1975 0.2126 0.2286 0.2464 0.2667 0.2735 0.2653 0.2531 0.2351 0.1866;
            0.0019 0.0049 0.0097 0.0158 0.0232 0.0323 0.0430 0.0553 0.0691 0.0844 0.1012 0.1196 0.1395 0.1610 0.1749 0.1891 0.2038 0.2191 0.2354 0.2534 0.2740 0.2810 0.2726 0.2603 0.2420 0.1927;
        ]'
    )

    ############################# C_Y data #####################################

    # Sideslip angle and flap deflection combined effect
    C_Y_β_δf = (
        β = deg2rad.([-20, 0, 20]),
        δf = deg2rad.([0, 40]),
        data = [0.1370 0.0957;
                0.0000 0.0000;
                -0.1370 -0.0957]
    )

    # Roll rate derivative with angle of attack and flap deflection
    C_Y_p = (
        α = [0.0, 0.094],
        δf = deg2rad.([0, 30]),
        data = [-0.0750 -0.1897;
                -0.1450 -0.2597]
    )

    # Yaw rate derivative with angle of attack and flap deflection
    C_Y_r = (
        α = [0.0, 0.094],
        δf = deg2rad.([0, 30]),
        data = [0.2140 0.1447;
                0.2670 0.1977]
    )

    ############################# C_L data #####################################

    # Flap deflection contribution to lift
    C_L_δf = (
        δf = deg2rad.(range(0, 40, step = 10)),
        data = [0.0000, 0.2, 0.3, 0.35, 0.375]
    )

    ################################## C_l data ################################

    # Yaw rate derivative with angle of attack and flap deflection for rolling moment
    C_l_r = (
        α = [0.0, 0.094],
        δf = deg2rad.([0, 40]),
        data = [0.0798 0.1395;
                0.1869 0.2466]
    )

    ################################## C_m data ################################

    # Flap deflection contribution to pitching moment
    C_m_δf = (
        δf = deg2rad.(range(0, 40, step = 10)),
        data = [0.0000, -0.0654, -0.0981, -0.1140, -0.122]
    )

    ############################## Interpolations ##############################

    # Create the eight relevant interpolations
    C_D_δf_interp = linear_interpolation(C_D_δf.δf, C_D_δf.data, extrapolation_bc = Line())
    C_D_α_δf_interp = linear_interpolation((C_D_α_δf.α, C_D_α_δf.δf), C_D_α_δf.data, extrapolation_bc = Line())
    C_Y_β_δf_interp = linear_interpolation((C_Y_β_δf.β, C_Y_β_δf.δf), C_Y_β_δf.data, extrapolation_bc = Line())
    C_Y_p_interp = linear_interpolation((C_Y_p.α, C_Y_p.δf), C_Y_p.data, extrapolation_bc = Line())
    C_Y_r_interp = linear_interpolation((C_Y_r.α, C_Y_r.δf), C_Y_r.data, extrapolation_bc = Line())
    C_L_δf_interp = linear_interpolation(C_L_δf.δf, C_L_δf.data, extrapolation_bc = Line())
    C_l_r_interp = linear_interpolation((C_l_r.α, C_l_r.δf), C_l_r.data, extrapolation_bc = Line())
    C_m_δf_interp = linear_interpolation(C_m_δf.δf, C_m_δf.data, extrapolation_bc = Line())

    # Return only the interpolations with their original data
    return (
        C_D_δf = (interp = C_D_δf_interp, orig_x = C_D_δf.δf, orig_y = C_D_δf.data),
        C_D_α_δf = (interp = C_D_α_δf_interp, orig_x1 = C_D_α_δf.α, orig_x2 = C_D_α_δf.δf, orig_y = C_D_α_δf.data),
        C_Y_β_δf = (interp = C_Y_β_δf_interp, orig_x1 = C_Y_β_δf.β, orig_x2 = C_Y_β_δf.δf, orig_y = C_Y_β_δf.data),
        C_Y_p = (interp = C_Y_p_interp, orig_x1 = C_Y_p.α, orig_x2 = C_Y_p.δf, orig_y = C_Y_p.data),
        C_Y_r = (interp = C_Y_r_interp, orig_x1 = C_Y_r.α, orig_x2 = C_Y_r.δf, orig_y = C_Y_r.data),
        C_L_δf = (interp = C_L_δf_interp, orig_x = C_L_δf.δf, orig_y = C_L_δf.data),
        C_l_r = (interp = C_l_r_interp, orig_x1 = C_l_r.α, orig_x2 = C_l_r.δf, orig_y = C_l_r.data),
        C_m_δf = (interp = C_m_δf_interp, orig_x = C_m_δf.δf, orig_y = C_m_δf.data)
    )
end


"""
    visualize_interpolations(interp)

Evaluates the current interpolations from 0 to 40 degrees, generates plots,
and saves them in the aero_update subfolder. This helps visualize how
the current extrapolation methods behave beyond the original data range (0-30 degrees).

Parameters:
- interp: The interpolations from aero_update()
"""
function visualize_interpolations(interp)
    # Automatically determine output directory based on script location
    output_dir = joinpath(@__DIR__, "aero_update")

    # Create a range of flap deflection values from 0 to 40 degrees
    δf_deg = 0:0.5:40
    δf_rad = deg2rad.(δf_deg)

    # Vertical line at 30 degrees to show where extrapolation begins
    extrapolation_boundary = 30.0

    # Set plot defaults
    default(size=(800, 600), fontfamily="Computer Modern", linewidth=2,
            legend=:outertopright, grid=true, framestyle=:box)

    # 1. Plot C_D_δf - Direct flap contribution to drag
    p1 = plot(δf_deg, [interp.C_D_δf.interp(δf) for δf in δf_rad],
              label="C_D_δf", xlabel="Flap Deflection (degrees)", ylabel="C_D",
              title="Drag Coefficient vs Flap Deflection")
    # Add vertical line at 30 degrees
    vline!([extrapolation_boundary], label="Extrapolation Boundary", linestyle=:dash, color=:red)
    # Add points for original data
    scatter!(rad2deg.(interp.C_D_δf.orig_x), interp.C_D_δf.orig_y, label="Original Data", markersize=6)
    savefig(p1, joinpath(output_dir, "CD_flap.png"))

    # 2. Plot C_L_δf - Direct flap contribution to lift
    p2 = plot(δf_deg, [interp.C_L_δf.interp(δf) for δf in δf_rad],
              label="C_L_δf", xlabel="Flap Deflection (degrees)", ylabel="C_L",
              title="Lift Coefficient vs Flap Deflection")
    vline!([extrapolation_boundary], label="Extrapolation Boundary", linestyle=:dash, color=:red)
    scatter!(rad2deg.(interp.C_L_δf.orig_x), interp.C_L_δf.orig_y, label="Original Data", markersize=6)
    savefig(p2, joinpath(output_dir, "CL_flap.png"))

    # 3. Plot C_m_δf - Direct flap contribution to pitching moment
    p3 = plot(δf_deg, [interp.C_m_δf.interp(δf) for δf in δf_rad],
              label="C_m_δf", xlabel="Flap Deflection (degrees)", ylabel="C_m",
              title="Pitching Moment Coefficient vs Flap Deflection")
    vline!([extrapolation_boundary], label="Extrapolation Boundary", linestyle=:dash, color=:red)
    scatter!(rad2deg.(interp.C_m_δf.orig_x), interp.C_m_δf.orig_y, label="Original Data", markersize=6)
    savefig(p3, joinpath(output_dir, "Cm_flap.png"))

    # 4. 2D interpolation visualizations
    # For C_D_α_δf - create a slice at a few different alpha values
    alphas = [-0.05, 0.0, 0.05, 0.1, 0.15, 0.2]
    p4 = plot(title="Drag Coefficient vs Flap Deflection at Various Angles of Attack")
    for α in alphas
        plot!(δf_deg, [interp.C_D_α_δf.interp(α, δf) for δf in δf_rad],
              label=@sprintf("α = %.2f rad", α))
    end
    xlabel!("Flap Deflection (degrees)")
    ylabel!("C_D")
    vline!([extrapolation_boundary], label="Extrapolation Boundary", linestyle=:dash, color=:red)
    savefig(p4, joinpath(output_dir, "CD_alpha_flap.png"))

    # 5. For C_Y_β_δf - create a slice at a few different beta values
    betas = [-0.3, -0.15, 0.0, 0.15, 0.3]
    p5 = plot(title="Side Force Coefficient vs Flap Deflection at Various Sideslip Angles")
    for β in betas
        plot!(δf_deg, [interp.C_Y_β_δf.interp(β, δf) for δf in δf_rad],
              label=@sprintf("β = %.2f rad", β))
    end
    xlabel!("Flap Deflection (degrees)")
    ylabel!("C_Y")
    vline!([extrapolation_boundary], label="Extrapolation Boundary", linestyle=:dash, color=:red)
    savefig(p5, joinpath(output_dir, "CY_beta_flap.png"))

    # 6. For C_Y_p - create a slice at different alpha values
    p6 = plot(title="Side Force Roll Rate Derivative vs Flap Deflection")
    for α in [0.0, 0.094]
        plot!(δf_deg, [interp.C_Y_p.interp(α, δf) for δf in δf_rad],
              label=@sprintf("α = %.3f rad", α))
    end
    xlabel!("Flap Deflection (degrees)")
    ylabel!("C_Y_p")
    vline!([extrapolation_boundary], label="Extrapolation Boundary", linestyle=:dash, color=:red)
    # Add points for original data
    for (i, α) in enumerate(interp.C_Y_p.orig_x1)
        scatter!(rad2deg.(interp.C_Y_p.orig_x2), interp.C_Y_p.orig_y[i,:], label=false, markersize=6)
    end
    savefig(p6, joinpath(output_dir, "CY_p_flap.png"))

    # 7. For C_Y_r - create a slice at different alpha values
    p7 = plot(title="Side Force Yaw Rate Derivative vs Flap Deflection")
    for α in [0.0, 0.094]
        plot!(δf_deg, [interp.C_Y_r.interp(α, δf) for δf in δf_rad],
              label=@sprintf("α = %.3f rad", α))
    end
    xlabel!("Flap Deflection (degrees)")
    ylabel!("C_Y_r")
    vline!([extrapolation_boundary], label="Extrapolation Boundary", linestyle=:dash, color=:red)
    # Add points for original data
    for (i, α) in enumerate(interp.C_Y_r.orig_x1)
        scatter!(rad2deg.(interp.C_Y_r.orig_x2), interp.C_Y_r.orig_y[i,:], label=false, markersize=6)
    end
    savefig(p7, joinpath(output_dir, "CY_r_flap.png"))

    # 8. For C_l_r - create a slice at different alpha values
    p8 = plot(title="Rolling Moment Yaw Rate Derivative vs Flap Deflection")
    for α in [0.0, 0.094]
        plot!(δf_deg, [interp.C_l_r.interp(α, δf) for δf in δf_rad],
              label=@sprintf("α = %.3f rad", α))
    end
    xlabel!("Flap Deflection (degrees)")
    ylabel!("C_l_r")
    vline!([extrapolation_boundary], label="Extrapolation Boundary", linestyle=:dash, color=:red)
    # Add points for original data
    for (i, α) in enumerate(interp.C_l_r.orig_x1)
        scatter!(rad2deg.(interp.C_l_r.orig_x2), interp.C_l_r.orig_y[i,:], label=false, markersize=6)
    end
    savefig(p8, joinpath(output_dir, "Cl_r_flap.png"))

    println("Plots saved to: $output_dir")

    # Return a summary of all plots as a dictionary for reference
    return Dict(
        "CD_flap" => p1,
        "CL_flap" => p2,
        "Cm_flap" => p3,
        "CD_alpha_flap" => p4,
        "CY_beta_flap" => p5,
        "CY_p_flap" => p6,
        "CY_r_flap" => p7,
        "Cl_r_flap" => p8
    )
end


function generate_aero_data_tuple_ext()

    ################################ C_D data ##################################

    C_D = (
        zero = 0.027,
        δe = (
            x = [-1.0, 0.0, 1.0],
            y = [0.06, 0.0, 0.06]
        ),
        β = (
            x = [-1.0, 0.0, 1.0],
            y = [0.17, 0.0, 0.17]
        ),
        ge = (
            x = [0.0000, 0.1000, 0.1500, 0.2000, 0.3000, 0.4000, 0.5000, 0.6000, 0.7000, 0.8000, 0.9000, 1.0000, 1.1000], #non-dimensional height
            y = [0.4800, 0.5150, 0.6290, 0.7090, 0.8150, 0.8820, 0.9280, 0.9620, 0.9880, 1.0000, 1.0000, 1.0000, 1.0000]
        ),
        δf = (
            x = deg2rad.(range(0, 40, step = 10)),
            y = [0.0000, 0.0070, 0.0120, 0.0180, 0.0240]
        ),
        α_δf = (
            x_α = deg2rad.(range(-5, 20, step = 1)),
            x_δf = deg2rad.(range(0, 40, step = 10)),
            y = [
                0.0041 0.0013 0.0001 0.0003 0.0020 0.0052 0.0099 0.0162 0.0240 0.0334 0.0442 0.0566 0.0706 0.0860 0.0962 0.1069 0.1180 0.1298 0.1424 0.1565 0.1727 0.1782 0.1716 0.1618 0.1475 0.1097;
                0.0000 0.0004 0.0023 0.0057 0.0105 0.0168 0.0248 0.0342 0.0452 0.0577 0.0718 0.0874 0.1045 0.1232 0.1353 0.1479 0.1610 0.1746 0.1892 0.2054 0.2240 0.2302 0.2227 0.2115 0.1951 0.1512;
                0.0005 0.0025 0.0059 0.0108 0.0172 0.0251 0.0346 0.0457 0.0583 0.0724 0.0881 0.1053 0.1240 0.1442 0.1573 0.1708 0.1849 0.1995 0.2151 0.2323 0.2521 0.2587 0.2507 0.2388 0.2214 0.1744;
                0.0014 0.0041 0.0084 0.0141 0.0212 0.0299 0.0402 0.0521 0.0655 0.0804 0.0968 0.1148 0.1343 0.1554 0.1690 0.1830 0.1975 0.2126 0.2286 0.2464 0.2667 0.2735 0.2653 0.2531 0.2351 0.1866;
                0.0019 0.0049 0.0097 0.0158 0.0232 0.0323 0.0430 0.0553 0.0691 0.0844 0.1012 0.1196 0.1395 0.1610 0.1749 0.1891 0.2038 0.2191 0.2354 0.2534 0.2740 0.2810 0.2726 0.2603 0.2420 0.1927
        ]'
    )
    )

    ############################# C_Y data #####################################

    C_Y = (
        δr = 0.1870,
        δa = 0.0,
        β_δf = (
            x_β = deg2rad.([-20, 0, 20]),
            x_δf = deg2rad.([0, 40]),
            y = [0.1370 0.0957;
                    0.0000 0.0000;
                    -0.1370 -0.0957]
        ),
        p = (
            x_α = [0.0, 0.094],
            x_δf = deg2rad.([0, 30]),
            y = [-0.0750 -0.1897;
                    -0.1450 -0.2597]
        ),
        r = (
            x_α = [0.0, 0.094],
            x_δf = deg2rad.([0, 30]),
            y = [0.2140 0.1447;
                    0.2670 0.1977]
        )
    )

    ############################# C_L data #####################################

    C_L = (
        δe = 0.4300,
        q = 3.900,
        α_dot = 1.700,
        ge = (
            x = [0.0000, 0.1000, 0.1500, 0.2000, 0.3000, 0.4000, 0.5000, 0.6000, 0.7000, 0.8000, 0.9000, 1.0000, 1.1000],
            y = [1.2030, 1.1270, 1.0900, 1.0730, 1.0460, 1.0550, 1.0190, 1.0130, 1.0080, 1.0060, 1.0030, 1.0020, 1.0000]
        ),
        α = (
            x_α = [-0.0900, 0.0000, 0.0900, 0.1000, 0.1200, 0.1400, 0.1600, 0.1700, 0.1900, 0.2100, 0.2400, 0.2600, 0.2800, 0.3000, 0.3200, 0.3400, 0.3600],
            x_stall = [0.0, 1.0],
            y = [
                -0.2200  0.2500  0.7300  0.8300  0.9200  1.0200  1.0800  1.1300  1.1900  1.2500  1.3500  1.4400  1.4700  1.4300  1.3800  1.3000  1.1500;
                -0.2200  0.2500  0.7300  0.7800  0.7900  0.8100  0.8200  0.8300  0.8500  0.8600  0.8800  0.9000  0.9200  0.9500  0.9900  1.0500  1.1500
            ]'
        ),
        δf = (
            x = deg2rad.(range(0, 40, step = 10)),
            y = [0.0000, 0.2, 0.3, 0.35, 0.375]
        )
    )

    ################################## C_l data ################################

    C_l = (
        δa = 0.229,
        δr = 0.0147,
        β = -0.09226,
        p = -0.4840,
        r = (
            x_α = [0.0, 0.094],
            x_δf = deg2rad.([0, 40]),
            y = [0.0798 0.1395;
                    0.1869 0.2466]
        )
    )

    ################################## C_m data ################################

    C_m = (
        zero = 0.100,
        δe = -1.1220,
        α = -1.8000,
        q = -12.400,
        α_dot = -7.2700,
        δf = (
            x = deg2rad.(range(0, 40, step = 10)),
            y = [0.0000, -0.0654, -0.0981, -0.1140, -0.122]
        )
    )

    ################################## C_n data ################################
    C_n = (
        δr = -0.0430,
        δa = -0.0053,
        β = 0.05874,
        p = -0.0278,
        r = -0.0937
    )

    ############################## Interpolations ##############################
    C_D_interp = (
        z = C_D.zero,
        β = linear_interpolation(C_D.β.x, C_D.β.y, extrapolation_bc = Flat()),
        δe = linear_interpolation(C_D.δe.x, C_D.δe.y, extrapolation_bc = Flat()),
        δf = linear_interpolation(C_D.δf.x, C_D.δf.y, extrapolation_bc = Flat()),
        α_δf = linear_interpolation((C_D.α_δf.x_α, C_D.α_δf.x_δf), C_D.α_δf.y, extrapolation_bc = Flat()),
        ge = linear_interpolation(C_D.ge.x, C_D.ge.y, extrapolation_bc = Flat())
    )

    C_Y_interp = (
        δr = C_Y.δr, δa = C_Y.δa,
        β_δf = linear_interpolation((C_Y.β_δf.x_β, C_Y.β_δf.x_δf), C_Y.β_δf.y, extrapolation_bc = Flat()),
        p = linear_interpolation((C_Y.p.x_α, C_Y.p.x_δf), C_Y.p.y, extrapolation_bc = Flat()),
        r = linear_interpolation((C_Y.r.x_α, C_Y.r.x_δf), C_Y.r.y, extrapolation_bc = Flat())
    )

    C_L_interp = (
        δe = C_L.δe, q = C_L.q, α_dot = C_L.α_dot,
        α = linear_interpolation((C_L.α.x_α, C_L.α.x_stall), C_L.α.y, extrapolation_bc = Flat()),
        δf = linear_interpolation(C_L.δf.x, C_L.δf.y, extrapolation_bc = Flat()),
        ge = linear_interpolation(C_L.ge.x, C_L.ge.y, extrapolation_bc = Flat())
    )

    C_l_interp = (
        δa = C_l.δa, δr = C_l.δr, β = C_l.β, p = C_l.p,
        r = linear_interpolation((C_l.r.x_α, C_l.r.x_δf), C_l.r.y, extrapolation_bc = Flat())
    )

    C_m_interp = (
        z = C_m.zero, δe = C_m.δe, α = C_m.α, q = C_m.q, α_dot = C_m.α_dot,
        δf = linear_interpolation(C_m.δf.x, C_m.δf.y, extrapolation_bc = Flat())
    )

    C_n_interp = ( δr = C_n.δr, δa = C_n.δa, β = C_n.β, p = C_n.p, r = C_n.r
    )

    return (
        C_D = C_D_interp, C_Y = C_Y_interp, C_L = C_L_interp,
        C_l = C_l_interp, C_m = C_m_interp, C_n = C_n_interp
    )
end

# Example of how to use the visualization function
# aero_data = aero_update()
# visualize_interpolations(aero_data)