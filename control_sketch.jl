using ControlSystems
using Flight

function test(save = true)

    ac = Cessna172R(LTF()) |> System
    env = SimpleEnvironment() |> System
    trim_params = C172R.Trim.Parameters(
        loc = LatLon(),
        h = HOrth(1000),
        ψ_nb = 0.0,
        TAS = 50.0,
        γ_wOb_n = 0.0,
        ψ_lb_dot = 0.0,
        θ_lb_dot = 0.0,
        β_a = 0.0,
        fuel = 0.5,
        mixture = 0.5,
        flaps = 0.0)

    C172R.Trim.trim!(ac, env, trim_params)

    sim = Simulation(ac; args_ode = (env, ), t_end = 150, adaptive = true)
    Sim.run!(sim, verbose = true)
    # plots = make_plots(sim; Plotting.defaults...)
    plots = make_plots(TimeHistory(sim).kinematics; Plotting.defaults...)
    save ? save_plots(plots, save_folder = joinpath("tmp", "trim_sim_test")) : nothing

    #recreate the aircraft with NED kinematics, suitable for linearization
    ac = Cessna172R(NED()) |> System
    lm = Generic.StateSpaceModel(ac; env, trim_params) #retrims

    #the couplings of a specific state S into the time derivative of any other
    #state can be examined in the column of A corresponding to S. from this, we
    #can confirm that, for a wings-level, horizontal flight condition,
    #longitudinal states (v_x, v_z, θ, q, α_filt) are indeed uncoupled from
    #lateral-directional states (ψ, φ, β_filt, v_y, p, r)

    #other states, such as position (ϕ, λ, h) or fuel content are so weakly
    #coupled from both sets of states that can and should be ignored, in
    #order to avoid unnecessary complexity in the resulting system's transfer
    #functions, in which they would yield additional poles nearly cancelled by
    #close zeros.

    #extract a submodel retaining only inputs, states and outputs relevant to
    #longitudinal dynamics on a wings-level steady-state flight condition.
    #states that do not couple or couple very weakly into the longitudinal
    #dynamics are excluded. some of them (ψ, φ, β, v_y, p, r) have a place in a
    #lateral-directional dynamics submodel.
    long_dyn = filter(lm;
        u = (:elevator, :throttle),
        x = (:v_x, :v_z, :θ, :q, :α_filt, :ω_eng),
        y = (:TAS, :α, :θ, :q, :f_x, :f_z, :ω_eng))

    ss_long = ss(long_dyn)
    tf_long = tf(ss_long)

    # return lm
    # return long_dyn

    #let's get the transfer function from elevator input to q. we have chosen q
    #as the fourth output and elevator as the first input

    e2q = tf_long[4,1]
    @show poles(e2q)
    rlocusplot(e2q, 1)
    bodeplot(e2q)
    step(e2q)
    # nyquistplot(e2q)
    # nicholsplot(e2q)


end