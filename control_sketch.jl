using ControlSystems
using Flight

function test()
    ac = Cessna172R(NED()) |> System
    env = SimpleEnvironment() |> System
    trim_params = C172R.Trim.Parameters()
    lm = C172R.Linear.Model(; ac, env, trim_params)

    #create a submodel retaining only inputs, states and outputs relevant to
    #longitudinal dynamics
    long_dyn = lm
    long_dyn = C172R.Linear.filter_u(long_dyn, (:yoke_y, :throttle))
    long_dyn = C172R.Linear.filter_x(long_dyn, (:v_x, :v_z, :θ, :q, :α_filt, :ω_eng))
    long_dyn = C172R.Linear.filter_y(long_dyn, (:TAS, :α, :θ, :q, :f_x, :f_z, :ω_eng))

    ss_long = ss(long_dyn)
    tf_long = tf(ss_long)

    #let's get the transfer function from yoke_y (elevator control input) to q.
    #we have chosen q as the fourth output and yoke_y as the first input

    e2q = tf_long[4,1]
    poles(e2q)


end