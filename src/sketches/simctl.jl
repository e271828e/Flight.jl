using DifferentialEquations

# problema: sleep hace basicamente lo que le sale de los cojones. tiene una
# resolucion de mierda. por mucho que se pueda comandar en incrementos de 1ms,
# en realidad los quanta son de 0.015. por tanto, para hacer la simulacion real
# time, en vez de marcar a priori el dtsim y despues dormir el tiempo restante
# entre tstops, la solucion puede ser la contraria: intentar dormir un valor
# nominal, que no se cumplira. a continuacion ver cuanto tiempo ha dormido
# realmente, y avanzar la simulacion ese mismo dt. aqui sleep, por su
# impredecibilidad, hace el papel del render (ver GameProgrammingPatterns)

#we sleep before, then update simulation. therefore, we need to sleep at the
#first call. but we cannot set tsim_callback = 0 so that the affect! function is
#called at tsim = 0, because the integrator is already at tsim = 0 and
#understands that the dt it needs to take to hit the next callback is dt =
#tsim_callback - tsim_current = 0. to bypass this we can set initial_affect =
#true, and then in that very first call to affect!, we can already compute
#tsim_callback after t=0. drawback: this is called as soon as the integrator is
#created, so we need to use solve(), not init() and then solve!(integrator)

#an alternative is to initialize tsim_callback = dtwall_ref, set initial_affect
#false, and then allow the integrator to take that first leap to dtwall_ref
#without sleeping first. in practice, this shouldn't be noticeable.

#CAUTION: if we choose an inaccurate method (such as Heun) but we don't set a
#fixed dt, it may happen that the tstops chosen by the integrator are so close
#together than when we hit the tstop for one of our time_ctl_callbacks and we go
#to sleep, we sleep for longer than the dtsim to the next tstop chosen by the
#integrator. then we have overflow.

#when we solve without real time constraints, we can leave the adaptive time
#step on, so that the only tstops are those chosen by the integrator. in this
#scenario Tsit5 is much, much faster than Heun for the same absolute and
#relative tolerance, because although the method is more expensive, the tstops
#are more spaced out. however, when we have real time constraints, and therefore
#use the iterative callback, the tstops set by the callback (which in turn are
#the outcome of the sleep() calls) will be much closer together, so the Tsit5
#algorithm doesn't have a chance to "shine".

Base.@kwdef mutable struct TimeControl
    dtwall_ref::Float64 = 0.1
    tsim_callback::Union{Nothing, Float64} = 0
    twall_last::Union{Nothing, Float64} = nothing
end

Base.@kwdef struct ModelInputs
    f_forcing::Function = t->0
end

Base.@kwdef struct MyParams
    inputs::ModelInputs = ModelInputs()
    time_ctl::TimeControl = TimeControl()
end

function make_problem()

    #equations
    l = 1.0
    m = 1.0
    g = 9.81

    function pendulum!(du,u,p,t)
        du[1] = u[2] #dθ = ω
        du[2] = -3g/(2l) * sin(u[1]) + 3/(m * l^2) * p.inputs.f_forcing(t)
    end

    θ₀ = 0.01
    ω₀ = 0.0
    u₀ = [θ₀, ω₀]
    tspan = (0.0, 1.0)

    params = MyParams()

    function affect!(integrator)
        println("Simulation time: t = $(integrator.t)")
        time_ctl = integrator.p.time_ctl
        if time_ctl.twall_last === nothing
            time_ctl.twall_last = time()
        end
        println("Going to sleep for $(time_ctl.dtwall_ref)")
        sleep(time_ctl.dtwall_ref)
        twall_current = time()
        dtwall_actual = twall_current - time_ctl.twall_last
        println("Slept for $dtwall_actual")
        time_ctl.twall_last = twall_current
        time_ctl.tsim_callback += dtwall_actual #tsim for next callback
        println("Next callback: tsim = $(time_ctl.tsim_callback)")
        u_modified!(integrator, false)
        # now the callback returns and the integrator steps to the next
        # tsim_callback, taking the appropriate intermediate steps
    end

    function tsim_next(integrator)
        #this adds a new tstop in which affect! will be called
        integrator.p.time_ctl.tsim_callback
    end

    time_ctl_callback = IterativeCallback(tsim_next, affect!, initial_affect = true, save_positions = (false, false))
    ODEProblem(pendulum!, u₀, tspan, params, callback = time_ctl_callback)
    #then call:
    # solve(make_problem(), Heun(), adaptive = true)
    #... letting the integrator choose the intermediate integration steps


    # ODEProblem(pendulum!, u₀, tspan, params)

end


########## ALTERNATIVE IMPLEMENTATION ###############

#we use a callback that fires on every integration step (IterativeCallback with
#initial_affect = true. on that callback, we go to sleep for dtwall_ref. when we
#wake-up, we compute dtwall_actual, save twall_current and
#set_proposed_dt!(dtwall_actual). the time_choice(integrator) function is simply
#integrator.t + get_proposed_dt(integrator). in this approach, the integrator
#does not take any intermediate steps between callback times. in fact, we can
#set adaptive = false and the initial dt to an arbitrary value (doesn't matter),
#because the iterative callback, which sets the dt for the next step, will be
#called before the first step (since initial_affect = true)

function make_problem2()

    #equations
    l = 1.0
    m = 1.0
    g = 9.81

    function pendulum!(du,u,p,t)
        du[1] = u[2] #dθ = ω
        du[2] = -3g/(2l) * sin(u[1]) + 3/(m * l^2) * p.inputs.f_forcing(t)
    end

    θ₀ = 0.01
    ω₀ = 0.0
    u₀ = [θ₀, ω₀]
    tspan = (0.0, 10.0)

    params = MyParams()

    function affect!(integrator)
        println("Simulation time: t = $(integrator.t)")
        time_ctl = integrator.p.time_ctl
        if time_ctl.twall_last === nothing
            time_ctl.twall_last = time()
        end
        println("Going to sleep for $(time_ctl.dtwall_ref)")
        sleep(time_ctl.dtwall_ref)
        twall_current = time()
        dtwall_actual = twall_current - time_ctl.twall_last
        println("Slept for $dtwall_actual")
        time_ctl.twall_last = twall_current
        set_proposed_dt!(integrator, dtwall_actual)
        time_ctl.tsim_callback = integrator.t + dtwall_actual
        println("Next tsim = $(time_ctl.tsim_callback)")
        u_modified!(integrator, false)
    end

    function tsim_next(integrator)
        #this adds a new tstop in which affect! will be called
        integrator.p.time_ctl.tsim_callback
    end

    time_ctl_callback = IterativeCallback(tsim_next, affect!, initial_affect = true, save_positions = (false, false))
    ODEProblem(pendulum!, u₀, tspan, params, callback = time_ctl_callback)

    #now do:
    # solve(make_problem2(), Heun(), adaptive = false, dt = 0.02)
    # we need to disable adaptive, otherwise the proposed dt set in affect! is
    # overwritten later. here, dt doesn't do anything, because it gets
    # overwritten on each callback

end