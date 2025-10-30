# Flight.jl

[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://e271828e.github.io/Flight.jl/dev/)
[![CI](https://github.com/e271828e/Flight.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/e271828e/Flight.jl/actions/workflows/CI.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


*A high-performance, extensible aircraft GNC framework for Julia.*

![Flight.jl GUI with X-Plane 12 Visualization](docs/src/assets/front.png?raw=true)

## Overview

`Flight.jl` offers a powerful and versatile toolkit for aircraft GNC modeling, design and simulation
tasks. It fully leverages Julia's expressiveness, extensibility and performance.

It is organized in three layers:

- A lightweight, domain-agnostic engine for causal modeling and simulation of
  complex systems with hybrid dynamics (`FlightCore`).

- A library of high-fidelity, reusable physics and engineering models (`FlightLib`).

- A collection of application examples (`FlightExamples`).

Key features:

*   **Hierarchical Modeling:** Enables building complex systems from simpler, reusable components,
    leveraging `ComponentArrays.jl` for clarity and convenience.

*   **High Performance:** Its core simulation loop is built on `DifferentialEquations.jl` and
    designed from the ground up to be allocation-free. This enables extremely fast headless
    execution and smooth performance on interactive runs.

*   **Interactive GUI:** Integrates an extensible GUI based on `CImGui.jl` for live model
    inspection and manipulation.

*   **External Visualization & I/O:** Offers out-of-the-box integration with [X-Plane
    12](https://www.x-plane.com/desktop/try-it/) for high-fidelity 3D visualization, joystick
    support via SDL2, and a generic interface layer for custom I/O functionality.

*   **Solid Physics Foundation:** Features built-in modules for attitude representation, geodesy,
    kinematics and rigid body dynamics, providing efficient and ergonomic types and operations.

*   **Pre-Built Aircraft Components:** Includes high-fidelity, customizable models for propellers,
    piston engines and landing gear.

*   **Case Study:** A custom fly-by-wire Cessna 172 model is provided as an example demonstrating a
    complete design workflow, from vehicle systems modeling to the implementation and testing of a
    complex, gain-scheduled autopilot.


## Installation

```julia
using Pkg
Pkg.add("Flight")
```

## Documentation

Documentation is still in its infancy. Please check out the [tutorials](https://e271828e.github.io/Flight.jl/dev/tutorials/tutorial01/tutorial01/)
for a first glance at the package's capabilities. If you're in a hurry, you can try the self-contained examples below.


## Examples
Automated turning climb under constant wind conditions:

```julia
using Flight

#1. Set up simulation

    #custom fly-by-wire Cessna172 variant with default environment
    world = SimpleWorld(; aircraft = Cessna172Xv2()) |> Model
    sim = Simulation(world; t_end = 600)

    #initialize using default trim conditions
    init!(sim, C172.TrimParameters())

#3. Set wind conditions

    world.atmosphere.wind.u.N = 1.0 #1 m/s North
    world.atmosphere.wind.u.E = 0.5 #0.5 m/s East

#3. Configure autopilot for turning climb

    #extract control laws submodel
    ctl = world.aircraft.avionics.ctl

    #set longitudinal control laws to track airspeed and climb rate
    ctl.lon.u.mode_req = C172XControl.ModeControlLon.EAS_clm

    #set lateral control laws to track bank angle and sideslip angle
    ctl.lat.u.mode_req = C172XControl.ModeControlLat.φ_β

    #set climb rate reference to 2.0 m/s, keep trim airspeed
    ctl.lon.u.clm_ref = 2.0

    #set bank angle reference to 30 degrees, keep trim sideslip angle
    ctl.lat.u.φ_ref = deg2rad(30)

#4. Run Simulation and extract results

    run!(sim)
    ts = TimeSeries(sim)

#5. Plot 3D trajectory

    kin_plots = make_plots(ts.aircraft.vehicle.kinematics; size = (900, 600))
    display(kin_plots[:Ob_t3d])
```
![Turning climb 3D trajectory](docs/src/assets/turning_climb_3d.png?raw=true)


Comparing elevator step response between nonlinear and linearized Cessna172S models:
```julia
using Flight
using ControlSystems, RobustAndOptimalControl, Plots, LaTeXStrings

#1. Set up and run nonlinear simulation

        #instantiate the aircraft with NED kinematics (required for linearization)
        world = SimpleWorld(aircraft = Cessna172Sv0(NED())) |> Model
        sim = Simulation(world; t_end = 10)

        #define trim conditions and initialize Simulation
        trim_params = C172.TrimParameters()
        init!(sim, trim_params)

        #advance 1 second from trim condition
        step!(sim, 1, true)

        #apply 10% elevator increment
        world.aircraft.vehicle.systems.act.u.elevator += 0.1

        #run to completion
        run!(sim)

        #extract pitch angle TimeSeries from simulation results
        ts = TimeSeries(sim)
        θ_nonlinear = ts.aircraft.vehicle.kinematics.e_nb.θ

    #2. Obtain linear SISO system

        #extract aircraft submodel and linearize it around the trim condition
        lss = linearize(world.aircraft, trim_params)

        #convert to NamedStateSpace
        nss = named_ss(lss)

        #extract elevator-to-pitch angle SISO system
        e2θ = nss[:θ, :elevator]

    #3. Compute linear response to elevator step input

        #simulate a 0.1 step input applied at t=1
        y, t, _, _ = lsim(e2θ, (x, t)->[0.1]*(t>=1), 0:0.01:10)

        #get perturbation Δθ around trim condition
        Δθ_linear = vec(y)

        #retrieve trim θ value from linearized aircraft model
        θ_trim = lss.y0[:θ]

        #compute total linear θ response
        θ_linear = θ_trim .+ Δθ_linear

    #4. Compare responses

        plot(θ_nonlinear; plot_title = "Pitch Angle", label = "Nonlinear", ylabel=L"$\theta \ (rad)$") |> display
        plot!(t, θ_linear; label = "Linear")
```
![Elevator step responses](docs/src/assets/elevator_step_response.png?raw=true)

## License

`Flight.jl` is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
