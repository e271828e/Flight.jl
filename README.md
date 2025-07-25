# Flight.jl

[![Build Status](https://github.com/e271828e/Flight.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/e271828e/Flight.jl/actions/workflows/CI.yml)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://e271828e.github.io/Flight.jl/dev/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

*A high-performance, extensible framework for flight dynamics modeling, analysis, design and simulation.*

![Flight.jl GUI with X-Plane 12 Visualization](https://github.com/e271828e/Flight.jl/blob/main/docs/src/assets/images/cessna172x_showcase.png?raw=true)

## Overview

It provides a complete ecosystem for creating, analyzing and flying detailed aircraft models,
designed for high performance and interactivity.



The framework is architected in distinct layers: a domain-agnostic core simulation engine, a rich
library of physics and engineering models, and a top-level aircraft implementation layer. This
modular design allows for unparalleled flexibility and extensibility.

Flight.jl is a powerful tool for aerospace engineers, researchers, and enthusiasts who need a framework for rigorous, end-to-end simulation. It is built on modern software engineering principles and leverages the full performance potential of the Julia language.

A library of high-fidelity of physics and engineering models.

Integrates with OrdinaryDiffEq

## Core Features

*   **Hierarchical Compositional Modeling:** Build complex systems by assembling simpler, reusable components. The framework's core is a lightweight, domain-agnostic engine for simulating hybrid (continuous-discrete) dynamic systems. It allows you to define a complete model as a tree of nested, stateful submodels.
*   **High Performance:** The simulation engine is designed from the ground up to be allocation-free in its core update loop. This enables extremely fast execution speeds—often hundreds of times faster than real-time—and smooth, consistent performance for interactive applications.
*   **Integrated Control Design Workflow:** Go from a nonlinear aircraft model to a fully implemented flight controller. The framework provides tools to:
    1.  **Trim** the aircraft at any point in the flight envelope.
    2.  **Linearize** the model to obtain a linear state-space representation.
    3.  **Design** controllers using standard techniques (PID, LQR) and analyze them with the `ControlSystems.jl` ecosystem.
    4.  **Implement** the controllers with gain scheduling directly within the simulation.
*   **Interactive GUI:** A built-in GUI based on `CImGui.jl` allows for live inspection and
    modification of simulated model variables during interactive runs.
*   **External Visualization & I/O:** Includes out-of-the-box support for the [X-Plane 12 Flight Simulator](https://www.x-plane.com/desktop/try-it/) as a high-fidelity 3D visualization tool, as well as a generic interface for connecting hardware like joysticks.
*   **Extensible Library:** Comes with a rich library of pre-built components, including detailed models for aerodynamics, piston engines, propellers, landing gear, atmospheric conditions, and terrain.

## Quick Start: Interactive Simulation

Get a feel for `Flight.jl` by running an interactive simulation of the Cessna 172X, a custom fly-by-wire variant of the classic aircraft.

!!! warning "Multithreading Required"
    Interactive simulations use multithreading to handle the GUI and I/O devices concurrently with the simulation loop. You will need to start Julia with multiple threads enabled. In VS Code, you can do this by adding `"julia.additionalArgs": ["--threads=auto"]` to your `settings.json`.

1.  **Add the Package:**
    ```julia
    using Pkg
    Pkg.add("Flight")
    ```

2.  **Run the Simulation:**
    ```julia
    using Flight
    using Sockets #for IPv4

    #Create a world with a C172X aircraft, atmosphere, and terrain
    world = SimpleWorld(Cessna172Xv1(), SimpleAtmosphere(), HorizontalTerrain())

    #Set up the simulation with a 1/60s integration step
    sim = Simulation(world; dt = 1/60, t_end = 3600)

    #Define an initial condition (trimmed, in-air)
    initializer = C172.TrimParameters(EAS = 55, h = 1000)
    Sim.init!(sim, initializer)

    #Attach an X-Plane 12 interface (optional)
    xp = XPlane12Control(address = IPv4("127.0.0.1"), port = 49000)
    Sim.attach!(sim, xp)

    #Attach any connected joysticks
    for joystick in update_connected_joysticks()
        Sim.attach!(sim, joystick)
    end

    #Run the interactive simulation (GUI runs on the main thread)
    Sim.run_interactive!(sim)
    ```

This will launch the simulation and open the GUI, allowing you to control the aircraft and monitor its state.

## Documentation

For a complete guide, including tutorials, usage examples, and the API reference, please see our [documentation](https://e271828e.github.io/Flight.jl/dev/).

## Contributing

We welcome contributions of all kinds! If you've found a bug, have a feature request, or would like to contribute code, please feel free to open an issue or a pull request on our [GitHub repository](https://github.com/e271828e/Flight.jl).

## License

`Flight.jl` is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

### Justification for the Latest Changes:

*   **Image Placement:** The screenshot is now the first major element. This immediately establishes the package's capability and visual appeal. I've also used a direct GitHub raw link for the image, which is a common and robust way to embed images in a `README.md`. You'll need to make sure that image path is correct in your repository.
*   **"General Core" as a Feature:** The "A General Modeling Core" section has been merged into the first feature bullet point:
    > **Hierarchical Compositional Modeling:** Build complex systems by composing simpler, reusable components. The framework's core is a lightweight, domain-agnostic engine for simulating hybrid (continuous-discrete) dynamic systems. It allows you to define a complete model as a tree of nested, stateful subsystems.
    This is a much more integrated and powerful way to present it. It starts with the benefit ("Build complex systems...") and then explains the underlying technology ("The framework's core is..."). It also naturally includes the key descriptive terms we discussed: "lightweight," "domain-agnostic," "compositional," and "stateful."
*   **Minor Wording Tweak:** I changed "External Visualization" to "External Visualization & I/O" to explicitly mention the joystick support, which is a key part of the interactive experience.

This version feels very strong. It's visually appealing, leads with its most exciting features, accurately describes its capabilities without overstating them, and provides a clear path for a new user to get started. It strikes an excellent balance between being a showcase for a flight simulator and an introduction to a powerful underlying framework.