## Installation

```julia
using Pkg
Pkg.add("Flight")
```

!!! warning "Configuring Julia for Multithreading"

    During simulation, ```Flight.jl``` uses multithreading to concurrently handle the simulation loop, the built-in GUI and any attached I/O devices. Therefore, if you intend to use its interactive or I/O capabilities, you will need to start Julia with multiple threads enabled. You can find out how to do this in the [manual](https://docs.julialang.org/en/v1/manual/multi-threading/). However, if you are using the Julia extension for VS Code, the easiest way is to add the following entry to your
    ```settings.json```:

    ```json
    "julia.additionalArgs": [
        "--threads=auto",
    ],
    ```

    This should work well for most CPUs and use cases. If you run into issues, you can set the number of threads manually to cover your specific needs.

## Overview

`Flight.jl` offers a powerful and versatile framework for aircraft modeling, analysis and
simulation. Its design fully leverages Julia's expressiveness, extensibility and performance.

It is organized in three distinct layers:

- A lightweight, domain-agnostic engine for modeling and simulation of
  complex systems with hybrid dynamics (`FlightCore`).

- A library of high-fidelity, reusable physics and engineering models (`FlightLib`).

- A collection of specific aircraft implementations (`FlightAircraft`).

Key features:

*   **Hierarchical Modeling:** Enables building complex systems from simpler, reusable components,
    leveraging `ComponentArrays.jl` for clarity and convenience.

*   **High Performance:** Its core simulation loop is built on `DifferentialEquations.jl` and
    designed from the ground up to be allocation-free. This enables extremely fast headless execution and smooth performance on interactive runs.

*   **Interactive GUI:** Offers an extensible GUI based on `CImGui.jl` for live model
    inspection and manipulation.

*   **External Visualization & I/O:** Features out-of-the-box integration with [X-Plane
    12](https://www.x-plane.com/desktop/try-it/) for high-fidelity 3D visualization,
    joystick support via `SDL2_jll`, and a generic interface layer for custom I/O functionality.

*   **Solid Physics Foundation:** Includes built-in modules for attitude representation, geodesy, kinematics
    and rigid body dynamics, providing fast, accurate and ergonomic types and operations.

*   **Pre-Built Aircraft Components:** Comes with high-fidelity, efficient and customizable models for
    propellers, piston engines and landing gear.

*   **Integrated Control Design Workflow:** Provides general-purpose trimming and linearization
    functions. Seamlessly import linearized models into the `ControlSystems.jl` ecosystem for
    controller synthesis. Then, realize your design in a practical, discrete-time implementation, and
    validate it via nonlinear simulation.
