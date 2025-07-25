# Flight.jl

[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://e271828e.github.io/Flight.jl/dev/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

*A high-performance, extensible flight dynamics framework for Julia.*

## Overview

`Flight.jl` offers a powerful and versatile framework for aircraft modeling, analysis and
simulation. Its design fully leverages Julia's expressiveness, extensibility and performance.

It is architected in three distinct layers:

- A lightweight, domain-agnostic engine for modeling and interactive simulation of
  complex systems with hybrid dynamics (`FlightCore`).

- A library of high-fidelity, reusable physics and engineering models (`FlightLib`).

- A collection of specific aircraft implementations (`FlightAircraft`).

## Feature Summary

Here are some of `Flight.jl`'s highlights:

*   **Hierarchical Modeling:** The modeling framework leverages `ComponentArrays.jl` to construct
    complex systems from simpler, reusable components.

*   **High Performance:** The core simulation loop is built on `DifferentialEquations.jl`. It is
    designed from the ground up to be allocation-free, enabling extremely fast execution of headless
    simulations and smooth performance for interactive applications.

*   **Interactive GUI:** Offers an extensible GUI based on `CImGui.jl` for real-time inspection and
    manipulation during interactive simulation.

*   **External Visualization & I/O:** Features out-of-the-box integration with [X-Plane
    12](https://www.x-plane.com/desktop/try-it/) for high-fidelity 3D visualization,
    joystick support via `SDL2_jll`, and a generic interface layer for custom I/O functionality.

*   **Solid Physics Foundation:** Built-in modules for attitude representation, geodesy, kinematics
    and rigid body dynamics, providing fast, accurate, robust and ergonomic types and operations.

*   **Pre-Built Components:** High-fidelity, efficient and customizable models for propellers,
    piston engines and landing gear.

*   **Integrated Control Design Workflow:** Provides general-purpose trimming and linearization
    functions. Seamlessly import linearized models into the `ControlSystems.jl` ecosystem for
    controller synthesis. Then, realize your design in a practical, discrete-time implementation, and
    validate it via nonlinear simulation.


## Installation
`Flight.jl` is not registered yet, so it must be installed from its URL:

```julia
using Pkg
Pkg.add(url="https://github.com/e271828e/Flight.jl.git")
```

## Documentation

[Documentation](https://e271828e.github.io/Flight.jl/dev/) is still in its infancy, but you can
check out the [Interactive Simulation](https://e271828e.github.io/Flight.jl/dev/examples/ex01/ex01/)
showcase. Tutorials and usage examples are in the works.

![Flight.jl GUI with X-Plane 12 Visualization](docs/src/showcase/ex01/github.png?raw=true)




## License

`Flight.jl` is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
