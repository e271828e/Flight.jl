# Flight.jl

[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://e271828e.github.io/Flight.jl/dev/)
[![CI](https://github.com/e271828e/Flight.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/e271828e/Flight.jl/actions/workflows/CI.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


*A high-performance, extensible flight dynamics framework for Julia.*

![Flight.jl GUI with X-Plane 12 Visualization](docs/src/assets/github.png?raw=true)

## Documentation

Documentation is still in its infancy, but you can
check out the [tutorials](https://e271828e.github.io/Flight.jl/dev/tutorials/tutorial01/tutorial01/)
for a first glance at the package's capabilities.


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


## License

`Flight.jl` is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
