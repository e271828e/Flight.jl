# Modeling

## Overview

The `Modeling` module lies at the very heart of `Flight.jl`. It provides a lightweight, flexible,
and high-performance framework for building complex dynamical systems from simple, reusable
components. This framework was developed specifically for `Flight.jl` to meet the demanding
requirements of real-time flight simulation and advanced control system design.

`Flight.jl` employs a **causal** (or signal-flow) modeling paradigm. This means the flow of
computation is explicitly defined by the developer. Every `Model` must adhere to a standardized
interface of update functions that define its behavior:

While Julia's ecosystem offers the powerful `ModelingToolkit.jl` (MTK) for acausal modeling,
`Flight.jl` uses a custom, causal framework for several deliberate reasons that are critical to its
design goals:

1.  **Guaranteed, Hand-Tuned Performance:** The causal approach gives the developer absolute control
    over the computational path in the performance-critical `f_ode!` function. This allows for
    manual optimization to guarantee a type-stable, zero-allocation simulation loop, which is
    essential for real-time applications and large-scale analysis.

2.  **Intuitive for Control System Design:** The signal-flow paradigm is the natural mental model
    for control engineers. The explicit definition of inputs, outputs, and discrete update logic
    (`f_periodic!`) maps directly to the implementation of digital flight control and guidance
    systems.

3.  **Simplicity and Transparency:** The framework's abstractions are simple and direct. A `Model`
    is just a Julia `struct`, and its dynamics are just Julia functions. This makes the entire
    system highly transparent and easy to debug with standard Julia tools.

In essence, `Flight.jl`'s framework is not an attempt to reinvent a general-purpose modeling
language. It is a specialized, domain-specific tool forged for the express purpose of building the
most performant and flexible flight simulation systems possible in Julia.

#### The `ModelDefinition` / `Model` Paradigm

The framework is built upon two core abstractions:

1.  **`ModelDefinition` (The Blueprint):** An abstract type for user-defined `struct`s that serve as
    a declarative blueprint for a system component. A `ModelDefinition` specifies a component's
    constant parameters and its hierarchical structure by including other `ModelDefinition`s as
    fields.

    ```julia
    using Flight.FlightCore # Assuming this is where Modeling is exported

    struct MyActuator <: ModelDefinition
        τ::Float64 #A constant parameter (time constant)
    end

    struct MyAircraft <: ModelDefinition
        actuator::MyActuator #A sub-component
        mass::Float64 #Another constant parameter
    end
    ```

2.  **`Model` (The Instance):** The concrete, stateful object used by the simulation engine. When a
    `ModelDefinition` is passed to the `Model` constructor, it recursively builds a tree of `Model`
    instances and allocates the necessary memory for its state, inputs, and outputs.

    ```julia
    aircraft_def = MyAircraft(actuator = MyActuator(τ = 0.1), mass = 1000.0)
    aircraft_model = Model(aircraft_def)
    ```

Every `Model` instance in the simulation holds five key properties:

*   `x`: The **continuous state** vector, which is managed by the ODE solver.
*   `ẋ`: The **continuous state derivative**, which is computed by the model's dynamics.
*   `s`: The **discrete state**, a mutable `struct` for states that change at discrete time
    instances.
*   `u`: The **inputs**, a mutable `struct` for variables that are assigned externally.
*   `y`: The **outputs**, an immutable `struct` holding the results of the `Model`'s computations.


#### Hierarchical Composition with `ComponentArrays`

The framework's power comes from its hierarchical nature. The state vector `x` of a parent `Model`
is a `ComponentVector` that is recursively built from the state vectors of its children. This
provides two key benefits:

1.  **Human-Readable State:** You can access the state of any component with an intuitive dot
    notation (e.g., `aircraft.vehicle.systems.pwp.engine.x.ω`).
2.  **Performance:** The entire state vector remains a single, contiguous block of memory, which is
    essential for the performance of the ODE solver.

[**NEW**] The `Model`'s output `y` is also constructed hierarchically. After a component's update
function is called, its `y` field is updated. Parent models can then access the outputs of their
children via `child.y` to compute their own `y`. This bottom-up assembly ensures that all
information is current. While this means some output data may be present at multiple levels of the
hierarchy, the performance impact is negligible because the output `struct`s are required to be
`isbits` types (immutable and stack-allocated), making these operations extremely fast. It is
ultimately the developer's choice which child outputs to propagate upwards. However, since the
`Simulation` log only records the root `Model`'s output, any information that is not included in it
will not be saved.

#### The Update Interface: Causal, Programmatic Dynamics

*   `f_ode!(model, ...)`: Defines the continuous dynamics (`ẋ = f(x, u, t)`). This is called by the
    ODE solver at each integration step and is the performance-critical "hot loop" of the
    simulation.
*   `f_step!(model, ...)`: Defines logic that executes *after* each successful integration step. It
    is typically used for discrete state transitions, such as a "weight-on-wheels" flag changing
    from `false` to `true`.
*   `f_periodic!(model, ...)`: Defines logic that executes at fixed, discrete time intervals (`Δt`).
    This is the ideal place for implementing digital controllers, guidance laws, and other
    discrete-time algorithms.

[**EXPANDED**] These three core methods can be defined with any number of arguments beyond the
`Model` itself, allowing parent models to pass down necessary data (like atmospheric conditions or
kinematic data) to their children. However, a **root `Model`**—one that is passed directly to the
`Simulation` constructor—must define single-argument versions of these methods. If a root `Model`
needs to interact with the outside world, it should use the package's I/O capabilities (for
joysticks, network, etc.) or a user-defined callback, not method arguments.

[**NEW**] When extending the periodic update function, developers must define their method for the
`NoScheduling` dispatch type: `f_periodic!(::NoScheduling, mdl::Model, args...)`. The framework uses
this dispatch internally to distinguish between user-defined logic and the scheduling mechanism that
determines *if* the logic should run on a given step. Calling `f_periodic!(mdl, args...)` without
the `NoScheduling` dispatcher will correctly invoke the scheduling logic.

[**EXPANDED**] To enable efficient simulation of systems with components operating at different
frequencies (e.g., a high-frequency inner-loop controller and a low-frequency guidance system), the
framework provides the `Subsampled` wrapper. By wrapping a `ModelDefinition` in
`Subsampled(definition, K)`, you specify that its `f_periodic!` function should be called every
`K`-th period of its **parent's** periodic execution. The absolute sampling period of any component
is therefore the product of the root `Simulation`'s `Δt` and the `K` multipliers of all `Subsampled`
models in its branch of the hierarchy.


### Automatic Linkage for U and S

### Dispatching on Parametric `ModelDefinition`s