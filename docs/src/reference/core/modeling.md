# Modeling

## Overview

### Philosophy: Performance, Control, and Transparency

The `Modeling` module is the cornerstone of `Flight.jl`. It provides a lightweight, flexible and
high-performance framework, which enables creating complex systems from simpler,
reusable components in a hierarchical fashion.

This modeling framework is **causal** and **imperative**.
It is domain-agnostic by design, but it is tailored to / heavily oriented towards its primary goals:
GNC algorithm development and real-time simulation.

It can be understood as a generalization of the classic state-space representation used in
control engineering and other disciplines. Shares the same input, state, output, but generalized to
accomodate nonlinear and discrete dynamics
state-space oriented, which is a good fit for most tasks within the GNC domain

In essence, modeling.jl is not a symbolic modeling system; it is a framework for organizing
imperative, state-space code into a composable hierarchy.

This framework was developed specifically to meet the requirements of guidance,
navigation & control applications.
Its minimalistic design is a deliberate choice prioritizing **ultimate performance, fine-grained
control, and architectural transparency**. absolute control

The user explicitly:
- Specifies the properties of each component: input, output, continuous state, etc.
- Writes the Julia code within the update functions that define the system's dynamics.

This approach places a greater responsibility on the user to achieve performant, efficient code. In
return, it grants absolute control and predictability, which are desirable for real-time simulation
and deep systems analysis.

This is different from the **acausal**, symbolic and declarative paradigm of modeling frameworks
like `ModelingToolkit.jl`.

Its capabilities are less general

For the specific domain of real-time, high-fidelity aircraft GNC simulation, where the underlying
physics are well-described by causal ODEs and performance is paramount, this is an excellent and
highly justifiable engineering choice

favors

Do focus on what ModelingToolkit.jl is or is not. Focus on describing what Flight.jl does.

While Julia's ecosystem offers the powerful `ModelingToolkit.jl` (MTK) for acausal modeling,
`Flight.jl` uses a custom, causal framework for several deliberate reasons that are critical to its
design goals:

1.  **Guaranteed, Hand-Tuned Performance:** The causal approach gives the developer absolute control
    over the computational path in the performance-critical model update functions. This allows for
    manual optimization to guarantee a type-stable, zero-allocation simulation loop, which is
    essential for real-time applications and large-scale analysis.

2.  **Intuitive for GNC Algorithm Design:** The signal-flow paradigm is the natural mental model
    for control engineers. The explicit definition of inputs, outputs, and discrete update logic
    (`f_periodic!`) maps directly to the implementation of digital flight control and guidance
    systems.

3.  **Simplicity and Transparency:** The framework's abstractions are simple and direct. A `Model`
    is just a Julia `struct`, and its dynamics are just Julia functions. This makes the entire
    system highly transparent and easy to debug with standard Julia tools.


In summary, `Flight.jl`'s modeling framework is a specialized tool designed for a specific job. It
trades the generality and automatic code generation of symbolic systems for the raw performance,
transparency, and control of a carefully crafted imperative structure.

State-space

Flight dynamics modeling tasks fit very well

This guide explains the core concepts of the framework and the "performance contract" you implicitly
agree to when creating new components.

`Flight.jl` employs a **causal** (or signal-flow) modeling paradigm. This means the flow of
computation is explicitly defined by the developer. Every `Model` must adhere to a standardized
interface of update functions that define its behavior:


In essence, `Flight.jl`'s framework is not an attempt to reinvent a general-purpose modeling
language. It is a specialized, domain-specific tool forged for the express purpose of building the
most performant and flexible flight simulation systems possible in Julia.

### Core Concepts: Blueprints and Instances

### The `ModelDefinition` / `Model` Paradigm

The framework is built upon two core abstractions:

#### 1. `ModelDefinition`: The Static Blueprint

A `ModelDefinition` is a `struct` that acts as a **blueprint** for a simulation component. Its
purpose is to:
1.  Define the static parameters and constants of the component.
2.  Establish the system's hierarchy by including other `ModelDefinition`s as fields.

    For example, a `Vehicle` is *defined* as having `Systems`, `Kinematics`, and `Dynamics` components:

    ```julia
    # A simplified conceptual example
    struct Vehicle <: ModelDefinition
        systems::Systems   #<-- Child ModelDefinition
        kinematics::WA     #<-- Child ModelDefinition
        dynamics::VehicleDynamics
    end
    ```

    `ModelDefinition`s are lightweight, compile-time constructs that describe the structure of your
    simulation.


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


A `Model` is the **stateful instance** of a component used during the simulation. It is created by
passing a `ModelDefinition` to the `Model` constructor: `mdl = Model(my_definition)`. The `Model`
struct holds the "live" data:
Every `Model` instance in the simulation holds five key properties:

*   `x`: The **continuous state** vector, which is managed by the ODE solver.
*   `ẋ`: The **continuous state derivative**, which is computed by the model's dynamics.
*   `s`: The **discrete state**, a mutable `struct` for states that change at discrete time
    instances.
*   `u`: The **input**, a mutable data type for variables that are assigned externally.
*   `y`: The **output**, an immutable data type holding the results of the `Model`'s computations.



When a `Model` is created, the framework recursively constructs `Model` instances for all its
children, creating a tree of interconnected objects. `ComponentArrays` are used under the hood to
ensure that the state vector `x` of a parent `Model` is a single, contiguous block of memory
composed of the state vectors of its children. This is so that the ODE solver sees the
continuous state vector and its derivative as a single, contiguous blocks of memory.

## How to define a Model's X, Xdot, S, U, Y

The `ModelDefinition` struct itself contains fields that are either `ModelDefinition` subtypes or
other types. These end up as `parameters` and `submodels` respectively. But, how do you specify the
`Model`'s `x`, `xdot`, etc?

By default, all of these `Model` properties will be `nothing`. The only exception is x and x_dot for a node (non-leaf) `Model`. To specify what you want them to be, you
define constructor methods for each of these `ModelDescriptor`s that dispatch
on your `ModelDefinition` subtype. This sounds more complicated than it is. An example:
Let's say we want

The only exception is x and x_dot for a node (non-leaf) Model. You typically should not define
Modeling.X for such models. Their x and x dot will be automatically assembled from those of its
children, and they will share memory with them.

```
struct MyDef <: ModelDefinition end

Modeling.X(::MyDef) = zeros(3) #must be an AbstractVector
Modeling.Y(::MyDef) = zeros(SVector{3}) #must be an isbits type
Modeling.U(::MyDef) = Ref(0.0) #can be any mutable type

mdl = MyDef() |> Model
@show mdl.x
@show mdl.y
@show mdl.u
@show mdl.s

Modeling.X(::MyDef) = ComponentVector(a = 0.0, b = zeros(2))
```


Modeling.X: for a leaf component, it must return an AbstractVector. If this method is not defined,
the leaf component will have no continuous state by default. For a node component, you
typically don't need to extend this method; its continuous state vector will be automatically
assembled from those of its children, and it will share memory with them. Children

Modeling
Very simple example of composition. Show that

Question: couldn't the model's `parameters` and `submodels` also be specified via method extension?
Given that the `Model` constructor extracts the fields and processes it depending on what they are,
why not simply leave the `ModelDefinition` as an empty struct and define the following functions:
`Modeling.Constants(::MyDef)`
`Modeling.Submodels(::MyDef)`

The answer is that this would prevent the user from creating different instances of a `Model`, each
one with its own constant parameter values or submodels. The high-level explanation comes from
asking the question: what characterizes a specific instance of a `Model`? Its specific constant
parameters, and the specific instances of the submodels it holds. Therefore, these must be left to
be specified when creating each instance, so they should not be universally set by method
extensions.

In contrast, for example, we know that for a leaf model, its x will always have the same structure.
And for a node model, it will be assembled from its children's x.

---



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

### The Update Cycle: The API Contract

Every component you create must adhere to a simple but strict API contract by extending a set of
standard functions. The simulation engine guarantees when and how these functions are called.

-   `f_ode!(model, ...)` This is the heart of the continuous simulation. It is called multiple times
    per integration step by the ODE solver. Its single purpose is to compute the state derivative
    `model.ẋ` based on the current state `model.x`, inputs `model.u`, and time `model.t`. **This is
    the performance-critical hot loop of your simulation.**

-   `f_periodic!(model, ...)` This function is called by the simulation at fixed time intervals
    (`Δt`). It is used to model discrete-time dynamics. This is the correct place for:
    *   Digital control logic (PIDs, LQRs).
    *   State machines (e.g., `engine_state` transitions).
    *   Updating gain schedules.

-   `f_step!(model, ...)` This function is called once at the end of every successful ODE
    integration step. It is ideal for tasks that must happen after a state update but are not part
    of the continuous dynamics, such as re-normalizing a quaternion to prevent numerical drift.

-   `init!(model, ...)` This is called to set the initial conditions for the model's states (`x` and
    `s`) and inputs (`u`) before a simulation run begins.

The `@sm_updates` macro can be used on a parent `ModelDefinition` to automatically generate
recursive update methods that call the corresponding methods on all its submodels.

---

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


### Dispatching on Parametric `ModelDefinition`s

A potential source of perplexity and frustration occurs when working with parametric
ModelDefinitions. If MyModelDef is parametric, Julia will NOT dispatch on a method
f(::Model{MyModelDef}). One needs to define f(::Model{<:MyModelDef}) instead. This can be a source
of perplexity and frustration to anyone not too familiar with the Julia type system.

### Parent-Child Linkage Options for U and S
```
using Flight
import .Modeling: U

struct Leaf <: ModelDefinition end

Modeling.U(::Leaf) = [0.0]

@kwdef struct Node1 <: ModelDefinition
a::Leaf = Leaf()
b::Leaf = Leaf()
end

@kwdef struct Node2 <: ModelDefinition
a::Leaf = Leaf()
b::Leaf = Leaf()
end

U(node::Node1) = (a = U(node.a), b = U(node.b))
U(node::Node2) = ComponentVector(a = U(node.a), b = U(node.b))

mdl1 = Node1() |> Model
mdl2 = Node2() |> Model

@assert mdl1.u isa NamedTuple
@assert mdl1.u.a === mdl1.a.u

@assert mdl2.u isa ComponentVector
@assert mdl2.u.a === mdl2.a.u #Component mdl.u returns a view into its a block
```

What's happening is that upon construction, `mdl.a` input `u` is assigned a view of `mdl.u`'s `a`
block.



### The Performance Contract: Your Responsibility

*If you care about performance!*

The exceptional performance of `Flight.jl` is not magic; it comes from a disciplined implementation
approach that every developer must follow for their own components. By extending the framework, you
agree to this "performance contract":

!!! warning "The Performance Contract" To ensure the high-performance, real-time capabilities of the
    simulation, your component implementations must adhere to the following rules.

1.  **Be Type-Stable:** All functions in the update cycle, especially `f_ode!`, must be type-stable.
    The types of variables should not change within the function. Use tools like `@code_warntype` to
    verify.
2.  **`f_ode!` Must Not Allocate:** The continuous dynamics function is the hot loop. It **must not
    allocate any memory on the heap**.
    *   Use `StaticArrays` (`SVector`, `SMatrix`) for all small, fixed-size vector and matrix
        operations.
    *   Do not create standard `Vector`s, `Matrix`es, or other heap-allocated objects inside
        `f_ode!`.
    *   If you need a temporary buffer, pre-allocate it in your `ModelDefinition` and store it in
        `model.parameters`.
3.  **Use the Right Data Structures:**
    *   Define your outputs (`Y`) as immutable (`isbits`) `struct`s. This avoids allocations when
        returning outputs.
    *   Use `ComponentArrays` for your state vector (`X`) to combine performance with readability.
4.  **Test Your Performance:** Use the `@ballocated` macro from `BenchmarkTools.jl` in your unit
    tests to verify that your `f_ode!` implementation is allocation-free. The `Flight.jl` source
    code provides numerous examples of this.

Adhering to this contract ensures that your custom components integrate seamlessly into the
framework without degrading the performance of the entire simulation.

### Detailed explanation for Model state vector assembly

Excellent question. This is arguably the most elegant and crucial part of the entire framework, and it's powered by the magic of `ComponentArrays.jl`. Let's walk through the process step-by-step.

The core idea is a **two-pass process**:
1.  **First Pass (Assembly):** Recursively walk the `ModelDefinition` hierarchy from the bottom up to determine the *shape* and *size* of the complete state vector for the root model. This creates a single, contiguous `ComponentVector`.
2.  **Second Pass (Assignment):** Recursively walk the hierarchy from the top down to build the `Model` instances. Each `Model` instance is given a *view* (or more accurately, a property accessor that behaves like a view) into the relevant slice of the single `ComponentVector` created in the first pass.

Let's illustrate with a simple example hierarchy:

```julia
# 1. Define the blueprints (ModelDefinitions)

# The simplest component, a "leaf" in our hierarchy
struct Leaf <: ModelDefinition end
# A leaf's state is a simple Vector. This is the base case for our recursion.
Modeling.X(::Leaf) = [0.0]

# A composite component, a "node"
struct Node <: ModelDefinition
    leaf_C::Leaf = Leaf()
    leaf_D::Leaf = Leaf()
end

# The top-level component, the "root"
struct Root <: ModelDefinition
    leaf_A::Leaf = Leaf()
    node_B::Node = Node()
end
```
Our hierarchy looks like this:
```
Root
├── leaf_A (Leaf)
└── node_B (Node)
    ├── leaf_C (Leaf)
    └── leaf_D (Leaf)
```
Now, let's trace what happens when we call `Model(Root())`.

---

### Pass 1: Assembling the State Vector (Bottom-Up)

The entry point is `Model(Root())`. The very first thing this constructor needs to know is the complete state vector for a `Root` model. It finds this by calling `X(Root())`.

**Step 1.1: `X(Root())` is called.**

The default implementation of `X` for a composite `ModelDefinition` does the following:
*   It finds all fields that are themselves `ModelDefinition`s (here, `:leaf_A` and `:node_B`).
*   It recursively calls `X()` on each of them.
*   It assembles the results into a `ComponentVector`.

So, `X(Root())` triggers two calls:
1.  `X(mdl.leaf_A)` which is `X(Leaf())`
2.  `X(mdl.node_B)` which is `X(Node())`

**Step 1.2: Base Case - `X(Leaf())`**

This is the **base case** of our recursion. We have explicitly defined `Modeling.X(::Leaf) = [0.0]`. So, this call simply returns `[0.0]`. The recursion for this branch stops.

**Step 1.3: Recursive Call - `X(Node())`**

Just like `X(Root())`, the `X(Node())` method finds its own children `ModelDefinition`s (`:leaf_C` and `:leaf_D`) and calls `X()` on them.
*   `X(mdl.leaf_C)` calls `X(Leaf())` and gets `[0.0]`.
*   `X(mdl.leaf_D)` calls `X(Leaf())` and gets `[0.0]`.

Now, `X(Node())` assembles these results into a `ComponentVector`, using the field names as keys:
`ComponentVector(leaf_C = [0.0], leaf_D = [0.0])`. This object has a named structure but stores its data in a flat, contiguous vector: `[0.0, 0.0]`.

**Step 1.4: Final Assembly in `X(Root())`**

The execution now returns to the top-level call, `X(Root())`. It has received the results from its children:
*   From `leaf_A`: `[0.0]`
*   From `node_B`: `ComponentVector(leaf_C = [0.0], leaf_D = [0.0])`

It assembles these into the final, complete `ComponentVector` for the entire hierarchy:
`x_root = ComponentVector(leaf_A = [0.0], node_B = (leaf_C = [0.0], leaf_D = [0.0]))`

**Crucially, `x_root` is a single object whose data is stored in a single, contiguous block of memory:** `[0.0, 0.0, 0.0]`. The nested naming (`x_root.node_B.leaf_C`) is just a convenient and efficient way to access slices of this underlying flat vector. This answers your "contiguous in memory" question. The ODE solver will only ever see this single, flat `Vector{Float64}`.

---

### Pass 2: Assigning Views to `Model` Instances (Top-Down)

Now that the master state vector `x_root` is created, the `Model` constructor can create the instances.

**Step 2.1: The Root `Model` is Created**

A `Model{Root}` instance is created. Its `x` field is assigned the complete `x_root` `ComponentVector` we just built.

**Step 2.2: The Children `Model`s are Created**

The `Model(Root())` constructor now iterates through its children `ModelDefinition`s (`leaf_A` and `node_B`) to create their `Model` instances.

*   **For `leaf_A`:** It creates a `Model{Leaf}`. For the `x` field of this child model, it does not create a new vector. Instead, it passes `x_root.leaf_A`. This is the key step. `x_root.leaf_A` is not a copy; it's a "view" or a "property accessor" into the first element of `x_root`'s underlying data. The `leaf_A` model now holds a direct reference to its slice of the master state vector.

*   **For `node_B`:** It creates a `Model{Node}`. For its `x` field, it passes `x_root.node_B`. This is a view into the last two elements of `x_root`'s data.

**Step 2.3: The Grandchildren `Model`s are Created**

The process continues recursively. When the `Model{Node}` for `node_B` is being constructed, it in turn looks at its own children:

*   **For `leaf_C`:** It creates a `Model{Leaf}`. It receives the view `x_root.node_B` from its parent. It then passes `(x_root.node_B).leaf_C` to its child `leaf_C`. This is a view into the second element of the master `x_root` vector.

*   **For `leaf_D`:** It does the same, passing `(x_root.node_B).leaf_D` to its child `leaf_D`, which is a view into the third element of `x_root`.

### The Final Result

After the constructor finishes, we have a tree of `Model` objects.

*   The `root_model` holds the **one and only** complete, contiguous state vector: `x_root`.
*   Every child, grandchild, and so on, holds an `x` field that is simply a **lightweight view** pointing to its designated section within that single master vector.

When the ODE solver updates the state, it modifies the elements of `x_root`. Because all submodels hold views into this same block of memory, they all instantly see the updated state without any copying or communication required.

Conversely, when a leaf model like `leaf_C` computes its derivative in `f_ode!` and writes to its `ẋ` field, it is actually modifying the specific slice of the root model's `ẋ` vector that corresponds to its state.

This design is brilliant because it combines:
1.  **Readability:** You can access any state from the top down with a clean, hierarchical syntax (`root_model.node_B.leaf_C.x`).
2.  **Performance:** There is only one block of memory for the entire state, making it extremely friendly to the ODE solver and the CPU cache.
3.  **Composability:** The logic for defining the shape of the state (`X(::MyDef)`) is local to each component, but the framework automatically assembles it into a global, performant structure.
