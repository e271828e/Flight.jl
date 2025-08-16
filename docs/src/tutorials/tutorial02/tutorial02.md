# Headless Simulation

In this tutorial, we revisit our [interactive simulation](@ref "Interactive Simulation") setup. This
time, we will learn how to run the `Simulation` programmatically and extract results for inspection
and plotting.

### Peeking Into the `Model`

Let's start from a new `SimpleWorld` instance. Here, keeping consistency with X-Plane's visuals
is no longer a concern, so we can stick to the default constructors:
```@example tutorial02
using Flight
world = SimpleWorld(Cessna172Xv1(), SimpleAtmosphere(), HorizontalTerrain())
```

Inspecting its type hierarchy reveals that `SimpleWorld` is a concrete subtype of the abstract type
`ModelDefinition`:
```@example tutorial02
using InteractiveUtils
supertypes(SimpleWorld)
```

An instance of a `ModelDefinition` subtype can be thought of as a blueprint specifying how
a particular `Model` should be built. To instantiate this `Model`, we pass the `ModelDefinition`
object to the `Model` constructor:
```@example tutorial02
mdl = Model(world)
```

Let's take a look at this `Model`'s properties:
```@example tutorial02
propertynames(mdl)
```

Here is a brief description of each one:
- `x` (*continuous state*): A vector containing `Model` states that evolve continuously over time.
  It is updated by the `Simulation`'s ODE integrator.
- `ẋ` (*continuous state derivative*): A vector containing the time derivative of `x`.
- `s` (*discrete state*): A mutable type holding `Model` states that change only at specific
  moments in time.
- `u` (*input*): A mutable type holding the `Model`'s input variables. These are assigned either by
  a parent `Model` or an external source.
- `y` (*output*): An immutable type holding all the relevant results from the `Model`'s
  computations. These may be used by a parent `Model` or saved in the `Simulation`'s log.
- `t` (*time*): A reference to the current simulation time.
- `Δt` (*sampling period*): The time interval for the `Model`'s periodic updates.
- `constants`: A `NamedTuple` containing `Model` parameters that do not change
  during simulation.
- `submodels`: A `NamedTuple` of child `Model`s that represent components of the parent `Model`.

The remaining properties are entries from `constants` and `submodels`. These may be accessed
directly via dot notation:
```@repl tutorial02
keys(mdl.constants)
keys(mdl.submodels)
mdl.atmosphere
```

A submodel can have submodels of its own:
```@repl tutorial02
keys(mdl.atmosphere.submodels)
mdl.atmosphere.sl
```

Thus, a complex `Model` can be made up of multiple, hierarchically arranged components, each one of
them itself a `Model`. To visualize a `Model`'s hierarchy, you can use the `print_tree` function.
For example, let's see what's underneath the `mdl.aircraft.vehicle.ctl` node:
```@example tutorial02
using AbstractTrees
print_tree(mdl.aircraft.avionics.ctl; maxdepth = 10)
```

You can also inspect a specific property across a `Model`'s hierarchy. For instance, to view the
discrete state of `mdl.aircraft.avionics.ctl` and every node underneath:
```@example tutorial02
print_tree(mdl.aircraft.avionics.ctl, :s; maxdepth = 10);
```

### Simulating an Elevator Doublet

Our plan is as follows:
1. Create a `Simulation`.
2. Initialize the aircraft in a trimmed state.
3. Apply an elevator doublet by setting the appropriate inputs and stepping the `Simulation` manually.
4. Let the `Simulation` run to completion.
5. Extract and plot some logged variables to observe the aircraft's response.

Let's begin by creating a `Simulation` from our `Model`:

```@example tutorial02
sim = Simulation(mdl; dt = 0.02, t_end = 60)
```

!!! tip "From a ModelDefinition"
    You can also pass a `ModelDefinition` object directly to the `Simulation` constructor. In that
    case, the `Model` is instantiated automatically under the hood:
    ```@example tutorial02
    Simulation(world; dt = 0.02, t_end = 60)
    nothing #hide
    ```

The `Model` is now stored within the `Simulation`, and it can be accesed at any time for inspection
or manipulation:

```@example tutorial02
@assert mdl === sim.mdl
```

Next, let's define a trim condition and use it to initialize the `Simulation`:

```@example tutorial02
init_air = C172.TrimParameters(); #straight and level, default airspeed and altitude
Sim.init!(sim, init_air)
```

You may recall how during [interactive simulation](@ref "Interactive Simulation") we used the
*Aircraft > Avionics > Flight Control* GUI panel to control the aircraft. That panel belongs to the
`mdl.aircraft.avionics.ctl` node, which implements the flight control laws for the `Cessna172Xv1`
aircraft. Let's retrieve its input `struct`:

```@example tutorial02
u_ctl = mdl.aircraft.avionics.ctl.u;
nothing #hide
```

It is this `struct` that the control inputs in the *Aircraft > Avionics > Flight Control* panel
mapped to. Here, we will be writing to it directly. To inspect its fields, we can use the `shf`
helper function:

```@example tutorial02
shf(u_ctl)
```

Note how the trim function, invoked by the `Sim.init!` call, has already set the primary control
inputs (`throttle_axis`, `elevator_axis`, `aileron_axis` and `rudder_axis`) to the values required
by the trimmed flight condition. Also, the longitudinal and lateral control mode inputs
(`lon_ctl_mode_req` and `lat_ctl_mode_req`) are set to their direct values. This bypasses all
control loops, preserving the aircraft's natural dynamic response, which is what we want in this
case.

Let's first advance the `Simulation` a few seconds without perturbing the trim equilibrium. To do
so, we use the `Sim.step!` method, which essentially wraps the one from `DifferentialEquations.jl`'s
[integrator
interface](https://docs.sciml.ai/DiffEqDocs/stable/basics/integrator/#CommonSolve.step!).

```@example tutorial02
Sim.step!(sim, 5) #advance the simulation 5 seconds
```

To apply the elevator doublet, we can either modify the `elevator_axis` input or, more conveniently,
use `elevator_offset` instead. Here's how to do it:

```@example tutorial02
u_ctl.elevator_offset = 0.1 #10 percent positive offset
Sim.step!(sim, 2) #advance 2 seconds
u_ctl.elevator_offset = -0.1 #10 percent negative offset
Sim.step!(sim, 2) #advance 2 seconds
u_ctl.elevator_offset = 0.0 #return to trim position
nothing #hide
```

After this, we don't need to interact with the `Simulation` any further, so we can let it run to
completion:

```@example tutorial02
Sim.run!(sim)
```

Now let's get some results.

A `Simulation`'s log consists of the timestamped values of its root `Model`'s output. The easiest
way to retrieve and handle this data is through the `TimeSeries` type:
```@example tutorial02
ts = TimeSeries(sim)
```

!!! tip "Controlling the logging behavior"

    By default, the `Model`'s outputs are saved at every integration step. You can use the
    `saveat` keyword argument to control logging behavior:
    ```@example tutorial02
    Simulation(mdl; dt = 0.02, t_end = 60, saveat = 0.1) #every 0.1s
    Simulation(mdl; dt = 0.02, t_end = 60, saveat = [0:0.1:10..., 11:1:60...]) #specific instants
    nothing #hide
    ```

A `TimeSeries` object lets you inspect the properties of its underlying data type, and generate
another `TimeSeries` object for any of these properties. This is particularly convenient when
dealing with large, deeply nested types, as `Model` outputs often are. Let's see a few examples:

```@repl tutorial02
propertynames(ts)
propertynames(ts.aircraft)
propertynames(ts.aircraft.vehicle)
ts_kin = ts.aircraft.vehicle.kinematics
ts_sys = ts.aircraft.vehicle.systems
```

```@repl tutorial02
ts_ω = ts_kin.ω_wb_b #angular velocity, Wander frame to Body frame (rad/s)
ts_θ = ts_kin.e_nb.θ #pitch angle, NED frame to Body frame (rad)
ts_α = ts_sys.aero.α #angle of attack (rad)
_, ts_q, _ = get_components(ts_ω); #pitch rate (rad/s)
ts_el_cmd = ts_sys.act.elevator.cmd[4 .<= get_time(ts) .< 10] #elevator command
ts_el_pos = ts_sys.act.elevator.pos[4 .<= get_time(ts) .< 10] #elevator position
```

`Plots.jl` [recipes](https://docs.juliaplots.org/stable/recipes/) are available for many common
`TimeSeries` subtypes, so plotting is usually straightforward:

```@example tutorial02
import Plots
using LaTeXStrings
Plots.default(:size, (900, 600))
Plots.default(:left_margin, 16Plots.px)
nothing #hide
```

```@example tutorial02
#TimeSeries{<:Ranged} recipe
Plots.plot(ts_el_cmd; plot_title="Elevator Response", label = "Command")
Plots.plot!(ts_el_pos; label = "Position")
```

```@example tutorial02
#TimeSeries{<:AbstractVector{<:Real}} recipe
Plots.plot(ts_ω; plot_title="Angular Velocity", ylabel=L"$\omega \ (rad/s)$")
```

```@example tutorial02
#TimeSeries{Real} recipe
Plots.plot(ts_α; plot_title="AoA vs Pitch Angle", ylabel=L"$\alpha, \ \theta \ (rad)$", label="AoA")
Plots.plot!(ts_θ; label = "Pitch Angle")
```

Additionally, the `make_plots` function generates a batch of typically useful plots from a specific
`TimeSeries` subtype. It is available out of the box for `TimeSeries{KinData}`,
`TimeSeries{DynData}` and `TimeSeries{AirData}`. Here are some examples:

```@example tutorial02
plots_kin = make_plots(ts.aircraft.vehicle.kinematics)
```

```@example tutorial02
plots_kin[:Ob_t3d]
```

```@example tutorial02
plots_dyn = make_plots(ts.aircraft.vehicle.dynamics)
```

```@example tutorial02
plots_dyn[:f_c_c]
```

```@example tutorial02
plots_air = make_plots(ts.aircraft.vehicle.airflow)
```

```@example tutorial02
plots_air[:airspeed_M_q]
```

You can batch save them with save_plots

Proceeds recursively down the `Model` hierarchy? Nope


-----------------------------------------

### Automating Model Control With User Callbacks

Introduce callback for doublet.

----------------------------------------

### A More Complex Example

Then move on to complex multi-phase flight. Ground init, startup, take-off, spiral climb until h
condition


world.ac.avionics.ctl.u.lon_ctl_mode_req |> typeof
using .C172X.C172XControl: lon_direct, lon_sas


--------------------------------------------

These values come from fields in the `ModelDefinition` subtype that *are not*
  themselves `ModelDefinition`s.
    These `Model`s are built from those fields in the `ModelDefinition` subtype that *are* themselves
  `ModelDefinition`s.


When the `Model` constructor is called on a `ModelDefinition`, those fields that are themselves
`ModelDefinition`s are also turned into `Model`s, and they become children in the parent `Model`'s
hierarchy. Therefore, within `submodels` we will find `Model`s created from the `Cessna172Xv1`,
`SimpleAtmosphere` and `HorizontalTerrain` objects we passed to the `SimpleWorld` constructor.
```@repl tutorial02
mdl.submodels
```


The `ModelDefinition` subtype is preserved as the `Model`'s first type parameter. It serves both as
the `Model`'s primary identifier and as a dispatch mechanism when calling `Model` update functions.
