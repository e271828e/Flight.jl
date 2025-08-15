# Headless Simulation

In this tutorial, we revisit our [interactive simulation](@ref "Interactive Simulation") setup. This
time, we will learn how to run the `Simulation` programmatically and extract results for inspection
and plotting.

### Peeking Into the `Model`

Let's start from a new `SimpleWorld` instance. Here, keeping consistency with X-Plane's visuals
is no longer a concern, so we can stick to the default constructors:
```@repl tutorial02
using Flight

world = SimpleWorld(Cessna172Xv1(), SimpleAtmosphere(), HorizontalTerrain())
```

Inspecting its type hierarchy reveals that `SimpleWorld` is a concrete subtype of the abstract type
`ModelDefinition`:
```@repl tutorial02
using InteractiveUtils
supertypes(SimpleWorld)
```

An instance of a `ModelDefinition` subtype can be thought of as a blueprint specifying how
a particular `Model` should be built. To instantiate this `Model`, we pass the `ModelDefinition`
object to the `Model` constructor:
```@repl tutorial02
mdl = Model(world)
```

Let's take a look at our `Model`'s properties:
```@repl tutorial02
propertynames(mdl)
```

Here's a brief description of each one:
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
```@repl tutorial02
using AbstractTrees
print_tree(mdl.aircraft.avionics.ctl; maxdepth = 10)
```

You can also inspect a specific property across a `Model`'s hierarchy. For instance, to view the
discrete state of `mdl.aircraft.avionics.ctl` and every node underneath:
```@repl tutorial02
print_tree(mdl.aircraft.avionics.ctl, :s; maxdepth = 10);
```

### Simulating an Elevator Doublet

Our plan in this section is as follows:
1. Initialize the aircraft in a trimmed state.
2. Inject an elevator doublet by setting the appropriate inputs and stepping the `Simulation` manually.
3. Let the `Simulation` run to completion.
4. Extract and plot some useful variables from the `Simulation`'s log to observe the aircraft's
   response.

Let's first create a `Simulation` from our `Model`:

```@repl tutorial02
sim = Simulation(mdl; dt = 0.02, t_end = 60)
```

The `Model` is now stored within the `Simulation`. It can be accesed at any time for inspection
or manipulation:

```@repl tutorial02
:mdl ∈ propertynames(sim)
sim.mdl === mdl
```

We could also pass a `ModelDefinition` object directly to the `Simulation` constructor. In that
case, a new `Model` would be instantiated automatically.

The initialization step is straightforward:

```@repl tutorial02
init_air = C172.TrimParameters(); #straight and level, default airspeed and altitude
Sim.init!(sim, init_air)
```

What next? You may recall how during [interactive simulation](@ref "Interactive Simulation") we used
the *Aircraft > Avionics > Flight Control* GUI panel to control the aircraft. That panel belongs to
the `mdl.aircraft.avionics.ctl` node, which implements the flight control laws for the
`Cessna172Xv1` aircraft. Let's retrieve its input `struct`:

```@repl tutorial02
u_ctl = mdl.aircraft.avionics.ctl.u;
```

It is this `struct` that the control inputs in *Aircraft > Avionics > Flight Control* map to. Here,
we will be writing to it directly. Let's inspect its fields:

```@repl tutorial02
shf(u_ctl)
```

Note how the trim function, invoked by the `Sim.init!` call, has already set the primary control
inputs (`throttle_axis`, `elevator_axis`, `aileron_axis` and `rudder_axis`) to the values required
by the trimmed flight condition. Also, the longitudinal and lateral control mode inputs
(`lon_ctl_mode_req` and `lat_ctl_mode_req`) are set to their direct values. This bypasses all
control loops, preserving the aircraft's natural dynamic response, which is what we want here.

Let's first advance the `Simulation` a few seconds without perturbing the trim equilibrium. To do
so, we use the `Sim.step!` method, which essentially wraps the one from `DifferentialEquations.jl`'s
[integrator
interface](https://docs.sciml.ai/DiffEqDocs/stable/basics/integrator/#CommonSolve.step!).

```@repl tutorial02
Sim.step!(sim, 5) #advance the simulation 5 seconds
```

Now, to apply the elevator doublet, we can either modify the `elevator_axis` input directly or, more
conveniently, use `elevator_offset` instead. Here's how to do it:

```@repl tutorial02
u_ctl.elevator_offset = 0.1; #10 percent positive offset
Sim.step!(sim, 2) #advance 2 seconds
u_ctl.elevator_offset = -0.1; #10 percent negative offset
Sim.step!(sim, 2) #advance 2 seconds
u_ctl.elevator_offset = 0.0; #return to trim position
Sim.run!(sim) #run the Simulation to completion (t_end)
```

The data logged by a `Simulation` are the values of its root `Model`'s output. By default, these are
saved at every integration step, but different saving intervals can be specified. The easiest way to
retrieve and handle this data is through the `TimeSeries` type:
```@repl tutorial02
ts = TimeSeries(sim)
```

A `TimeSeries` object lets you inspect the properties of its underlying data type and generate
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

`Plots.jl` recipes are available for many common `TimeSeries` subtypes, so plotting is usually
straightforward:

```@example tutorial02
import Plots
using LaTeXStrings

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

Function `make_plots` generates a batch of useful plots from a specific `TimeSeries` subtype. Out of
the box, it is implemented for `TimeSeries{KinData}`, `TimeSeries{DynData}` and
`TimeSeries{AirData}`, so for example you can do:

Plots.default(:size, (900, 600))
display(t3d)

You can batch save them with save_plots

Don't mention recursive...
Proceeds recursively down the `Model` hierarchy.


-----------------------------------------

### Automating Simulation Control With User Callbacks

Introduce callback for doublet.

Then move on to complex multi-phase flight. Ground init, startup, take-off, spiral climb until h
condition


world.ac.avionics.ctl.u.lon_ctl_mode_req |> typeof
using .C172X.C172XControl: lon_direct, lon_sas



----------------------------------------

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
