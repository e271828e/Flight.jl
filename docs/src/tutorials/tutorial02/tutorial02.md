# Headless Simulation

In this tutorial, we revisit our [interactive simulation](@ref "Interactive Simulation") setup. This
time, we will learn how to run the `Simulation` programmatically and extract results for inspection
and plotting.

### Peeking Into the Simulation `Model`

Let's begin by recreating our `SimpleWorld` instance:
```@example tutorial02
using Flight

loc_LOWS15 = LatLon(ϕ = deg2rad(47.80433), λ = deg2rad(12.997))
h_LOWS15 = HOrth(427.2)
ψ_LOWS15 = deg2rad(157)
world = SimpleWorld(Cessna172Xv1(), SimpleAtmosphere(), HorizontalTerrain(h_LOWS15))
nothing #hide
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

The `ModelDefinition` subtype is preserved as the `Model`'s first type parameter. It serves as the
`Model`'s primary identifier, and it is used for dispatch when calling `Model` update functions.

Let's take a look at our `Model`'s properties:
```@repl tutorial02
propertynames(mdl)
```

Here's a brief description for each one:
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
them itself a `Model`. To visualize a `Model`'s hierarchy, you can do the following:
```@repl tutorial02
using AbstractTrees
print_tree(mdl; maxdepth = 10)
```

You can also inspect a specific property across a `Model`'s hierarchy. For example, to view the
discrete state `s` of the `avionics` node and every node below, you would do:
```@repl tutorial02
print_tree(mdl.aircraft.avionics, :s; maxdepth = 10);
```

You may recognize some of the `Model`s in the hierarchy from the GUI panels you saw in the
[interactive simulation](@ref "Interactive Simulation") tutorial.

Recall that we used the *Aicraft > Avionics > Flight Control* GUI
panel to control the aircraft. This panel belongs to the flight control `Model`, which we can
retrieve as.

We are particularly interested in the `Controller`'s input struct
shf(). This is precisely the structure that is written by the `Controller`'s GUI panel

It is the fields of this structure that the GUI panel's inputs actually mapped to. Here, we will be
assigning them directly.
shf(mdl.aircraft.avionics.ctl)

-----------------------------------------

### Simulation

Next now create a `Simulation` for our `Model`:
```@repl tutorial02
sim = Simulation(mdl; dt = 0.02)
```

The `Model` is stored within the `Simulation` and can be retrieved at any moment for inspection or
manipulation:
```@repl tutorial02
propertynames(sim)
mdl === sim.mdl
```

If an instance of a `ModelDefinition` subtype is passed directly to the `Simulation` constructor (as
we did in the [interactive simulation](@ref "Interactive Simulation") tutorial), the `Model` is
automatically instantiated under the hood.

--------------------------------------

world.ac.avionics.ctl.u.lon_ctl_mode_req |> typeof
using .C172X.C172XControl: lon_direct, lon_sas

------------------------------------

Furthermore, all of its fields are themselves `ModelDefinition` subtypes:
```@repl tutorial02
fieldtypes(SimpleWorld)
supertype.(fieldtypes(SimpleWorld))
```

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

```@example tutorial02
init_gnd = KinInit(;
    location = loc_LOWS15,
    h = h_LOWS15 + C172.Δh_to_gnd,
    q_nb = REuler(ψ_LOWS15, 0, 0),
    ) |> C172.Init

init_air = C172.TrimParameters(;
    Ob = Geographic(loc_LOWS15, h_LOWS15 + 500),
    EAS = 50.0,
    ψ_nb = ψ_LOWS15,
)
nothing #hide
```
