# Headless Simulation

In this tutorial, we revisit our setup from the [interactive simulation](@ref "Interactive
Simulation") walkthrough. This time, we will learn how to run the `Simulation` programmatically and
extract results for inspection and plotting.

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

A `ModelDefinition` subtype can be thought of as the blueprint for a specific `Model` instance. More
precisely, it defines a `Model`'s states, inputs, outputs, constants and subcomponents. To
instantiate a `Model`, we call the `Model` constructor on a `ModelDefinition` subtype:
```@repl tutorial02
mdl = Model(world)
```

A complex `Model` can contain multiple, hierarchically arranged subcomponents, each one of them
itself a `Model`. To inspect a `Model`'s hierarchy, you can do the following:
```@repl tutorial02
using AbstractTrees
print_tree(mdl, maxdepth = 10)
```

Recall that during interactive simulation, we used the *Aicraft > Avionics > Flight Control* GUI
panel to control the aircraft. This panel corresponds to the flight control `Model`, which we can
retrieve as.

We are particularly
interested in this model's input structure, which we can access as:

It is the fields of this structure that the GUI panel's inputs actually mapped to. Here, we will be
assigning them directly.
shf(mdl.aircraft.avionics.ctl)

You may recognize some of these components from their GUI panels.

Any node in a `Model`'s hierarchy can be easily accessed via dot notation. For example, to retrieve
the left landing gear...

You may recognize some of these `Model`s from their GUI elements. In particular...
We are particularly interested in the `Controller`'s input struct
shf(). This is precisely the structure that is written by the `Controller`'s GUI panel

The resulting `Model` is ready for simulation:
```@repl tutorial02
sim = Simulation(mdl; dt = 0.02)
```

If a `ModelDefinition` is passed directly to the `Simulation` constructor, `Model` instantiation
occurs implicitly as an intermediate step. The result is the same.

Once the `Simulation` is created, the underlying model can be easily retrieved:
```@repl tutorial02
propertynames(sim)
sim.mdl
```

---
Not necessary


Furthermore, all of its fields are themselves `ModelDefinition` subtypes:
```@repl tutorial02
fieldtypes(SimpleWorld)
supertype.(fieldtypes(SimpleWorld))
```

Let's see what's inside:
```@repl tutorial02
propertynames(mdl)
```

Brief explanation of properties

No need for this:

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