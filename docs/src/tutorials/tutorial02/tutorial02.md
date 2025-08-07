# Headless Simulation

In this tutorial, we revisit the `Simulation` we set up in the [previous one](@ref "Interactive
Simulation"). This time, we will learn how run it programmatically and extract results for
inspection and plotting.

This tutorial shows how run a `Simulation` programmatically and extract results for inspection and
plotting.

## Peeking Into the `Model`

Our setup will be similar to that from the [interactive simulation tutorial](@ref "Interactive
Simulation"). However, we no longer need to worry about keeping consistency with X-Plane visuals, so
we can make things a bit simpler:
```@repl tutorial02
using Flight
h_trn = HOrth() #zero MSL
world = SimpleWorld(Cessna172Xv1(), SimpleAtmosphere(), HorizontalTerrain(h_trn));
init_gnd = C172.Init(KinInit(h = h_trn + C172.Δh_to_gnd))
```

Our setup will be the same as in the [interactive simulation tutorial](@ref "Interactive
Simulation"):
```@repl tutorial02
init_gnd = KinInit(;
    h = h_trn + C172.Δh_to_gnd, #altitude, as an offset with respect to terrain elevation
    ) |> C172.Init
```

Recall that passing the `SimpleWorld`...

We can retrieve the underlying model `Simulation` underlying model as

Also, We can extract the underlying Model from a Simulation

The resulting `Model` is made up of multiple, hierarchically arranged components, each of them
itself a `Model`. The complete hierarchy of a `Model` can be inspected as:
```@repl tutorial02
using AbstractTrees
print_tree(Model(world))
```

You may recognize some of these `Model`s from their GUI panels.
Any node in the `Model` hierarchy can be easily accessed as... For example, to retrieve the left
landing gear...

You may recognize some of these `Model`s from their GUI elements. In particular...
We are particularly interested in the `Controller`'s input struct
shf(). This is precisely the structure that is written by the `Controller`'s GUI panel

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

First, note that `SimpleWorld` is a concrete subtype of the abstract type `ModelDefinition`:
```@repl tutorial02
using InteractiveUtils
supertypes(SimpleWorld)
```

A `ModelDefinition` can be understood as a blueprint defining how a specific `Model` should be
built. In the [interactive simulation](@ref "Interactive Simulation") tutorial, this `Model`
instantiation step happened implicitly when we passed our `SimpleWorld` object to the `Simulation`
constructor. We can also get a standalone `Model` by calling its constructor explicitly:
```@repl tutorial02
mdl = Model(world)
```


When the `Model` constructor is called on a `ModelDefinition`, those fields that are themselves
`ModelDefinition`s are also turned into `Model`s, and they become children in the parent `Model`'s
hierarchy. Therefore, within `submodels` we will find `Model`s created from the `Cessna172Xv1`,
`SimpleAtmosphere` and `HorizontalTerrain` objects we passed to the `SimpleWorld` constructor.
```@repl tutorial02
mdl.submodels
```