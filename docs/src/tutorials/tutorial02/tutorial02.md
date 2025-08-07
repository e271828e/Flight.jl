# Headless Simulation

In this tutorial we revisit the Simulation we set up in the [interactive simulation tutorial]. This
this time we will learn how to run it programmatically and extract the results for inspection and
plotting

## Peeking Into the `Model`

Let's begin by initializing the package and creating a `SimpleWorld` object:
```@repl tutorial02
using Flight
world = SimpleWorld(Cessna172Xv1(), SimpleAtmosphere(), HorizontalTerrain());
```

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

Also, We can extract the underlying Model from a Simulation

The resulting `Model` is made up of multiple, hierarchically arranged components, each of them
itself a `Model`. The complete hierarchy of a `Model` can be inspected as:
```@repl tutorial02
using AbstractTrees
print_tree(mdl)
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


When the `Model` constructor is called on a `ModelDefinition`, those fields that are themselves
`ModelDefinition`s are also turned into `Model`s, and they become children in the parent `Model`'s
hierarchy. Therefore, within `submodels` we will find `Model`s created from the `Cessna172Xv1`,
`SimpleAtmosphere` and `HorizontalTerrain` objects we passed to the `SimpleWorld` constructor.
```@repl tutorial02
mdl.submodels
```