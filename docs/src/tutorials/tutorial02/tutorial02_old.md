# Headless Simulation

In this tutorial we revisit the `Model` from the previous [interactive simulation](@ref "Interactive
Simulation") tutorial. We will take a look at its internals and use what we learn to set up and run
a simulation programmatically.

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

Furthermore, all of its fields are themselves `ModelDefinition` subtypes:
```@repl tutorial02
fieldtypes(SimpleWorld)
supertype.(fieldtypes(SimpleWorld))
```

A `ModelDefinition` can be understood as a blueprint defining how a specific `Model` should be
built. In the [interactive simulation](@ref "Interactive Simulation") tutorial, this `Model`
instantiation step happened implicitly when we passed our `SimpleWorld` object to the `Simulation`
constructor. However, we can get a standalone `Model` by calling its constructor instead:
```@repl tutorial02
mdl = Model(world)
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