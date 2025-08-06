# Headless Simulation

In this tutorial we revisit the `Model` from the previous [interactive simulation](@ref "Interactive
Simulation") tutorial. We will take a look at its internals and use what we learn to set up and run
a simulation programmatically.

## Peeking Into the `Model`

Let's begin by initializing the package and creating a `SimpleWorld` object. Since we are no longer
concerned with keeping consistency with X-Plane, we can simply use the default constructors:
```@repl tutorial02
using Flight
world = SimpleWorld(Cessna172Xv1(), SimpleAtmosphere(), HorizontalTerrain());
```

<!-- First, we note that `Cessna172Xv1`, `SimpleAtmosphere`, `HorizontalTerrain` and `SimpleWorld` is a concrete subtype of the abstract type `ModelDefinition`: -->
First, we note that `Cessna172Xv1`, `SimpleAtmosphere`, `HorizontalTerrain` and `SimpleWorld` are
all concrete subtypes of the abstract type `ModelDefinition`:
```@repl tutorial02
using InteractiveUtils
supertype(Cessna172Xv1)
supertypes(SimpleAtmosphere)
supertypes(HorizontalTerrain)
supertypes(SimpleWorld)
```

A `ModelDefinition` can be understood as a blueprint defining how a specific `Model` should be
built. In the [interactive simulation](@ref "Interactive Simulation") tutorial, this `Model`
instantiation step happened implicitly when we passed our `SimpleWorld` object to the `Simulation`
constructor. However, we can also get a standalone `Model` by calling its constructor instead:
```@repl tutorial02
mdl = Model(world)
```

When the `Model` constructor is called on a `ModelDefinition`, those fields that are themselves
`ModelDefinition`s are also converted into `Model`s, and they become submodels in the parent
`Model`'s hierarchy.

Let's see how this looks:
```@repl tutorial02
propertynames(mdl)
keys(world.submodels)
```

Inside the `submodels` field, we will find

Since these were also `ModelDefinition`s,
