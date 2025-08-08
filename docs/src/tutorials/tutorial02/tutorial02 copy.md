# Headless Simulation

In this tutorial, we revisit the `Simulation` we set up in the [previous one](@ref "Interactive
Simulation"). This time, we will learn how run it programmatically and extract results for
inspection and plotting.

This tutorial shows how run a `Simulation` programmatically and extract results for inspection and
plotting.

## Peeking Into the `Model`

We start from the same setup we had in the [interactive simulation](@ref "Interactive Simulation")
tutorial:
```@repl tutorial02
using Flight
```

```@repl tutorial02
loc_LOWS15 = LatLon(ϕ = deg2rad(47.80433), λ = deg2rad(12.997));
h_LOWS15 = HOrth(427.2);
ψ_LOWS15 = deg2rad(157);
world = SimpleWorld(Cessna172Xv1(), SimpleAtmosphere(), HorizontalTerrain(h_LOWS15));
```

```@repl tutorial02
init_gnd = KinInit(;
    location = loc_LOWS15, #2D location
    h = h_LOWS15 + C172.Δh_to_gnd, #altitude, as an offset with respect to terrain elevation
    q_nb = REuler(ψ_LOWS15, 0, 0), #attitude with respect to NED frame, as Euler angles
    ω_wb_b = zeros(3), #angular velocity with respect to local tangent frame, aircraft frame coordinates
    v_eb_n = zeros(3), #Earth-relative velocity, NED frame coordinates
    ) |> C172.Init;
```

```@repl tutorial02
init_air = C172.TrimParameters(;
    Ob = Geographic(loc_LOWS15, h_LOWS15 + 500), #500 m above runway 15
    EAS = 50.0, #equivalent airspeed (m/s)
    ψ_nb = ψ_LOWS15, #runway heading (rad)
    γ_wb_n = 0.0, #wind-relative flight path angle (rad)
    ψ_wb_dot = 0.0, #turn rate (rad/s)
    flaps = 0.0, #flaps position (0 to 1)
    fuel_load = 0.5, #normalized fuel load (0 to 1)
);
```

Recall that `SimpleWorld`
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

and then pass it to the `Simulation`

We can retrieve the underlying model from a `Simulation` underlying model as follows:

Recall that a complex model `Model` is made up of multiple, hierarchically arranged components, each of them
itself a `Model`. You can inspect the complete hierarchy of a `Model` as follows:
```@repl tutorial02
using AbstractTrees
print_tree(mdl, maxdepth = 10)
```

You may recognize some of these `Model`s from their GUI panels.
Any node in the `Model` hierarchy can be easily accessed as... For example, to retrieve the left
landing gear...

You may recognize some of these `Model`s from their GUI elements. In particular...
We are particularly interested in the `Controller`'s input struct
shf(). This is precisely the structure that is written by the `Controller`'s GUI panel

Our setup will be similar to that from the [interactive simulation tutorial](@ref "Interactive
Simulation"). However, we no longer need to worry about keeping consistency with X-Plane visuals, so
we can make things a bit simpler:
```@repl tutorial02
using Flight
h_trn = HOrth() #zero MSL
world = SimpleWorld(Cessna172Xv1(), SimpleAtmosphere(), HorizontalTerrain(h_trn));
init_gnd = C172.Init(KinInit(h = h_trn + C172.Δh_to_gnd))
```



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


When the `Model` constructor is called on a `ModelDefinition`, those fields that are themselves
`ModelDefinition`s are also turned into `Model`s, and they become children in the parent `Model`'s
hierarchy. Therefore, within `submodels` we will find `Model`s created from the `Cessna172Xv1`,
`SimpleAtmosphere` and `HorizontalTerrain` objects we passed to the `SimpleWorld` constructor.
```@repl tutorial02
mdl.submodels
```