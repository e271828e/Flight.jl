# API

## FlightPhysics

### Dynamics

```@docs
FrameTransform
```

```@docs
Dynamics.translate(t_bc::FrameTransform, r_cP_c::AbstractVector{<:Real})
```

```@docs
Base.:adjoint(t_bc::FrameTransform)
```

```@docs
Base.:∘(t_bc::FrameTransform, t_cd::FrameTransform)
```

```@docs
Wrench
```

```@index
```
