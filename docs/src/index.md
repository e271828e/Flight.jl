## Installation

```julia
using Pkg
Pkg.add("Flight")
```

!!! warning "Configuring Julia for Multithreading"

    During simulation, ```Flight.jl``` uses multithreading to concurrently handle the simulation loop, the built-in GUI and any attached I/O devices. Therefore, if you intend to use its interactive or I/O capabilities, you will need to start Julia with multiple threads enabled. You can find out how to do this in the [manual](https://docs.julialang.org/en/v1/manual/multi-threading/). However, if you are using the Julia extension for VS Code, the easiest way is to add the following entry to your
    ```settings.json```:

    ```json
    "julia.additionalArgs": [
        "--threads=auto",
    ],
    ```

    This should work well for most CPUs and use cases. If you run into issues, you can set the number of threads manually to cover your specific needs.