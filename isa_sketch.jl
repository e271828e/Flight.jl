using StructArrays

# layers = [
#     (β = -6.5e-3, h_ceil = 11000.0)
#     (β = 0.0, h_ceil = 20000.0)
#     (β = 1e-3, h_ceil = 32000.0)
#     (β = 2.8e-3, h_ceil = 47000.0)
#     (β = 0.0, h_ceil = 51000.0)
#     (β = -2.8e-3, h_ceil = 71000.0)
#     (β = -2e-3, h_ceil = 80000.0)
# ]

sa_layers = StructArray(
    β =      Vector{Float64}([-6.5e-3, 0, 1e-3, 2.8e-3, 0, -2.8e-3, -2e-3]),
    h_ceil = Vector{Int}([11000, 20000, 32000, 47000, 51000, 71000, 80000])
)

@inline function get_layer_data(h::Real)
    for i in 1:length(sa_layers)
        h < sa_layers.h_ceil[i] && return sa_layers[i]
    end
    throw(ArgumentError("Altitude out of bounds"))
end