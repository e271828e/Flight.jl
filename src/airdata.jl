module Airdata

export AirY

const p0 = 101325
const T0 = 288.15

Base.@kwdef struct AirY
    ps::Float64 = p0
    pt::Float64 = p0
    Tt::Float64 = T0
    α::Float64 = 0.0
    β::Float64 = 0.0
end

end