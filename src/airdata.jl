module Airdata

using ComponentArrays

import Flight.System: Y
export AirData, AirDataY

const p0 = 101325
const T0 = 288.15

struct AirData end

const AirDataYTemplate = ComponentVector(ps = p0, pt = p0, Tt = T0, α = 0.0, β = 0.0)
const AirDataY{D} = ComponentVector{Float64, D, typeof(getaxes(AirDataYTemplate))} where {D<:AbstractVector{Float64}}
Y(::AirData) = copy(AirDataYTemplate)

end