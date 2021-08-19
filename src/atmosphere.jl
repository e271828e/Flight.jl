module Atmosphere

export AbstractAtmosphericModel, DummyAtmosphericModel

abstract type AbstractAtmosphericModel end
struct DummyAtmosphericModel <: AbstractAtmosphericModel end

end