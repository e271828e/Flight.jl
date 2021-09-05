module Airdata

using StaticArrays

using Flight.Kinematics
using Flight.Atmosphere

export AirY

Base.@kwdef struct AirY
    v::SVector{3,Float64} = zeros(SVector{3}) #v_wOb_b
    α::Float64 = 0.0
    β::Float64 = 0.0
    ps::Float64 = 0.0
    pt::Float64 = 0.0
    Tt::Float64 = 0.0
    M::Float64 = 0.0
    TAS::Float64 = 0.0
    EAS::Float64 = 0.0
    CAS::Float64 = 0.0
end

# estas magnitudes van en y_air::AirY, parte de AircraftY, y por tanto estan
# disponibles para cualquier Component. A partir de ahi un Component (sobre todo
# Aerodynamics, pero no exclusivamente) tiene la posibilidad de usar su propio
# frame para expresar v_w, alpha y beta. Para ello simplemente convierte v =
# v_wOb_b a sus propios ejes v_wOb_c, y luego llama a get_airflow_angles para
# recalcular alpha_a y beta_a. Incluso puede tener en cuenta la velocidad
# angular y calcular v_wOa_a != v_wOb_b:

# v_wOb_b = v_eOb_b - v_ew_b
# v_eOa_b = v_eOb_b + ω_eb_b × r_ObOa_b
# v_wOa_b = v_eOa_b - v_ew_b = v_eOb_b + ω_eb_b × r_ObOa_b - v_ew_b
# v_wOa_b = v_wOb_b + ω_eb_b × r_ObOa_b
# v_wOa_a = q_ba'(v_wOa_b)

#en general, v_wOa ≈ v_wOb y no sera necesaria esta conversion. como mucho, si
#la base de datos aerodinamica esta definida en otros ejes, si puede ser
#necesaria la conversion a v_wOb_a para obtener α_a y β_a, que son los que usara
#la aerodynamic database.

get_airflow_angles() = nothing

function get_air_data(::KinY, ::AbstractAtmosphericModel)

    #get position from KinY, query the AtmosphericModel for the atmospheric
    #properties
    #get velocity from KinY, compute pt, Tt, M, alpha, beta, TAS
    # compute air data at the airframe reference frame

end

end