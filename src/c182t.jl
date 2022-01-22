using Flight
### mass properties

const skew = Attitude.skew

"""
la obligacion del airframe completo para con las ecuaciones dinamicas es
proporcionar m, r_ObG_b y J_Ob_b. ahora, pensemos que pasa si del airframe
cuelgan dos componentes c1 y c2. que tiene que pedirle el airframe a cada uno de
ellos para poder cumplir con su obligacion con las ecuaciones dinamicas?
necesita m_c1 y m_c2 para calcular m = m_c1 + m_c2, r_ObGC1_b y r_ObGC2 para
calcular r_ObG_b = 1/m * (m_c1 * r_ObG1_b + m_c2 * r_ObG2_b) y J1_Ob_b y
J2_Ob_b para calcular J_Ob_b = J1_Ob_b + J2_Ob_b
y que necesita a su vez
"""

"""
al airframe de la c182 le anado un subsystem que sea payload. este a su vez
tendra HasMass, HasNoWrench, HasNoAngularMomentum. pilot, copilot, passengerp
left, passenger right, baggage
airframe lo creo como SystemGroup, de manera que

"""

#define the pilot in its own reference frame as a point mass with m = 75,
#r_OcG_c = zeros(3), J_Oc_c = 0. its reference frame is simply a translation of
#fb.

Base.@kwdef struct Payload <: SystemDescriptor
    pilot_mass::MassProperties = translate(FrameTransform(r = SA[-0.914, -0.356, -0.610]), MassProperties(m = 75))
    copilot_mass::MassProperties = translate(FrameTransform(r = SA[-0.914, 0.356, -0.610]), MassProperties(m = 75))
end
#DISCRETE STATES: pilot present, pilot_not present

#al crear el systemdescriptor, las masas de piloto y copiloto ya aparecen
#expresadas en fb. por tanto, get_mass_properties(System{Payload}) es trivial,
#se limita a sumar las MassProperties de cada uno atendiendo a si estan
#presentes o no. pero ojo, que pasa si no hay nadie? no podemos sumar las
#contribuciones si son todas cero. si, deberiamos
mass = OEW(afr::System{Airframe})
function get_mass_properties(pl::System{Payload})
    mass = MassProperties()
    pl.d.pilot_present ? mass += payload.params.pilot : nothing
    pl.d.copilot_present ? mass += payload.params.copilot : nothing
    return mass
end

#we need to switch x and z signs

#let P be the origin of the JSBSim aircraft main reference frame and b the
#conventional flight physics axes. looks like JSBSim main reference frame is
#roughly at the center of the mean chord


#let 0 denote the empty aircraft (OEW)
m0 = 894
r_PG0_b = SA[-1.041, 0, -0.927]
J0_G0_b = SA[1285 0 0; 0 1825 0; 0 0 2667] #OEW

#ahora, la implementacion de get_mass_properties para airframe se limitara a
#llamar a la de payload para obtener su contribucion en fb, que a su vez
#dependera de




MTOW = 1406

#



#aerodynamics (SI)
S = 16.165
b = 10.912
c = 1.494

#aerodynamics frame origin:
r_POa_b = [-1.097, 0, -1.509]