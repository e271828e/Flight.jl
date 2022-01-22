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
tendra HasMass, HasNoWrench, HasNoAngularMomentum

"""

#let P be the origin of the JSBSim aircraft main reference frame and b the
#conventional flight physics axes. looks like JSBSim main reference frame is
#roughly at the center of the mean chord

#let 0 denote the empty aircraft (OEW)

m0 = 894
r_PG0_b = SA[-1.041, 0, -0.927]
J0_G0_b = SA[1285 0 0; 0 1825 0; 0 0 2667] #OEW
J0_P_b = J0_G0_b - m0 * skew(r_PG0_b)^2

#pilot
m1 = 75
r_PG1_b = SA[-0.914, -0.356, -0.610]
J1_G1_b = SA[0 0 0; 0 0 0; 0 0 0]
J1_P_b = J1_G1_b - m1 * skew(r_PG1_b)^2

#copilot
m2 = 75
r_PG2_b = SA[-0.914, 0.356, -0.610]
J2_G2_b = SA[0 0 0; 0 0 0; 0 0 0]
J2_P_b = J2_G2_b - m2 * skew(r_PG2_b)^2


#si queremos J_P_b:

J_P_b = J0_P_b + J1_P_b + J2_P_b



#one option is to make Ob = G0

MTOW = 1406

#

#we need to switch x and z signs


#aerodynamics (SI)
S = 16.165
b = 10.912
c = 1.494

#aerodynamics frame origin:
r_POa_b = [-1.097, 0, -1.509]