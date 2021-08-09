## Electric Powerplant

From Vehicle Dynamics, we have, for each rotating element, the following scalar angular momentum equation along its
x-axis (assumed to be the axis of rotation)

$J_k^{xx} \left( R^{b_0}_{s_k}\right)^T \dot{\omega}^{b_0}_{eb_0} + J_k^{xx} \dot{\omega}_{k} =$
$1_x^T M_{ext,B_k[O_k]}^{s_k} + M_{B_0,B_k}$

For most atmospheric vehicles, $\left| J_{B_k} \right| \ll \left| J_{sys}\right|$. This means that, the angular
acceleration experienced by the airframe due to the torque exchanged with the rotating element will be negligible
compared to that experienced by the rotating element. This justifies dropping the first term on the left hand side,
which uncouples the rotational dynamics of the rotating element and the airframe. We are then left with:

$J_k^{xx} \dot{\omega}_{k} = 1_x^T M_{ext,B_k[O_k]}^{s_k} + M_{B_0,B_k}$

Applying this to the motor shaft we can write:

$J_s \dot{\omega}_s = M_{gb,s} + M_{eng,s}$

Where:
- $J_s$: Axial moment of inertia of the motor shaft
- $\omega_s$: Angular velocity of the motor shaft relative to the airframe
- $M_{gb,s}$: Torque transmitted by the gearbox to the motor shaft
- $M_{eng,s}$: Torque output by the motor to the motor shaft

We now need to add an angular for the propeller:

$J_p \dot{\omega}_p = M_{air,p} + M_{gb,p}$

Where:
- $J_p$: Axial moment of inertia of the propeller
- $\omega_p$: Angular velocity of the propeller relative to the airframe
- $M_{gb,p}$: Torque transmitted by the gearbox to the propeller
- $M_{air,s}$: Torque exerted by the airflow on the propeller *along its rotating axis*

For the gearbox, neglecting its internal moment of inertia, we have:

$\omega_p = \omega_s / n$ \
$M_{gb,p} = \eta n M_{s,gb}$

Where:
- $n = N_p / N_s$: Gear ratio, where $N_p$ is the number of teeth of the propeller gear and $N_s$ is the number of teeth
  of the shaft gear. Typically, $N_p > N_s$, which reduces the propeller angular velocity with respect to that of the
  motor shaft and amplifies the torque transmitted through the gearbox. Setting $n<0$ represents an inverting gearbox,
  in which both the angular velocity and the transmitted torque have their sense switched by the gearbox.
- $M_{gb,p}$: Torque transmitted by the gearbox to the propeller
- $M_{s,gb}$: Torque transmitted by the motor shaft to the gearbox
- $\eta$: Gearbox efficiency

Finally, we have:\
$M_{s, gb} = -M_{gb, s}$

The torque at the motor shaft is given by the following first-order electric motor model [Drela] is:

$M_{eng,s} = \alpha \left[ \left(uV_b - \dfrac{\alpha\omega_s}{Rk_V} \right) - i_0 \right] \dfrac{1}{k_M}$

Where:
- $\alpha \in \{-1, 1\}$: Turning sense (-1: CCW, 1: CW)
- $u$: Throttle input
- $V_b$: Battery voltage
- $\omega_s$: Angular velocity of the motor shaft
- $R$: Internal resistance
- $i_0$: No-load current
- $k_V$: Velocity constant ($(rad/s)/V$)
- $k_M$: Torque constant ($A/(Nm)$), often assumed equal to $k_V$

An overly simplistic model for the aerodynamic drag torque on the propeller is:\
$M_{air,p} = -sgn(\omega_p) k_M \omega_p^2$

Combining the above equations we arrive at the complete system:

$M_{eng,s} + \dfrac{M_{air,p}}{\eta n}  = \left(J_s + \dfrac{J_p}{\eta n^2} \right)\dot{\omega}_s$

$M_{eng,s} = \alpha \left[ \left(uV_b - \dfrac{\alpha\omega_s}{Rk_V} \right) - i_0 \right] \dfrac{1}{k_M}$

$M_{air,p} = -sgn(\omega_p) k_M \omega_p^2$

In which the turn sense of the propeller is determined by the choice of $\alpha$ and $sgn(n)$.

The additional angular momentum contributed by the engine, gearbox and propeller assembly to the overall vehicle due to
their angular velocity with respect to the airframe is:

$h_{[G_c]}^{c} = \begin{pmatrix} J_p \omega_p + J_s \omega_s & 0 & 0 \end{pmatrix}$ \

Where $c$ are the local component axes.