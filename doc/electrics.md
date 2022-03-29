# Electric Powerplant

## Angular Momentum Equation

For each rotating element attached to the airframe, we have the following scalar angular momentum equation along its
axis of rotation (arbitrarily chosen as its x-axis):

$J_k^{xx} \left( R^{b_0}_{s_k}\right)^T \dot{\omega}^{b_0}_{eb_0} + J_k^{xx} \dot{\omega}_{k} =$
$1_x^T M_{ext,B_k[O_k]}^{s_k} + M_{B_0,B_k}$

For most atmospheric vehicles, $\left| J_{B_k} \right| \ll \left| J_{sys}\right|$. This means that, the angular
acceleration experienced by the airframe due to the torque exchanged with the rotating element will be negligible
compared to that experienced by the rotating element. This justifies dropping the first term on the left hand side,
which uncouples the rotational dynamics of the rotating element and the airframe. We are then left with:

$J_k^{xx} \dot{\omega}_{k} = 1_x^T M_{ext,B_k[O_k]}^{s_k} + M_{B_0,B_k}$

### Motor

Applying this to the motor shaft we can write:

$J_s \dot{\omega}_s = M_{gb,s} + M_{eng,s}$

Where:
- $J_s$: Axial moment of inertia of the motor shaft
- $\omega_s$: Angular velocity of the motor shaft relative to the airframe
- $M_{gb,s}$: Torque transmitted by the gearbox to the motor shaft
- $M_{eng,s}$: Torque output by the motor to the motor shaft

### Propeller

Similarly, for the propeller:

$J_p \dot{\omega}_p = M_{air,p} + M_{gb,p}$

Where:
- $J_p$: Axial moment of inertia of the propeller
- $\omega_p$: Angular velocity of the propeller relative to the airframe
- $M_{gb,p}$: Torque transmitted by the gearbox to the propeller
- $M_{air,s}$: Torque exerted by the airflow on the propeller *along its rotating axis*

An overly simplistic model for the aerodynamic drag torque on the propeller is:\
$M_{air,p} = -sgn(\omega_p) k_M \omega_p^2$

To avoid chattering around $\omega_p = 0$, we can replace $sgn(\omega_p)$ with a sigmoid-like function, such as $tanh$:

$sgn(\omega_p) \approx \tanh \left(\dfrac{\omega_p}{\omega_{ref}}\right)$

For $|\omega_{p}| / \omega_{ref} = 2$, $\tanh \left(\dfrac{\omega_p}{\omega_{ref}}\right) = 0.964 sgn(\omega_p)$

With this:

$M_{air,p} = -\tanh \left(\dfrac{\omega_p}{\omega_{ref}}\right) k_M \omega_p^2$

### Gearbox

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

## Motor Model

We use a simple first-order electric motor model [Drela], modified to account for both CW and CCW nominal turn senses.

The voltage across the motor's terminals is given by:

$e_m(\omega_s) = \alpha \omega_s / k_V$\
$v_m = e_m(\omega_s) + iR_m$

Where:
- $e_m (V)$: Back-EMF
- $i (A)$: Current
- $\alpha \in \{-1, 1\}$: Motor's nominal turn sense (-1: CCW, 1: CW)
- $k_V ((rad/s)/V)$: Motor's velocity constant
- $R_m (\Omega)$: Motor's internal resistance.
- $\omega_s (rad/s)$: Signed angular velocity of the motor shaft

The shaft torque is given by:

$M_{eng,s} = \alpha \dfrac{i}{k_M} + M_{f}$

Where:
- $k_M$: Motor's torque constant ($A/(Nm)$), often assumed equal to $k_V$
- $M_{f}$: Motor's friction torque

The motor's friction torque opposes the angular velocity. Substituting the sign function for $tanh$ we have:

$M_{f} = -sgn(\omega_s) \left| M_f \right| \approx -\tanh \left(\dfrac{\omega_s}{\omega_{ref}}\right) |M_f|$

If the friction torque is assumed constant, we can express it as:

$|M_f| \approx \dfrac{i_0}{k_M}$

Where $i_0$ is the motor's no-load current. Then:

$M_f \approx -\tanh \left(\dfrac{\omega_s}{\omega_{ref}}\right) \dfrac{i_0}{k_M}$

We then have:

$M_{eng,s} = \dfrac{1}{k_M} \left[ \alpha i -\tanh \left(\dfrac{\omega_s}{\omega_{ref}}\right) i_0 \right]$

## Battery Model

The effective voltage supplied by the battery can be modeled as:

$V_b(c) = n_c V_c f(c)$\
$R_b = n_c R_c$\
$v_b = uV_b(c) - iR_b$

Where:
- $V_b (V)$: No-load battery voltage
- $n_c$: Number of cells connected in series
- $V_c (V)$: No-load, maximum charge cell voltage
- $c = C/C_{max}\in [0,1]$: Battery charge ratio
- $f \in [0,1]$: Cell voltage curve
- $R_b$: Battery internal resistance
- $R_c$: Cell internal resistance
- $u$: Throttle

The battery charge ratio evolves as:

$\dot{c} = \dot{C} / C_max = -i / C_{max}$

## Electric Circuit

The voltage equation for the electric circuit is:

$v_b - v_m = uV_b(c) - iR_b - (e_m(\omega_s) + iR_m) =
uV_b(c) - e_m(\omega_s) - i(R_b +R_m) = 0$

From which:

$i = \dfrac{uV_b(c) - e_m(\omega_s)}{R_b +R_m}$


### System Dynamics

Combining the above equations we arrive at the complete system:

$M_{eng,s} + \dfrac{M_{air,p}}{\eta n}  = \left(J_s + \dfrac{J_p}{\eta n^2} \right)\dot{\omega}_s$

$i = uV_b - \dfrac{\alpha\omega_s}{Rk_V}$

$M_{eng,s} = \dfrac{1}{k_M} \left[ \alpha i -\tanh \left(\dfrac{\omega_s}{\omega_{ref}}\right) i_0 \right]$

$M_{air,p} = -\tanh \left(\dfrac{\omega_p}{\omega_{ref}}\right) k_M \omega_p^2$

In which the turn sense of the propeller is determined by the choice of $\alpha$ and $sgn(n)$.

The additional angular momentum contributed by the engine, gearbox and propeller assembly to the overall vehicle due to
their angular velocity with respect to the airframe is:

$h_{[G_c]}^{c} = \begin{pmatrix} J_p \omega_p + J_s \omega_s & 0 & 0 \end{pmatrix}$

Where $c$ are the local component axes.