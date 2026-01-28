## Problem Statement

The figure shows a simplified 2D model for a Segway-like self-balancing vehicle. The model consists
of two rigid bodies and a DC motor. The main body (1) represents the vehicle's chassis, to which the
motor's casing is attached. The rolling body (2) comprises the motor's rotor, the axle and the
wheels.

The system's motion with respect to an inertial reference frame (i) is described by the following kinematic variables:
- $\theta$: Angle of the main body's z-axis with respect to the vertical
- $\phi$: Angle of the rolling body's z-axis with respect to the vertical
- $\eta$: Horizontal position of the system's origin with respect to the inertial frame origin

Additionally, the following quantities are defined as known parameters:
- $L$: Distance from the main body's origin ($O_{b1}$) to its CoM ($G_{b1}$)
- $R$: Wheel radius
- $m_1$: Main body's mass
- $J_1$: Main body's moment of inertia with respect to its CoM
- $m_2$: Rolling body's mass
- $J_2$: Rolling body's moment of inertia with respect to its CoM

The wheels roll without slipping, from which we have the kinematic constraint
$$\dot{\phi} = \dot{\eta} / R$$

The net output torque $\tau_m$ produced by the DC motor is given by the following [model](https://thingsdaq.org/2022/07/05/dc-motor-characterization-1-of-2/):
$$
\tau_{m} = k_m u_m - b_m \omega_m - J_m \dot{\omega}_m
$$

Where $u_m \in [-1, 1]$ is the motor's control input, $\omega_m$ is the angular rate of the motor's
axis
with respect to its casing, and the following are known parameters:
- $k_m$: Motor's torque constant
- $J_m$: Motor's effective moment of inertia
- $b_m$: Motor's effective damping coefficient

The motor's output torque can also be written as:
$$
\tau_{m} = \tau_{ss} - J_m \dot{\omega}_m
$$

Where the steady-state torque $\tau_{ss} = k_m u_m - b_m \omega_m$ is that produced by the motor for
a constant angular rate $\omega_m$.


Notes:
- Positive angles are clockwise and $y$ axes point inwards.
- The rolling body has axial symmetry, so it's CoM is coincident with origins $O_{b1}$ and $O_{b2}$


## Equations of Motion

To obtain the system's equations of motion we apply Lagrangian mechanics, using $\eta$ and
$\theta$ as generalized coordinates, and treating $\tau_m$ as a non-conservative internal torque,
exerted by the main body on the rolling body through the DC motor.

The resulting equations of motion are:
$$ (m_1 L^2 + J_1)\ddot{\theta} + \left(m_1 L \cos\theta\right)\ddot{\eta} = m_1 L g \sin\theta - \tau_{m} $$

$$ \left(m_1 L \cos\theta \right)\ddot{\theta} + \left(m_1 + m_2 + \frac{J_2}{R^2}\right)\ddot{\eta} = m_1 L \omega^2 \sin\theta + \frac{\tau_{m}}{R} $$

Introducing the motor model:
$$ (m_1 L^2 + J_1 + J_m)\ddot{\theta} + \left(m_1 L \cos\theta - \frac{J_m}{R}\right)\ddot{\eta} = m_1 L g \sin\theta - \tau_{ss} $$

$$ \left(m_1 L \cos\theta - \frac{J_m}{R}\right)\ddot{\theta} + \left(m_1 + m_2 + \frac{J_2 +
J_m}{R^2}\right)\ddot{\eta} = m_1 L \omega^2 \sin\theta + \frac{\tau_{ss}}{R} $$

For convenience, we define:
$$ M_{11} = m_1 L^2 + J_1 + J_m $$
$$ M_{22} = m_1 + m_2 + \frac{J_2 + J_m}{R^2} $$
$$ M_{12}(\theta) = m_1 L \cos\theta - \frac{J_m}{R} $$

So that:
$$ M_{11} \ddot{\theta} + M_{12}(\theta) \ddot{\eta} = m_1 L g \sin\theta - \tau_{ss} $$

$$ M_{12}(\theta) \ddot{\theta} + M_{22} \ddot{\eta} = m_1 L \omega^2 \sin\theta + \frac{\tau_{ss}}{R} $$

This system can be rewriten as a first-order state-space system:

$$ M_{11} \dot{\omega} + M_{12}(\theta) \dot{v} = m_1 L g \sin\theta - \tau_{ss} $$

$$ M_{12}(\theta) \dot{\omega} + M_{22} \dot{v} = m_1 L \omega^2 \sin\theta + \frac{\tau_{ss}}{R} $$

$$ \dot{\theta} = \omega $$

$$ \dot{\eta} = v $$

With:
$$
\tau_{ss} = k_m u_m - b_m \omega_m = k_m u_m - b_m \left( v / R - \omega \right)
$$

In matrix form:
$$ \boldsymbol{x} = \begin{bmatrix} \omega \\ v \\ \theta \\ \eta \end{bmatrix} $$
$$ \boldsymbol{M}(\boldsymbol{x}) \dot{\boldsymbol{x}} = \boldsymbol{b}(\boldsymbol{x}) $$
$$ \boldsymbol{M}(\boldsymbol{x}) = \begin{pmatrix} M_{11} & M_{12}(\theta) & 0 & 0 \\ M_{12}(\theta) & M_{22} & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix} $$

$$ \boldsymbol{b}(\boldsymbol{x}) = \begin{bmatrix} m_1 L g \sin\theta - \tau_{ss} \\ m_1 L \omega^2 \sin\theta + \frac{\tau_{ss}}{R} \\ \omega \\ v \end{bmatrix} $$

## Initialization and Steady State Solution
We have 10 unknowns: $\dot{\omega}$, $\dot{v}$, $\omega$, $v$, $\dot{\theta}$, $\dot{\eta}$,
$\theta$, $\eta$, $u_m$, $\tau_{ss}$.

And 5 equations:
$$ M_{11} \dot{\omega} + M_{12}(\theta) \dot{v} = m_1 L g \sin\theta - \tau_{ss} $$

$$ M_{12}(\theta) \dot{\omega} + M_{22} \dot{v} = m_1 L \omega^2 \sin\theta + \frac{\tau_{ss}}{R} $$

$$ \dot{\theta} = \omega $$

$$ \dot{\eta} = v $$

$$
\tau_{ss} = k_m u_m - b_m \left( v / R - \omega \right)
$$


Thus, we need 5 constraints. We choose $u_m$, $\theta$, $\eta$, $\omega$ and $\dot{\omega}$.

Rearranging the two equations of motion:
$$ \tau_{ss} + M_{12}(\theta) \dot{v} = m_1 L g \sin\theta - M_{11} \dot{\omega} $$

$$ -\frac{1}{R} \tau_{ss} + M_{22} \dot{v} = m_1 L \omega^2 \sin\theta - M_{12}(\theta)
\dot{\omega}$$

With $\theta$, $\omega$ and $\dot{\omega}$ as parameters, we can solve for $\tau_{ss}$ and
$\dot{v}$.

Then, we find $v$ as:
$$ \omega_m = (k_m u_m - \tau_{ss}) / b_m $$
$$ v = (\omega + \omega_m ) * R $$

Finally, $\dot{\theta}$ and $\dot{\eta}$ are given by their kinematic equations.

Eliminating $\tau_{ss}$ from the equations of motion yields:
$$ (M_{11} +R M_{12}(\theta)) \dot{\omega} + (M_{12}(\theta) + R M_{22}) \dot{v} =
(m_1 L g + m_1 L R \omega^2) \sin\theta$$

This shows that, in order to have both zero angular acceleration ($\dot{\omega} = 0$) and zero linear acceleration
($\dot{v}=0$), the chassis must be perfectly vertical ($\theta = 0$ or $\theta = \pi$). Since
$\dot{\theta} = \omega$, for this steady state to be maintained we must also have $\omega = 0$.

In the steady state, $\tau_{ss} = 0$, and linear velocity is given by:
$$ v = \frac{k_m u_m}{b_m} R$$
