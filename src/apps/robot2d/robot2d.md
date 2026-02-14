## Problem Statement

The figure shows a simplified 2D model for a Segway-like self-balancing vehicle. The model consists
of two rigid bodies. The main body ($b$) represents the vehicle's chassis, and the rolling body ($r$) represents the axle and the
wheels. The rolling body is driven by a DC motor attached to the main body.

The rolling body has axial symmetry. Its origin is located at its center of mass ($O_r = G_r$), and
it is fixed to the origin of the main body ($O_r \equiv O_b$).

The system's motion is analyzed with respect to an inertial reference frame ($i$), whose axis $z_i$ is aligned with
the local vertical. Its kinematic state is described by the following variables:
- $\theta$: Angle from $z_i$ to $z_b$
- $\phi$: Angle from $z_i$ to $z_r$
- $\eta$: Horizontal position from $O_i$ to $O_b$

The following quantities are known parameters:
- $L$: Distance from $O_{b}$ to $G_{b}$ along $z_b$
- $R$: Wheel radius
- $m_b$: Main body's mass
- $J_b$: Main body's moment of inertia with respect to $G_b$
- $m_r$: Rolling body's mass
- $J_r$: Rolling body's moment of inertia with respect to $G_r$

The wheels roll without slipping, from which we have the kinematic constraint
$$\dot{\phi} = \dot{\eta} / R$$

The net output torque $\tau_m$ produced by the DC motor is given by the following [model](https://thingsdaq.org/2022/07/05/dc-motor-characterization-1-of-2/):
$$
\tau_{m} = k_m u - b_m \omega_m - J_m \dot{\omega}_m
$$

Where $u \in [-1, 1]$ is the motor's control input, $\omega_m$ is the angular rate of the motor's
axis
with respect to its casing, and the following are known parameters:
- $k_m$: Motor's torque constant
- $J_m$: Motor's effective moment of inertia
- $b_m$: Motor's effective damping coefficient

The motor's output torque can also be written as:
$$
\tau_{m} = \tau_{ss} - J_m \dot{\omega}_m
$$

Where the steady-state torque $\tau_{ss} = k_m u - b_m \omega_m$ is that produced by the motor for
a constant angular rate $\omega_m$.

Note: Clockwise angles are positive and $y$ axes point inwards.


## Equations of Motion

To obtain the system's equations of motion we apply Lagrangian mechanics, using $\eta$ and
$\theta$ as generalized coordinates, and treating $\tau_m$ as a non-conservative internal torque,
exerted by the main body on the rolling body through the DC motor.

The resulting equations of motion are:
$$ (m_b L^2 + J_b)\ddot{\theta} + \left(m_b L \cos\theta\right)\ddot{\eta} = m_b L g \sin\theta - \tau_{m} $$

$$ \left(m_b L \cos\theta \right)\ddot{\theta} + \left(m_b + m_r + \frac{J_r}{R^2}\right)\ddot{\eta} = m_b L \omega^2 \sin\theta + \frac{\tau_{m}}{R} $$

Introducing the motor model:
$$ (m_b L^2 + J_b + J_m)\ddot{\theta} + \left(m_b L \cos\theta - \frac{J_m}{R}\right)\ddot{\eta} = m_b L g \sin\theta - \tau_{ss} $$

$$ \left(m_b L \cos\theta - \frac{J_m}{R}\right)\ddot{\theta} + \left(m_b + m_r + \frac{J_r +
J_m}{R^2}\right)\ddot{\eta} = m_b L \omega^2 \sin\theta + \frac{\tau_{ss}}{R} $$

For convenience, we define:
$$ M_{11} = m_b L^2 + J_b + J_m $$
$$ M_{22} = m_b + m_r + \frac{J_r + J_m}{R^2} $$
$$ M_{12}(\theta) = m_b L \cos\theta - \frac{J_m}{R} $$

So that:
$$ M_{11} \ddot{\theta} + M_{12}(\theta) \ddot{\eta} = m_b L g \sin\theta - \tau_{ss} $$

$$ M_{12}(\theta) \ddot{\theta} + M_{22} \ddot{\eta} = m_b L \omega^2 \sin\theta + \frac{\tau_{ss}}{R} $$

This can be rewritten as a first-order state-space system:

$$ M_{11} \dot{\omega} + M_{12}(\theta) \dot{v} = m_b L g \sin\theta - \tau_{ss} \tag{1}$$

$$ M_{12}(\theta) \dot{\omega} + M_{22} \dot{v} = m_b L \omega^2 \sin\theta + \frac{\tau_{ss}}{R} $$

$$ \dot{\theta} = \omega $$

$$ \dot{\eta} = v $$

Or, in matrix form:
$$ \boldsymbol{M}(\boldsymbol{x}) \dot{\boldsymbol{x}} = \boldsymbol{b}(\boldsymbol{x}) $$

$$ \boldsymbol{x} = \begin{bmatrix} \omega \\ v \\ \theta \\ \eta \end{bmatrix} $$
$$ \boldsymbol{M}(\boldsymbol{x}) = \begin{pmatrix} M_{11} & M_{12}(\theta) & 0 & 0 \\ M_{12}(\theta) & M_{22} & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix} $$

$$ \boldsymbol{b}(\boldsymbol{x}) = \begin{bmatrix} m_b L g \sin\theta - \tau_{ss} \\ m_b L \omega^2 \sin\theta + \frac{\tau_{ss}}{R} \\ \omega \\ v \end{bmatrix} $$

With:
$$
\tau_{ss} = k_m u - b_m \omega_m = k_m u - b_m \left( v / R - \omega \right) \tag{5}
$$

## Initialization and Steady State Solution
To initialize the system, a total of 10 unknowns ($\dot{\omega}$, $\dot{v}$, $\omega$, $v$,
$\dot{\theta}$, $\dot{\eta}$, $\theta$, $\eta$, $u$, $\tau_{ss}$) must be determined for $t=0$. We have 5
equations, so we need to impose 5 constraints.

First, we'd like the vehicle to start with zero angular acceleration and zero linear acceleration. Eliminating $\tau_{ss}$ from $(1)$ and $(2)$ yields: $$ (M_{11} +R M_{12}(\theta))
\dot{\omega} + (M_{12}(\theta) + R M_{22}) \dot{v} = (m_b L g + m_b L R \omega^2) \sin\theta$$

This shows that $\dot{\omega}(0) = \dot{v}(0)=0$ requires the chassis to be perfectly vertical ($\theta(0) = 0$ or
$\theta(0) = \pi$). We'd like the vehicle to be initially upright, so we choose $\theta(0) = 0$.

With $\dot{\omega}(0) = \dot{v}(0) = \theta(0) = 0$, $(1)$ yields $\tau_{ss}(0) = 0$.  Then, if we choose an
arbitrary $u(0) = u_0$, linear velocity can be found from $(5)$ as:
$$v(0) = \frac{k_m u_0}{b_m} R$$

In principle, $\omega(0)$ can also be set arbitrarily to some $\omega_0$. Then, $\dot{\theta}(0)$ is
given by $(3)$. For the initial $\theta(0) = 0$ to be maintained, we must obviously choose $\omega_0 = 0$.

Finally, $\dot{\eta}(0)$ is given by $(4)$, while the horizontal position $\eta(0)$ is decoupled from all other
variables and can be set to an arbitrary value $\eta_0$.

To recap, our constraints are:
- $\dot{\omega}(0) = 0$
- $\dot{v}(0) = 0$
- $\omega(0)=\omega_0$
- $\eta(0) = \eta_0$
- $u(0) = u_0$

And, from the system's equations, we find:
- $\theta(0) = 0$
- $\tau_{ss}(0) = 0$
- $v(0) = {R k_m u_0}/{b_m}$
- $\dot{\theta}(0) = \omega_0$
- $\dot{\eta}(0) = v(0)$
