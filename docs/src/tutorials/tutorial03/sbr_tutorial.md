### Intro

This classic problem from dynamics and control theory will allow us to explore the core principles
of the modeling framework by building a new "plant" model from first principles.

Write it from scratch?



Let's consider the following variation on the classic 2D inverted pendulum problem. The system
depicted in the diagram consists consists of:
- Chassis:


Describe bodies forces, moments, etc

Applying the Newton-Euler equations to each body and the no-slipping kinematic constraint ($v =
\omega_2 R$) yields the following six scalar equations.

Applying the Newton-Euler rigid body equations to the chassis, we get

Doing the same for the wheel and axle yields

We also have the no-slipping kinematic constraint:

Finally, we will use the following simplified DC motor model


### Dynamics
$$
F_{21x} = m_1 \dot{v} + m_1 L \left(\dot{\omega}_1 \cos \theta - \omega_{1}^2 \sin \theta \right)
$$

$$
F_{21z} - m_1 g = m_1 L  \left( -\dot{\omega}_{1} \sin \theta -\omega_{1}^2 \cos \theta\right)
$$

$$
-\tau_{m} + L \left(F_{21z} \sin \theta - F_{21x} \cos \theta \right) = J_{1yy} \dot{\omega}_{1}
$$

$$
F_{i2x} - F_{21x} = m_2 \dot{v}
$$

$$
F_{i2z} - F_{21z} - m_2 g = 0
$$

$$
\tau_{m} - RF_{i2x} = J_{2yy} \dot{\omega}_{2}
$$

The no-slip kinematic constraint is

Additionally, let's consider the following DC motor model:


Substituting the no-slip constraint


State vector:
