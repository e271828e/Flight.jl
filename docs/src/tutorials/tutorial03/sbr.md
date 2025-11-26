Notes:

- Main body assembly includes the rigidly attached motor case
- Wheel assembly includes the rigidly attached motor shaft

Parameters:
- $L$: Distance from $O_{b1}$ to $O_{c1}$ (CoM) along $z_{b1}$

Caution: Need to revise J1

## Main Body Linear Momentum Equation
$$
F_{ext,b_1}^i + m_1 g^i = m_1 \dot{v}_{ic_1}^i
$$

Where
- $F_{ext,b_1}^i$: External force applied on $b_1$ expressed in $i$ coordinates
- $g^{i}$: Gravity vector expressed in $i$ coordinates
- ${v}_{ic_1}^i$: Velocity of frame $c_1$ with respect to $i$, expressed in $i$coordinates

Let's develop the terms:
$$
r_{ic_1}^{i} = \begin{bmatrix}
p + L \sin \theta &&
0 &&
R + L \cos \theta
\end{bmatrix} ^T
$$

$$
v_{ic_1}^{i} = \begin{bmatrix}
\dot{p} + L \dot{\theta} \cos \theta &&
0 &&
- L \dot{\theta} \sin \theta
\end{bmatrix} ^T
$$

$$
\dot{v}_{ic_1}^{i} = \begin{bmatrix}
\ddot{p} + L \left( \ddot{\theta} \cos \theta - \dot{\theta}^2 \sin \theta \right) &&
0 &&
- L \left( \ddot{\theta} \sin \theta + \dot{\theta}^2 \cos \theta\right)
\end{bmatrix} ^T
$$

$$
\dot{v}_{ic_1}^{i} = \begin{bmatrix}
\dot{v} + L \left( \dot{\omega}_{1} \cos \theta - \omega_{1}^2 \sin \theta \right) &&
0 &&
- L \left( \dot{\omega}_{1} \sin \theta + \omega_{1}^2 \cos \theta\right)
\end{bmatrix} ^T
$$

$$
m_1g^i = \begin{bmatrix} 0 && 0 && -m_1 g \end{bmatrix}^T
$$

$$
F_{ext,b_1}^i = \begin{bmatrix} F_{21x} && 0 && F_{21z} \end{bmatrix}^T
$$

Finally:

$$
F_{21x} = m_1 \dot{v} + m_1 L \left(\dot{\omega}_1 \cos \theta - \omega_{1}^2 \sin \theta \right)
$$

$$
F_{21z} - m_1 g = m_1 L  \left( -\dot{\omega}_{1} \sin \theta -\omega_{1}^2 \cos \theta\right)
$$


## Rod Assembly Angular Momentum Equation

$$
\tau_{ext,b_1[c_1]}^{c_1} = \dot{h}_{b_1[c_1]}^{c_1} + \omega_{ib_1}^{c_1}\times h_{b_1[c_1]}^{c_1}
$$
$$
h_{b_1[c_1]}^{c_1} = J_{b_1[c_1]}^{c_1} \omega_{ib_1}^{c_1}
$$
$$
\tau_{ext,b_1[c_1]}^{c_1} = J_{b_1[c_1]}^{c_1} \dot{\omega}_{ib_1}^{c_1} + \Omega_{ib_1}^{c_1}J_{b_1[c_1]}^{c_1} \omega_{ib_1}^{c_1}
$$

Where
- $\tau_{ext,b_1[c_1]}^{c_1}$: External torque applied on $b_1$ with respect to the origin of frame
  $c_1$ expressed in $c_1$ coordinates.
- $h_{b_1[c_1]}^{c_1}$: Angular momentum of $b_1$ with respect to the origin of $c_1$, expressed in
  $c_1$ coordinates.
- $J_{b_1[c_1]}^{c_1}$: Inertia tensor of $b_1$ with respect to the origin of $c_1$, expressed in
  $c_1$ coordinates.

$$
\tau_{ext,b_1[c_1]}^{c_1} = \begin{bmatrix}
0 &&
\tau_{21} + L \left(F_{21z} \sin \theta - F_{21x} \cos \theta \right)   &&
0
\end{bmatrix}^T
$$

$$
\omega_{ib_1}^{c_1} = \begin{bmatrix} 0 && \omega_{1} && 0 \end{bmatrix}^T
$$

$$
h_{b_1[c_1]}^{c_1} = J_{b_1[c_1]}^{c_1} \omega_{ib_1}^{c_1} = \begin{bmatrix}
J_{1xx} && 0 && 0 \\
0 && J_{1yy} && 0 \\
0 && 0 && J_{1zz} \\
\end{bmatrix}
\begin{bmatrix}
0 \\ \omega_{1} \\ 0
\end{bmatrix} =
\begin{bmatrix}
0 && J_{1yy} \omega_{1} && 0
\end{bmatrix}^T
$$

$$
\omega_{ib_1}^{c_1} \times h_{b_1[c_1]}^{c_1} = \begin{bmatrix} 0 && 0 && 0 \end{bmatrix}^T
$$

$$
\dot{h}_{b_1[c_1]}^{c_1} = \begin{bmatrix} 0 && J_{1yy} && 0
\end{bmatrix}
$$

$$
\tau_{21} + L \left(F_{21z} \sin \theta - F_{21x} \cos \theta \right) = J_{1yy} \dot{\omega}_{1}
$$

The torque exerted by $b_1$ on $b_2$ is the torque produced by the motor:
$$
\tau_{12} = \tau_{m}
$$

Therefore
$$
\tau_{21} = -\tau{12} = -\tau_{m}
$$

$$
-\tau_{m} + L \left(F_{21z} \sin \theta - F_{21x} \cos \theta \right) = J_{1yy} \dot{\omega}_{1}
$$

The rod assembly's moment of inertia about its center of mass might be approximated as that of a
thin rod:
$$ J_{1yy} = \dfrac{1}{12} m_1 (2L)^2 \omega_{1} $$

## Wheel Assembly Linear Momentum Equation

$$
F_{ext,b_2}^i + m_2 g^i = m_2 \dot{v}_{ic_2}^i
$$

Where
- $F_{ext,b_2}^i$: External force applied on $b_2$ expressed in $i$ coordinates
- $g^{i}$: Gravity vector expressed in $i$ coordinates
- ${v}_{ic_2}^i$: Velocity of frame $c_2$ with respect to $i$, expressed in $i$coordinates

$$
r_{ic_2}^{i} = \begin{bmatrix} p && 0 && R \end{bmatrix} ^T
$$

$$
v_{ic_2}^{i} = \begin{bmatrix} \dot{p} && 0 && 0 \end{bmatrix} ^T
$$

$$
\dot{v}_{ic_2}^{i} = \begin{bmatrix}
\dot{v} && 0 && 0\end{bmatrix} ^T
$$

$$
m_2g^i = \begin{bmatrix} 0 && 0 && -m_2 g \end{bmatrix}^T
$$

$$
F_{ext,b_2}^i = \begin{bmatrix} F_{i2x} + F_{12x} && 0 && F_{i2z} + F_{12z} \end{bmatrix}^T
$$
$$
F_{ext,b_2}^i = \begin{bmatrix} F_{i2x} - F_{21x} && 0 && F_{i2z} - F_{21z} \end{bmatrix}^T
$$

Equating components gives
$$
F_{i2x} - F_{21x} = m_2 \dot{v}
$$
$$
F_{i2z} - F_{21z} - m_2 g = 0
$$


## Wheel Assembly Angular Momentum Equation

$$
\tau_{ext,b_2[c_2]}^{c_2} = \dot{h}_{b_2[c_2]}^{c_2} + \omega_{ib_2}^{c_2}\times h_{b_2[c_2]}^{c_2}
$$
$$
h_{b_2[c_2]}^{c_2} = J_{b_2[c_2]}^{c_2} \omega_{ib_2}^{c_2}
$$
$$
\tau_{ext,b_2[c_2]}^{c_2} = J_{b_2[c_2]}^{c_2} \dot{\omega}_{ib_2}^{c_2} + \Omega_{ib_2}^{c_2}J_{b_2[c_2]}^{c_2} \omega_{ib_2}^{c_2}
$$

Where
- $\tau_{ext,b_2[c_2]}^{c_2}$: External torque applied on $b_2$ with respect to the origin of frame
  $c_2$ expressed in $c_2$ coordinates.
- $h_{b_2[c_2]}^{c_2}$: Angular momentum of $b_2$ with respect to the origin of $c_2$, expressed in
  $c_2$ coordinates.
- $J_{b_2[c_2]}^{c_2}$: Inertia tensor of $b_2$ with respect to the origin of $c_2$, expressed in
  $c_2$ coordinates.

$$
\tau_{ext,b_2[c_2]}^{c_2} = \begin{bmatrix}
0 &&
\tau_{12} - RF_{i2x}    &&
0
\end{bmatrix}^T
$$

$$
\omega_{ib_2}^{c_2} = \begin{bmatrix} 0 && \omega_{2} && 0 \end{bmatrix}^T
$$

$$
h_{b_2[c_1]}^{c_1} = J_{b_1[c_1]}^{c_1} \omega_{ib_1}^{c_1} =
\begin{bmatrix}
0 && J_{2yy} \omega_{2} && 0
\end{bmatrix}^T
$$

$$
\omega_{ib_2}^{c_2} \times h_{b_2[c_2]}^{c_2} = \begin{bmatrix} 0 && 0 && 0 \end{bmatrix}^T
$$

$$
\dot{h}_{b_2[c_2]}^{c_2} = \begin{bmatrix} 0 && J_{2yy}  \dot{\omega}_{2} && 0
\end{bmatrix}
$$

$$
\tau_{12} - RF_{i2x} =  \dot{\omega}_{2}
$$

$$
\tau_{m} - RF_{i2x} =  \dot{\omega}_{2}
$$

The moment of inertia of the wheel assembly is that of the wheel itself plus that of the motor:
$$
J_{2yy} = J_{m} + \dfrac{1}{2} m_2 R^2
$$

## DC Motor Model
https://thingsdaq.org/2022/07/05/dc-motor-characterization-1-of-2/
$$
\tau_{m} = K_V  u_{m} - b_{m} \omega_{m}
$$

$$
\tau_{m} = K_V  u_{m} - b_{m} (\omega_{2} - \omega_{1})
$$

Parameters:
- Effective moment of inertia at motor shaft: $J_{m} = 0.0014 \ kg\cdot m^2$
- Effective damping coefficient: $b_{m} = 0.0189 \ N \cdot m / (rad /s)$
- Torque constant: $K_V = 0.32 \ N \cdot m$

Variables:
- Angular velocity of motor shaft: $\omega_{m}$
- Voltage ratio: $u_{m} \in [-1, 1]$

## No-Slip Condition

$$
v_{iQ}^{i} = v_{ib_2}^{i} + \omega_{ib_2}^{i} \times r_{O_2Q}^{i} =
\begin{bmatrix}
v - \omega_{2} R &&
0 &&
0
\end{bmatrix}^T = 0
$$

$$
\dot{v} = \dot{\omega}_{2} R
$$

## Kinematics

$$
\dot{\theta} = \omega_{1}
$$
$$
\dot{p} = v = \omega_{2} R
$$

## System Equations
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

### Motor Model
$$
\tau_{m} = K_V u_{m} - b_{m} (\omega_{2} - \omega_{1})
$$

### No-Slip Condition

$$
v = \omega_{2} R
$$

### Kinematics

$$
\dot{\theta} = \omega_{1}
$$

$$
\dot{p} = v
$$

---

Setting the decoupled kinematic equations aside, we have the
following unknowns:

$\dot{\omega}_{1}$, $\dot{\omega}_{2}$, $\dot{v}$, $F_{21x}$, $F_{21z}$, $F_{i2x}$, $F_{i2z}$, $\tau_{m}$

---
Reorder terms:

### Dynamics
$$
\left(m_1 L \cos \theta \right) \dot{\omega}_1+ m_1 \dot{v} - F_{21x} = m_1 L \omega_{1}^2 \sin \theta
$$

$$
\left(m_1 L \sin \theta \right) \dot{\omega}_{1} + F_{21z} = m_1 g -m_1 L \omega_{1}^2 \cos \theta
$$

$$
J_{1yy} \dot{\omega}_{1} + \left(L \cos \theta
\right)F_{21x} + \left( -L \theta \right)F_{21z} + \tau_{m} = 0
$$

$$
m_2 \dot{v} + F_{21x} - F_{i2x} = 0
$$

$$
- F_{21z} + F_{i2z} = m_2 g
$$

$$
J_{2yy} \dot{\omega}_{2} + RF_{i2x} - \tau_{m} = 0
$$

### Motor Model
$$
\tau_{m} = K_V u_{m} - b_{m} (\omega_{2} - \omega_{1})
$$

### No-Slip Condition
$$
R \omega_{2} - v = 0
$$

---
Substitute the no-slip condition
Applying the linear and angular momentum equations to rod + motor case assembly and wheel + rotor
shaft assembly, and using the
no-slip rolling constraint on the wheel:

### Dynamics
$$
\left(m_1 L \cos \theta \right) \dot{\omega}_1+ m_1 R \dot{\omega}_{2} - F_{21x} = m_1 L \omega_{1}^2 \sin \theta
$$

$$
\left(m_1 L \sin \theta \right) \dot{\omega}_{1} + F_{21z} = m_1 g -m_1 L \omega_{1}^2 \cos \theta
$$

$$
J_{1yy} \dot{\omega}_{1} + \left(L \cos \theta
\right)F_{21x} + \left( -L \sin \theta \right)F_{21z} + \tau_{m} = 0
$$

$$
m_2 R \dot{\omega}_{2} + F_{21x} - F_{i2x} = 0
$$

$$
- F_{21z} + F_{i2z} = m_2 g
$$

$$
J_{2yy} \dot{\omega}_{2} + RF_{i2x} - \tau_{m} = 0
$$

### Motor Model
$$
\tau_{m} = K_V u_{m} - b_{m} (\omega_{2} - \omega_{1})
$$

---

Matrix form:
$$

\begin{bmatrix}
m_1 L \cos \theta && m_1 R && -1 && 0 && 0 && 0 && 0 \\
m_1 L \sin \theta && 0 && 0 && 1 && 0 && 0 && 0 \\
J_{1yy} && 0 && L \cos \theta && -L \sin \theta && 0 && 0 && 1 \\
0 && m_2R && 1 && 0 && -1 && 0 && 0 \\
0 && 0 && 0 && -1 && 0 && 1 && 0 \\
0 && J_{2yy} && 0 && 0 && R && 0 && -1 \\
0 && 0 && 0 && 0 && 0 && 0 && 1 \\
\end{bmatrix}

\begin{bmatrix}
\dot{\omega}_{1} \\
\dot{\omega}_{2} \\
F_{21x} \\
F_{21z} \\
F_{i2x} \\
F_{i2z} \\
\tau_{m} \\
\end{bmatrix}
=
\begin{bmatrix}
m_1 L \omega_{1}^2 \sin \theta \\
m_1 g -m_1 L \omega_{1}^2 \cos \theta \\
0 \\
0 \\
m_2 g \\
0 \\
K_V u_{m} - b_{m} (\omega_{2} - \omega_{1})
\end{bmatrix}