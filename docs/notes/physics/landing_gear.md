## Landing Gear Model

### Assumptions

1. The inertia of the landing gear assembly is neglected, so compression forces are transmitted directly to the airframe

2. The landing gear strut (whether it is an oleo or a spring steel strut) compresses along a straight line. Its force depends on compression length and compression rate.

3. No detailed tire dynamics

### Definitions

Strut frame $s$:
- $O_s$: Strut-airframe attachment point. Its position vector with respect to the vehicle origin $r_{O_bO_s}^b$ is known and constant.
- $z_s$: Parallel to the landing gear compression line, positive along strut elongation.
- $y_s$: Parallel to the wheel axle when the steering angle $\psi_{sw}$ is zero.
The attitude of the strut frame axes with respect to the vehicle axes is constant and known.

Wheel frame $w$:
- $O_w$: Wheel's endpoint along $z_s$.
- $z_w$: Parallel to $z_s$.
- $y_w$: Parallel to the wheel axle. The angle from $y_s$ to $y_w$ is the steering angle $\psi_{sw}$.

Contact frame $c$:
- $O_c$: Intersection of the strut line (origin $O_s$, direction $z_s$) with the terrain surface. When the wheel is in contact with the ground, $O_c = O_w$.
- $z_c$: Parallel to the terrain's inward surface normal at $O_c$, given by its ECEF components $k_t^e$.
- $x_c$: Parallel to the projection of the $x_w$ axis onto the terrain tangent plane at $O_c$.

Distance parameters:
- $l$: Signed distance from $O_s$ to $O_c$ along $z_s$.
- $l_0$: Strut natural length.
- $\xi = l - l_0$: Strut deformation. Ground contact requires $\xi \leq 0$. When the wheel is in the air, the strut is at its natural length and $\xi = 0$.

The force exerted by the strut outwards along the positive $z_s$ is a function $F(\xi, \dot{\xi})$. Therefore, $\xi < 0$ does not necessarily mean that the strut is actually pushing on the ground. In particular, if the strut is compressed ($\xi<0$) but expanding ($\dot{\xi} > 0$), the function $F(\xi, \dot{\xi})$ could become negative. Of course, this is not actually possible, because in that case the wheel would lift off the ground momentarily. Therefore, the actual force exerted on the ground is $\max\{0,F(\xi, \dot{\xi})\}$.

We define a maximum force $F_{max}$ such that if $F(\xi, \dot{\xi}) > F_{max}$, we consider the mechanical load has exceeded its admissible limit and we abort the simulation.

### Contact Determination

We start by expressing the strut line in ECEF coordinates:

$$R^e_s = R^e_b R^b_s$$

$$k^e_s = R^e_s e_3$$

$$r_{O_eO_s}^e = r_{O_eO_b}^e + R^e_b r_{O_bO_s}^b$$

Contact determination is delegated to the terrain model through a bounded ray query. Given the ray origin $O_s$, its unit direction $k_s^e$ and the bound $l_{max} = l_0$, the terrain model returns the nearest surface point $O_c$ along the ray with parameter $l \in [0, l_0]$, together with the inward unit surface normal $k_t^e$ at $O_c$ and the surface type. If no such point exists, there is no contact and we are done. Otherwise:

$$l = (k_s^e)^T (r_{O_eO_c}^e - r_{O_eO_s}^e)$$

$$\xi = l - l_0 \leq 0$$

The bound on the ray parameter is what turns the query into a contact test: the natural-length wheel endpoint $O_{w0}$ lies at $l = l_0$ along the strut line, so a hit at $l \leq l_0$ means that $O_{w0}$ is at or below the surface. Any cheap altitude-based pre-test to skip the intersection when the aircraft is far from the ground is also the terrain model's responsibility. How the intersection is computed depends on the terrain model; see [Terrain](terrain.md).

A hit does not by itself imply a physically meaningful contact. If the aircraft is inverted or steeply banked, the strut line may still pierce the surface at $l \leq l_0$, with the strut nearly tangent to it. We therefore compute the angle from the contact normal to the strut axis:

$$\cos \alpha_{cs} = (k_t^e)^T k_s^e$$

and abort the simulation with a ground crash when $\alpha_{cs}$ exceeds a threshold (60°). Note that the terrain model guarantees $\cos \alpha_{cs} > 0$ for any hit, since a ray parallel to the surface or pointing away from it cannot intersect it.


### Contact Frame Construction

We start by constructing the wheel axes $w$ as a rotation of the strut axes $s$ by an angle
$\psi_{sw}$ (the steering angle) around $z_s$.

$$R^s_w = R_z(\psi_{sw})$$

$$R^e_w = R^e_s R^s_w$$

$$i^e_w = R^e_w i^w_w = R^e_w e_1$$

The direction of the $x_c$ axis is given by the projection of the $x_w$ axis onto the terrain
tangent plane. Thus, $i_c$ can be computed by subtracting from $i_w$ its projection along the terrain
normal $k_t$, and then normalizing:

$$i_c^e = \dfrac{i_w^e - (i_w^e \cdot k_t^e) k_t^e}{\left|i_w^e - (i_w^e \cdot k_t^e) k_t^e\right|}$$

The $z_c$ axis is parallel to $k_t$, and $y_c$ can be constructed from $z_c$ and $x_c$:

$$k_c^e = k_t^e$$

$$j_c^e = k_c^e \times i_c^e$$

Then:

$$R^e_c = {\begin{bmatrix} i_c^e & j_c^e & k_c^e \end{bmatrix}}$$

$$R^s_c = (R^e_s)^T R^e_c$$

$$R^b_c = R^b_s R^s_c$$

The position of $O_c$ is given by:

$${r}_{O_sO_c}^s = {\begin{pmatrix} 0 & 0 & l \end{pmatrix}}^T$$

$$r_{O_bO_c}^b = r_{O_bO_s}^b + R^b_s r_{O_sO_c}^s$$

### Computing $v_{eO_c}^c$

$$r_{O_eO_c}^e = r_{O_eO_b}^e + r_{O_bO_c}^e = r_{O_eO_b}^e + R^e_b r_{O_bO_c}^b$$

$$\dot{r}_{O_eO_c}^e = \dot{r}_{O_eO_b}^e + R^e_b \Omega_{eb}^b r_{O_bO_c}^b + R^e_b \dot{r}_{O_bO_c}^b$$

$${v}_{eO_c}^b = {v}_{eO_b}^b + \Omega_{eb}^b r_{O_bO_c}^b + \dot{r}_{O_bO_c}^b$$

Now:

$$r_{O_bO_c}^b = r_{O_bO_s}^b + R^b_s r_{O_sO_c}^s$$

$$\dot{r}_{O_bO_c}^b = \dot{r}_{O_bO_s}^b + R^b_s \dot{r}_{O_sO_c}^s = R^b_s \dot{r}_{O_sO_c}^s$$

$$\dot{r}_{O_sO_c}^s = {\begin{pmatrix} 0 & 0 & \dot{\xi} \end{pmatrix}}^T = e_3 \dot{\xi}$$

$$\dot{r}_{O_bO_c}^b = R^b_s e_3 \dot{\xi}$$

Then:

$${v}_{eO_c}^b = {v}_{eO_b}^b + \Omega_{eb}^b r_{O_bO_c}^b + R^b_s e_3 \dot{\xi}$$

$${v}_{eO_c}^b = {v}_{eO_c(b)}^b + R^b_s e_3 \dot{\xi}$$

Where we have defined:

$${v}_{eO_c(b)}^b = {v}_{eO_b}^b + \Omega_{eb}^b r_{O_bO_c}^b$$

The term ${v}_{eO_c(b)}^b$ represents the velocity of the contact point due to vehicle motion, and
$R^b_s e_3 \dot{\xi}$ is its velocity due to strut deformation.

Now:

$$R^c_s = R^c_b R^b_s$$

$${v}_{eO_c(b)}^c = R^c_b {v}_{eO_c(b)}^b$$

$${v}_{eO_c}^c = {v}_{eO_c(b)}^c + R^c_s e_3 \dot{\xi} = {v}_{eO_c(b)}^c + k^c_s \dot{\xi}$$

The non-penetration constraint requires that:

$$e_3^T {v}_{eO_c}^c = 0 = e_3^T {v}_{eO_c(b)}^c + e_3^T k^c_s \dot{\xi} = {v}_{eO_c(b)}^c[3] + \dot{\xi} k^c_s[3]$$

Since $k^c_s[3] = e_3^T R^c_s e_3 = (k_c^e)^T k_s^e = \cos \alpha_{cs}$:

$$\dot{\xi} = \dfrac{-{v}_{eO_c(b)}^c [3]}{k^c_s[3]} = \dfrac{-{v}_{eO_c(b)}^c [3]}{\cos \alpha_{cs}}$$

Once $\dot{\xi}$ is known, we can compute:

$${v}_{eO_c}^c = {v}_{eO_c(b)}^c + k^c_s \dot{\xi}$$


### Maximum Friction Coefficients

We now determine the maximum lateral and longitudinal friction coefficients. These are only valid when the contact point
velocity is non-zero. When this is not the case, the friction coefficients take values between zero and their maximum
values as required to keep the contact point velocity at zero.

We start by computing the tire slip angle from the horizontal components of $v_{eO_c}^c = {R^e_c}^T v_{eO_c}^e$ as:

$\psi_{cv} = atan2(v_{eO_c}^{y_c}, v_{eO_c}^{x_c}) \in [-\pi, \pi)$ (assertion here)

The cornering friction coefficient $\mu_y$ is a function of the following parameters:
- $\mu_{skid}$: Skidding friction coefficient. Typical values: $\mu_{skid} = 0.5$ for dry surfaces, $\mu_{skid} = 0.15$
  for wet surfaces, $\mu_{skid} = 0.05$ for icy surfaces.

- $|\psi_{cv}|$: Absolute value of the tire slip angle. Note that $|\psi_{cv}| \leq \pi$.

- $\psi_{cv(skid)}$: Tire slip angle at which full skidding friction is reached. Typical value: $\psi_{cv(skid)} = 10°$.


The cornering friction coefficient is given by:

- $\mu_{y(max)} = \mu_{skid} \dfrac{|\psi_{cv}|}{\psi_{cv(skid)}} ,$
$\forall |\psi_{cv}| \in [0, \psi_{cv(skid)})$

- $\mu_{y(max)} = \mu_{skid},$
$\forall |\psi_{cv}| \in [\psi_{cv(skid)}, \pi - \psi_{cv(skid)}]$

- $\mu_{y(max)} =\mu_{skid} \left( 1 - \dfrac{ \psi_{cv(skid)}-\left(\pi-|\psi_{cv}| \right)}{\psi_{cv(skid)}} \right) ,$
  $\forall |\psi_{cv}| \in (\pi - \psi_{cv(skid)}, \pi]$


The (maximum) longitudinal friction coefficient $\mu_{x(max)}$ is a function of three parameters:

- $\mu_{roll}$: Free (no braking) rolling friction coefficient for zero slip angle. Typical value: $\mu_{roll} = 0.02$.

- $\mu_{max}$: Maximum braking friction coefficient. Assuming anti-skid is available, its maximum theoretically
  achievable value is the skidding friction coefficient $\mu_{skid}$.

- $u_{brk}$: Braking intensity, with $u_{brk} \in [0, 1]$

Then, $\mu_{x(max)}$ can be computed as:

$$\mu_{x(max)} = \mu_{roll} + \left( \mu_{skid} - \mu_{roll} \right) u_{brk}$$

From the above, we can see that both $\mu_{x(max)}$ and $\mu_{y(max)}$ have $\mu_{skid}$ as their maximum value.
However, the magnitude of the friction coefficient, both components considered, cannot exceed $\mu_{skid}$. This makes
intuitive sense. For example, if the tire slip angle is large enough that the tire is already fully skidding, applying
the brakes should make no difference to the overall friction coefficient. And, if the tire is only slipping, applying
the brakes should not allow the overall friction coefficient to surpass that of pure skidding.

Let:

$$\mu_{mag} = \sqrt{\mu_{x(max)}^2 + \mu_{y(max)}^2}$$

To keep $\mu_{mag} \leq \mu_{skid}$, the friction coefficients are scaled as follows:

$$\mu_{x(max)}:= \mu_{x(max)} \, \min \left( 1, \mu_{skid} / \mu_{mag} \right)$$

$$\mu_{y(max)}:= \mu_{y(max)} \, \min \left( 1, \mu_{skid} / \mu_{mag} \right)$$

### $\mu$ Scaling

In our friction model, the force exerted by the ground on the tire at the contact point is given by:

$$F_{gnd,P}^c = F_N \begin{pmatrix} \mu_x & \mu_y & -1 \end{pmatrix}$$

Where $F_N \ge 0$ is the magnitude of the force along the outward terrain normal, and $\mu_x$ and $\mu_y$ are signed friction coefficients.

As described above, the magnitude and sign of $\mu_x$ and $\mu_y$ must be such that the contact point velocity is kept at zero. This constraint is enforced as long as the friction coefficients remain below their maximum values, computed above. Once a friction coefficient reaches its maximum, it will saturate, and the contact point will start moving along that coefficient's direction. In this new condition, the friction force along that direction will have a constant value, given by the maximum value of that friction coefficient, and will oppose the direction of movement.

This behavior can be simulated by defining the friction coefficients as follows:

$$\mu_{x} = \alpha_x \mu_{x(max)}$$

$$\mu_{y} = \alpha_y \mu_{y(max)}$$

Each $\alpha$ is a signed scaling coefficient, whose value is determined by a PI regulator acting on the contact velocity
along the corresponding axis. The value of $\alpha$ is given by:

$$\alpha_x = -k_p v_{eO_c}^{x_c} - k_i \int{1_A(\alpha_x) v_{eO_c}^{x_c} dt}$$

$$\alpha_y = -k_p v_{eO_c}^{y_c} - k_i \int{1_A(\alpha_y) v_{eO_c}^{y_c} dt}$$

Where $k_p$ and $k_i$ are constant gains, and $1_A(x)$ is the indicator function, with $A = \{-1, 1\}$. Its purpose is
to stop the integration whenever the overall output value of $\alpha$ is saturated at $-1$ or $1$ in order to avoid
windup.

For a moving contact point, the constant (saturated) values will be given by:

$$\alpha_x = -sgn(v_{eO_c}^{x_c})$$

$$\alpha_y = -sgn(v_{eO_c}^{y_c})$$

Since these friction coefficients are signed, the friction force components are given simply by:

$$F_{O_c}^{x_c} = \mu_x F_N$$

$$F_{O_c}^{y_c} = \mu_y F_N$$


### Computing ground forces

For the purposes of the linear momentum equation, let $P$ denote the location of the center of mass of the piston rod
assembly (which includes the tire). This point moves along the strut compression line, that is, along the $z_s$ axis.
The linear momentum equation applied to $P$ is:

$$\dot{v}_{iP}^i = \dfrac{1}{m} F_{ext,P}^i + G^i$$

Expressed in the ECEF frame, the above takes the form:

$$\dot{v}_{eP}^e = \dfrac{1}{m} F_{ext,P}^e + g^e - 2 \Omega_{ie}^e v_{eP}^e$$

Also:

$$\dot{v}_{eP}^e = a_{eO_s}^e + R^e_s (\ddot{r}_{O_sP}^s + (\Omega_{es}^s \Omega_{es}^s  +
\dot{\Omega}_{es}^s )r_{O_sP}^s + \Omega_{es}^s \dot{r}_{sP}^s)$$

Substituting:

$$a_{eO_s}^s + \ddot{r}_{O_sP}^s + (\Omega_{es}^s \Omega_{es}^s +
\dot{\Omega}_{es}^s) r_{O_sP}^s + \Omega_{es}^s \dot{r}_{sP}^s =
\dfrac{1}{m} F_{ext,P}^s + g^s - 2 \Omega_{ie}^s v_{eP}^s$$

We now neglect:
- Inertia terms due to the Earth rotation
- Centrifugal and angular acceleration terms due to the rotation of the strut frame

The result is:

$$a_{eO_s}^s + \ddot{r}_{O_sP}^s = \dfrac{1}{m} F_{ext,P}^s + g^s$$

Now, the position of $P$ can be expressed as:

$$r_{O_sP}^s = r_{O_sO_w}^s + r_{O_wP}^s = {\begin{pmatrix} 0 & 0 & l + d_{O_wP} \end{pmatrix}}^T$$

Where $d_{O_wP}$ is a constant signed distance that locates $P$ with respect to $O_w$. Then:

$$\ddot{r}_{O_sP}^s = {\begin{pmatrix} 0 & 0 & \ddot{\xi} \end{pmatrix}}^T$$

Projecting onto the $z_s$ axis yields:

$$\ddot{\xi} = e_3^T \left(\dfrac{1}{m} F_{ext,P}^s + g^s - a_{eO_s}^s \right)$$

The external forces acting on the piston rod assembly are:
- The force $F_{oleo,P}^s(\xi, \dot{\xi})$ along $z_s$ due to the shock absorber:

  $$F_{oleo,P}^s(\xi, \dot{\xi}) = {\begin{pmatrix} 0 & 0 & F_{oleo,P}^{z_s}(\xi, \dot{\xi})  \end{pmatrix}}^T$$

- The constraint forces along $x_s$ and $y_s$ transmitted by the strut casing:

  $$F_{c,P}^s = {\begin{pmatrix} F_{c,P}^{x_s} & F_{c,P}^{y_s} & 0 \end{pmatrix}}^T$$

- Ground force $F_{gnd,P}^s$ (normal and friction), transmitted through the tire:

  $$F_{gnd,P}^s = R^s_c F_{gnd,P}^c(F_N) = R^s_c F_N {\begin{pmatrix} \mu_x \alpha_x & \mu_y \alpha_y & -1 \end{pmatrix}}^T$$

Substituting:

  $$\ddot{\xi} = \dfrac{1}{m} \left( e_3^T R^s_c F_{gnd,P}^c(F_N) + F_{oleo,P}^{z_s}(\xi, \dot{\xi})\right) + e_3^T \left(g^s - a_{eO_s}^s \right)$$

When there is no ground contact, the ground force is zero and we can solve for $\ddot{\xi}$.

However, when there is ground contact we have a scalar equation with two unknowns: $\ddot{\xi}$ and the normal force
$F_N$. In this case, the first can be obtained from the non-penetration constraint. Formally, we must require that
the acceleration of the contact point $O_c$ along $z_s$ be zero. This will provide an additional scalar equation
which we can solve for $\ddot{\xi}$.

To avoid having to handle the non-penetration constraint explicitly, we can neglect the inertia of the piston assembly
by setting $m \approx 0$. This causes the ground forces to be transmitted instantaneously and directly to
the shock absorber, which is often a reasonable approximation. Multiplying the linear momentum equation by $m$ we have:

$$e_3^T R^s_c F_{gnd,P}^c(F_N) + F_{oleo,P}^{z_s}(\xi, \dot{\xi}) = m \ddot{\xi} - m e_3^T \left(g^s - a_{eO_s}^s \right) = 0$$

Then we have the following scalar equation:

$$
e_3^T R^s_c
\begin{pmatrix}
\mu_x \alpha_x \\
\mu_y \alpha_y \\
-1
\end{pmatrix}
F_N + F_{oleo,P}^{z_s}(\xi, \dot{\xi}) = 0
$$

Let us define the non-dimensional contact force as:

$$
f_{gnd,P}^c =
\begin{pmatrix}
\mu_x \alpha_x \\
\mu_y \alpha_y \\
-1
\end{pmatrix}
$$

We can then write the previous equation as:

$$e_3^T f_{gnd,P}^s F_N + F_{oleo,P}^{z_s}(\xi, \dot{\xi}) = 0$$

Solving for $F_N$:

$$F_N = \dfrac{-F_{oleo,P}^{z_s}(\xi, \dot{\xi})}{f_{gnd,P}^{z_s}}$$

Since the contact constraint is unilateral, we must have $F_N \ge 0$, so we must bound the result accordingly. $F_N > 0$
means that the vertical ground force is negative along $z_c$, as it should.

With $F_N$, we can now compute $F_{gnd,P}^c(F_N)$, which is the ground force transmitted to the vehicle through the tire:

$$F_{gnd,P}^s = R^s_c F_N {\begin{pmatrix} \mu_x \alpha_x & \mu_y \alpha_y & -1 \end{pmatrix}}^T$$
