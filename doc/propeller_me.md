Define:

- Propeller frame ($\varepsilon_p$): Origin $O_p$ at the intersection of the propeller axis with the propeller reference plane. The $x_p$ axis is coincident with the propeller axis, and pointing forward. Axes $y_p$ and $z_p$ are arbitrarily defined with respect to the vehicle axes. Both $O_p$ and $p$ must be either fixed to the vehicle or parametrically (slowly) changing with respect to it.
- Blade section frame ($\varepsilon_s$): This frame is attached to a propeller blade section. The $x_s$ axis is parallel to $x_p$, and the $y_s$ axis points radially outward along the blade airfoil center line. Thus, axes $s$ are defined by an elemental rotation of $p$ around $x_p$ by a (time-varying) angle we will denote $\phi$, which defines the azimuthal location of the blade. The  origin $O_s$ is at a blade section located by a (fixed) radial coordinate $r$ measured from the propeller axis along the blade, so that $r_{O_pO_s}^s = [0,r,0]^T$. The definition of the section frame is the same for a CW and a CCW propeller.
- Airfoil frame ($\varepsilon_a)$: Origin $O_a$ coincident with $O_s$. The $x_a$ axis points along the circumferential direction within the propeller plane and towards the nominal advance direction of the airfoil. The $z_a$ axis is perpendicular to the propeller plane, pointing downwards towards the pressure surface of the airfoil. And $y_a$ as required. This means:

    a) For a CW propeller, a rotation $\theta = -\pi/2$ around $y_s$.
$$ R^s_a = R_x(-\pi/2) = \begin{pmatrix} 0 & 0 & -1 \\
                                       0 & 1 & 0 \\
                                        1 & 0 & 0 \end{pmatrix}$$
b) For a CCW propeller, a rotation $\psi = \pi$ around $z_s$ to invert $y_s$, then a $\theta = \pi/2$ rotation around $y_s$.
$$ R^s_a = R_z(\pi)R_x(-\pi/2) =
    \begin{pmatrix} -1 & 0 & 0 \\
                    0 & -1 & 0 \\
                    0 & 0 & 1 \end{pmatrix}
    \begin{pmatrix} 0 & 0 & 1 \\
                    0 & 1 & 0 \\
                    -1 & 0 & 0 \end{pmatrix} =
    \begin{pmatrix} 0 & 0 & -1 \\
                    0 & -1 & 0 \\
                    -1 & 0 & 0 \end{pmatrix}$$

Our goal is to compute the average aerodynamic wrench produced by the propeller at $O_j$. As it sweeps its complete azimuth range $\phi \in [0,2\pi]$, the propeller blade section located at a given $r$ will produce be two 1D distributions of average force and moment per unit length, $f_{[O_p]}^p(\phi)$ and $m_{[O_p]}^p(\phi)$. When considering all blade sections along the blade radius, we get two 2D distributions, $f_{[O_p]}^p(r, \phi)$ and $m_{[O_p]}^p(r, \phi)$. These are then integrated over $r\in[r_h, R]$ and averaged over $\phi \in [0,2\pi]$ to yield the average aerodynamic wrench at $O_p$.
for a single blade. This is finally multiplied by the number of blades to yield the total
aerodynamic wrench.

$$dM_{[O_p]}^p=m_{[O_p]}^p(r,\phi) dr d\phi$$
$$dF_{[O_p]}^p=f_{[O_p]}^p(r,\phi) dr d\phi$$

The first step is to compute the aerodynamic velocity at the origin $O_s$ of an arbitrary section determined by its coordinates $(r, \phi)$, resolved in section axes $s$.

We start from:
$$v_{eO_p}^e = \dot{r}_{O_eO_p}^e = \dot{r}_{O_eO_b}^e + \dot{r}_{O_bO_p}^e =
    v_{eO_b}^e + R^e_b \left(\dot{r}_{O_bO_p}^b + \Omega_{eb}^b r_{O_bO_p}^b\right) = \\
    v_{eO_b}^e + R^e_b \Omega_{eb}^b r_{O_bO_p}^b $$

$$ v_{eO_p}^b = v_{eO_b}^b + \Omega_{eb}^b r_{O_bO_p}^b $$

The aerodynamic velocity is computed by subtracting the local wind velocity with respect to the ECEF frame:
$$ v_{wO_p}^b = v_{eO_p}^b - v_{ew}^b(O_p) = v_{eO_b}^b + \Omega_{eb}^b r_{O_bO_p}^b - v_{ew}^b(O_p) \approx \\
v_{eO_b}^b + \Omega_{eb}^b r_{O_bO_p}^b - v_{ew}^b(O_b) = v_{wO_b}^b + \Omega_{eb}^b r_{O_bO_p}^b $$

Now:
$$v_{eO_s}^e = \dot{r}_{O_eO_s}^e = \dot{r}_{O_eO_p}^e + \dot{r}_{O_pO_s}^e =
    v_{eO_p}^e +  R^e_s \left(\dot{r}_{O_pO_s}^s + \Omega_{es}^s r_{O_pO_s}^s\right) = \\
    v_{eO_p}^e +  R^e_s \Omega_{es}^s r_{O_pO_s}^s
    $$

$$v_{eO_s}^s = v_{eO_p}^s + \Omega_{es}^s r_{O_pO_s}^s $$

$$v_{eO_s}^s = \left(R^b_s\right)^T v_{eO_p}^b + \Omega_{es}^s r_{O_pO_s}^s $$

$$v_{wO_s}^s = v_{eO_s}^s - v_{ew}^s(O_s) = \left(R^b_s\right)^T v_{eO_p}^b + \Omega_{es}^s
r_{O_pO_s}^s - v_{ew}^s(O_s) =\\
\left(R^b_s\right)^T \left(v_{eO_p}^b - v_{ew}^b(O_s) \right) + \Omega_{es}^s r_{O_pO_s}^s \approx \\
\left(R^b_s\right)^T \left(v_{eO_p}^b - v_{ew}^b(O_p) \right) + \Omega_{es}^s r_{O_pO_s}^s = \\
\left(R^b_s\right)^T v_{wO_p}^b + \Omega_{es}^s r_{O_pO_s}^s$$


Now:
$$\omega_{es}^s = \omega_{eb}^s + \omega_{bc}^s+ \omega_{cs}^s \approx \omega_{cs}^s $$

Thus:
$$v_{wO_s}^s \approx \left(R^b_s\right)^T v_{wO_p}^b + \Omega_{cs}^s r_{O_pO_s}^s = \left(R^p_s\right)^T v_{wO_p}^p + \Omega_{cs}^s r_{O_pO_s}^s $$
$$v_{wO_a}^a = v_{wO_s}^a = \left(R^s_a\right)^T \left( \left(R^p_s\right)^T v_{wO_p}^p + \Omega_{cs}^s r_{O_pO_s}^s \right)$$

Note that:
$$\omega_{cs}^s = \begin{bmatrix} \omega & 0 & 0\end{bmatrix}^T$$
$$r_{O_pO_s}^s = \begin{bmatrix} 0 & r & 0\end{bmatrix}^T$$

Where $\omega > 0$ for a CW propeller and $\omega <0$ for a CCW propeller.

With this:
$$\Omega_{cs}^s r_{O_pO_s}^s = [0,0,r\omega]$$
$$v_{wO_a}^a = \left(R^s_a\right)^T \left( \left(R^p_s\right)^T v_{wO_p}^p + \omega r \ e_3 \right)$$

Now we define a non-dimensional radial coordinate:

$$\zeta = \dfrac{r}{R}$$

Where $R$ is the propeller radius.

And the following non-dimensional velocity vectors:

$$\bar{v}_{wO_p}^{p} = \dfrac{1 }{|\omega| R} v_{wO_p}^{p}$$
$$\bar{v}_{wO_a}^{a} = \dfrac{1 }{|\omega| R}v_{wO_a}^{a}$$

Using these definitions in the previous equation:

$$\bar{v}_{wO_a}^{a} = \left(R^s_a\right)^T \left( \left(R^p_s\right)^T \bar{v}_{wO_p}^{p} + sgn(\omega) \ \zeta \ e_3 \right)$$

Now, for a given azimuth $\phi$ we have:
$$ R^p_s = R_x(\phi) = \begin{pmatrix} 1 & 0 & 0 \\
                                       0 & \cos\phi & -\sin\phi \\
                                        0 & \sin\phi & \cos\phi \end{pmatrix}$$

So:

$$\left(R^p_s\right)^T \bar{v}_{wO_b}^p + sgn(\omega) \ \zeta \ e_3 =
\begin{pmatrix} 1 & 0 & 0 \\
                0 & \cos\phi & \sin\phi \\
                0 & -\sin\phi & \cos\phi \end{pmatrix}


\begin{pmatrix} \bar{v}_{wO_p}^{x_p} \\
                \bar{v}_{wO_p}^{y_p} \\
                \bar{v}_{wO_p}^{z_p} \end{pmatrix} +

\begin{pmatrix} 0 \\
                0 \\
                sgn(\omega) \ \zeta \end{pmatrix} =

\begin{pmatrix} \bar{v}_{wO_p}^{x_p} \\
                \bar{v}_{wO_p}^{y_p} \cos\phi + \bar{v}_{wO_p}^{z_p} \sin \phi \\
                -\bar{v}_{wO_p}^{y_p} \sin\phi + \bar{v}_{wO_p}^{z_p} \cos\phi + sgn(\omega) \ \zeta \end{pmatrix}$$


<!-- And the vector advance ratio, whose components represent distances travelled by the propeller along each of its axes in one revolution, divided by the propeller diameter:

$$\bar{v} = \dfrac{v_{wO_p}^{p} T}{2R} = \dfrac{v_{wO_p}^{p}}{2 R} \dfrac{2\pi}{\omega} = \dfrac{\pi }{\omega R}v_{wO_p}^{p}$$

With this we can write:

$$ v_{wO_s}^s = \dfrac{\omega R}{\pi}
\begin{pmatrix} \bar{v}_x \\
                \bar{v}_y \cos\phi + \bar{v}_z \sin \phi \\
                -\bar{v}_y \sin\phi + \bar{v}_z \cos\phi + \pi \zeta\end{pmatrix} =
\dfrac{\omega R}{\pi}f(\bar{v},\phi,\zeta)$$ -->

For a CW propeller, $sgn(\omega) = 1$ and:
$$ (R^s_a)^T = \begin{pmatrix} 0 & 0 & 1 \\
                                0 & 1 & 0 \\
                                -1 & 0 & 0 \end{pmatrix}$$

Using these values yields:

$$
\bar{v}_{wO_a}^{a} =
\begin{pmatrix}
    -\bar{v}_{wO_p}^{y_p} \sin\phi + \bar{v}_{wO_p}^{z_p} \cos\phi + sgn(\omega) \ \zeta \\
    \bar{v}_{wO_p}^{y_p} \cos\phi + \bar{v}_{wO_p}^{z_p} \sin \phi \\
    -\bar{v}_{wO_p}^{x_p} \\
\end{pmatrix} =
\begin{pmatrix}
    -\bar{v}_{wO_p}^{y_p} \sin\phi + \bar{v}_{wO_p}^{z_p} \cos\phi + \zeta \\
    \bar{v}_{wO_p}^{y_p} \cos\phi + \bar{v}_{wO_p}^{z_p} \sin \phi \\
    -\bar{v}_{wO_p}^{x_p} \\
\end{pmatrix}
$$

For a CCW propeller, $sgn(\omega) = -1$ and:
$$ (R^s_a)^T =
    \begin{pmatrix} 0 & 0 & -1 \\
                    0 & -1 & 0 \\
                    -1 & 0 & 0 \end{pmatrix}$$

Using these values yields:

$$
\bar{v}_{wO_a}^{a} = \begin{pmatrix}
    \bar{v}_{wO_p}^{y_p} \sin\phi - \bar{v}_{wO_p}^{z_p} \cos\phi - sgn(\omega) \ \zeta \\
    -\bar{v}_{wO_p}^{y_p} \cos\phi - \bar{v}_{wO_p}^{z_p} \sin \phi \\
    -\bar{v}_{wO_p}^{x_p} \\
\end{pmatrix} =
\begin{pmatrix}
    \bar{v}_{wO_p}^{y_p} \sin\phi - \bar{v}_{wO_p}^{z_p} \cos\phi + \zeta \\
    -\bar{v}_{wO_p}^{y_p} \cos\phi - \bar{v}_{wO_p}^{z_p} \sin \phi \\
    -\bar{v}_{wO_p}^{x_p} \\
\end{pmatrix}
$$

Note that:
- $\zeta \in [\zeta_h, 1]$, where $\zeta_h$ corresponds to the propeller hub radius and is typically larger than $0.1$. Therefore, $\zeta \sim 1$, while generally $\bar{v}_{wO_p}^{y_p} \ll 1$ and $\bar{v}_{wO_p}^{z_p} \ll 1$. So we will generally have, $\bar{v}_{wO_a}^{x_a} \gt 1$, as we would expect. However, very slow propeller rotation and large propeller off-axis velocities may result in $\bar{v}_{wO_a}^{x_a} \le 1$, particularly near the hub. This corresponds to a locally receding airflow around the airfoil and should not be allowed in the numerical computations.
- For a forward moving propeller, ${v}_{wO_p}^{x_p} \ge 0$

Pero al final como esto lo vamos a usar para generar tablas, es decir, no hay que resolver en tiempo real. Lo que se puede hacer es limitar los valores de vbar que vamos a tabular, de manera que su suma nunca exceda el valor de zeta_h de nuestra helice. Por ejemplo, limitar ambos valores a 0.25zeta_h. Si zeta_h = 0.1, tenemos limitadas las componentes de v_bar a 0.025. A 500 rpm tenemos 52 rad/s. v_y = omega R * 0.025. Si R=1, queda que por encima de 1.3 m/s estamos fuera de tablas (y hay que saturar). A 1500 rpm que es un valor bajo pero posible, tendriamos unos 4m/s. Si estamos volando a 40m/s, eso es unos 5.5 deg de AoA. Es bastante restrictivo.

Otra opcion es simplemente lidiar con las secciones problematicas. Primero, imponer que vza sea menor o igual que cero (es decir, que sea negativa en za, es decir, que vaya hacia arriba, es decir, que tienda a reducir el angulo de ataque, como corresponde). Por ejemplo, limitar inferiormente \bar{v}_x_a a un valor positivo y estrictamente mayor que 5 \bar{v}

OJO: EN PHI PROMEDIAMOS, NO INTEGRAMOS SIN MAS. NO OLVIDAR DIVIDIR POR 2PI

OJO: SER CONSISTENTE CON MIS SIGNOS EN LOS DIBUJOS QUE HAGA Y EN LOS DESARROLLOS. Phillips usa velocidades de la corriente respecto del perfil, yo hasta aqui he trabajado con velocidades del perfil respecto de la corriente. Lo mejor va a ser volver a deducirme todo por mi cuenta, con mi criterio de signos, es decir, invirtiendo el sentido de todos los vectores del dibujo

With $v_{sO_a}^a$, the 2D airfoil velocity x and z components ($\omega r$ and $V_{\infty}$ respectively for Phillips) are known. We can now switch on the aerodynamic mindset (fixed airfoil, incident airflow), and naturally introduce the induced velocity following Phillips.

At the end of this process the result will be a force per unit length applied at the airfoil frame origin, $f_{[O_a]}^a$. We switch the components to project it back to section axes, obtaining $f_{[O_a]}^s = f_{[O_s]}^s$. From it we define a Wrench, construct the frame transform from $\varepsilon_p$ to $\varepsilon_s$ and transfer the Wrench back to $\varepsilon_p$. This could also be done manually as $m_{[O_p]}^s = r_{O_pO_s}^s \times f_{[O_s]}^s$. And then we rotate back to $c$ axes using the 2D rotation $R^p_s$.

We should pre-allocate six 2D $n_r \times n_\phi$ arrays to hold the distribution corresponding to each component of $f_{[O_p]}^p$ and $m_{[O_p]}^p$.

The inner loop corresponds to $r$, the second to $\phi$. Since Julia has column-major arrays, $r$ should grow along columns and $\phi$ along rows. Our 2D arrays should be accessed as $f_x[r,\phi]$

<!-- A final note: we may want to define a non-dimensional radial coordinate $\zeta \in [0,1]$ such that:
$$r = r_h + (R - r_h)\zeta  $$

With this, the integral along the radial direction is from 0 to 1. -->