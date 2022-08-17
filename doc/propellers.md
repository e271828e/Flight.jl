Reference:
Mechanics of Flight (Second Edition) - Warren F. Phillips

The implementation follows Phillips almost verbatim. A few noteworthy considerations follow.

The book implicitly defines a propeller frame ($\varepsilon_p$) with origin $O_p$ at the intersection of the propeller axis with the propeller reference plane. The $x_p$ axis is coincident with the propeller axis, and pointing forward. The orientation of $y_p$ and $z_p$ within the propeller plane are arbitrary.

In the book, the axial aerodynamic torque is by convention positive along the *negative* $x_p$ axis. Therefore:

$$C_{Mx} = C_l = -\pi^2/8 \int ... > 0$$

In the book, the in-plane aerodynamic force $N$ due to $\alpha_p$ is by convention positive along the *negative* $z_p$ axis as defined above. Therefore:

$$C_{Fz, \alpha} = -C_{N, \alpha} = -\pi^2/8 \alpha_p \int ... < 0$$

In contrast, the book considers the yawing moment $n$ positive along $z_p$, so:

$$C_{Mz, \alpha} = C_{n, \alpha} = -\pi^2 / 16 \ \alpha_p \int ... <0$$

By symmetry, it is straightforward to see that $C_{Fy,\beta} = C_{Fz,\alpha}$ and $C_{My,\beta} = C_{Mz,\alpha}$.

All of the above is derived for a CW propeller. By symmetry, we can reason that for a CCW propeller, force coefficients must have the same sign. Traction is still traction. And in-plane forces are caused by one half of the propeller disc seeing larger airflow velocity. In a CCW propeller, this half will be the opposite to that in a CW propeller. However, the net resultant force is the same.

The opposite is true for moment coefficients: they all switch signs when $\omega$ changes from positive to negative. This is obvious for $Mx$. For in-plane moments, if the propeller disc half seeing higher airflow velocities changes from one side of an axis to the other, the moment with respect to that axis will change signs.

Because we need to consider both the CW and CCW cases, we define the power coefficient as:

$$C_P = \dfrac{P}{\rho |\omega / 2\pi|^3 d^5} = \dfrac{M_x \omega}{\rho |\omega / 2\pi|^3 d^5} = \dfrac{\rho (\omega / 2\pi)^2 d^5 C_{Mx} \omega}{\rho |\omega / 2\pi|^3 d^5} = 2\pi \ sgn(\omega)C_{M_x}$$

With our definitions:
- For a CW propeller, $C_{Mx} < 0$ and $\omega > 0$, so $C_{P} <0$
- For a CCW propeller, $C_{Mx} > 0$ and $\omega < 0$, so $C_{P} <0$

That is, $C_P$ is always negative, which is consistent with the fact that the aerodynamic torque produced by the propeller opposes the angular velocity and thus generates negative power ($P = Mx \omega$).

## Compressibility

For axial airflow, the aerodynamic velocity vector at an any blade section is parallel to the section's plane, and its magnitude is:

$$V_b = \sqrt{\omega ^2 r^2 + V_\infty ^2} \cos \varepsilon_i = \omega R \sqrt{\zeta^2 + \left(\dfrac{J}{\pi}\right)^2} \cos \varepsilon_i$$

The Mach number is then:

$$M = \dfrac{\omega R}{a} \sqrt{\zeta^2 + \left(\dfrac{J}{\pi}\right)^2} \cos \varepsilon_i$$

And the Mach number at the blade tip:

$$M_{tip} = \dfrac{\omega R}{a} \sqrt{1 + \left(\dfrac{J}{\pi}\right)^2} \cos \varepsilon_i$$

With these we can write:

$$M = M_{tip} \sqrt{\dfrac{\zeta^2 + \left(\dfrac{J}{\pi}\right)^2}{1 + \left(\dfrac{J}{\pi}\right)^2}} \cos \varepsilon_i = M(\zeta,J,M_{tip},\varepsilon_i)$$

When the airflow is non-axial, the airspeed at a given section (and therefore the Mach number) will also depend on the airflow angle $\alpha_p$ and the section's azimuth $\theta$ within the propeller disc. However, as long as the airflow angle is small, these contributions will be minor and the above expression remains a good approximation.

The propeller's force and torque coefficients for non-axial incompressible airflow were:

$$C_{F_x} = \dfrac{\pi^2}{4} \int_{\zeta_h}^{1} \zeta^2 k \tilde{c} \dfrac{\cos^2 \varepsilon_i}{\cos^2 \varepsilon_\infty} \left[c_{L}(\alpha) \cos(\varepsilon_\infty + \varepsilon_i) - c_{D}(\alpha) \sin(\varepsilon_\infty + \varepsilon_i)\right] d\zeta$$

$$C_{M_x} = -\dfrac{\pi^2}{8} \int_{\zeta_h}^{1} \zeta^3 k \tilde{c} \dfrac{\cos^2 \varepsilon_i}{\cos^2 \varepsilon_\infty}
            \left[c_{D}(\alpha) \cos(\varepsilon_\infty + \varepsilon_i) + c_{L}(\alpha) \sin(\varepsilon_\infty + \varepsilon_i)\right] d\zeta$$

$$C_{M_x} = -\dfrac{\pi^2}{8} \int_{\zeta_h}^{1} \zeta^3 k \tilde{c} \dfrac{\cos^2 \varepsilon_i}{\cos^2 \varepsilon_\infty}
            \left[c_{D}(\alpha) \cos(\varepsilon_\infty + \varepsilon_i) + c_{L}(\alpha) \sin(\varepsilon_\infty + \varepsilon_i)\right] d\zeta$$

etc.

Where $c_{L}(\alpha)$, $c_{D}(\alpha)$ and $c_{L,\alpha}(\alpha)$ were the airfoil's coefficients for the incompressible limit $M \rightarrow 0$.

To account for compressibility, we can replace the incompressible airfoil coefficients with suitable functions $c_L(M, \alpha)$, $c_D(M, \alpha)$ and $c_{L,\alpha}(M, \alpha)$.  The same can be done in the equation for $\varepsilon_i$.

Now the propeller coefficients will be functions of $J$ and $M_{tip}$, so they must be computed and tabulated as such.
