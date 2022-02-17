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