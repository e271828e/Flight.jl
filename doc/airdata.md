The stability axes are defined by a rotation $-\alpha$ around $y_b$. Thus:

$$R^{b}_{s} = \begin{pmatrix} \cos\alpha & 0 & -\sin\alpha \\
                                0        & 1 &          0 \\
                                \sin\alpha & 0 & \cos\alpha \end{pmatrix}$$

The wind axes are defined by a rotation $\beta$ around $z_s$. Thus:

$$R^{s}_{w} = \begin{pmatrix} \cos\beta & -\sin\beta & 0 \\
                                \sin\beta & \cos\beta & 0 \\
                                0 & 0 & 1 \end{pmatrix}$$

With this:

$$R^{b}_{w} = \begin{pmatrix} \cos\alpha \cos\beta & -\cos\alpha\sin\beta & -\sin\alpha \\
                                \sin\beta & \cos\beta & 0 \\
                                \sin\alpha \cos\beta & -\sin\alpha\sin\beta & \cos\alpha
                                \end{pmatrix}$$

The aerodynamic velocity components in body axes are then:

$$v_{wOb}^{b} = R^{b}_{w} v_{wOb}^{w} = R^{b}_{w} \begin{pmatrix} 1\\0\\0\end{pmatrix} =
    \begin{pmatrix} \cos\alpha \cos\beta \\ \sin\beta \\ \sin\alpha \cos\beta \end{pmatrix}$$

And the airflow angles with respect to the body axes can be computed as:

$$\alpha_b = atan2(v_{wOb}^{b}(3), v_{wOb}^{b}(1))$$
$$\beta_b = atan2\left(v_{wOb}^{b}(2), \sqrt{v_{wOb}^{b}(1)^2 + v_{wOb}^{b}(3)^2}\right)$$