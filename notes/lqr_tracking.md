# LQR Tracking

## Linearization and Set Point

We have the nonlinear system:
$$
\dot{x} = f(x, u)\\
z = h(x, u)
$$

Here, $z$ represents a vector of output variables to be commanded, rather than the system's complete measurement vector, and our goal is to drive it to a reference value $z^*$.

We start by finding an equilibrium point $(x_{trim}, u_{trim})$, so that:
$$ \dot{x}_{trim} = f(x_{trim}, u_{trim}) = 0 $$

We linearize around it:
$$
\dot{x} - \dot{x}_{trim}  \approx \left. \frac{\partial f}{\partial x} \right \rvert_{trim}(x-x_{trim}) + \left. \frac{\partial f}{\partial u} \right \rvert_{trim}(u-u_{trim})
$$
$$
z_{trim}  \approx \left. \frac{\partial h}{\partial x} \right \rvert_{trim}(x-x_{trim}) + \left. \frac{\partial h}{\partial u} \right \rvert_{trim}(u-u_{trim})
$$
$$ \Delta \dot{x} = F \Delta x + G \Delta u $$
$$ \Delta z = H_x \Delta x + H_u \Delta u $$

To achieve $z = z^*$, we need to drive $\Delta z$ to $\Delta z^* = z^* - z_{trim}$ in the linearized system.

Now we seek an equilibrium reference $(\Delta x^*, \Delta u^*)$ such that:
$$ 0 = F \Delta x^* + G \Delta u^* $$
$$ \Delta z^* = H_x \Delta x^* + H_u \Delta u^* $$

In matrix form:
$$
A \begin{pmatrix} \Delta x^* \\ \Delta u^* \end{pmatrix} = \begin{pmatrix} 0 \\ \Delta z^* \end{pmatrix}
$$
$$
\begin{pmatrix} \Delta x^* \\ \Delta u^* \end{pmatrix} = B \begin{pmatrix} 0 \\ \Delta z^* \end{pmatrix}
$$

Where:
$$
A \triangleq \begin{pmatrix} F & G\\ H_x & H_u \end{pmatrix}
$$
$$
B \triangleq A^{-1} = \begin{pmatrix} B_{11} & B_{12}\\ B_{21} & B_{22} \end{pmatrix}
$$
$$ B_{22} = (-H_x F^{-1}G + H_u)^{-1}$$
$$ B_{21} = -B_{22} H_x F^{-1}$$
$$ B_{11} = F^{-1}(-GB_{21} + I)$$
$$ B_{12} = F^{-1}GB_{22}$$

For $A$ to be invertible and $B$ to exist, both $F$ and $(-H_x F^{-1}G + H_u)$ must be invertible. If this is so:
$$ \Delta x^* = B_{12} \Delta z^*$$
$$ \Delta u^* = B_{22} \Delta z^*$$

Subtracting equations:
$$ \Delta \dot{x} = F (\Delta x - \Delta x^*) + G (\Delta u - \Delta u^*)$$
$$ \Delta z - \Delta z^{*} = H_x (\Delta x - \Delta x^*) + H_u (\Delta u - \Delta u^*)$$

If we define:
$$\Delta \tilde{x} \triangleq \Delta x - \Delta x^*$$
$$\Delta \tilde{u} \triangleq \Delta u - \Delta u^*$$
$$\Delta \tilde{z} \triangleq \Delta z - \Delta z^*$$

Then we have:
$$ \Delta \dot{\tilde{x}} = F \Delta \tilde{x} + G \Delta \tilde{u}$$
$$ \Delta \tilde{z} = H_x \Delta \tilde{x} + H_u \Delta \tilde{u}$$

## LQR Tracker

If we design a LQR for this dynamic system, our optimal control law will be:
$$ \Delta \tilde{u} = - C \Delta \tilde{x} $$

The resulting closed-loop is:
$$ \Delta \dot{\tilde{x}} = (F - G C) \Delta \tilde{x}$$

Since this homogeneous system is guaranteed to be asymptotically stable by the LQR, in the steady-state:
$$\Delta \dot{\tilde{x}} = 0$$
$$\Delta \tilde{x} = 0$$
$$\Delta \tilde{u} = - C \Delta \tilde{x} = 0$$

Which means:
$$\dot{x} = 0$$
$$\Delta x = \Delta x^*$$
$$\Delta u = \Delta u^*$$
$$\Delta z = H_x \Delta x^* + H_u \Delta u^* = \Delta z^*$$

The optimal control law is implemented as:
$$\Delta u = \Delta u^* - C\Delta x + \Delta x^*= (B_{22} + CB_{12}) \Delta z^* - C\Delta x = C_{fwd}\Delta z^* - C_{fbk} \Delta x$$
$$u = u_{trim} + C_{fwd}(z^* - z_{trim}) - C_{fbk} (x - x_{trim})$$

## LQR Tracker With Integral Compensation
We choose the subset $z_{i}$ of the command variables for which we want to include integral compensation.
We augment the system dynamics with an integral state vector $\xi$ such that:
$$ \dot{\xi} = \Delta \tilde{z}_{i} = H_{x,i} \Delta \tilde{x} + H_{u,i} \Delta \tilde{u}$$

Where $H_{x,i}$ and $H_{u,i}$ are respectively the row blocks of $H_x$ and $H_u$ corresponding to $z_i$.

This integral state is computed as:
$$ \xi = \xi(0) + \int^t_{0} \Delta \tilde{z}_{i} d\tau = \xi(0) + \int_0^t (\Delta z_i - \Delta z_i^*) d\tau = \xi(0) + \int_0^t (z_i - z_i^*) d\tau $$

The augmented system is
$$
\begin{pmatrix} \Delta \dot{\tilde{x}} \\ \dot{\xi} \end{pmatrix} = \begin{pmatrix} F & 0\\ H_{x,i} & 0 \end{pmatrix} \begin{pmatrix} \Delta \tilde{x} \\ \xi\end{pmatrix} + \begin{pmatrix} G \\ H_{u,i} \end{pmatrix} \Delta \tilde{u}
$$
$$
 \Delta {\tilde{z}} = \begin{pmatrix} H_x & 0 \end{pmatrix} \begin{pmatrix} \Delta \tilde{x} \\ \xi\end{pmatrix} + H_u \Delta \tilde{u}
$$

Defining:
$$ F_{aug} \triangleq \begin{pmatrix} F & 0\\ H_{x,i} & 0 \end{pmatrix}$$
$$ G_{aug} \triangleq \begin{pmatrix} G \\ H_{u,i} \end{pmatrix} $$
$$ x_{aug} \triangleq \begin{pmatrix} \Delta \tilde{x} \\ \xi \end{pmatrix} $$

We can write:
$$\dot{x}_{aug} = F_{aug} x_{aug} + G_{aug} \Delta \tilde{u}$$
$$
 \Delta {\tilde{z}} = \begin{pmatrix} H_x & 0 \end{pmatrix} x_{aug} + H_u \Delta \tilde{u}
$$

Now, if we design a LQR for the augmented system:
$$\Delta \tilde{u} = -C_{aug} x_{aug} = \begin{pmatrix} C_x & C_{\xi}\end{pmatrix} x_{aug}$$
$$\dot{x}_{aug} = (F_{aug} - G_{aug}C)x_{aug}$$

This system is guaranteed to be stable by the LQR, so in the steady-state we have:
$$\dot{x}_{aug}=0$$

Which means:
$$\Delta \dot{\tilde{x}} = 0$$
$$\dot{\xi} = \Delta \tilde{z}_i =0 = z_i - z_i^*$$

Let's assume we have an external disturbance such that the actual system dynamics are instead:
$$
\begin{pmatrix} \Delta \dot{\tilde{x}} \\ \dot{\xi} \end{pmatrix} = \begin{pmatrix} F & 0\\ H_{x,i} & 0 \end{pmatrix} \begin{pmatrix} \Delta \tilde{x} \\ \xi\end{pmatrix} + \begin{pmatrix} G \\ H_{u,i} \end{pmatrix} \Delta \tilde{u} + \begin{pmatrix} L w \\ 0 \end{pmatrix}
$$

The closed-loop system is now:
$$\dot{x}_{aug} = (F_{aug} - G_{aug}C_{aug})x_{aug} + \begin{pmatrix} L w \\ 0 \end{pmatrix}$$

Stability depends only on the dynamics matrix $F_{aug} - G_{aug}C_{aug}$, so this system will still be asymptotically stable. This means that in the steady state $\dot{x}_{aug} = 0$, so $\dot{\xi} = 0$, and therefore necessarily $z = z*$. Because the system is no longer homogeneous, in general we will have $x_{aug} \neq 0$. In particular, $\xi$ will converge to the values required by $\Delta \tilde{z}_i = 0$.

The optimal control law is implemented as:
$$\Delta u = \Delta u^* - C_{aug}x_{aug} = \Delta u^* - C_x \Delta \tilde{x} - C_{\xi} \xi = \Delta u^* - C_x \Delta x + C_x \Delta x^* - C_{\xi} \xi = (B_{22} + C_x B_{12}) \Delta z^* - C_x \Delta x - C_{\xi} \xi$$
$$\Delta u = C_{fwd} \Delta z^* - C_{fbk} \Delta x - C_{\xi} \xi$$
$$u = u_{trim} + C_{fwd}(z^* - z_{trim}) - C_{fbk} (x - x_{trim}) - C_{\xi} \xi$$

With:
$$C_{fwd} = B_{22} + C_x B_{12}$$
$$C_{fbk} = C_x$$
$$ \xi = \xi(0) + \int_0^t (z_i - z_i^*) d\tau $$

For a practical implementation, we write:
$$u = u_{trim} + C_{fwd}(z^* - z_{trim}) - C_{fbk} (x - x_{trim}) + u_{int}$$

Where:
$$u_{int} = -C_{\xi} \xi(0) - \int_0^t C_{\xi} (z_i - z_i^*) d\tau = u_{int}(0)  - \int_0^t C_{\xi} (z_i - z_i^*) d\tau = u_{int}(0) - \int_0^t C_{int} (z - z^*) d\tau$$

That is, we move $C_\xi$ inside the integrator, which will now have dimension $n_u$ instead of $n_i$. The practical advantage is that having an integration path per control input allows to handle saturation of each control input independently. We also define a $C_{int}$ is a $n_u \times n_z$ matrix with all rows set to zero except for those $n_i$ rows corresponding to the integrated command variables, which will be taken from $C_{\xi}$. This allows for a more general LQR tracker implementation in which the integral gain matrix always has size $n_u \times n_z$ and always multiplies the complete command vector.