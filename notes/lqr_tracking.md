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
$$ \Delta \dot{x} = A \Delta x + B \Delta u $$
$$ \Delta z = C \Delta x + D \Delta u $$

To achieve $z = z^*$, we need to drive $\Delta z$ to $\Delta z^* = z^* - z_{trim}$ in the linearized system.

Now we seek an equilibrium reference $(\Delta x^*, \Delta u^*)$ such that:
$$ 0 = A \Delta x^* + B \Delta u^* $$
$$ \Delta z^* = C \Delta x^* + D \Delta u^* $$

In matrix form:
$$
L \begin{pmatrix} \Delta x^* \\ \Delta u^* \end{pmatrix} = \begin{pmatrix} 0 \\ \Delta z^* \end{pmatrix}
$$
$$
\begin{pmatrix} \Delta x^* \\ \Delta u^* \end{pmatrix} = M \begin{pmatrix} 0 \\ \Delta z^* \end{pmatrix}
$$

Where:
$$
L \triangleq \begin{pmatrix} A & B\\ C & D \end{pmatrix}
$$
$$
M \triangleq L^{-1} = \begin{pmatrix} M_{11} & M_{12}\\ M_{21} & M_{22} \end{pmatrix}
$$
$$ M_{22} = (-C A^{-1}B + D)^{-1}$$
$$ M_{21} = -M_{22} C A^{-1}$$
$$ M_{11} = A^{-1}(-BM_{21} + I)$$
$$ M_{12} = A^{-1}BM_{22}$$

For $L$ to be invertible and $M$ to exist, both $A$ and $(-C A^{-1}B + D)$ must be invertible. If this is so:
$$ \Delta x^* = M_{12} \Delta z^*$$
$$ \Delta u^* = M_{22} \Delta z^*$$

Subtracting equations:
$$ \Delta \dot{x} = A (\Delta x - \Delta x^*) + B (\Delta u - \Delta u^*)$$
$$ \Delta z - \Delta z^{*} = C (\Delta x - \Delta x^*) + D (\Delta u - \Delta u^*)$$

If we define:
$$\Delta \tilde{x} \triangleq \Delta x - \Delta x^*$$
$$\Delta \tilde{u} \triangleq \Delta u - \Delta u^*$$
$$\Delta \tilde{z} \triangleq \Delta z - \Delta z^*$$

Then we have:
$$ \Delta \dot{\tilde{x}} = A \Delta \tilde{x} + B \Delta \tilde{u}$$
$$ \Delta \tilde{z} = C \Delta \tilde{x} + D \Delta \tilde{u}$$

## LQR Tracker

If we design a LQR for this dynamic system, our optimal control law will be:
$$ \Delta \tilde{u} = - K \Delta \tilde{x} $$

The resulting closed-loop is:
$$ \Delta \dot{\tilde{x}} = (A - B K) \Delta \tilde{x}$$

Since this homogeneous system is guaranteed to be asymptotically stable by the LQR, in the steady-state:
$$\Delta \dot{\tilde{x}} = 0$$
$$\Delta \tilde{x} = 0$$
$$\Delta \tilde{u} = - K \Delta \tilde{x} = 0$$

Which means:
$$\dot{x} = 0$$
$$\Delta x = \Delta x^*$$
$$\Delta u = \Delta u^*$$
$$\Delta z = C \Delta x^* + D \Delta u^* = \Delta z^*$$

The optimal control law is implemented as:
$$\Delta u = \Delta u^* + \Delta \tilde{u} = \Delta u^* - K\Delta \tilde{x} = \Delta u^* - K \Delta
x + K \Delta x^* = M_{22} \Delta z^* - K \Delta x + K M_{12} \Delta z^* = (M_{22} + KM_{12}) \Delta z^* - K\Delta x$$
$$\Delta u = (M_{22} + KM_{12}) \Delta z^* - K\Delta x
= K_{fwd}\Delta z^* - K_{fbk} \Delta x$$
$$u = u_{trim} + K_{fwd}(z^* - z_{trim}) - K_{fbk} (x - x_{trim})$$

## LQR Tracker With Integral Compensation
We choose the subset $z_{i}$ of command variables for which we want to include integral compensation.
We augment the system dynamics with an integral state vector $\xi$ such that:
$$ \dot{\xi} = \Delta \tilde{z}_{i} = C_{i} \Delta \tilde{x} + D_{i} \Delta \tilde{u}$$

Where $C_{i}$ and $D_{i}$ are respectively the row blocks of $C$ and $D$ corresponding to $z_i$.

This integral state is computed as:
$$ \xi = \xi(0) + \int^t_{0} \Delta \tilde{z}_{i} d\tau = \xi(0) + \int_0^t (\Delta z_i - \Delta z_i^*) d\tau = \xi(0) + \int_0^t (z_i - z_i^*) d\tau $$

The augmented system is
$$
\begin{pmatrix} \Delta \dot{\tilde{x}} \\ \dot{\xi} \end{pmatrix} = \begin{pmatrix} A & 0\\ C_{i} & 0 \end{pmatrix} \begin{pmatrix} \Delta \tilde{x} \\ \xi\end{pmatrix} + \begin{pmatrix} B \\ D_{i} \end{pmatrix} \Delta \tilde{u}
$$
$$
 \Delta {\tilde{z}} = \begin{pmatrix} C & 0 \end{pmatrix} \begin{pmatrix} \Delta \tilde{x} \\ \xi\end{pmatrix} + D \Delta \tilde{u}
$$

Defining:
$$ A_{aug} \triangleq \begin{pmatrix} A & 0\\ C_{i} & 0 \end{pmatrix}$$
$$ B_{aug} \triangleq \begin{pmatrix} B \\ D_{i} \end{pmatrix} $$
$$ \Delta \tilde{x}_{aug} \triangleq \begin{pmatrix} \Delta \tilde{x} \\ \xi \end{pmatrix} $$

We can write:
$$\Delta \dot{\tilde{x}}_{aug} = A_{aug} \Delta{\tilde{x}}_{aug} + B_{aug} \Delta \tilde{u}$$
$$
 \Delta {\tilde{z}} = \begin{pmatrix} C & 0 \end{pmatrix} \Delta{\tilde{x}}_{aug} + D \Delta \tilde{u}
$$

Now, if we design a LQR for the augmented system, we have the control law:
$$\Delta \tilde{u} = -K_{aug} \Delta{\tilde{x}}_{aug} = \begin{pmatrix} K_x & K_{\xi}\end{pmatrix} \Delta{\tilde{x}}_{aug}$$

And the resulting closed-loop system is:
$$\Delta \dot{\tilde{x}}_{aug} = (A_{aug} - B_{aug}K)\Delta{\tilde{x}}_{aug}$$

This system is guaranteed to be stable by the LQR, so in the steady-state we have:
$$\Delta \dot{\tilde{x}}_{aug}=0$$

Which means:
$$\Delta \dot{\tilde{x}} = 0$$
$$\dot{\xi} = \Delta \tilde{z}_i =0 = z_i - z_i^*$$

Let's assume we have an external disturbance such that the actual system dynamics are instead:
$$
\begin{pmatrix} \Delta \dot{\tilde{x}} \\ \dot{\xi} \end{pmatrix} = \begin{pmatrix} A & 0\\ C_{i} & 0 \end{pmatrix} \begin{pmatrix} \Delta \tilde{x} \\ \xi\end{pmatrix} + \begin{pmatrix} B \\ D_{i} \end{pmatrix} \Delta \tilde{u} + \begin{pmatrix} L w \\ 0 \end{pmatrix}
$$

The closed-loop system is now:
$$\Delta \dot{\tilde{x}}_{aug} = (A_{aug} - B_{aug}K_{aug}) \Delta \tilde{x}_{aug} + \begin{pmatrix} L w \\ 0 \end{pmatrix}$$

Stability depends only on the dynamics matrix $A_{aug} - B_{aug}K_{aug}$, so this system will still
be asymptotically stable. This means that in the steady state $\Delta \dot{\tilde{x}}_{aug} = 0$, so $\dot{\xi} =
0$, and therefore necessarily $z_i = z_i^*$. Because the system is no longer homogeneous, in general we
will have $\Delta \tilde{x}_{aug} \neq 0$. In particular, $\xi$ will converge to the values required by $\Delta
\tilde{z}_i = 0$.

The optimal control law is implemented as:
$$\Delta u = \Delta u^* + \Delta \tilde{u}= \Delta u^* - K_{aug}\Delta \tilde{x}_{aug} = \Delta u^*
- K_x \Delta \tilde{x} - K_{\xi} \xi = \Delta u^* - K_x \Delta x + K_x \Delta x^* - K_{\xi} \xi $$
$$\Delta u = M_{22} \Delta z^* - K_x \Delta x + K_x M_{12} \Delta z^* - K_{\xi} \xi= (M_{22} + K_x M_{12}) \Delta z^* - K_x \Delta x - K_{\xi} \xi$$$
$$\Delta u = K_{fwd} \Delta z^* - K_{fbk} \Delta x - K_{\xi} \xi$$
$$u = u_{trim} + K_{fwd}(z^* - z_{trim}) - K_{fbk} (x - x_{trim}) - K_{\xi} \xi$$

With:
$$K_{fwd} = M_{22} + K_x M_{12}$$
$$K_{fbk} = K_x$$
$$ \xi = \xi(0) + \int_0^t (z_i - z_i^*) d\tau $$

For a practical implementation, we write:
$$u = u_{trim} + K_{fwd}(z^* - z_{trim}) - K_{fbk} (x - x_{trim}) + u_{int}$$

Where:
$$u_{int} = -K_{\xi} \xi(0) - \int_0^t K_{\xi} (z_i - z_i^*) d\tau = u_{int}(0)  - \int_0^t K_{\xi} (z_i - z_i^*) d\tau = u_{int}(0) - \int_0^t K_{int} (z - z^*) d\tau$$

That is, we have moved the $n_u \times n_i$ gain matrix $K_\xi$ inside the integrator, which will
now have dimension $n_u$ instead of $n_i$. Having an integration path per control input allows us to
handle saturation independently on each control input by halting its specific integrator. Finally,
we have defined $K_{int}$ as a $n_u \times n_z$ matrix with all columns set to zero, except for
those $n_i$ columns corresponding to the integrated command variables, which are taken from
$K_{\xi}$. This allows for a more general implementation wherein the integral gain matrix always has
size $n_u \times n_z$ and multiplies the complete command vector.