## Linearization and Set Point

We have the nonlinear system:
$$
\dot{x} = f(x, u)\\
y = h(x, u)
$$

Here, $y$ represents a command vector rather than the complete measurement vector, and our goal is to drive it to a reference value $y^*$.

We start by finding an equilibrium point $(x_{trim}, u_{trim})$, so that:
$$ \dot{x}_{trim} = f(x_{trim}, u_{trim}) = 0 $$

We linearize around it:
$$
\dot{x} - \dot{x}_{trim}  \approx \left. \frac{\partial f}{\partial x} \right \rvert_{trim}(x-x_{trim}) + \left. \frac{\partial f}{\partial u} \right \rvert_{trim}(u-u_{trim})
$$
$$
y_{trim}  \approx \left. \frac{\partial h}{\partial x} \right \rvert_{trim}(x-x_{trim}) + \left. \frac{\partial h}{\partial u} \right \rvert_{trim}(u-u_{trim})
$$
$$ \Delta \dot{x} = F \Delta x + G \Delta u $$
$$ \Delta y = H_x \Delta x + H_u \Delta u $$

To achieve $y = y^*$, we need to drive $\Delta y$ to $\Delta y^* = y^* - y_{trim}$ in the linearized system.

Now we seek an equilibrium setpoint $(\Delta x^*, \Delta u^*)$ such that:
$$ 0 = F \Delta x^* + G \Delta u^* $$
$$ \Delta y^* = H_x \Delta x^* + H_u \Delta u^* $$

In matrix form:
$$
A \begin{pmatrix} \Delta x^* \\ \Delta u^* \end{pmatrix} = \begin{pmatrix} 0 \\ \Delta y^* \end{pmatrix}
$$
$$
\begin{pmatrix} \Delta x^* \\ \Delta u^* \end{pmatrix} = B \begin{pmatrix} 0 \\ \Delta y^* \end{pmatrix}
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
$$ \Delta x^* = B_{12} \Delta y^*$$
$$ \Delta u^* = B_{22} \Delta y^*$$

Subtracting equations:
$$ \Delta \dot{x} = F (\Delta x - \Delta x^*) + G (\Delta u - \Delta u^*)$$
$$ \Delta y - \Delta y^{*} = H_x (\Delta x - \Delta x^*) + H_u (\Delta u - \Delta u^*)$$

If we define:
$$\Delta \tilde{x} \triangleq \Delta x - \Delta x^*$$
$$\Delta \tilde{u} \triangleq \Delta u - \Delta u^*$$
$$\Delta \tilde{y} \triangleq \Delta y - \Delta y^*$$

Then we have:
$$ \Delta \dot{\tilde{x}} = F \Delta \tilde{x} + G \Delta \tilde{u}$$
$$ \Delta \tilde{y} = H_x \Delta \tilde{x} + H_u \Delta \tilde{u}$$

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
$$\Delta y = H_x \Delta x^* + H_u \Delta u^* = \Delta y^*$$

The optimal control law is implemented as:
$$\Delta u = \Delta u^* - C\Delta x + \Delta x^*= (B_{22} + CB_{12}) \Delta y^* - C\Delta x = C_{fwd}\Delta y^* - C_{fbk} \Delta x$$
$$u = u_{trim} + C_{fwd}(y^* - y_{trim}) - C_{fbk} (x - x_{trim})$$

## LQR Tracker With Integral Compensation
We augment the system dynamics with an integral state vector $\xi$ such that:
$$ \dot{\xi} = \Delta \tilde{y} = H_x \Delta \tilde{x} + H_u \Delta \tilde{u}$$

This integral state is computed as:
$$ \xi = \xi(0) + \int^t_{0} \Delta \tilde{y} d\tau = \xi(0) + \int_0^t (\Delta y - \Delta y^*) d\tau = \xi(0) + \int_0^t (y - y^*) d\tau $$

The augmented system is
$$
\begin{pmatrix} \Delta \dot{\tilde{x}} \\ \dot{\xi} \end{pmatrix} = \begin{pmatrix} F & 0\\ H_x & 0 \end{pmatrix} \begin{pmatrix} \Delta \tilde{x} \\ \xi\end{pmatrix} + \begin{pmatrix} G \\ H_u \end{pmatrix} \Delta \tilde{u}
$$

Defining:
$$ F_{aug} \triangleq \begin{pmatrix} F & 0\\ H_x & 0 \end{pmatrix}$$
$$ G_{aug} \triangleq \begin{pmatrix} G \\ H_u \end{pmatrix} $$
$$ x_{aug} \triangleq \begin{pmatrix} \Delta \tilde{x} \\ \xi \end{pmatrix} $$

We can write:
$$\dot{x}_{aug} = F_{aug} x_{aug} + G_{aug} \Delta \tilde{u}$$

Now, if we design a LQR for the augmented system:
$$\Delta \tilde{u} = -C x_{aug} = \begin{pmatrix} C_x & C_{\xi}\end{pmatrix} x_{aug}$$
$$\dot{x}_{aug} = (F_{aug} - GC)x_{aug}$$

This system is guaranteed to be stable by the LQR, so in the steady-state we have:
$$\dot{x}_{aug}=0$$

Which means:
$$\Delta \dot{\tilde{x}} = 0$$
$$\dot{\xi} = \Delta \tilde{y}=0 = y - y^*$$

Let's assume we have an external disturbance such that the actual system dynamics are instead:
$$
\begin{pmatrix} \Delta \dot{\tilde{x}} \\ \dot{\xi} \end{pmatrix} = \begin{pmatrix} F & 0\\ H_x & 0 \end{pmatrix} \begin{pmatrix} \Delta \tilde{x} \\ \xi\end{pmatrix} + \begin{pmatrix} G \\ H_u \end{pmatrix} \Delta \tilde{u} + \begin{pmatrix} L w \\ 0 \end{pmatrix}
$$

The closed-loop system is now:
$$\dot{x}_{aug} = (F_{aug} - GC)x_{aug} + \begin{pmatrix} L w \\ 0 \end{pmatrix}$$

Stability depends only on the dynamics matrix $F_{aug} - GC$, so this system will still be asymptotically stable. This means that in the steady state $\dot{x}_{aug} = 0$, so $\dot{\xi} = 0$, and therefore necessarily $y = y*$. Because the system is no longer homogeneous, in general we will have $x_{aug} \neq 0$. In particular, $\xi$ will converge to the values required by $\Delta \tilde{y} = 0$.

The optimal control law is implemented as:
$$\Delta u = \Delta u^* - Cx_{aug} = \Delta u^* - C_x \Delta \tilde{x} - C_{\xi} \xi = \Delta u^* - C_x \Delta x + C_x \Delta x^* - C_{\xi} \xi = (B_{22} + C_x B_{12}) \Delta y^* - C_x \Delta x - C_{\xi} \xi$$
$$\Delta u = C_{fwd} \Delta y^* - C_{fbk} \Delta x - C_{\xi} \xi$$
$$u = u_{trim} + C_{fwd}(y^* - y_{trim}) - C_{fbk} (x - x_{trim}) - C_{\xi} \xi$$

With:
$$C_{fwd} = B_{22} + C_x B_{12}$$
$$C_{fbk} = C_x$$
$$ \xi = \xi(0) + \int_0^t (y - y^*) d\tau $$
