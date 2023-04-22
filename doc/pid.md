# Discretized PID

Parameter definitions:
- $k_p$: Proportional gain
- $k_i$: Integral gain
- $k_d$: Derivative gain
- $\tau_d$: Derivative filter time constant
- $\beta_p$: Proportional path setpoint weighting factor
- $\beta_d$: Derivative path setpoint weighting factor

Variable definitions:
- $r$: Setpoint
- $v$: Feedback
- $u_p$: Proportional path input
- $u_i$: Integral path input
- $u_d$: Derivative path input
- $y_p$: Proportional path output
- $y_i$: Integral path output
- $y_d$: Derivative path output
- $y$: Total output

Path inputs are computed as:
$$u_p = \beta_p r - v$$
$$u_i = r - v$$
$$u_d = \beta_d r - v$$

Proportional path output is simply:
$$y_p(t) = k_p u_p(t)$$

Integral path output is given by:
$$y_i(t) = y_i(t_0) + \int^{t}_{t_0}{k_pu(t)dt}$$

The corresponding state-space model is simply:
$$\dot{x}_i = k_p u_i$$
$${y}_i = x_{i}$$

We now discretize the differential equation using backward differences:
$$\frac{x(t_k) - x(t_{k-1})}{T} = k_p u_p(t_k)$$

The discretized model is:
$$x(t_k) = x(t_{k-1}) + T k_p u_p(t_k)$$
$$y(t_k) = x(t_k)$$

Derivative path output is given by its transfer function:
$$\frac{Y_d(s)}{U_d(s)}= F_d(s) = \frac{k_d s}{1 + \tau_d s}$$
$$(1 + \tau_d s){Y_d(s)}= k_d s U_d(s) $$

Inverting the Laplace transform:
$$\dot{y} = -\frac{1}{\tau_d} y + \frac{k_d}{\tau_d}\dot{u}$$

An equivalent state-space representation is:
$$\dot{x} = -\frac{1}{\tau_d} x + \frac{k_d}{\tau_d} u$$
$$y = -\frac{1}{\tau_d} x + \frac{k_d}{\tau_d} u$$

This can be easily verified by taking the derivative of the output equation and using $\dot{x}=y$:
$$\dot{y} = -\frac{1}{\tau_d} \dot{x} + \frac{k_d}{\tau_d} \dot{u} = -\frac{1}{\tau_d} y + \frac{k_d}{\tau_d} \dot{u}$$

Applying backward differences:
$$\frac{x(t_k) - x(t_{k-1})}{T} = -\frac{1}{\tau_d}x(t_k)+\frac{k_d}{\tau_d}u(t_k)$$

After some manipulation:
$$x(t_k) = \frac{\tau_d}{T + \tau_d} x(t_{k-1}) + k_d \frac{T}{T + \tau_d}u(t_k)$$
$$y(t_k) = -\frac{1}{\tau_d}x(t_k) + \frac{k_d}{\tau_d}u(t_k)$$

The output equation is ill-conditioned for $\tau_d \rightarrow 0$ (unfiltered derivative). We can rewrite it as:
$$y(t_k) =
-\frac{1}{\tau_d}\left(\frac{\tau_d}{T + \tau_d} x(t_{k-1}) + k_d \frac{T}{T + \tau_d}u(t_k)\right) + \frac{k_d}{\tau_d}u(t_k)$$
$$=-\frac{1}{T + \tau_d} x(t_{k-1}) + \frac{k_d}{\tau_d} \left(-\frac{T}{T + \tau_d} + 1\right) u(t_k)$$
$$=\frac{1}{T + \tau_d} \left(x(t_{k-1}) + k_d u(t_k)\right)$$