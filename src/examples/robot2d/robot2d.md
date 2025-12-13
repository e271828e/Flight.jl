Trim system:
$$
J \dot{\omega} + (K \cos\theta - \frac{J_m}{R}) \dot{v} = K g \sin\theta - \tau_{ss}
$$

$$
(K \cos\theta - \frac{J_m}{R}) \dot{\omega} + M \dot{v} = K \omega^2 \sin\theta + \frac{1}{R}\tau_{ss}
$$

$$ \tau_{ss} = k_m u_m - b_m \left(\frac{v}{R} - \omega\right) $$

$$ \tau_{m} = \tau_{ss} - J_m \left(\frac{\dot{v}}{R} - \dot{\omega}\right)$$

Steps:
- Solve for $\dot{v}$ and $\tau_{ss}$ in (1) and (2)
- With $\tau_{ss}$ and constraints $u_m$ and $\omega$, solve for $v$ in (3)
- With $\dot{v}$, $\tau_{ss}$ and constraint $\dot{\omega}$, solve for $\tau_m$ in (4)