Variables `q` and `υ` are defined as:

$$q(t_k) = q^{c_0}_{c_k}$$
$$\upsilon(t_k) = \int_{t_0}^{t_{k}} R^{c_0}_{c} f^{c} dt$$

At each $t_k$, we need to provide
$$\int_{t_{k-1}}^{t_{k}} R^{c_{k-1}}_{c} f^{c} dt$$

Which we can expand as
$$\int_{t_{k-1}}^{t_{k}} R^{c_{k-1}}_{c} f^{c} dt = (R^{c_{0}}_{c_{k-1}})^T \int_{t_{k-1}}^{t_{k}} R^{c_{0}}_{c} f^{c} dt =
(R^{c_{0}}_{c_{k-1}})^T \left( \int_{t_{0}}^{t_{k}} R^{c_{0}}_{c} f^{c} dt -
\int_{t_{0}}^{t_{k-1}} R^{c_{0}}_{c} f^{c} dt \right) = (R^{c_{0}}_{c_{k-1}})^T
(\upsilon(t_k) - \upsilon(t_{k-1}))$$

In code this corresponds to the line:

`υ_c_sc = z.q'(u.V - z.V)`