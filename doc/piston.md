# Power and SFC Computation for Lycoming IO 360

References:
Lycoming IO-360 Operator's Manual, Figures 3-21 & 3-5
Internal Combustion Engine Fundamentals (Second Edition)
Kevin Horton's Spreadsheets: https://vansairforce.net/community/showthread.php?t=61330

## ISA Relations

The ISA model relates ambient pressure and temperature at a given geopotential altitude $h$ with their corresponding values at sea level, $p_{sl}$ and $T_{sl}$:
$$\dfrac{T}{T_{sl}} = 1 + \dfrac{\beta}{T_{sl}} h$$
$$\dfrac{p}{p_{sl}} = \left( 1 + \dfrac{\beta}{T_{sl}} h\right)^{-\dfrac{g_{std}}{\beta R}} = \left(\dfrac{T}{T_{sl}} \right)^{-\dfrac{g_{std}}{\beta R}}$$

For the first ISA layer ($h<11000 m$), when the sea level conditions are ISA standard conditions:
$$p = p_{std}\left( 1 + \dfrac{\beta}{T_{std}} h\right)^{-\dfrac{g_{std}}{\beta R}}$$
$$T = T_{std}\left( \dfrac{p}{p_{std}} \right)^{-\dfrac{\beta R}{g_{std}}}$$

From these we can define the following functions:
$$p_{ISA}(h) := p_{std}\left( 1 + \dfrac{\beta}{T_{std}} h\right)^{-\dfrac{g_{std}}{\beta R}}$$
$$T_{ISA}(p) = T_{std}\left( \dfrac{p}{p_{std}} \right)^{-\dfrac{\beta R}{g_{std}}}$$

The first one gives the static ambient pressure one would encounter at an altitude $h$ if the sea level conditions were standard ISA conditions. This function can be used to compute the static pressure $p$ corresponding to the pressure altitudes used by Figure 3-21.

The second one gives the temperature one would encounter at an altitude at which the static pressure is $p$ if the sea level conditions were standard conditions.


## Engine Brake Power Output

From theory, we know that engine power is a function of engine speed ($\omega$), manifold air pressure ($MAP$), ambient pressure $p$ and ambient temperature $T$:

$$P = P(\omega, MAP, p, T)$$

More specifically, engine power can be expressed as follows (see *Internal Combustion Engine Fundamentals*, Section 2.12. E):
$$P = \eta_{f,b} Q_{HV} \chi_{f} \dot{m}_a$$

Where $\eta_{f,b}$ is the brake fuel conversion efficiency, $Q_{HV}$ is the fuel's heating value, $\chi_f$ is the fuel/air ratio and $\dot{m}_a$ is the air mass flow rate entering the engine.

The air flow through the inlet duct can be modeled as an isentropic flow. Therefore, air mass flow ratio can be expressed as:
$$\dot{m}_a =
A_e \dfrac{p}{\sqrt{RT}}\sqrt{\dfrac{2\gamma}{\gamma -1} \left[ \left(\dfrac{p}{p_s}\right)^{-2/\gamma} - \left(\dfrac{p}{p_s}\right)^{-(\gamma+1)/\gamma}\right]}
$$

Where $p$ and $T$ are the total (stagnation) pressure and temperature in the duct (respectively equal to the ambient static pressure and temperature), $p_s$ is the static pressure measured at some section inside the duct, and $A_e$ is an equivalent area corresponding to that section.

For a given throttle setting and engine speed, the ratio $p/p_s$ can be assumed to be constant. Under these conditions, we have:
$$\dot{m}_a \propto \dfrac{p}{\sqrt{T}}$$

And therefore:
$$P \propto \dfrac{p}{\sqrt{T}} = \dfrac{p}{\sqrt{T_{ISA}}}\sqrt{\dfrac{T_{ISA}}{T}} \propto
\dfrac{p/p_{std}}{\sqrt{T_{ISA}/T_{std}}}\sqrt{\dfrac{T_{ISA}}{T}}$$

For convenience, we define the following non-dimensional quantity:
$$\delta =\dfrac{p/p_{std}}{\sqrt{T_{ISA}/T_{std}}}$$

Then:
$$P \propto  \delta \sqrt{\dfrac{T_{ISA}}{T}}$$

Now, introducing $T_{ISA}(p)$ we have the following relation:
$$\delta(p) =\dfrac{p/p_{std}}{\sqrt{T_{ISA}(p)/T_{std}}} =
\dfrac{p}{p_{std}} \left(\dfrac{T_{ISA}(p)}{T_{std}}\right) ^{-\dfrac{1}{2}} =
 \dfrac{p}{p_{std}} \left( \dfrac{p}{p_{std}} \right)^{\dfrac{\beta R}{2g_{std}}}=
 \left( \dfrac{p}{p_{std}} \right)^{1+\dfrac{\beta R}{2g_{std}}}$$

Since $\delta = \delta(p)$, we can write, without loss of generality:
$$P(n,MAP,p,T) =P(n,MAP,\delta(p),T) = P(n, MAP, \delta, T)$$

$$\dfrac{P(n,MAP,\delta,T)}{P(n,MAP,\delta,T_{ISA})} \propto
\dfrac{\delta \sqrt{\dfrac{T_{ISA}}{T}}}{\delta \sqrt{\dfrac{T_{ISA}}{T_{ISA}}}} = \sqrt{\dfrac{T_{ISA}}{T}}$$

All power values in Figure 3-21 are valid for $T_{ISA}$, and therefore so will be the result of any calculation using it. The actual $P(n,MAP,\delta,T)$ for an arbitrary $T$ is then computed as:
$$P(\omega,MAP,\delta,T) = P_{ISA}(\omega,MAP,\delta)\sqrt{\dfrac{T_{ISA}}{T}}$$

Where:
$$P_{ISA}(\omega, MAP, \delta) := P(\omega, MAP, \delta, T_{ISA})$$

## Wide-Open Throttle Power at Altitude
The right graph in Figure 3-21 provides $P_{ISA,wot}(\omega,h)$, that is, the brake power produced by the engine at wide-open throttle as a function of altitude $h$ for each engine speed $n$, under the assumption $p = p_{ISA}(h)$ and $T = T_{ISA}(h)$. Although this function appears linear in $h$ for each $n$ it actually isn't, because the x-axis is not uniformly scaled. This non-linearity is to be expected. We know from theory that, for a given $n$ and throttle setting (wide-open throttle in this case), $P_{ISA}(\delta)$ is a linear function. Since $\delta = \delta(p) = \delta(p(h))$ is non-linear, so must be $P_{ISA}(h)$.

For convenience, we would like to replace $h$ with $\delta$ as the x-axis variable. To do this, for each $h$ value we simply evaluate:
$$\delta(h) = \delta(p(h)) = \delta(p_{ISA}(h)) = \dfrac{p_{ISA}(h)/p_{std}}{\sqrt{T_{ISA}(h)/T_{std}}}$$

When using the engine model, we will get $p$ as an input rather than $h$. The corresponding $\delta$ to be entered into the performance tables is given simply by the previous expression:
$$\delta({p}) = \left( \dfrac{p}{p_{std}} \right)^{1+\dfrac{\beta R}{2g_{std}}}$$

The linearity of $P_{ISA,wot}(\omega, \delta)$ with $\delta$ allows us to construct a linear interpolation for it using only two values $P_{ISA,wot}(\omega, \delta_1)$ and $P_{ISA,wot}(\omega, \delta_2)$ for each $n$, where $\delta_1 = \delta_{std}$.

From the right graph in Figure 3-21  we can see that, for a given engine speed $\omega$, the MAP at wide-open throttle is fully determined by altitude $h$ or, equivalently, by $\delta$. Let:
$$MAP_{wot} = MAP_{wot}(\omega, \delta)$$

This relation between wide-open throttle $MAP$ and $\delta$ is one-to-one. Thus, for each $\omega$ and $MAP$ combination, there is a $\delta$ value $\delta_{wot}(\omega, MAP)$ for which $MAP$ is the wide-open throttle $MAP$. This inverse function is:
$$\delta_{wot} = \delta_{wot}(\omega, MAP)$$

The power output of the engine at a given $\omega$ and $\delta$ combination when operated at wide-open throttle is a function $P_{ISA,wot}(\omega, \delta)$, which can also be constructed from the right graph.

## Part Throttle Power at Standard Conditions

At standard sea level conditions, for an arbitrary throttle setting, we have:
$$P_{ISA,std}(\omega, MAP) = P_{ISA}(\omega, MAP, \delta_{std})$$

This function is provided by the left graph in Figure 3-21, with $MAP$ on the x axis and $P$ on the y axis, for different values of $n$.

## Part Throttle Power at Altitude

Our goal is to compute $P_{ISA}(\omega, MAP, \delta)$ for the general case. For this, we can write:
$$P_{ISA}(\omega, MAP, \delta) \approx P_{ISA}(\omega, MAP, \delta_{std}) + \left.\dfrac{\partial P_{ISA}}{\partial \delta}\right|_{\omega,MAP,\delta_{std}}(\delta - \delta_{std})$$

The first term is directly available from the left graph.

For the second term, we can estimate the derivative at $\delta_{std}$ as:
$$\left.\dfrac{\partial P_{ISA}}{\partial \delta}\right|_{\omega,MAP,\delta_{std}}\approx \frac{P_{ISA}(\omega,MAP,\delta_{ref}) - P_{ISA}(\omega,MAP,\delta_{std})}{\delta_{ref} - \delta_{std}}$$

Where $\delta_{ref}$ denotes some arbitrary reference $\delta$ value and $P_{ISA}(\omega,MAP,\delta_{ref})$ is the corresponding power output.

We can choose $\delta_{ref} = \delta_{wot}(\omega, MAP)$, that is, a $\delta$ value such that $MAP$ is the wide-open throttle $MAP$. In that case, the power output $P_{ISA}(\omega,MAP,\delta_{ref})$ will be:

$$P_{ISA}(\omega, MAP, \delta_{ref}) = P_{ISA,wot}(\omega,\delta_{wot}(\omega, MAP))$$

Once we have $P_{ISA}(\omega,MAP,\delta)$, we compute:
$$P(\omega,MAP,\delta,T) = P_{ISA}(\omega,MAP,\delta)\sqrt{\dfrac{T_{ISA}}{T}}$$

Where:
$$T_{ISA}(p) = T_{std}\left( \dfrac{p}{p_{std}} \right)^{-\dfrac{\beta R}{g_{std}}}$$

## Throttle to MAP

Our actual input to the engine is not a $MAP$ value, but a throttle setting $thr$. What should the relation between the engine throttle input and the actual $MAP$ be in each operating condition?

For each $(\omega, \delta)$ combination, the right graph gives us $MAP_{wot}$. Then we simply compute:
$$MAP(thr, \delta, \omega) = MAP_{wot} \left(k + thr (1 - k)\right)$$

Where $k$ is a suitably chosen $MAP_{idle}/MAP_{wot}$ ratio that achieves the desired target idle engine RPMs at each operating condition.

With this relation, engine power becomes a function of the form:
$$P = P(\omega, MAP, \delta, T) = P(thr, \omega, \delta, T)$$

## Continuous Operation Limits

Both graphs in Fig. 3-21 show a line labelled "limiting manifold pressure for continuous operation". From looking at both graphs, it is apparent that this line actually enforces a maximum power for each $\omega$. That is, it defines a function $P_{ISA,lim}(\omega)$

So far we have ignored this limit when computing $P_{ISA}$, so we may obtain values $P_{ISA} > P_{ISA, lim}(\omega)$. In principle, this power limitation is not enforced by the engine itself; they are under the pilot's responsibility. Therefore, we can choose to ignore them when modeling the engine. We could easily enforce these limits by setting $P_{ISA} = min\{P_{ISA}, P_{ISA,lim}(\omega)\}$.

## Effect of Mixture

Figure 3-21 assumes maximum power mixture ($mix = 1$). We can use the ratios in figure 3-1 to scale down the power output for $mix \in [0,1]$.

## Fuel consumption

Figure 3-5 provides fuel consumption as:
$$\dot{m} = \dot{m}(\omega,P,mix)$$

Where $P$ is the actual engine power output.


To model any fuel injected piston engine in the Lycoming lineup, for example IO-540, non-dimensionalize:
- Power values in Figure 3-21 with rated power
- Fuel consumption in figure 3-5 with fuel consumption at rated power with maximum power mixture.

## Piston Thruster

We define a piston thruster as the combination of a piston engine, a propeller and a transmission mechanism connecting them. The transmission imposes a kinematic constraint on the angular velocities of the engine and the propeller.

For the engine the angular momentum equation along the crankshaft can be written as:
$$M_{tr, cs} + M_{eng,cs} = J_{cs} \dot{\omega}_{cs}$$

Where:
- $M_{tr,cs}$: Torque exerted on the crankshaft by the transmission.
- $M_{eng, cs}$: Output torque exerted by the engine on the crankshaft.
- $J_{cs}$: Crankshaft moment of inertia.
- $\omega_{cs}$: Crankshaft angular velocity.

For the propeller:
$$M_{tr, p} + M_{air,p} = J_{p} \dot{\omega}_{p}$$

Where:
- $M_{tr,p}$: Torque exerted on the propeller by the transmission.
- $M_{air, p}$: Aerodynamic torque on the propeller.
- $J_{p}$: Propeller moment of inertia.
- $\omega_{p}$: Propeller angular velocity.

The transmission transfers mechanical power from the engine crankshaft to the propeller  with a mechanical efficiency $\eta$:
$$P_{out} = \eta P_{in}$$.

Where $P_{out}$ is the power exerted on the propeller and $P_{in}$ is the power provided by the engine. These are given by:
$$P_{out} = M_{tr,p} \omega_{p}$$
$$P_{in} = M_{cs,tr} \omega_{cs}$$

The angular velocity of the crankshaft and the propeller are related by the transmission's (signed) gear ratio $n$:
$$\omega_{p} = n \omega_{cs}$$

All kinematic and dynamic magnitudes defined so far are positive clockwise.

From the above equations we have:
$$M_{tr, p} \omega_{p} = \eta M_{cs, tr} \omega_{cs} $$
$$M_{tr, cs} = -M_{cs,tr} = -n \dfrac{1}{\eta} M_{tr, p}$$

Using this equality we can combine the engine and propeller angular momentum equations to yield:
$$-n \dfrac{1}{\eta} M_{tr, p} + M_{eng,cs} = J_{cs} \dot{\omega}_{cs}$$

$$M_{eng, cs} + n \dfrac{1}{\eta} M_{air, p} = \left(J_{cs} + \dfrac{n^2}{\eta}
J_{p} \right) \dot{\omega}_{cs}$$

Or, equivalently:
$$M_{air, p} + \dfrac{\eta}{n} M_{eng, cs} = \left(J_{p} + \dfrac{\eta}{n^2}
J_{cs} \right) \dot{\omega}_{p}$$