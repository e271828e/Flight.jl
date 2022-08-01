See "Manual of the ICAO Standard Atmosphere", "Principles of GNSS...", 2.138

A couple of comments...


### Latitude Dependence

Applying fluidostatic balance to a differencial slice of atmosphere gives:

$dp = -\rho g dh = -\rho g(\phi, h) dh$

Where $h$ denotes MSL geometric altitude and $\phi$ is geodetic latitude.

Define the geopotential altitude $z$ as:

$dz = \dfrac{g(\phi, h)}{g(\phi, 0)}dh$

$z = \dfrac{1}{g(\phi, 0)} \int_{0}^{h}{g(\phi, h) dh}$

So that:

$dp = -\rho g(\phi, 0) dh$

Thus, strictly speaking, since gravity depends on latitude, the relationship $z=z(h)$ altitude does
too.

However, in this context, a sufficiently accurate approximation for $g(\phi,h)$ is:

$g(\phi,h) \approx g(\phi,0) \left({\dfrac{r_s(\phi)}{r_s(\phi) + h}}\right)^2$

Where $r_s$ denotes the length of the ECEF position vector to the geoid surface.

Inserting this in the definition of $z$:

<!--
$z = \int_{0}^{h}{\left({\dfrac{r_s(\phi)}{r_s(\phi) + h}}\right)^2 dh}$
$z = \displaystyle \int_{0}^{h}{\left({\dfrac{r_s(\phi)}{r_s(\phi) + h}}\right)^2 dh}$
-->

$$z = \int_{0}^{h}{\left({\dfrac{r_s(\phi)}{r_s(\phi) + h}}\right)^2 dh}$$

Integrating gives:

$z = \dfrac{r_s(\phi) h}{r_s(\phi) + h}$

A further, and also sufficiently accurate approximation is:

$z = \dfrac{ah}{a + h}$

Where $a$ is the Earth's equatorial radius. This eliminates the latitude dependency altogether.

### Non-Standard Base Temperatures and Pressures

The temperature variation with geopotential altitude $T = T(z)$ is assumed linear and it is known
for each atmospheric layer.

The temperature at the base of each atmospheric layer can be precomputed for a standard day, in
which $T(0) = T_{0} = 288.15\;K$). Then, since its dependence with geopotential height is linear, if
$T(0) = T_0 + \Delta T$, we can simply apply the same offset $\Delta T$ to the standard base
temperature for each layer.