## Terrain Model

### Queries

A terrain model answers two queries:

1. Surface data at a 2D location: the surface point $P$ at that location, its inward unit normal $k_t^e$ and its surface type.

2. Surface intersection with a bounded ray. Given the ray origin $O$, its unit direction $u^e$ and a bound $l_{max}$, the model returns the nearest surface point $P$ with ray parameter $l \in [0, l_{max}]$, together with $k_t^e$ at $P$ and the surface type, or reports that there is none. $k_t^e$ must be unit-norm and inward-pointing, so that callers can adopt it directly as a frame axis. This is the query used for [landing gear contact determination](landing_gear.md).

Any cheap pre-test to skip the intersection when the ray origin is far from the surface belongs inside the terrain model. For a mesh terrain, an altitude-based cull is not safe as a strict test on steep slopes; for a uniform terrain, it fuses with the intersection itself.

### Uniform Terrain

A uniform terrain has constant elevation and surface type. The elevation may be given as orthometric (constant altitude above the geoid) or ellipsoidal (constant altitude above the WGS84 ellipsoid). The latter is cheaper to query and analytically exact, since no geoid interpolation is involved.

The surface is approximated near the ray origin $O$ by its local tangent plane, anchored at the terrain point $P_0$ that shares the 2D location of $O$:

$$n^e(P_0) = n^e(O)$$

$$h(P_0) = h_{trn}$$

The plane's inward normal is the reversed local vertical at $P_0$, that is, the reversed n-vector:

$$k_t^e = -n^e(P_0)$$

For orthometric elevation the geoid slope is neglected, consistent with the 2D query. Within this approximation the normal field is uniform over the ray span, so its value at the intersection point is that at $P_0$.

A point $P$ lies on the plane when:

$$(k_t^e)^T (r_{O_eP}^e - r_{O_eP_0}^e) = 0$$

Imposing this on the ray $r_{O_eP}^e = r_{O_eO}^e + l u^e$:

$$\cos \alpha = (k_t^e)^T u^e$$

$$l = \dfrac{(k_t^e)^T (r_{O_eP_0}^e - r_{O_eO}^e)}{\cos \alpha}$$

If $\cos \alpha \leq 0$ the ray is parallel to the surface or points away from it, and there is no intersection. Otherwise the intersection is valid when $0 \leq l \leq l_{max}$. Since $r_{O_eO}^e - r_{O_eP_0}^e$ is along $n^e(P_0)$, the numerator is the height of $O$ above the surface, so the bound check doubles as the altitude cull.

The approximation error grows with the intersection distance, roughly $l / R_e$ rad in normal direction and $l^2 / 2 R_e$ in surface height, so results are only meaningful for $l \ll R_e$. Should longer rays ever require it, the solution can be refined by re-anchoring the tangent plane at the 2D location of the candidate intersection and intersecting again.
