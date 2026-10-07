"""Route score: did the visitor come from, or is it heading to, a nearby star?

Everything here is 3D and far-field: the visitor moves on straight lines (its inbound and outbound asymptotes, through the Sun), stars move linearly. Units: pc, km/s, Myr.

For star k at position s (now, t = 0) with velocity w, and the visitor's asymptotic velocity V, the relative motion is r(t) = s + (w - V) t. Inbound (t < 0) and outbound (t > 0) closest approaches give a miss distance d_k and the time t_k of the encounter.

Null hypothesis: the visitor's direction is isotropic and unrelated to the stars. Then u = w - V is uniform on a sphere of radius v_inf about w, and a miss within angle g of star k needs u inside a thin cone around +-s_hat:

	p_k = sum over cone/sphere crossings  lam^2 g^2 / (4 v_inf |lam - s_hat.w|)

for small g (the cone's cross-section projected onto the sphere, over its area); `p_values` evaluates the exact finite-cone fraction numerically.

Look-elsewhere uses a weighted Bonferroni over stars with weights falling as 1/rank by distance among stars reachable at v_inf, encoding a planner's preference for short hops: score = min_k p_k / w_k, so the nearest reachable stars are cheap to "try".
"""
import math
import os
from dataclasses import dataclass

import numpy as np
from scipy.stats import rice

from .constants import V_UNIT_KMS

KMS_PC_MYR = 1.0227  # 1 km/s in pc/Myr
SCREEN = 0.2         # exact p only for stars passed within this fraction of D
N_NOISE = 24         # grid points over the true miss for the noise convolution


@dataclass
class Errors:
	star_v: float = 0.3        # km/s per axis, star space-velocity error
	parallax_mas: float = 0.03
	visitor_dir: float = 1e-5  # rad, asymptote direction error
	visitor_v: float = 0.02    # km/s, v_inf error

	def miss_sigma(self, D, t):
		"""1-sigma error on the miss distance (pc) for distance D, time t."""
		tt = np.abs(t) * KMS_PC_MYR
		sd = D * D * self.parallax_mas * 1e-3
		return np.sqrt((tt * self.star_v) ** 2 + (tt * self.visitor_v) ** 2
					   + (D * self.visitor_dir) ** 2 + sd ** 2)


@dataclass
class Field:
	pos: np.ndarray   # (N, 3) pc, relative to the Sun now
	vel: np.ndarray   # (N, 3) km/s, relative to the Sun
	D: np.ndarray
	ids: np.ndarray = None

	@classmethod
	def from_gaia(cls, path, radius_pc=20.0):
		"""Gaia DR3 stars with radial velocities, ICRS barycentric frame.

		Downloads the table from the Gaia archive if `path` is missing. Returns (field, errors, n_without_rv). Per-star errors: velocity error per axis from parallax, proper-motion and RV errors.
		"""
		import csv
		if not os.path.exists(path):
			fetch_gaia(path, radius_pc)
		rows = list(csv.DictReader(open(path)))
		n_all = len(rows)
		rows = [r for r in rows if r["radial_velocity"]]
		f = lambda k: np.array([float(r[k]) for r in rows])
		ra, dec, plx = np.radians(f("ra")), np.radians(f("dec")), f("parallax")
		D = 1000.0 / plx
		rhat = np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], 1)
		e_ra = np.stack([-np.sin(ra), np.cos(ra), np.zeros_like(ra)], 1)
		e_dec = np.stack([-np.sin(dec) * np.cos(ra), -np.sin(dec) * np.sin(ra), np.cos(dec)], 1)
		k = 4.740470446
		vra, vdec, vr = k * f("pmra") / plx, k * f("pmdec") / plx, f("radial_velocity")
		vel = vra[:, None] * e_ra + vdec[:, None] * e_dec + vr[:, None] * rhat
		s_plx = f("parallax_error")
		s_ra = k * np.hypot(f("pmra_error") / plx, f("pmra") * s_plx / plx ** 2)
		s_dec = k * np.hypot(f("pmdec_error") / plx, f("pmdec") * s_plx / plx ** 2)
		s_v = np.sqrt((s_ra ** 2 + s_dec ** 2 + f("radial_velocity_error") ** 2) / 3)
		field = cls(D[:, None] * rhat, vel, D, np.array([r["source_id"] for r in rows]))
		errors = Errors(star_v=s_v, parallax_mas=s_plx)
		return field, errors, n_all - len(rows)

	@classmethod
	def synthetic(cls, rng, density=0.08, radius=20.0,
				  v_mean=(-11.1, -12.2, -7.3), v_sigma=(35.0, 25.0, 18.0)):
		"""Uniform stars, Gaussian velocities relative to the Sun.

		Defaults: ~0.08 systems/pc^3 (the 10 pc census) and the local disk velocity ellipsoid, offset by the Sun's peculiar motion.
		"""
		n = rng.poisson(density * 4 / 3 * math.pi * radius ** 3)
		r = radius * rng.random(n) ** (1 / 3)
		d = rng.normal(size=(n, 3))
		pos = d / np.linalg.norm(d, axis=1)[:, None] * r[:, None]
		vel = rng.normal(v_mean, v_sigma, size=(n, 3))
		return cls(pos, vel, np.linalg.norm(pos, axis=1))


GAIA_TAP = "https://gea.esac.esa.int/tap-server/tap/sync"


def fetch_gaia(path, radius_pc=20.0):
	"""Gaia DR3 sources within radius_pc with parallax S/N > 10, as CSV."""
	from urllib.parse import urlencode
	from urllib.request import urlopen
	query = ("SELECT source_id, ra, dec, parallax, parallax_error, pmra, pmra_error, "
			 "pmdec, pmdec_error, radial_velocity, radial_velocity_error, phot_g_mean_mag "
			 f"FROM gaiadr3.gaia_source WHERE parallax > {1000.0 / radius_pc} "
			 "AND parallax_over_error > 10")
	data = urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv",
					  "QUERY": query}).encode()
	with urlopen(GAIA_TAP, data=data, timeout=600) as r:
		body = r.read()
	if not body.startswith(b"source_id"):
		raise RuntimeError("unexpected Gaia archive response")
	os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
	with open(path, "wb") as f:
		f.write(body)


def end_weights(field, v_inf, sign, tol=0.5):
	"""Look-elsewhere weights for one end: 1/rank by distance, reachable only.

	A star can be met at speed v_inf only if some direction puts the relative velocity u = w - V along +-s_hat. Directions to the sphere of possible u (center w, radius v_inf) form a cone about w_hat of half-angle asin(v_inf/|w|) (all directions if |w| <= v_inf), so the smallest achievable miss is D sin(max(0, angle(axis, w) - half-angle)). Stars that cannot be approached within `tol` pc get no weight. This depends only on v_inf, never on the visitor's direction, so the test stays valid while the budget is not wasted on impossible stars.
	"""
	axis = field.pos / field.D[:, None] * (-sign)
	wn = np.linalg.norm(field.vel, axis=1)
	cosang = np.einsum("ij,ij->i", axis, field.vel) / wn
	ang = np.arccos(np.clip(cosang, -1, 1))
	half = np.arcsin(np.minimum(1.0, v_inf / wn))
	gmin = np.where(wn <= v_inf, 0.0, np.maximum(0.0, ang - half))
	ok = field.D * np.sin(np.minimum(gmin, 0.5 * math.pi)) < tol
	rank = np.empty(len(field.D))
	order = np.argsort(np.where(ok, field.D, np.inf))
	rank[order] = np.arange(1, len(order) + 1)
	w = np.where(ok, 1 / rank, 0.0)
	return w / w.sum()


def encounters(field, V, sign):
	"""Miss distance and time for every star; sign -1 inbound, +1 outbound.

	Stars whose closest approach falls on the wrong side of t = 0 get the t = 0 separation (the visitor never came nearer).
	"""
	u = field.vel - V
	t = -np.einsum("ij,ij->i", field.pos, u) / (KMS_PC_MYR * np.einsum("ij,ij->i", u, u))
	t = np.where(sign * t > 0, t, 0.0)
	r = field.pos + u * (t * KMS_PC_MYR)[:, None]
	return np.linalg.norm(r, axis=1), t


def p_values(field, v_inf, d_eff, sign, n_grid=160):
	"""Chance that an isotropic visitor passes star k within d_eff.

	Exact for any cone angle: u = w + v n_hat with n_hat uniform; write n_hat in polar angle theta from the cone axis a (a = +-s_hat) and azimuth phi from w's transverse part (size rho). The cone condition |u_perp| <= tan(g) u_par is, for each theta, cos(phi) <= c(theta), a closed-form fraction of phi; integrate over cos(theta)/2. The grid is concentrated near the angles where the sphere comes closest to the axis.
	"""
	a = field.pos / field.D[:, None] * (-sign)
	w_par = np.einsum("ij,ij->i", a, field.vel)
	rho = np.sqrt(np.maximum(np.einsum("ij,ij->i", field.vel, field.vel) - w_par ** 2, 0.0))
	v = v_inf
	g = np.minimum(np.arcsin(np.minimum(d_eff / field.D, 1.0)), 0.5 * math.pi - 1e-6)
	tg = np.tan(g)
	th0 = np.arcsin(np.minimum(1.0, rho / v))
	T_max = tg * (np.abs(w_par) + v)
	h = np.minimum(0.5 * math.pi, 4 * np.sqrt(T_max / v) + 4 * T_max / v + 1e-9)
	x = np.linspace(-1, 1, n_grid)
	grids = [np.broadcast_to(np.linspace(0, math.pi, n_grid), (len(rho), n_grid))]
	for c in (th0, math.pi - th0):
		grids.append(np.clip(c[:, None] + h[:, None] * x, 0, math.pi))
	th = np.sort(np.concatenate(grids, axis=1), axis=1)
	sin, cos = np.sin(th), np.cos(th)
	u_par = w_par[:, None] + v * cos
	T2 = np.where(u_par > 0, (tg[:, None] * u_par) ** 2, -1.0)
	num = T2 - rho[:, None] ** 2 - (v * sin) ** 2
	den = 2 * rho[:, None] * v * sin
	with np.errstate(divide="ignore", invalid="ignore"):
		c = np.where(den > 0, num / den, np.where(num >= 0, 1.0, -1.0))
	frac = np.where(T2 < 0, 0.0, 1 - np.arccos(np.clip(c, -1, 1)) / math.pi)
	# integrate frac over d(cos theta) / 2
	dmu = -np.diff(cos, axis=1)
	p = 0.5 * np.sum(0.5 * (frac[:, 1:] + frac[:, :-1]) * dmu, axis=1)
	return np.clip(p, 0.0, 1.0)


def end_score(field, V, sign, errors, rng=None, weights=None):
	"""Score one end of the route. Returns (p_end, best star index, info).

	rng: if given, simulate measurement noise on the misses (for the null). The observed miss carries a 2D Gaussian measurement error of size sigma, so p is the null chance of an *observed* miss this small: the null distribution of the true miss r (from p_geom) convolved with the Rice distribution of |r + noise|, p = integral P(Rice(r, sigma) <= d_obs) dp_geom(r). Where the null density is locally uniform this equals p_geom(d_obs); it matters near the hard floor of barely reachable stars, where noise below the floor would otherwise score as impossible.
	"""
	v_inf = float(np.linalg.norm(V))
	d, t = encounters(field, V, sign)
	sig = errors.miss_sigma(field.D, t)
	if rng is not None:
		# the transverse miss is a 2D vector; add 2D Gaussian error
		d = np.hypot(d + rng.normal(size=d.shape) * sig, rng.normal(size=d.shape) * sig)
	# Only stars passed within a fraction of their distance can be rare; others (including stars never approached, d = D) get p = 1, which is exact for d >= D and conservative otherwise.
	p = np.ones_like(d)
	near = np.flatnonzero(d + 2 * sig < SCREEN * field.D)
	if len(near):
		dn, sn = d[near], sig[near]
		r = np.linspace(0, 1, N_NOISE + 1) * (dn + 8 * sn)[:, None]   # true miss grid
		rows = np.repeat(near, N_NOISE + 1)
		sub = Field(field.pos[rows], field.vel[rows], field.D[rows])
		pg = p_values(sub, v_inf, r.ravel(), sign).reshape(len(near), N_NOISE + 1)
		rm = 0.5 * (r[:, 1:] + r[:, :-1])
		kern = rice.cdf(dn[:, None] / sn[:, None], rm / sn[:, None])
		p[near] = np.sum(np.diff(pg, axis=1) * kern, axis=1)
	if weights is None:
		weights = end_weights(field, v_inf, sign)
	with np.errstate(divide="ignore"):
		ratio = np.where(weights > 0, p / weights, np.inf)
	k = int(np.argmin(ratio))
	return min(1.0, float(ratio[k])), k, {"d": d, "sigma": sig, "t": t, "p": p,
										  "weights": weights}


def combine(p1, p2):
	"""Exact tail of the product of two independent uniforms."""
	q = p1 * p2
	return q * (1 - math.log(q)) if q > 0 else 0.0


def sun_turn(V_in, q_au, phi):
	"""Outbound asymptote after the Sun bends an arrival V_in (km/s).

	q_au: perihelion; phi: azimuth of the orbit plane about V_in.
	"""
	v = np.linalg.norm(V_in)
	e = 1 + q_au * (v / V_UNIT_KMS) ** 2
	delta = 2 * math.asin(1 / e)
	a = V_in / v
	tmp = np.array([1.0, 0, 0]) if abs(a[0]) < 0.9 else np.array([0, 1.0, 0])
	n1 = np.cross(a, tmp)
	n1 /= np.linalg.norm(n1)
	n2 = np.cross(a, n1)
	n = math.cos(phi) * n1 + math.sin(phi) * n2
	return v * (math.cos(delta) * a + math.sin(delta) * n)


def random_unit(rng):
	x = rng.normal(size=3)
	return x / np.linalg.norm(x)


def null_route(field, v_inf, n, errors, seed=0, q_max=30.0):
	"""Route scores of n visitors with random arrival, perihelion, plane.

	Impact parameters are uniform in area out to perihelion q_max.
	"""
	rng = np.random.default_rng(seed)
	g = 2 / (v_inf / V_UNIT_KMS) ** 2        # b^2 = q^2 + g q (AU)
	b2max = q_max ** 2 + g * q_max
	w_in, w_out = end_weights(field, v_inf, -1), end_weights(field, v_inf, +1)
	out = np.empty((n, 3))
	for k in range(n):
		V_in = v_inf * random_unit(rng)
		b2 = b2max * rng.random()
		q = 0.5 * (-g + math.sqrt(g * g + 4 * b2))
		V_out = sun_turn(V_in, q, rng.uniform(0, 2 * math.pi))
		p_in, _, _ = end_score(field, V_in, -1, errors, rng, weights=w_in)
		p_out, _, _ = end_score(field, V_out, +1, errors, rng, weights=w_out)
		out[k] = (p_in, p_out, combine(p_in, p_out))
	return out


def aim_at(star_pos, star_vel, v_inf, sign):
	"""Asymptotic velocity that meets the star exactly (|V| = v_inf).

	Solves s + (w - V) t = 0 for t with sign(t) = sign: V = w + s / t. Returns (V, t) or (None, None) if the star is unreachable at v_inf.
	"""
	s, w = star_pos, star_vel
	# |w + s/t|^2 = v^2 with x = 1/t:  |s|^2 x^2 + 2 (s.w) x + |w|^2 - v^2 = 0
	A, B, C = s @ s / KMS_PC_MYR ** 2, 2 * (s @ w) / KMS_PC_MYR, w @ w - v_inf ** 2
	disc = B * B - 4 * A * C
	if disc < 0:
		return None, None
	for x in sorted(((-B + math.sqrt(disc)) / (2 * A), (-B - math.sqrt(disc)) / (2 * A)),
					key=abs):
		if x != 0 and sign * x > 0:
			t = 1 / x
			return w + s / (t * KMS_PC_MYR), t
	return None, None


def plan_route(field, i_from, i_to, v_inf):
	"""Inbound from star i_from, outbound to i_to, via a Sun flyby.

	Returns (V_in, V_out, q_au) or None if the turn needs q below the Sun's radius. The outbound speed is v_inf (Sun-only, energy conserved).
	"""
	V_in, _ = aim_at(field.pos[i_from], field.vel[i_from], v_inf, -1)
	V_out, _ = aim_at(field.pos[i_to], field.vel[i_to], v_inf, +1)
	if V_in is None or V_out is None:
		return None
	cosd = V_in @ V_out / v_inf ** 2
	delta = math.acos(max(-1.0, min(1.0, cosd)))
	if delta < 1e-9:
		return None
	e = 1 / math.sin(delta / 2)
	q = (e - 1) / (v_inf / V_UNIT_KMS) ** 2
	if q < 0.00465:
		return None
	return V_in, V_out, q


def _null_chunk(args):
	return null_route(*args)


def null_route_parallel(field, v_inf, n, errors, pool, chunks=72, seed=0):
	per = n // chunks
	parts = pool.map(_null_chunk, [(field, v_inf, per, errors, seed + c) for c in range(chunks)])
	return np.vstack(parts)
