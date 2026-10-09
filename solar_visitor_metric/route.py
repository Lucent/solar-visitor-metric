"""Route score: did the visitor come from, or is it heading to, a nearby star?

Everything here is 3D and far-field: the visitor moves on straight lines (its inbound and outbound asymptotes, through the Sun), stars move linearly. Units: pc, km/s, Myr.

For star k at position s (now, t = 0) with velocity w, and the visitor's asymptotic velocity V, the relative motion is r(t) = s + (w - V) t. Inbound (t < 0) and outbound (t > 0) closest approaches give a miss distance d_k and the time t_k of the encounter.

Null hypothesis: the visitor's direction is isotropic and unrelated to the stars. Then u = w - V is uniform on a sphere of radius v_inf about w, and a miss within angle g of star k needs u inside a thin cone around +-s_hat:

	p_k = sum over cone/sphere crossings  lam^2 g^2 / (4 v_inf |lam - s_hat.w|)

for small g (the cone's cross-section projected onto the sphere, over its area); `p_values` evaluates the exact finite-cone fraction numerically.

Candidates are the stars a visitor at v_inf could meet within HOP of the Sun: the hop is where the meeting happens, not where the star is today, so a star 60 pc away now that crossed 10 pc from the Sun when the visitor passed is a short hop. A field is complete for visitors at least as fast as its v_floor: a star met at hop L after time L / v_inf is now at most L (1 + |w| / v_inf) away, so it is enough to hold every star within HOP (1 + |w| / v_floor) of the Sun today.

Look-elsewhere uses a weighted Bonferroni over the candidates with weights falling as 1/rank by hop, encoding a planner's preference for short hops: score = min_k p_k / w_k, so the nearest reachable stars are cheap to "try".

The route is scored as one itinerary, origin, Sun, destination: the arrival end against an isotropic null (`end_score`), the departure end given the arrival (`turn_score`), where the only freedom left is how the Sun's turn is oriented. The two p's are then independent by construction and combine exactly; scored separately, a pair of stars on one line through the Sun would count twice.
"""
import math
import os
from dataclasses import dataclass, fields
from types import SimpleNamespace

import numpy as np
from scipy.spatial import cKDTree
from scipy.special import ndtr
from scipy.stats import rice

from .constants import V_UNIT_KMS

KMS_PC_MYR = 1.0227  # 1 km/s in pc/Myr
SCREEN = 0.2         # exact p only for stars passed within this fraction of D
N_NOISE = 24         # grid points over the true miss for the noise convolution
N_TURN = 128         # orientations of the Sun's turn on the coarse circle
N_FINE = 128         # points across one star's window on the circle
HOP = 20.0           # pc from the Sun: the longest hop counted at either end
V_FLOOR = 20.0       # km/s: fields are complete for visitors at least this fast
R_MAX = 200.0        # pc: completeness cap, missing only stars faster than 180 km/s
TOL = 0.5            # pc: stars a visitor can pass this close count as reachable
MATCH = 2.0          # arcsec at the Gaia epoch: a Hipparcos star this close to a usable Gaia star is that star
GAIA_YEAR = 2016.0
HIPPARCOS_YEAR = 1991.25
RV_SNR = 5.0         # Gaia DR3 radial velocities need rv_expected_sig_to_noise >= this
V_ESC = 600.0        # km/s in the Galactic rest frame: faster stars are unbound, i.e. mismeasured
ICRS_TO_GAL = np.array([[-0.0548755604162154, -0.8734370902348850, -0.4838350155487132],
						[0.4941094278755837, -0.4448296299600112, 0.7469822444972189],
						[-0.8676661490190047, -0.1980763734312015, 0.4559837761750669]])
V_SUN_GAL = np.array([11.1, 245.6, 7.25])  # the Sun's velocity about the Galactic center, km/s


def reach_radius(w):
	"""Present distance within which a field must hold stars of speed w (km/s)."""
	return np.minimum(HOP * (1 + w / V_FLOOR) + TOL, R_MAX)


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

	def take(self, idx):
		"""Errors of a subset of stars (per-star arrays indexed, scalars kept)."""
		return Errors(*(x[idx] if np.ndim(x) else x
						for x in (getattr(self, f.name) for f in fields(self))))


@dataclass
class Field:
	pos: np.ndarray   # (N, 3) pc, relative to the Sun now
	vel: np.ndarray   # (N, 3) km/s, relative to the Sun
	D: np.ndarray
	ids: np.ndarray = None

	def take(self, idx):
		return Field(self.pos[idx], self.vel[idx], self.D[idx],
					 None if self.ids is None else self.ids[idx])

	@classmethod
	def from_catalogs(cls, gaia_path, hipparcos_path):
		"""Stars with radial velocities that could meet a visitor within HOP, ICRS barycentric at the Gaia epoch: Gaia DR3 wherever it has a usable radial velocity, and Hipparcos (van Leeuwen 2007) with Pulkovo radial velocities (Gontcharov 2006) for the stars it lacks.

		Downloads either table if its path is missing. Returns (field, errors, n_untraced), the last counting stars within HOP today, in either catalog, with no usable radial velocity, which therefore cannot be traced. Per-star errors: velocity error per axis from parallax, proper-motion and RV errors.

		Gaia saturates on the brightest stars and measures radial velocities only for stars bright and cool enough, so Sirius, Vega and α Cen A, the neighbors everyone knows, come from Hipparcos, less precisely. A Hipparcos star counts only where no usable Gaia star lies within MATCH of it at the Gaia epoch; Gaia's own cross-match misses about 800 such pairs, most of them fast movers.

		A Gaia radial velocity is usable at rv_expected_sig_to_noise >= RV_SNR, and any star only with a space velocity the Galaxy can hold. Low-S/N spectra leave spurious velocities of hundreds of km/s in DR3 (Katz et al. 2023, Sect. 9, who recommend the S/N cut for high-velocity work), and this selection is dominated by fast stars, since they can arrive from farther away. A spurious radial velocity is worse than noise here: real proper motion plus a huge radial speed makes a star's track run almost through the Sun, a fake close encounter. Stars faster than V_ESC about the Galactic center (escape is ~530 km/s, Deason et al. 2019) are such errors.
		"""
		import csv
		if not os.path.exists(gaia_path):
			fetch_gaia(gaia_path)
		if not os.path.exists(hipparcos_path):
			fetch_hipparcos(hipparcos_path)
		gaia = list(csv.DictReader(open(gaia_path)))
		hip = list(csv.DictReader(open(hipparcos_path)))
		usable = np.array([bool(r["radial_velocity"]) and float(r["rv_expected_sig_to_noise"]) >= RV_SNR for r in gaia])
		g = _kinematics(gaia, "ra dec parallax parallax_error pmra pmra_error pmdec pmdec_error radial_velocity radial_velocity_error", 0.0)
		h = _kinematics(hip, "RArad DErad Plx e_Plx pmRA e_pmRA pmDE e_pmDE RV e_RV", GAIA_YEAR - HIPPARCOS_YEAR)
		has_rv = np.isfinite(h.vel[:, 0])
		fill = has_rv & ~_near(h.unit, g.unit[usable]) & (h.D <= reach_radius(np.linalg.norm(np.nan_to_num(h.vel), axis=1)))
		near = lambda x: x.D < HOP
		n_untraced = ((near(g) & ~usable & ~_near(g.unit, h.unit[fill])).sum()
					  + (near(h) & ~has_rv & ~_near(h.unit, g.unit)).sum())
		ids = np.concatenate([[f"Gaia DR3 {r['source_id']}" for r in gaia], [f"HIP {r['HIP']}" for r in hip]])
		pick = np.concatenate([usable, fill])
		pos, vel, D = (np.concatenate([getattr(g, k), getattr(h, k)])[pick] for k in ("pos", "vel", "D"))
		s_v, s_plx = (np.concatenate([getattr(g, k), getattr(h, k)])[pick] for k in ("s_v", "s_plx"))
		bound = np.flatnonzero(np.linalg.norm(vel @ ICRS_TO_GAL.T + V_SUN_GAL, axis=1) < V_ESC)
		return (cls(pos, vel, D, ids[pick]).take(bound), Errors(star_v=s_v, parallax_mas=s_plx).take(bound),
				int(n_untraced))

	@classmethod
	def synthetic(cls, rng, density=0.08,
				  v_mean=(-11.1, -12.2, -7.3), v_sigma=(35.0, 25.0, 18.0)):
		"""Uniform stars, Gaussian velocities relative to the Sun, complete for visitors faster than V_FLOOR.

		Defaults: ~0.08 systems/pc^3 (the 10 pc census) and the local disk velocity ellipsoid, offset by the Sun's peculiar motion. A Poisson field out to R_MAX thinned to each star's reach_radius: draw velocities for the full ball, keep a star with the chance its radius falls inside its reach.
		"""
		n = rng.poisson(density * 4 / 3 * math.pi * R_MAX ** 3)
		vel = rng.normal(v_mean, v_sigma, size=(n, 3))
		reach = reach_radius(np.linalg.norm(vel, axis=1))
		keep = rng.random(n) < (reach / R_MAX) ** 3
		vel, reach = vel[keep], reach[keep]
		r = reach * rng.random(len(vel)) ** (1 / 3)
		d = rng.normal(size=(len(vel), 3))
		pos = d / np.linalg.norm(d, axis=1)[:, None] * r[:, None]
		return cls(pos, vel, r)


def _kinematics(rows, columns, years):
	"""Unit vectors moved on by `years` of proper motion, distances, positions, space velocities (NaN without a radial velocity), and per-axis velocity and parallax errors, from catalog columns named in the order ra dec parallax and errors, pmra (times cos dec) and error, pmdec and error, radial velocity and error."""
	ra, dec, plx, s_plx, pmra, s_pmra, pmdec, s_pmdec, rv, s_rv = (
		np.array([float(r[c] or "nan") for r in rows]) for c in columns.split())
	mas = math.radians(1 / 3.6e6)
	dec = np.radians(dec)
	ra = np.radians(ra) + pmra * years * mas / np.cos(dec)
	dec = dec + pmdec * years * mas
	unit = np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], 1)
	e_ra = np.stack([-np.sin(ra), np.cos(ra), np.zeros_like(ra)], 1)
	e_dec = np.stack([-np.sin(dec) * np.cos(ra), -np.sin(dec) * np.sin(ra), np.cos(dec)], 1)
	k = 4.740470446
	D = 1000.0 / plx
	vel = (k * pmra / plx)[:, None] * e_ra + (k * pmdec / plx)[:, None] * e_dec + rv[:, None] * unit
	s_ra = k * np.hypot(s_pmra / plx, pmra * s_plx / plx ** 2)
	s_dec = k * np.hypot(s_pmdec / plx, pmdec * s_plx / plx ** 2)
	s_v = np.sqrt((s_ra ** 2 + s_dec ** 2 + s_rv ** 2) / 3)
	return SimpleNamespace(unit=unit, D=D, pos=D[:, None] * unit, vel=vel, s_v=s_v, s_plx=s_plx)


def _near(a, b):
	"""Which unit vectors in a have one in b within MATCH."""
	if not len(a) or not len(b):
		return np.zeros(len(a), bool)
	d, _ = cKDTree(b).query(a, distance_upper_bound=math.radians(MATCH / 3600))
	return np.isfinite(d)


GAIA_TAP = "https://gea.esac.esa.int/tap-server/tap/sync"


def fetch_gaia(path):
	"""Gaia DR3 sources with parallax S/N > 10 that are within reach_radius of their speed, plus all within HOP (to count those without radial velocity), as CSV."""
	from urllib.parse import urlencode
	from urllib.request import urlopen
	k = 4.740470446
	w = (f"SQRT(POWER({k}*pmra/parallax, 2) + POWER({k}*pmdec/parallax, 2) "
		 "+ POWER(radial_velocity, 2))")
	query = ("SELECT source_id, ra, dec, parallax, parallax_error, pmra, pmra_error, "
			 "pmdec, pmdec_error, radial_velocity, radial_velocity_error, rv_expected_sig_to_noise, "
			 "phot_g_mean_mag "
			 f"FROM gaiadr3.gaia_source WHERE parallax > {1000.0 / R_MAX} "
			 f"AND parallax_over_error > 10 AND (parallax > {1000.0 / HOP} OR "
			 f"(radial_velocity IS NOT NULL AND "
			 f"1000/parallax < {HOP} * (1 + {w} / {V_FLOOR}) + {TOL}))")
	data = urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv",
					  "QUERY": query}).encode()
	with urlopen(GAIA_TAP, data=data, timeout=600) as r:
		body = r.read()
	if not body.startswith(b"source_id"):
		raise RuntimeError("unexpected Gaia archive response")
	os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
	with open(path, "wb") as f:
		f.write(body)


VIZIER_TAP = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"


def fetch_hipparcos(path):
	"""Hipparcos (new reduction, I/311) stars with parallax S/N > 10 within R_MAX, with Pulkovo radial velocities (III/252) where they exist, as CSV."""
	from urllib.parse import urlencode
	from urllib.request import urlopen
	query = ('SELECT h.HIP, h.RArad, h.DErad, h.Plx, h.e_Plx, h.pmRA, h.e_pmRA, h.pmDE, h.e_pmDE, p.RV, p.e_RV '
			 'FROM "I/311/hip2" AS h LEFT JOIN "III/252/table8" AS p ON h.HIP = p.HIP '
			 f"WHERE h.Plx > {1000.0 / R_MAX} AND h.Plx > 10 * h.e_Plx")
	data = urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "MAXREC": 200000,
					  "QUERY": query}).encode()
	with urlopen(VIZIER_TAP, data=data, timeout=600) as r:
		body = r.read()
	if not body.startswith(b"HIP,RArad"):
		raise RuntimeError("unexpected VizieR response")
	os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
	with open(path, "wb") as f:
		f.write(body)


def hops(field, v_inf, sign):
	"""Shortest hop (pc from the Sun at the meeting) at which a visitor at v_inf meets each star; inf where it cannot pass within TOL.

	A star can be met at speed v_inf only if some direction puts the relative velocity u = w - V along +-s_hat. Directions to the sphere of possible u (center w, radius v_inf) form a cone about w_hat of half-angle asin(v_inf/|w|) (all directions if |w| <= v_inf), so the smallest achievable miss is D sin(max(0, angle(axis, w) - half-angle)). An exact meeting at time t has |s + w t| = v_inf |t|, a quadratic in x = 1/t (as in aim_at); the hop is v_inf |t| for the earliest root on this side of t = 0, or at the grazing time -B / 2A for stars only reachable within TOL. All of this depends on v_inf alone, never on the visitor's direction.
	"""
	axis = field.pos / field.D[:, None] * (-sign)
	wn = np.linalg.norm(field.vel, axis=1)
	sw = np.einsum("ij,ij->i", field.pos, field.vel)
	ang = np.arccos(np.clip(-sign * sw / (field.D * wn), -1, 1))
	half = np.arcsin(np.minimum(1.0, v_inf / wn))
	gmin = np.where(wn <= v_inf, 0.0, np.maximum(0.0, ang - half))
	ok = field.D * np.sin(np.minimum(gmin, 0.5 * math.pi)) < TOL
	A, B, C = field.D ** 2 / KMS_PC_MYR ** 2, 2 * sw / KMS_PC_MYR, wn ** 2 - v_inf ** 2
	disc = B * B - 4 * A * C
	sq = np.sqrt(np.maximum(disc, 0.0))
	x = np.stack([(-B + sq) / (2 * A), (-B - sq) / (2 * A)])
	x = np.where(sign * x > 0, np.abs(x), 0.0).max(axis=0)   # earliest meeting: largest |1/t|
	graze = np.where(sign * B < 0, np.abs(B) / (2 * A), 0.0)
	x = np.where(disc >= 0, x, graze)
	with np.errstate(divide="ignore"):
		return np.where(ok & (x > 0), v_inf * KMS_PC_MYR / x, np.inf)


def end_weights(field, v_inf, sign):
	"""Look-elsewhere weights for one end: 1/rank by hop among stars met within HOP.

	Depends only on v_inf, so the test stays valid while the budget is not wasted on impossible stars; the field must be complete at this speed.
	"""
	assert v_inf >= V_FLOOR, f"field is complete only for v_inf >= {V_FLOOR} km/s"
	L = hops(field, v_inf, sign)
	ok = L <= HOP
	rank = np.empty(len(L))
	rank[np.argsort(L)] = np.arange(1, len(L) + 1)
	w = np.where(ok, 1 / rank, 0.0)
	return w / w.sum()


def miss(pos, vel, V, sign):
	"""Miss distance and time of stars (pos, vel) against visitor velocities V, broadcast over leading axes; sign -1 inbound, +1 outbound.

	Stars whose closest approach falls on the wrong side of t = 0 get the t = 0 separation (the visitor never came nearer).
	"""
	u = vel - V
	t = -np.sum(pos * u, axis=-1) / (KMS_PC_MYR * np.sum(u * u, axis=-1))
	t = np.where(sign * t > 0, t, 0.0)
	return np.linalg.norm(pos + u * (t * KMS_PC_MYR)[..., None], axis=-1), t


def rice_cdf(x, b):
	"""P(|b + unit 2D Gaussian| <= x), elementwise. SciPy's series slows without bound as b grows; beyond b = 100 the Gaussian limit about sqrt(b^2 + 1/2) is within 1.3% of it."""
	x, b = np.broadcast_arrays(x, b)
	far = b > 100
	out = np.empty(b.shape)
	out[~far] = rice.cdf(x[~far], b[~far])
	out[far] = ndtr(x[far] - np.sqrt(b[far] ** 2 + 0.5))
	return out


def encounters(field, V, sign):
	return miss(field.pos, field.vel, V, sign)


def observe(d, sig, rng):
	"""A measured miss: the transverse miss is a 2D vector, so add 2D Gaussian error."""
	return np.hypot(d + rng.normal(size=d.shape) * sig, rng.normal(size=d.shape) * sig)


def best(p, weights, info):
	"""Weighted Bonferroni over stars: (p_end, star index, info)."""
	with np.errstate(divide="ignore"):
		ratio = np.where(weights > 0, p / weights, np.inf)
	k = int(np.argmin(ratio))
	return min(1.0, float(ratio[k])), k, {**info, "p": p, "weights": weights}


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
		d = observe(d, sig, rng)
	# Only stars passed within a fraction of their distance can be rare; others (including stars never approached, d = D) get p = 1, which is exact for d >= D and conservative otherwise.
	p = np.ones_like(d)
	if weights is None:
		weights = end_weights(field, v_inf, sign)
	near = np.flatnonzero((d + 2 * sig < SCREEN * field.D) & (weights > 0))
	# A star sets the score only if p / w < 1. The kernel is at least P(|noise| <= 3 sigma) = 0.989 out to d - 3 sigma, so p >= 0.989 p_geom(d - 3 sigma): convolve only stars whose bound (with slack for the quadrature) is below their weight.
	lo = p_values(field.take(near), v_inf, np.maximum(d[near] - 3 * sig[near], 0.0), sign)
	near = near[0.9 * lo < weights[near]]
	if len(near):
		dn, sn = d[near], sig[near]
		r = np.linspace(0, 1, N_NOISE + 1) * (dn + 8 * sn)[:, None]   # true miss grid
		rows = np.repeat(near, N_NOISE + 1)
		pg = p_values(field.take(rows), v_inf, r.ravel(), sign).reshape(len(near), N_NOISE + 1)
		rm = 0.5 * (r[:, 1:] + r[:, :-1])
		kern = rice_cdf(dn[:, None] / sn[:, None], rm / sn[:, None])
		p[near] = np.sum(np.diff(pg, axis=1) * kern, axis=1)
	return best(p, weights, {"d": d, "sigma": sig, "t": t})


def turn_circle(V_in, V_out):
	"""Departures the Sun allows an arrival V_in: the circle at V_out's angle from it, as a function of the turn's orientation x (x = 0 is V_out)."""
	v = np.linalg.norm(V_out)
	a = V_in / np.linalg.norm(V_in)
	c = a @ V_out / v
	n1 = V_out / v - c * a
	s = np.linalg.norm(n1)
	n1 /= s
	n2 = np.cross(a, n1)
	return lambda x: v * (c * a + s * (np.cos(x)[..., None] * n1 + np.sin(x)[..., None] * n2))


def turn_score(field, V_in, V_out, errors, rng=None, weights=None):
	"""Score the departure end given the arrival. Returns (p_end, best star index, info) like end_score.

	The Sun turns an arrival through an angle set by the perihelion distance, and discovery, not design, favors close perihelia; so the angle is taken as observed, and the null spins the turn uniformly about the arrival direction. The departure is then a point on a circle of directions, and p_j is the fraction of that circle on which star j would be passed at least as closely as observed, convolved with the Rice distribution of the measured miss as in end_score. A star's closeness on the circle is integrated finely across a window about its closest orientation (from the local curvature where the window is narrower than a coarse step), and coarsely everywhere else, so a second close stretch is never dropped. A path the Sun barely turns has a circle that passes a star straight ahead all the way round: p = 1, and a destination collinear with the origin earns nothing.
	"""
	v_inf = float(np.linalg.norm(V_out))
	if weights is None:
		weights = end_weights(field, v_inf, +1)
	turn = turn_circle(V_in, V_out)
	d, t = encounters(field, V_out, +1)
	sig = errors.miss_sigma(field.D, t)
	if rng is not None:
		d = observe(d, sig, rng)
	p = np.ones_like(d)
	step = 2 * math.pi / N_TURN
	x = np.arange(N_TURN) * step
	j = np.flatnonzero(weights > 0)
	r = miss(field.pos[j, None], field.vel[j, None], turn(x), +1)[0]
	# A star sets the score only if p / w < 1, so bound p from below and drop the rest; dropping sets p to 1, so a loose bound costs power, never validity. Two bounds: where the circle passes within d - 3 sigma the kernel is at least 0.989 (count those coarse points, less one step at every edge of a stretch); and nowhere is the kernel below its value at the circle's farthest pass (padded by the largest change between neighboring points).
	inside = r < (d[j] - 3 * sig[j])[:, None]
	edges = (inside != np.roll(inside, 1, axis=1)).sum(axis=1)
	far = r.max(axis=1) + np.abs(np.diff(r, axis=1)).max(axis=1)
	lower = np.maximum(0.989 * np.maximum(inside.sum(axis=1) - edges, 0) * step / (2 * math.pi), rice_cdf(d[j] / sig[j], far / sig[j]))
	keep = lower < weights[j]
	j, r = j[keep], r[keep]
	pos, vel, dj, sj = field.pos[j], field.vel[j], d[j], sig[j]
	cut = (dj + 8 * sj)[:, None]   # beyond this true miss the Rice kernel is negligible

	def kernel(r):
		"""Rice CDF at the observed miss for true misses r, one row per star."""
		out = np.zeros_like(r)
		m = r < cut
		out[m] = rice_cdf(np.broadcast_to((dj / sj)[:, None], r.shape)[m], (r / sj[:, None])[m])
		return out

	x0 = x[np.argmin(r, axis=1)]
	f = lambda x: miss(pos, vel, turn(x), +1)[0] ** 2
	h = 1e-5
	for _ in range(8):
		fm, f0, fp = f(x0 - h), f(x0), f(x0 + h)
		curv = fm - 2 * f0 + fp
		x0 -= np.where(curv > 0, np.clip(h * (fp - fm) / (2 * np.where(curv > 0, curv, 1.0)), -step, step), 0.0)
	fm, f0, fp = f(x0 - h), f(x0), f(x0 + h)
	c = np.maximum((fm - 2 * f0 + fp) / (2 * h * h), 1e-300)
	narrow = 1.5 * np.sqrt(np.maximum(cut[:, 0] ** 2 - f0, 0.0) / c)
	W = np.minimum(math.pi, np.where(narrow < step, narrow, ((r < cut).sum(axis=1) + 1) * step))
	xs = x0[:, None] + W[:, None] * np.linspace(-1, 1, N_FINE)
	fine = kernel(miss(pos[:, None], vel[:, None], turn(xs), +1)[0])
	outside = np.abs((x - x0[:, None] + math.pi) % (2 * math.pi) - math.pi) > W[:, None]
	p[j] = np.minimum(1.0, (np.trapezoid(fine, xs, axis=1) + step * (kernel(r) * outside).sum(axis=1)) / (2 * math.pi))
	return best(p, weights, {"d": d, "sigma": sig, "t": t})


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


def null_route(field, v_inf, n, errors, seed=0, q_max=30.0, directions=None):
	"""Route scores of n visitors with random arrival, perihelion, plane.

	Impact parameters are uniform in area out to perihelion q_max. Arrival directions are isotropic, or drawn from `directions` (unit vectors, e.g. `kinematic_directions`).
	"""
	rng = np.random.default_rng(seed)
	g = 2 / (v_inf / V_UNIT_KMS) ** 2        # b^2 = q^2 + g q (AU)
	b2max = q_max ** 2 + g * q_max
	# only weighted stars can set the score, so each end works on its candidates alone
	ends = []
	for sign in (-1, +1):
		w = end_weights(field, v_inf, sign)
		idx = np.flatnonzero(w)
		ends.append((field.take(idx), errors.take(idx), w[idx]))
	(f_in, e_in, w_in), (f_out, e_out, w_out) = ends
	out = np.empty((n, 3))
	for k in range(n):
		V_in = v_inf * (random_unit(rng) if directions is None else directions[rng.integers(len(directions))])
		b2 = b2max * rng.random()
		q = 0.5 * (-g + math.sqrt(g * g + 4 * b2))
		V_out = sun_turn(V_in, q, rng.uniform(0, 2 * math.pi))
		p_in, _, _ = end_score(f_in, V_in, -1, e_in, rng, weights=w_in)
		p_out, _, _ = turn_score(f_out, V_in, V_out, e_out, rng, weights=w_out)
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
					key=lambda x: -abs(x)):   # earliest meeting, the shorter hop
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



def reachable(field, ranks, v_inf, sign, exclude=()):
	"""First star at or after each hop rank that v_inf can meet exactly: (index, rank)."""
	srt = np.argsort(hops(field, v_inf, sign))
	out = []
	for r in ranks:
		for kk in range(r - 1, len(srt)):
			if srt[kk] in exclude:
				continue
			V, t = aim_at(field.pos[srt[kk]], field.vel[srt[kk]], v_inf, sign)
			if V is not None:
				out.append((srt[kk], kk + 1))
				break
	return out


def planted_bits(V_in, V_out, field, errors, draws=31, seed=0):
	"""Median bits (origin, destination, combined) of a planned route over measurement-noise draws: the visitor passes exactly through both stars, the catalog sees it through its errors."""
	rng = np.random.default_rng(seed)
	rows = []
	for _ in range(draws):
		p_in, _, _ = end_score(field, V_in, -1, errors, rng)
		p_out, _, _ = turn_score(field, V_in, V_out, errors, rng)
		rows.append((p_in, p_out, combine(p_in, p_out)))
	return np.median(-np.log2(np.maximum(rows, 1e-300)), axis=0)

def kinematic_directions(field, v_inf, width=0.15):
	"""Directions of motion of the field's stars whose speed relative to the Sun is within `width` of v_inf: the arrivals interstellar objects make if, natural or not, they move as stars do. A field's reach depends on speed alone, so within the band the directions are unbiased."""
	w = np.linalg.norm(field.vel, axis=1)
	band = np.abs(w / v_inf - 1) < width
	return field.vel[band] / w[band, None]


def _null_chunk(args):
	return null_route(*args)


def null_route_parallel(field, v_inf, n, errors, pool, chunks=72, seed=0, directions=None):
	per = n // chunks
	parts = pool.map(_null_chunk, [(field, v_inf, per, errors, seed + c, 30.0, directions) for c in range(chunks)])
	return np.vstack(parts)
