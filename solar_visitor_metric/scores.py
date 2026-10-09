"""Visit and assist scores for planar paths, and their null distribution.

Visit score: for every ring crossing, p = chance a random planet phase puts the planet at least this close (in planetocentric impact parameter b_p). Crossings are folded into one p per planet (`Path.planet_ps`). The score asks "how surprising is it that k planets were this close", for the best k: min over k of P(Binomial(n, p_(k)) >= k), times n for trying every k. This rewards several close planets without one very close pass being diluted by the many distant ones (as Fisher's method would).

Assist score: how much the planets changed the outbound asymptote, as |v_out - v_out(Sun only)| / v_inf. Its null tail comes from Monte Carlo.

Paths are planar patched conics tracked as (E, h, omega, t_peri): energy, signed angular momentum, perihelion direction and time. A planet's phase enters only through its along-track offset at the crossing, which sets b_p. Two nulls: independent random phases per crossing (`fly`), or real J2000 phases on a common clock, swept over heading and epoch (`fly_timed`).
"""
import math
import random
from dataclasses import dataclass, field

import numpy as np
from scipy.stats import beta, binom

from .constants import PLANETS
from .flyby import Crossing, crossings, ecc, local_velocity, state_after

DEFLECT_SOI = 3.0   # apply a flyby's deflection when |b_p| < this many r_SOI


def true_anomaly(E, h, r, sign):
	e = ecc(E, h)
	c = max(-1.0, min(1.0, (h * h / r - 1) / e))
	return sign * math.acos(c)


def polar(h, omega, f):
	return omega + (f if h > 0 else -f)


def asymptote(E, h, omega):
	"""Outbound asymptotic velocity vector (vx, vy), or None if bound."""
	if E <= 0:
		return None
	e = ecc(E, h)
	f_inf = math.acos(-1 / e)
	v = math.sqrt(2 * E)
	th = polar(h, omega, f_inf)   # far away, motion is radial at f_inf
	return v * math.cos(th), v * math.sin(th)


@dataclass
class Path:
	v_inf: float
	b: float
	records: list = field(default_factory=list)   # (planet, sign, p, b_p, dv/v)
	segments: list = field(default_factory=list)  # (E, h, omega, f0, f1)
	marks: dict = field(default_factory=dict)     # planned (planet, sign) -> xy
	impact: bool = False
	captured: bool = False
	heading: float = 0.0  # perihelion direction of the Sun-only path
	E: float = 0.0
	h: float = 0.0
	omega: float = 0.0

	def planet_ps(self, records=None):
		"""One p per planet: chance that any of its crossings came this close.

		p = min(1, n_c p_min) over that planet's n_c crossings. On a common clock one phase sets the planet's offset at every crossing, so its windows are arcs of one circle: the chance of either is their union, n_c p_min when they are apart and less when they overlap. The independent-crossing fold 1 - (1 - p_min)^n_c falls short of the union and overstated common-clock nulls by up to 10%.
		"""
		best = {}
		for r in (self.records if records is None else records):
			pm, n = best.get(r[0], (1.0, 0))
			best[r[0]] = (min(pm, r[2]), n + 1)
		return [min(1.0, n * pm) for pm, n in best.values()]

	def visit_bits(self):
		return order_bits(self.planet_ps())

	def excess_visit_bits(self, assist_dv=0.01):
		"""Visit score over crossings that did not themselves bend the path.

		Planet phases are independent under the null, so this part is independent of the assist and its bits add to the assist's bits.
		"""
		return order_bits(self.planet_ps([r for r in self.records if not r[4] >= assist_dv]))

	def assist(self):
		"""|v_out - v_out(Sun only)| / v_inf; inf if captured or impacted."""
		if self.impact or self.captured:
			return math.inf
		a = asymptote(self.E, self.h, self.omega)
		a0 = asymptote(0.5 * self.v_inf ** 2, self.b * self.v_inf, self.heading)
		return math.hypot(a[0] - a0[0], a[1] - a0[1]) / self.v_inf


def order_bits(ps):
	"""-log2 of the exact null chance of a planet pattern this strong: the statistic t = min_k P(at least k of n uniforms <= k-th smallest p) counts one very close pass and several fairly close ones alike, and choosing the best k is paid for by t's own null distribution (`order_null`), not by a Bonferroni factor of n."""
	n = len(ps)
	if n == 0:
		return 0.0
	p = np.sort(np.asarray(ps))
	t = binom.sf(np.arange(n), n, p).min()
	return max(0.0, -math.log2(order_null(t, n)))


def order_null(t, n):
	"""P(min_k P(at least k of n <= U_(k)) <= t) for n independent uniforms, exactly.

	The minimum is at most t when some order statistic U_(k) falls below c_k, the level-t quantile of its own Beta(k, n - k + 1) law; c_k rises with k, so the complement is the chance that at most k - 1 points lie below every c_k, summed over how many fall in each interval between the c's (multinomial weights). Planet p's at least as large as uniform under the null only make this conservative. Below t = 1e-10 the complement cancels in floating point, and the Bonferroni bound n t is used.
	"""
	if t >= 1:
		return 1.0
	if t < 1e-10:
		return n * t
	edges = np.concatenate([[0.0], beta.ppf(t, np.arange(1, n + 1), np.arange(n, 0, -1)), [1.0]])
	ways = np.zeros(n + 1)
	ways[0] = 1.0
	for i in range(1, n + 2):
		width = edges[i] - edges[i - 1]
		cap = min(i - 1, n) if i <= n else n
		ways = np.array([sum(ways[j] * width ** (m - j) / math.factorial(m - j) for j in range(m + 1)) if m <= cap else 0.0
						 for m in range(n + 1)])
	return 1.0 - math.factorial(n) * ways[n]


# J2000 mean longitudes (deg), JPL approximate elements; time unit yr/(2 pi)
LAMBDA_J2000 = (252.25, 181.98, 100.46, 355.45, 34.40, 49.94, 313.23, 304.88)


def _t_from_peri(E, h, r, sign):
	"""Time from perihelion to radius r on the conic (sign -1: before)."""
	e = ecc(E, h)
	if E > 0:
		A = 1 / (2 * E)
		F = math.acosh(max(1.0, (1 + r / A) / e))
		return sign * (e * math.sinh(F) - F) * A ** 1.5
	A = -1 / (2 * E)
	Ec = math.acos(max(-1.0, min(1.0, (1 - r / A) / e)))
	return sign * (Ec - e * math.sin(Ec)) * A ** 1.5


def _walk(v_inf, b, offset, heading=0.0, t_peri=0.0):
	"""Walk one planar patched-conic path through the rings.

	offset(i, sign, crossing, theta, t) -> (p, b_p) gives the null p and planetocentric impact parameter at each crossing of ring i, where theta is the visitor's polar angle there and t the time.
	"""
	E, h, omega = 0.5 * v_inf ** 2, b * v_inf, heading
	path = Path(v_inf, b, heading=heading)
	f_cur = -math.acos(-1 / ecc(E, h))
	tp = t_peri
	seq, bound = crossings(E, h, math.inf, -1)
	k = 0
	while k < len(seq):
		i, sg = seq[k]
		k += 1
		p = PLANETS[i]
		cr = Crossing(i, *local_velocity(E, h, p.a, sg))
		if cr.v_r == 0:
			continue
		f_here = true_anomaly(E, h, p.a, sg)
		theta = polar(h, omega, f_here)
		t_c = tp + _t_from_peri(E, h, p.a, sg)
		pv, b_p = offset(i, sg, cr, theta, t_c)
		if abs(b_p) < cr.b_surface():
			path.records.append((i, sg, pv, b_p, math.nan))
			path.segments.append((E, h, omega, f_cur, f_here))
			path.impact = True
			path.E, path.h, path.omega = E, h, omega
			return path
		dvv = 0.0
		if abs(b_p) < DEFLECT_SOI * p.soi:
			E2, h2, sg2, dv = state_after(cr, b_p)
			dvv = dv / cr.v_mag
			path.segments.append((E, h, omega, f_cur, f_here))
			f2 = true_anomaly(E2, h2, p.a, sg2)
			E, h, omega = E2, h2, theta - (f2 if h2 > 0 else -f2)
			f_cur = f2
			tp = t_c - _t_from_peri(E, h, p.a, sg2)
			seq, bound = crossings(E, h, p.a, sg2)
			if sg2 < 0 and not bound:
				# the return crossing of this planet's ring, after perihelion
				at = next((n for n, (j, s) in enumerate(seq)
						   if s > 0 and PLANETS[j].a > p.a), len(seq))
				seq.insert(at, (i, +1))
			k = 0
		path.records.append((i, sg, pv, b_p, dvv))
	path.E, path.h, path.omega = E, h, omega
	if E <= 0:
		path.captured = True
		path.segments.append((E, h, omega, f_cur, f_cur + 2 * math.pi))
	else:
		path.segments.append((E, h, omega, f_cur, math.acos(-1 / ecc(E, h))))
	return path


def fly(v_inf, b, plan=None, rng=None, default_p=None):
	"""Path under the independent random-phase null.

	plan: {(planet, sign): b_p} forces those crossings (sign -1 inbound, +1 outbound). Other crossings draw p ~ U(0,1) from `rng`, or use `default_p` if given (a fixed, typical miss for illustrations).
	"""
	plan = plan or {}
	marks = {}

	def offset(i, sg, cr, theta, t):
		dens = cr.density
		if (i, sg) in plan:
			a = PLANETS[i].a
			marks[(i, sg)] = (a * math.cos(theta), a * math.sin(theta))
			b_p = plan[(i, sg)]
			return min(1.0, 2 * abs(b_p) * dens), b_p
		pv = default_p if default_p is not None else rng.random()
		side = 1 if (default_p is not None or rng.random() < 0.5) else -1
		return pv, side * pv / (2 * dens)

	path = _walk(v_inf, b, offset)
	path.marks = marks
	return path


def fly_timed(v_inf, b, heading, t_peri, lambdas=LAMBDA_J2000):
	"""Path with planets at real phases on a common clock.

	The arrival is rotated by `heading` (perihelion direction) and reaches perihelion at `t_peri`. The planet's along-track offset s from the visitor sets b_p = s v_r / |w|; the null p is |s| / (pi a), the same quantity `fly` draws at random.
	"""
	def offset(i, sg, cr, theta, t):
		a = PLANETS[i].a
		ds = (math.radians(lambdas[i]) + t * a ** -1.5 - theta + math.pi) % (2 * math.pi) - math.pi
		return min(1.0, abs(ds) / math.pi), ds * a * cr.v_r / cr.w_mag

	return _walk(v_inf, b, offset, heading, t_peri)


def b_max(v_inf, q_max=None):
	q = q_max or PLANETS[-1].a
	return math.sqrt(q * q + 2 * q / (v_inf * v_inf))


def null_sample(v_inf, n, seed=0, visits=False, timed=False, span_yr=1e5):
	"""Monte Carlo of the null: b uniform over q < Neptune, both senses.

	timed=True uses real phases on a common clock, with heading and epoch (over span_yr years) drawn at random. Returns columns (assist, captured, visit bits or nan); visit bits are optional as they cost ~30x more.
	"""
	rng = random.Random(seed)
	B = b_max(v_inf)
	out = np.empty((n, 3))
	for k in range(n):
		b = rng.uniform(-B, B)
		if timed:
			path = fly_timed(v_inf, b, rng.uniform(0, 2 * math.pi),
							 rng.uniform(0, span_yr * 2 * math.pi))
		else:
			path = fly(v_inf, b, rng=rng)
		out[k] = (path.assist(), path.captured,
				  path.visit_bits() if visits else math.nan)
	return out


def _null_chunk(args):
	return null_sample(*args)


def null_parallel(v_inf, n, pool, chunks=72, seed=0, visits=False, timed=False):
	per = n // chunks
	parts = pool.map(_null_chunk, [(v_inf, per, seed + c, visits, timed) for c in range(chunks)])
	return np.vstack(parts)
