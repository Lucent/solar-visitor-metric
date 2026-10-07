"""Semi-analytic encounter-chain probabilities under the random-phase null.

Null hypothesis: planets on circular coplanar orbits with independent, uniformly random phases. A visitor arrives in-plane with speed v_inf and signed impact parameter b. For each ring crossing in time order we know the probability of a qualifying encounter in closed form; after an encounter we integrate over where in the window it happened (which sets the deflection) and recurse on the deflected conic. No sampling noise, no timing needed.
"""
import math
from dataclasses import dataclass, field

import numpy as np

from .constants import PLANETS
from .flyby import Crossing, crossings, local_velocity, state_after


# --- encounter criteria: Crossing -> largest qualifying |b_p| --------------

class dv_criterion:
	"""Flyby changes heliocentric velocity by >= frac of its magnitude."""
	def __init__(self, frac):
		self.frac = frac
		self.label = f"|dv|/v >= {frac:g}"

	def __call__(self, cr):
		return cr.b_for_dv(self.frac)


class soi_criterion:
	"""Closest approach inside scale * Laplace sphere of influence."""
	def __init__(self, scale=1.0):
		self.scale = scale
		self.label = f"r_p < {scale:g} r_SOI"

	def __call__(self, cr):
		return cr.b_for_rp(self.scale * cr.p.soi)


# --- recursion ------------------------------------------------------------

@dataclass
class Outcome:
	"""Probability distribution over the number of qualifying encounters."""
	dist: np.ndarray                 # dist[n] = P(n encounters), last bin is ">= K"
	impact: float = 0.0              # P(hits a planet surface)
	capture: float = 0.0             # P(ends bound to the Sun)
	first: np.ndarray = field(default=None)   # P(first encounter is planet i)
	pairs: dict = field(default_factory=dict) # (i, j) -> P(first two are i then j)

	def at_least(self, n):
		return float(self.dist[n:].sum())


def _gauss_log(lo, hi, n):
	"""Nodes/weights for integrating over b in [lo, hi] with log spacing."""
	x, w = np.polynomial.legendre.leggauss(n)
	u0, u1 = math.log(lo), math.log(hi)
	u = 0.5 * (u1 - u0) * x + 0.5 * (u1 + u0)
	b = np.exp(u)
	return b, w * 0.5 * (u1 - u0) * b


def walk(E, h, r0, sign, crit, K=3, nodes=(16, 8), skip=frozenset(), level=0):
	seq, bound = crossings(E, h, r0, sign, skip)
	n_pl = len(PLANETS)
	out = Outcome(np.zeros(K + 1), first=np.zeros(n_pl))
	survive = 1.0
	for i, sg in seq:
		cr = Crossing(i, *local_velocity(E, h, PLANETS[i].a, sg))
		if cr.v_r == 0:
			continue
		dens = cr.density
		lo, hi = cr.b_surface(), crit(cr)
		p_imp = min(2 * lo * dens, 1.0)
		p_hit = min(2 * max(hi - lo, 0.0) * dens, 1.0 - p_imp)
		if p_hit > 0:
			out.first[i] += survive * p_hit
			if level + 1 >= K or level >= len(nodes):
				out.dist[K if level + 1 >= K else level + 1] += survive * p_hit
			else:
				bs, ws = _gauss_log(lo, hi, nodes[level])
				norm = (p_hit / 2) / ws.sum()   # renormalize after clipping
				for b, wt in zip(bs, ws):
					for b_p in (b, -b):
						E2, h2, sg2, _ = state_after(cr, b_p)
						sub = walk(E2, h2, PLANETS[i].a, sg2, crit, K, nodes,
								   skip | {i}, level + 1)
						m = survive * wt * norm
						out.dist[1:] += m * sub.dist[:-1]
						out.dist[K] += m * sub.dist[-1]
						out.impact += m * sub.impact
						out.capture += m * sub.capture
						if level == 0:
							for j in range(n_pl):
								if sub.first[j]:
									out.pairs[(i, j)] = out.pairs.get((i, j), 0) + m * sub.first[j]
		out.impact += survive * p_imp
		survive *= 1 - p_hit - p_imp
	out.dist[0] += survive
	if bound:
		out.capture += survive
	return out


def from_infinity(v_inf, b, crit, **kw):
	"""Outcome for an in-plane arrival; b > 0 prograde, b < 0 retrograde."""
	return walk(0.5 * v_inf * v_inf, b * v_inf, math.inf, -1, crit, **kw)


# --- sweep over impact parameter ------------------------------------------

def q_to_b(q, v_inf):
	return math.sqrt(q * q + 2 * q / (v_inf * v_inf))


def ring_edges():
	return [0.0] + [p.a for p in PLANETS]


def _one(args):
	v_inf, b, crit, kw = args
	return from_infinity(v_inf, b, crit, **kw)


def sweep(v_inf, crit, n_per_ring=24, mapper=map, **kw):
	"""Integrate over all in-plane arrivals with perihelion inside Neptune.

	Perihelion bins are the rings: (0, Mercury), (Mercury, Venus), ... Within each bin q = top - width * x^2, which removes the integrable 1/sqrt singularity where the conic grazes the outer ring tangentially. `mapper` may be a process pool's map. Returns one dict per ring.
	"""
	x, wx = np.polynomial.legendre.leggauss(n_per_ring)
	x, wx = 0.5 * (x + 1), 0.5 * wx
	edges = ring_edges()
	jobs, meta = [], []
	for k in range(len(PLANETS)):
		lo, hi = edges[k], edges[k + 1]
		for direction in (+1, -1):
			for xi, wi in zip(x, wx):
				q = hi - (hi - lo) * xi * xi
				b = q_to_b(q, v_inf)
				db = wi * 2 * (hi - lo) * xi * (q + 1 / v_inf ** 2) / b
				jobs.append((v_inf, direction * b, crit, kw))
				meta.append((k, db))
	rows = [{"ring": k, "q_lo": edges[k], "q_hi": edges[k + 1], "db": 0.0, "dist": 0,
			 "impact": 0.0, "capture": 0.0, "first": 0, "pairs": {}}
			for k in range(len(PLANETS))]
	for (k, db), o in zip(meta, mapper(_one, jobs)):
		acc = rows[k]
		acc["db"] += db
		acc["dist"] = acc["dist"] + db * o.dist
		acc["impact"] += db * o.impact
		acc["capture"] += db * o.capture
		acc["first"] = acc["first"] + db * o.first
		for key, val in o.pairs.items():
			acc["pairs"][key] = acc["pairs"].get(key, 0) + db * val
	return rows
