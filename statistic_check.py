"""How to combine the planets into one score: power of candidate statistics on random planned tours, and the visitors under each statistic and under alternative nulls.

	.venv/bin/python statistic_check.py > out/statistic.txt      # about 40 minutes on 18 cores

Every statistic is calibrated by Monte Carlo on the same planet null, so each is a valid score; they differ only in which planned paths they catch. Random planned tours (one to four random planets, each at a random distance from 0.03 to 0.5 AU, at two speeds and perihelia) choose the statistic, never the visitors' scores. The visitors are then scored under every statistic, and under nulls that randomize less (orientation as observed, timing only) or differently (arrivals moving as the catalog's stars do, a shorter span of perihelion times).
"""
import math
import sys
from multiprocessing import Pool

import numpy as np
from scipy.stats import binom

from solar_visitor_metric import planets, real, route

N_CURVE = 1_000_000
N_STAT = 2_000_000
CHUNK = 50_000
TOURS = 30
SETTINGS = ((26.3, 0.5), (58.0, 1.0))   # km/s, perihelion AU
A = np.array([0.387, 0.723, 1.0, 1.524, 5.203, 9.537, 19.19, 30.07])
RADIUS_AU = np.array(planets.RADIUS_KM) / 1.495978707e8


def candidates(p, d):
	"""Candidate statistics of passes, larger is more surprising."""
	s = np.sort(p, axis=1)
	return {
		"order statistic over k": -np.log(binom.sf(np.arange(8), 8, s).min(axis=1)),
		"Fisher (all eight)": -np.log(p).sum(1),
		"smallest p": -np.log(s[:, 0]),
		"two smallest p": -np.log(s[:, :2]).sum(1),
		"three smallest p (adopted)": -np.log(s[:, :3]).sum(1),
		"four smallest p": -np.log(s[:, :4]).sum(1),
		"planets with p < 0.05": (p < 0.05).sum(1) - 1e-3 * np.log(s[:, 0]),
		"nearest in planet radii": -np.log((d / RADIUS_AU).min(1)),
		"nearest relative to orbit": -np.log((d / A).min(1)),
		"three nearest relative to orbit": -np.log(np.sort(d / A, axis=1)[:, :3]).sum(1),
	}


def arrivals(rng, u, q, v):
	"""Rotations taking the hyperbola's arrival direction onto each of u (ecliptic), with a uniform roll about it."""
	e = 1 + q * v * v / planets.MU
	a = np.array([1 / e, math.sqrt(1 - 1 / e ** 2), 0.0])
	b = np.cross(a, [0, 0, 1.0]) / np.linalg.norm(np.cross(a, [0, 0, 1.0]))
	Bp = np.stack([a, b, np.cross(a, b)], -1)
	u = u / np.linalg.norm(u, axis=1)[:, None]
	x = np.cross(u, rng.normal(size=u.shape))
	x /= np.linalg.norm(x, axis=1)[:, None]
	return np.einsum("nab,cb->nac", np.stack([u, x, np.cross(u, x)], -1), Bp)


def _chunk(args):
	q, v, orbits, m, seed, span, directions, R = args
	rng = np.random.default_rng(seed)
	tp = orbits.t0 + (rng.random(m) - 0.5) * span
	if R is not None:
		Rs = np.broadcast_to(R, (m, 3, 3))
	elif directions is not None:
		Rs = arrivals(rng, directions[rng.integers(len(directions), size=m)], q, v)
	else:
		Rs = planets.random_rotations(rng, m)
	return planets.passes(q, v, Rs, tp, orbits)[0]


class Calibrated:
	"""Each candidate statistic's null distribution for one null."""

	def __init__(self, pool, q, v, orbits, span=planets.SPAN, directions=None, R=None):
		draw = lambda n, s: np.vstack(pool.map(_chunk, [(q, v, orbits, min(CHUNK, n - k), (s, k), span, directions, R) for k in range(0, n, CHUNK)]))
		self.curve = np.sort(draw(N_CURVE, 1), axis=0)
		d = draw(N_STAT, 2)
		self.T = {name: np.sort(x) for name, x in candidates(self.p(d), d).items()}

	def p(self, d):
		d = np.atleast_2d(d)
		out = np.empty(d.shape)
		for j in range(8):
			c = self.curve[:, j]
			count = np.searchsorted(c, d[:, j], side="right")
			out[:, j] = np.where(count >= planets.FIT, count / len(c), planets.FIT / len(c) * (d[:, j] / c[planets.FIT - 1]) ** 2)
		return out

	def bits(self, d):
		d = np.atleast_2d(d)
		return {name: -math.log2((len(self.T[name]) - np.searchsorted(self.T[name], x[0], side="left") + 1) / (len(self.T[name]) + 1))
				for name, x in candidates(self.p(d), d).items()}


def visitor(tag):
	el = real.elements(tag)
	orbits = planets.Orbits.at(real.planets_at(el["epoch"]), el["epoch"])
	v = math.sqrt(planets.MU / abs(el["a"]))
	R = planets.rotation(*(np.array([math.radians(el[k])]) for k in ("i", "om", "w")))[0]
	return el, orbits, v, R, planets.passes(el["q"], v, R[None], np.array([el["tp"]]), orbits)[0][0]


def main():
	field, _, _ = route.Field.from_catalogs("data/gaia_reach.csv", "data/hipparcos_reach.csv")
	with Pool() as pool:
		_, orbits, _, _, _ = visitor("3I")
		rng = np.random.default_rng(2026)
		tours = {}
		for v_kms, q in SETTINGS:
			v = v_kms / real.AU_DAY_KMS
			cal = Calibrated(pool, q, v, orbits)
			for n in range(TOURS):
				targets = {int(j): float(10 ** rng.uniform(math.log10(0.03), math.log10(0.5))) for j in rng.choice(8, size=rng.integers(1, 5), replace=False)}
				_, _, d = planets.plan(q, v, orbits, targets, seed=n, n=40_000, starts=5)
				for name, b in cal.bits(d).items():
					tours.setdefault(name, []).append(b)
		ceiling = math.log2(N_STAT + 1)
		print(f"random planned tours ({len(SETTINGS) * TOURS}; scores cap at {ceiling:.1f} bits)")
		print(f"  {'statistic':<34}{'mean':>7}{'median':>8}{'>= 10 bits':>12}")
		for name, b in tours.items():
			b = np.array(b)
			print(f"  {name:<34}{b.mean():7.2f}{np.median(b):8.2f}{(b >= 10).mean():11.0%}")
		for tag in ("1I", "2I", "3I"):
			el, orbits, v, R, d = visitor(tag)
			dirs = real.icrs_to_ecl(route.kinematic_directions(field, v * real.AU_DAY_KMS))
			print(f"\n{tag}, v_inf {v * real.AU_DAY_KMS:.1f} km/s, perihelion {el['q']:.2f} AU")
			for label, kw in (("random orientation, 10,000 years (adopted)", {}), ("random orientation, 1,000 years", {"span": 1000 * 365.25}),
							  ("arrivals moving as the catalog's stars do", {"directions": dirs}), ("orientation as observed, timing only", {"R": R})):
				res = Calibrated(pool, el["q"], v, orbits, **kw).bits(d)
				print(f"  {label:<44}" + "  ".join(f"{name}: {b:.1f}" for name, b in res.items()), flush=True)


if __name__ == "__main__":
	sys.exit(main())
