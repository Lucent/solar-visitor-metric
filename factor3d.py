"""Chance a random-phase circular Jupiter comes within B of a hyperbola: 3D isotropic orientation vs in-plane. Units AU, yr/2pi (GM_sun = 1).

	.venv/bin/python factor3d.py
"""
import math
from multiprocessing import Pool

import numpy as np

from solar_visitor_metric.constants import PLANETS, kms

A = PLANETS[4].a
N_J = A ** -1.5
B = 0.016


def hit_fraction(v, q, pole, omega):
	"""Measure of Jupiter phases within B of the path (union over time)."""
	E = 0.5 * v * v
	a_h = 1 / (2 * E)
	e = 1 + q / a_h
	# sample hyperbolic anomaly finely where r is within B of A
	F_A = math.acosh((1 + A / a_h) / e)
	out = []
	for sgn in (-1, 1):
		F = sgn * F_A + np.linspace(-1, 1, 4001) * 0.02 * F_A
		x = a_h * (e - np.cosh(F))
		y = a_h * math.sqrt(e * e - 1) * np.sinh(F)
		t = a_h ** 1.5 * (e * np.sinh(F) - F)
		# orbit-plane frame -> 3D: perihelion along p, motion along q, pole n
		n = pole
		tmp = np.array([1.0, 0, 0]) if abs(n[0]) < 0.9 else np.array([0, 1.0, 0])
		p = np.cross(n, tmp); p /= np.linalg.norm(p)
		qv = np.cross(n, p)
		p, qv = math.cos(omega) * p + math.sin(omega) * qv, -math.sin(omega) * p + math.cos(omega) * qv
		P = x[:, None] * p + y[:, None] * qv
		rxy = np.hypot(P[:, 0], P[:, 1])
		z = P[:, 2]
		th = np.arctan2(P[:, 1], P[:, 0])
		# Jupiter at angle phi0 + N_J t: distance^2 = (rxy-A)^2 + z^2 + 2 A rxy (1-cos d)
		rem = B * B - (rxy - A) ** 2 - z * z
		ok = rem > 0
		if not ok.any():
			continue
		half = np.arccos(np.clip(1 - rem[ok] / (2 * A * rxy[ok]), -1, 1))
		center = th[ok] - N_J * t[ok]
		out += list(zip(center - half, center + half))
	if not out:
		return 0.0
	# union of angular intervals (mod 2 pi)
	iv = sorted(((lo % (2 * math.pi), (lo % (2 * math.pi)) + (hi - lo)) for lo, hi in out))
	tot, cur_lo, cur_hi = 0.0, None, None
	for lo, hi in iv:
		if cur_hi is None or lo > cur_hi:
			if cur_hi is not None:
				tot += cur_hi - cur_lo
			cur_lo, cur_hi = lo, hi
		else:
			cur_hi = max(cur_hi, hi)
	tot += cur_hi - cur_lo
	return min(1.0, tot / (2 * math.pi))


def chunk(seed):
	rng = np.random.default_rng(seed)
	v = kms(26.3)
	res = []
	for _ in range(400):
		q = A * rng.random()   # both cases use the same q distribution
		om = rng.uniform(0, 2 * math.pi)
		n3 = rng.normal(size=3); n3 /= np.linalg.norm(n3)
		n2 = np.array([0, 0, 1.0 if rng.random() < 0.5 else -1.0])
		res.append((hit_fraction(v, q, n3, om), hit_fraction(v, q, n2, om)))
	return res

if __name__ == "__main__":
	with Pool() as pool:
		r = np.array([x for c in pool.map(chunk, range(72)) for x in c])
	p3, p2 = r[:, 0].mean(), r[:, 1].mean()
	print(f"mean P(within {B} AU of Jupiter): 3D {p3:.3e}, planar {p2:.3e}, ratio {p3 / p2:.2e}; B/a = {B / A:.2e}")
