"""Planets and route scored together: calibration under one null visitor for both.

	.venv/bin/python together_check.py > out/together.txt      # about an hour on 18 cores

The planet null and the route's origin end both randomize the arrival direction, so the two scores could move together under chance, and combining their p's as independent would then overstate. Here each null visitor, a random orientation and perihelion time at a real visitor's speed and perihelion distance, is scored for planets and for route on the same path, and P(bits >= x) * 2^x is checked for each score and for the two combined.
"""
import math
import sys
from multiprocessing import Pool

import numpy as np

from solar_visitor_metric import planets, real, route

N = 36_000
CHUNK = 500


def asymptotes(q, v_inf, R):
	"""Arriving and leaving directions of motion (ICRS unit vectors) of hyperbolas with orientations R (n, 3, 3)."""
	e = 1 + q * v_inf ** 2 / planets.MU
	s = math.sqrt(1 - 1 / e ** 2)
	u_in = np.einsum("nab,b->na", R, np.array([1 / e, s, 0.0]))
	u_out = np.einsum("nab,b->na", R, np.array([-1 / e, s, 0.0]))
	return np.array([real.ecl_to_icrs(u) for u in u_in]), np.array([real.ecl_to_icrs(u) for u in u_out])


def _route(args):
	ends, V_in, V_out, seed = args
	(f_in, e_in, w_in), (f_out, e_out, w_out) = ends
	rng = np.random.default_rng(seed)
	out = np.empty(len(V_in))
	for k in range(len(V_in)):
		p_in, _, _ = route.end_score(f_in, V_in[k], -1, e_in, rng, weights=w_in)
		p_out, _, _ = route.turn_score(f_out, V_in[k], V_out[k], e_out, rng, weights=w_out)
		out[k] = route.combine(p_in, p_out)
	return out


def calibration(bits):
	return "  ".join(f"{x}: {(bits >= x).mean() * 2 ** x:.2f}" for x in (1, 2, 4, 6, 8))


def main():
	field, errors, _ = route.Field.from_catalogs("data/gaia_reach.csv", "data/hipparcos_reach.csv")
	with Pool() as pool:
		for tag in ("1I", "2I", "3I"):
			el = real.elements(tag)
			orbits = planets.Orbits.at(real.planets_at(el["epoch"]), el["epoch"])
			v = math.sqrt(planets.MU / abs(el["a"]))
			null = planets.Null(el["q"], v, orbits, pool.map)
			rng = np.random.default_rng(11)
			R = planets.random_rotations(rng, N)
			tp = orbits.t0 + (rng.random(N) - 0.5) * planets.SPAN
			t = planets.statistic(null.p(planets.passes(el["q"], v, R, tp, orbits)[0]))
			p_planets = (np.searchsorted(null.T, t, side="right") + 1) / (len(null.T) + 1)
			v_kms = v * real.AU_DAY_KMS
			u_in, u_out = asymptotes(el["q"], v, R)
			ends = []
			for sign in (-1, 1):
				w = route.end_weights(field, v_kms, sign)
				idx = np.flatnonzero(w)
				ends.append((field.take(idx), errors.take(idx), w[idx]))
			jobs = [(ends, v_kms * u_in[k:k + CHUNK], v_kms * u_out[k:k + CHUNK], (12, k)) for k in range(0, N, CHUNK)]
			p_route = np.concatenate(pool.map(_route, jobs))
			bits = lambda p: -np.log2(np.minimum(p, 1.0))
			together = np.array([route.combine(a, b) for a, b in zip(p_planets, p_route)])
			print(f"{tag} ({N:,} null visitors at {v_kms:.1f} km/s, perihelion {el['q']:.2f} AU), P(bits >= x) * 2^x", flush=True)
			print(f"  planets   {calibration(bits(p_planets))}")
			print(f"  route     {calibration(bits(p_route))}")
			print(f"  together  {calibration(bits(together))}   correlation of log p's {np.corrcoef(np.log(p_planets), np.log(np.minimum(p_route, 1)))[0, 1]:+.3f}", flush=True)


if __name__ == "__main__":
	sys.exit(main())
