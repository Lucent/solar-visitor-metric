"""Null-model encounter-chain rates for in-plane interstellar visitors.

	.venv/bin/python run.py            # tables + figures into out/
"""
import math
import sys
from multiprocessing import Pool

import numpy as np

from solar_visitor_metric.chain import _one, dv_criterion, q_to_b, soi_criterion, sweep
from solar_visitor_metric.constants import PLANETS, VISITORS, kms
from solar_visitor_metric.flyby import Crossing, crossings, local_velocity, state_after
from solar_visitor_metric import nbody, svg

OUT = "out"


def bits(p):
	return -math.log2(p) if p > 0 else math.inf


def report(name, v_kms, crit, pool):
	rows = sweep(kms(v_kms), crit, mapper=pool.map)
	tot_db = sum(r["db"] for r in rows)
	tot = sum(r["dist"] for r in rows)
	print(f"\n{name}  v_inf = {v_kms} km/s   [{crit.label}]")
	print(f"  {'perihelion ring':<22}{'share':>8}{'P(>=1)':>11}{'P(>=2)':>11}{'P(>=3)':>11}"
		  f"{'capture':>10}")
	edges = ["Sun"] + [p.name for p in PLANETS]
	for r in rows:
		f = r["dist"] / r["db"]
		print(f"  {edges[r['ring']] + '-' + edges[r['ring'] + 1]:<22}"
			  f"{r['db'] / tot_db:>8.1%}{f[1:].sum():>11.2e}{f[2:].sum():>11.2e}"
			  f"{f[3:].sum():>11.2e}{r['capture'] / r['db']:>10.1e}")
	f = tot / tot_db
	p2 = f[2:].sum()
	print(f"  {'all q < Neptune':<22}{1:>8.0%}{f[1:].sum():>11.2e}{p2:>11.2e}"
		  f"{f[3:].sum():>11.2e}   = {bits(p2):.1f} bits for a chain")
	pairs = {}
	for r in rows:
		for k, val in r["pairs"].items():
			pairs[k] = pairs.get(k, 0) + val / tot_db
	top = sorted(pairs.items(), key=lambda kv: -kv[1])[:5]
	print("  top chains: " + ", ".join(f"{PLANETS[i].name[:3]}->{PLANETS[j].name[:3]} {p:.1e}"
									   for (i, j), p in top))
	return f


def profile(name, v_kms, crit, pool, n=90):
	"""P(>=1), P(>=2) vs perihelion, prograde then retrograde."""
	v = kms(v_kms)
	qs = np.geomspace(0.05, 35, n)
	cols = []
	for sgn in (+1, -1):
		outs = pool.map(_one, [(v, sgn * q_to_b(q, v), crit, {}) for q in qs])
		cols += [[o.at_least(1) for o in outs], [o.at_least(2) for o in outs]]
	return (name, v_kms, qs, *cols)


def find_chain(v, b, i, j, crit):
	"""Pick a qualifying i-flyby whose deflected path can qualify at j."""
	E, h = 0.5 * v * v, b * v
	seq, _ = crossings(E, h, math.inf, -1)
	for pi, sg in seq:
		if pi != i:
			continue
		cr = Crossing(i, *local_velocity(E, h, PLANETS[i].a, sg))
		lo, hi = cr.b_surface(), crit(cr)
		for bp in sorted(np.geomspace(lo * 1.05, hi * 0.95, 40), key=lambda x: -x):
			for s in (1, -1):
				E2, h2, sg2, _ = state_after(cr, s * bp)
				seq2, _ = crossings(E2, h2, PLANETS[i].a, sg2, {i})
				for pj, sgj in seq2:
					if pj != j:
						continue
					crj = Crossing(j, *local_velocity(E2, h2, PLANETS[j].a, sgj))
					bcj = crit(crj)
					if bcj > 2 * crj.b_surface():
						return sg, s * bp, 0.5 * bcj
	return None


def main():
	dv10 = dv_criterion(0.10)
	soi = soi_criterion(1.0)
	with Pool() as pool:
		for crit in (dv10, soi):
			for name, v in VISITORS.items():
				report(name, v, crit, pool)
		for crit, tag in ((dv10, "dv10"), (soi, "soi")):
			svg.profile_figure(f"{OUT}/profile_{tag}.svg",
							   [profile(n, v, crit, pool) for n, v in VISITORS.items()],
							   crit.label)

	# one targeted chain, flown in REBOUND with and without the first planet
	v, b, i, j = kms(26.3), 2.0, 4, 5
	sg, bp_i, bp_j = find_chain(v, b, i, j, dv10)
	tracks = {}
	res = nbody.check_chain(v, b, i, sg, bp_i, j, bp_j, dv10, tracks=tracks)
	ci, cj = res["chain"][i], res["chain"][j]
	print(f"\nChain check {PLANETS[i].name}->{PLANETS[j].name}, v=26.3 km/s, b={b} AU (REBOUND):")
	print(f"  {PLANETS[i].name}: b_p {ci['b_meas']:+.5f} AU (target {bp_i:+.5f}), "
		  f"|dv|/v {res['dv_i']:.3f}")
	print(f"  {PLANETS[j].name}: b_p {cj['b_meas']:+.5f} AU (target {bp_j:+.5f}, "
		  f"qualifies below {res['b_crit_j']:.5f}), |dv|/v {res['dv_j']:.3f}")
	print(f"  removal ({PLANETS[i].name} massless): {PLANETS[j].name} r_min "
		  f"{res['removed'][j]['rmin']:.3f} AU -> "
		  f"{'still qualifies' if res['j_survives'] else 'chain broken (dependent)'}")
	svg.chain_figure(f"{OUT}/chain_jupiter_saturn.svg", res, tracks, i, j, v, b,
					 "Jupiter → Saturn chain, 1I-like speed, with removal test")
	print(f"\nFigures in {OUT}/")


if __name__ == "__main__":
	sys.exit(main())
