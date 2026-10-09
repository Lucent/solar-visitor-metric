"""Visit and assist scores for example planar paths, against the null.

	.venv/bin/python score.py > out/scores.txt
"""
import json
import math
import random
import sys
from multiprocessing import Pool

import numpy as np

from solar_visitor_metric import nbody, scores, svg
from solar_visitor_metric.chain import dv_criterion
from solar_visitor_metric.constants import PLANETS, VISITORS, kms
from solar_visitor_metric.flyby import Crossing, local_velocity

OUT = "out"
FLYBYS = f"{OUT}/flybys.json"   # read by real_visitors.py for the site's missions
# Each planned flyby by its pass distances, as the site names it.
NAMES = {"Earth + Mars visits, Jupiter assist": "Earth and Mars at 1.5 million km, then a Jupiter assist",
			  "Jupiter assist only": "Jupiter assist alone",
			  "Earth + Mars visits only": "Earth and Mars at 1.5 million km, no assist",
			  "Close Earth + Mars, Jupiter assist": "Earth and Mars at 150,000 km, then a Jupiter assist"}
N_NULL = 4_000_000
N_FILL = 201
EARTH, MARS, JUPITER = 2, 3, 4


def jupiter_bp(v, b, sign=+1, frac=0.5):
	"""Impact parameter giving well over a 10% |dv|/v at Jupiter."""
	E, h = 0.5 * v * v, b * v
	cr = Crossing(JUPITER, *local_velocity(E, h, PLANETS[JUPITER].a, sign))
	return frac * dv_criterion(0.10)(cr)


def examples():
	out = []
	for v_kms in (26.3, 58.0):
		v, b = kms(v_kms), 0.8
		jb = jupiter_bp(v, b)
		tag = f"{v_kms:g} km/s"
		out += [
			(f"Earth + Mars visits, Jupiter assist ({tag})", v_kms, b,
			 {(EARTH, 1): 0.01, (MARS, 1): 0.01, (JUPITER, 1): jb}),
			(f"Jupiter assist only ({tag})", v_kms, b, {(JUPITER, 1): jb}),
		]
	v = kms(26.3)
	out += [
		("Earth + Mars visits only (26.3 km/s)", 26.3, 0.8, {(EARTH, 1): 0.01, (MARS, 1): 0.01}),
		("Close Earth + Mars, Jupiter assist (26.3)", 26.3, 0.8,
		 {(EARTH, 1): 0.001, (MARS, 1): 0.001, (JUPITER, 1): jupiter_bp(v, 0.8)}),
	]
	return out


def main():
	with Pool() as pool:
		nulls = {v: np.sort(scores.null_parallel(kms(v), N_NULL, pool)[:, 0])
				 for v in VISITORS.values()}
		cal = scores.null_parallel(kms(26.3), 720_000, pool, seed=10_000, visits=True)[:, 2]

	def assist_bits(v_kms, a):
		arr = nulls[v_kms] if v_kms in nulls else nulls[min(nulls, key=lambda u: abs(u - v_kms))]
		k = len(arr) - np.searchsorted(arr, a, side="left")
		if k == 0:
			return math.log2(len(arr)), True
		return -math.log2(k / len(arr)), False

	print(f"Null: {N_NULL:,} paths per speed. Unplanned crossings are random; "
		  f"scores are medians over {N_FILL} draws.\n")
	print(f"{'example':<46}{'visit':>7}{'assist':>9}{'a-bits':>8}{'excess':>8}{'total':>8}")
	panels, marks, flybys = [], [], []
	for title, v_kms, b, plan in examples():
		draws = [scores.fly(kms(v_kms), b, plan=plan, rng=random.Random(s))
				 for s in range(N_FILL)]
		path = draws[0]
		vb = float(np.median([d.visit_bits() for d in draws]))
		a = float(np.median([d.assist() for d in draws]))
		ab, floor = assist_bits(v_kms, a) if a > 0.005 else (0.0, False)
		xb = float(np.median([d.excess_visit_bits() for d in draws]))
		ge = "≥" if floor else ""
		print(f"{title:<46}{vb:>7.1f}{a:>9.3f}{ge + format(ab, '.1f'):>8}{xb:>8.1f}"
			  f"{ge + format(ab + xb, '.1f'):>8}")
		lines = [f"visit score: {vb:.1f} bits",
				 f"assist: {a:.3f} of v∞ → {ge}{ab:.1f} bits",
				 f"visits beyond the assist: {xb:.1f} bits",
				 f"combined: {ge}{ab + xb:.1f} bits"]
		for (i, sg), bp in plan.items():
			lines.append(f"{PLANETS[i].name}: b_p {bp * 1.495978707e8 / 1e3:,.0f} thousand km")
		panels.append((title, path, lines))
		flybys.append({"name": NAMES[title.split(" (")[0]], "speed": v_kms, "bits": ab + xb, "floor": floor})
		if "26.3" in title and a > 0.005:
			marks.append((title.split(" (")[0], a))

	with open(FLYBYS, "w") as f:
		json.dump(flybys, f, ensure_ascii=False, indent="\t")

	# random null paths for contrast
	rng = random.Random(7)
	v = kms(26.3)
	for n in range(2):
		b = rng.uniform(-1.5, 1.5)
		path = scores.fly(v, b, rng=rng)
		lines = [f"visit score: {path.visit_bits():.1f} bits", f"assist: {path.assist():.4f} of v∞"]
		panels.append((f"Chance path {n + 1} (26.3 km/s)", path, lines))

	bits = np.arange(0, 15.5, 0.5)
	ratio = [(cal >= x).mean() * 2 ** x for x in bits]
	svg.tails_figure(f"{OUT}/null_tails.svg",
					 [(n, v, nulls[v]) for n, v in VISITORS.items()], (bits, ratio), marks[:3])
	svg.paths_figure(f"{OUT}/scored_paths.svg", panels)

	# N-body confirmation of the first example
	title, v_kms, b, plan = examples()[0]
	order = [(i, s, bp) for (i, s), bp in plan.items()]
	enc, a_nb = nbody.fly_plan(kms(v_kms), b, order)
	a0 = scores.asymptote(0.5 * kms(v_kms) ** 2, b * kms(v_kms), 0.0)
	a_pc = scores.fly(kms(v_kms), b, plan=plan, default_p=0.99).assist()
	print(f"\nREBOUND check of '{title}':")
	for i, s, bp in order:
		print(f"  {PLANETS[i].name}: b_p target {bp:+.5f} AU, flown {enc[i]['b_meas']:+.5f} AU")
	print(f"  assist: patched conic {a_pc:.4f}, N-body {math.hypot(a_nb[0] - a0[0], a_nb[1] - a0[1]) / kms(v_kms):.4f}")


if __name__ == "__main__":
	sys.exit(main())
