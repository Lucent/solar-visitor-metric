"""Route score on a synthetic star field: null calibration and planted routes.

	.venv/bin/python route_check.py > out/route.txt
"""
import math
import sys
from dataclasses import replace
from multiprocessing import Pool

import numpy as np

from solar_visitor_metric import route, svg
from solar_visitor_metric.constants import VISITORS

OUT = "out"
N_NULL = 72_000
N_NOISE_DRAWS = 31


def planted(field, v, pair, errors):
	(i, ri), (j, rj) = pair
	plan = route.plan_route(field, i, j, v)
	if plan is None:
		return None
	V_in, V_out, q = plan
	return V_in, V_out, q, route.planted_bits(V_in, V_out, field, errors, N_NOISE_DRAWS)


def main():
	rng = np.random.default_rng(1)
	field = route.Field.synthetic(rng)
	err = route.Errors()
	print(f"Synthetic field: {len(field.D)} stars that could meet a visitor faster than "
		  f"{route.V_FLOOR:g} km/s within {route.HOP:g} pc of the Sun ({(field.D < route.HOP).sum()} "
		  f"within {route.HOP:g} pc today); errors: star v {err.star_v} km/s, "
		  f"parallax {err.parallax_mas} mas, visitor direction {err.visitor_dir:g} rad\n")

	bits = np.arange(0, 14.5, 0.5)
	panels = []
	nulls = {}
	with Pool() as pool:
		for name, v in VISITORS.items():
			nr = route.null_route_parallel(field, v, N_NULL, err, pool)
			b = -np.log2(np.maximum(nr, 1e-300))
			nulls[v] = np.sort(b[:, 2])
			series = [(lab, bits, [(b[:, c] >= x).mean() * 2 ** x for x in bits])
					  for c, lab in ((0, "from"), (1, "to"), (2, "combined"))]
			panels.append((f"{name} ({v:g} km/s)", series))
			worst = max(max(r for x, r in zip(bits, s[2]) if x <= 10) for s in series)
			print(f"{name:<12} null: median combined {np.median(b[:, 2]):.2f} bits, "
				  f"P(combined >= 10 bits) = {(b[:, 2] >= 10).mean():.1e} "
				  f"(exact would be {2 ** -10:.1e}); worst ratio to 10 bits {worst:.2f}x")
	svg.calibration_figure(f"{OUT}/route_null.svg", "Route score calibration (synthetic field)",
						   panels, n_samples=N_NULL)

	def calibrated(v, raw):
		"""-log2 of the null fraction scoring at least `raw` (MC, +1 smoothed)."""
		arr = nulls[v]
		k = len(arr) - np.searchsorted(arr, raw - 1e-9)
		return -math.log2((k + 1) / (len(arr) + 1)), k == 0

	print("\nPlanted routes (visitor passes exactly through both stars; median over noise).")
	print("'comb' is the conservative bound; 'calib' is calibrated on the null sample.")
	print(f"  {'case':<44}{'v':>6}{'from':>7}{'to':>7}{'comb':>7}{'calib':>8}{'q (AU)':>9}")
	map_done = False
	for name, v in VISITORS.items():
		first = route.reachable(field, (1,), v, -1)
		for sign, lab in ((-1, "from"), (+1, "to")):
			L = route.hops(field, v, sign)
			ok = L <= route.HOP
			print(f"  {name} {lab}: {ok.sum()} candidate stars, {(ok & (field.D > route.HOP)).sum()} "
				  f"of them beyond {route.HOP:g} pc today (farthest {field.D[ok].max():.0f} pc)")
		cases = {
			"near -> near": first + route.reachable(field, (2,), v, +1, {first[0][0]}),
			"near -> mid (rank ~100)": first + route.reachable(field, (100,), v, +1),
			"far -> far (rank ~1000, ~2000)": route.reachable(field, (1000,), v, -1)
			+ route.reachable(field, (2000,), v, +1),
		}
		for case, pair in cases.items():
			if len(pair) < 2 or pair[0][0] == pair[1][0]:
				print(f"  {case:<40}{v:>6g}   (no reachable pair)")
				continue
			res = planted(field, v, pair, err)
			if res is None:
				print(f"  {case:<40}{v:>6g}   (turn needs perihelion inside the Sun)")
				continue
			V_in, V_out, q, b = res
			label = f"{case} [ranks {pair[0][1]}, {pair[1][1]}]"
			cb, floor = calibrated(v, b[2])
			cal = ("≥" if floor else "") + f"{cb:.1f}"
			print(f"  {label:<44}{v:>6g}{b[0]:>7.1f}{b[1]:>7.1f}{b[2]:>7.1f}{cal:>8}{q:>9.2f}")
			if not map_done and case == "near -> near":
				_, t_in = route.aim_at(field.pos[pair[0][0]], field.vel[pair[0][0]], v, -1)
				_, t_out = route.aim_at(field.pos[pair[1][0]], field.vel[pair[1][0]], v, +1)
				svg.route_map(f"{OUT}/route_map.svg", field, V_in, V_out,
							  [(pair[0][0], t_in, f"origin (rank {pair[0][1]})"),
							   (pair[1][0], t_out, f"destination (rank {pair[1][1]})")],
							  f"Planted route at {v:g} km/s",
							  [f"{name}-speed visitor, perihelion {q:.2f} AU",
							   f"from: {b[0]:.1f} bits, {field.D[pair[0][0]]:.2f} pc away now, "
							   f"met {abs(t_in):.2f} Myr ago",
							   f"to: {b[1]:.1f} bits, {field.D[pair[1][0]]:.2f} pc away now, "
							   f"reached in {t_out:.2f} Myr",
							   f"combined: {b[2]:.1f} bits"])
				map_done = True

	print("\nPrecision sweep, near -> near at 26.3 km/s (star velocity error dominates):")
	v = 26.3
	first = route.reachable(field, (1,), v, -1)
	pair = first + route.reachable(field, (2,), v, +1, {first[0][0]})
	for sv in (3.0, 1.0, 0.3, 0.1, 0.03):
		res = planted(field, v, pair, replace(err, star_v=sv))
		print(f"  star velocity error {sv:>5g} km/s: combined {res[3][2]:.1f} bits")


if __name__ == "__main__":
	sys.exit(main())
