"""How steeply route look-elsewhere weights should fall with hop, and below which hop they stop rising: power on planted routes, and the visitors under every choice.

	.venv/bin/python falloff_check.py > out/falloff.txt      # about 15 minutes on 18 cores

Weights fall as hop^-k with hops under a floor counted as the floor, for each k in FALLOFFS and floor in FLOORS. Planted routes, from the nearest stars a visitor can meet out to the 300th, show which pair detects a traveler choosing near stops; the real visitors are scored under every pair, and their best over all of them is reported less log2 of the number tried, the strongest valid reading of the probe hypothesis. The pair is chosen by power on planted routes, never by the visitors' scores.
"""
import math
import sys
from multiprocessing import Pool

import numpy as np

from solar_visitor_metric import real, route

FALLOFFS = (0.0, 2.0, 3.0, 4.0, 5.0, 6.0)
FLOORS = (0.5, 1.0, 2.0, 3.0)
PAIRS = [(k, f) for k in FALLOFFS for f in FLOORS]
DRAWS = 31
SPEED = 26.3
RANKS = ((1, 2), (2, 3), (3, 3), (5, 5), (8, 8), (13, 13), (20, 20), (30, 30), (50, 50), (100, 100), (300, 300))


def weights(field, v, sign, k, floor):
	L = route.hops(field, v, sign)
	w = np.where(L <= route.HOP, np.maximum(L, floor) ** -k, 0.0)
	return w / w.sum()


def bits(field, errors, V_in, V_out, pair, rng=None):
	v = float(np.linalg.norm(V_in))
	p_in, _, _ = route.end_score(field, V_in, -1, errors, rng, weights=weights(field, v, -1, *pair))
	p_out, _, _ = route.turn_score(field, V_in, V_out, errors, rng, weights=weights(field, v, +1, *pair))
	return [max(0.0, -math.log2(min(1.0, p))) for p in (p_in, p_out, route.combine(p_in, p_out))]


def job(args):
	label, V_in, V_out, planted = args
	field, errors, _ = route.Field.from_catalogs("data/gaia_reach.csv", "data/hipparcos_reach.csv")
	rng = np.random.default_rng(0) if planted else None
	return label, [np.median([bits(field, errors, V_in, V_out, pair, rng) for _ in range(DRAWS)], axis=0) if planted
				   else bits(field, errors, V_in, V_out, pair) for pair in PAIRS]


def main():
	field, errors, _ = route.Field.from_catalogs("data/gaia_reach.csv", "data/hipparcos_reach.csv")
	jobs = [(f"{tag} (real path)", *real.asymptotes(tag), False) for tag in ("1I", "2I", "3I")]
	for a, b in RANKS:
		(i, ri), = route.reachable(field, (a,), SPEED, -1)
		# the first destination at or after rank b that a Sun flyby can turn the visitor to
		for rank in range(b, b + 200):
			(j, rj), = route.reachable(field, (rank,), SPEED, +1, {i})
			plan = route.plan_route(field, i, j, SPEED)
			if plan:
				break
		jobs.append((f"planted, hop ranks {ri} -> {rj} at {SPEED:g} km/s", plan[0], plan[1], True))
	ids = list(field.ids)
	plan = route.plan_route(field, ids.index("HIP 32349"), ids.index("HIP 71683"), 58.0)
	jobs.append(("mission: Sirius -> α Cen A at 58 km/s", plan[0], plan[1], True))
	with Pool(len(jobs)) as pool:
		results = dict(pool.imap(job, jobs))
	planted = [label for label in results if label.startswith("planted")]
	print("combined route bits (origin / destination) for hop^-k with hops under the floor counted as the floor")
	print(f"{'k':>4}{'floor':>7}{'mean planted':>14}{'ranks <= 30':>13}{'mission':>9}   " + "   ".join(f"{t} (real path)" for t in ("1I", "2I", "3I")))
	for n, (k, f) in enumerate(PAIRS):
		row = lambda label: results[label][n]
		visitors = "   ".join(f"{row(f'{t} (real path)')[2]:4.1f} ({row(f'{t} (real path)')[0]:.1f}/{row(f'{t} (real path)')[1]:.1f})" for t in ("1I", "2I", "3I"))
		print(f"{k:>4g}{f:>7g}{np.mean([row(l)[2] for l in planted]):>14.2f}{np.mean([row(l)[2] for l in planted[:8]]):>13.2f}"
			  f"{row('mission: Sirius -> α Cen A at 58 km/s')[2]:>9.1f}   {visitors}")
	for t in ("1I", "2I", "3I"):
		best = max(r[2] for r in results[f"{t} (real path)"])
		print(f"{t}: best combined over all {len(PAIRS)} pairs {best:.1f} bits, less log2 {len(PAIRS)}: {max(0.0, best - math.log2(len(PAIRS))):.1f}")


if __name__ == "__main__":
	sys.exit(main())
